import numpy as np

from grain_growth_model.core.geometry import compute_min_distances_between_lines
from grain_growth_model.core.meltpool import compute_normalized_gradient_at_point


def split_seeds_into_batches(
        size_limit:int, overlap_ratio:float, simulation_length:float, 
        coordinates_seeds:np.array, euler_angles:np.array, gradients:np.array, grain_ids:np.array
        ):
    """
    Split seeds into spatial batches to handle memory constraints. Split bead and assume there is 
    no conflicts when initial coordinates are far away from each other. A overlap area must be 
    considered to deal with conflicts at the boundaries of the sub-domain.

    Parameters
    ----------
        size_limits (int): limit of elements per batch.
        overlap_ratio (float): percentage of length for overlap checking.
        simulation_length (float): length of the simulation.
        coordinates_seeds (ndarray) (N, 3): Starting positions of grains.
        euler_angles (ndarray) (N, 3): Crystal orientations.
        gradients (ndarray) (N, 3): Initial thermal gradient vectors at seed points.
        grain_ids (ndarray) (N,): Unique grain identifiers.

    Returns
    -------
        batch_size (float): Length of each spatial batch.
        batch_dict (dict): Dictionary containing seed data for each batch.
    """

    total = len(coordinates_seeds)
    num_batches = int(np.ceil(total / size_limit))
    batch_size = simulation_length / num_batches

    results = {}
    for i in range(num_batches):
        in_zone = (
            (coordinates_seeds[:, 0] < batch_size * (i + 1 + overlap_ratio)) &
            (coordinates_seeds[:, 0] > batch_size * (i - overlap_ratio))
        )
        results[i] = {
            'id': grain_ids[in_zone],
            'coordinates': coordinates_seeds[in_zone],
            'orientations': euler_angles[in_zone],
            'gradient interface': gradients[in_zone],
        }

    return batch_size, results


def resolve_growth_conflicts(
        coordinates_seeds: np.ndarray,
        growth_directions: np.ndarray,
        global_id: np.ndarray,
        thermal_layers: dict, 
        layer_index: int, 
        meltpool_z:float,
        threshold: float,
        mean_grain_diameter:float,
        bd_increment: float,
        simulation_width: float,
        z_ultimate: float,
        ) -> np.ndarray:
    """
    Perform pairwise conflict resolution based on minimal distance and alignment to thermal gradients.

    Parameters
    ----------
        coordinates_seeds (ndarray) (N, 3): Starting positions of grains.
        growth_directions (ndarray) (N, 3): Easy growth direction vectors (OOI).
        global_id (ndarray) (N,): Global grain id of the simulation for grains involved in this batch.
        thermal_layers (dict): Full thermal history.
        layer_index (int): Index of the current layer.
        meltpool_z (float): Actual height of the z-coordinates of the ellipse center.
        threshold (float): Normalized distance for conflict detection among minimal distances between two trajectories.
        mean_grain_diameter (float): Mean grain size for substrate initialization.
        bd_increment (float): building direction increment (constant all along the wall)
        simulation_width (float): width of the simulation.
        z_ultimate (float): Top limit of the simulation.

    Returns
    -------
        final_lengths (ndarray): Length of each surviving dendrite.
        alive (ndarray): Boolean array of grains that did not lose during the competitive growth phase.
        grain_conflict_history (dict): History of grains conflict.
    """

    # Tracking of conflict history.
    grain_conflict_history = {
        int(id_):{
            "opponent":[],
            "coordinates conflict opponent": [],
            "direction conflict opponent": [],
            "gradient conflict opponent": [],
            "coordinates conflict": [],
            "gradient conflict": [],
        } for id_ in global_id      
    }

    alive = np.ones(len(coordinates_seeds), dtype=bool)
    dendrite_lengths = np.ones(len(coordinates_seeds)) * 1e4

    # Minimal distance between conflicts and \
    # parameter 'l' in 'vect(x) + l*vect(n)' to reach conflict

    dist_matrix, lambda_matrix = compute_min_distances_between_lines(
        coordinates_seeds=coordinates_seeds, 
        growth_directions=growth_directions, 
        mean_grain_diameter=mean_grain_diameter
    )

    lamb_matrix_copy = np.copy(lambda_matrix)
    # Valid conflicts when there are both inside the melt pool footprint \
    # and a minimal distance below the threshold
    valid = (dist_matrix < threshold) & (lambda_matrix > 0) & (lambda_matrix.T > 0)
    lambda_matrix[~valid] = np.nan

    while True:

        # Among the valid conflict, focus on conflict where both grains are involved for the first time.
        with np.errstate(invalid='ignore'):
            target = lambda_matrix / np.nanmin(lambda_matrix, axis=1)[:, None]
        i_indices, j_indices = np.where((target * target.T == 1))

        if len(i_indices) == 0:
            break

        for i, j in zip(i_indices, j_indices):

            # Check if they're both still in the competitive growth stage
            if alive[i] and alive[j]:

                # Coordinates of the conflicts.
                pi = coordinates_seeds[i] + lambda_matrix[i, j] * growth_directions[i]
                pj = coordinates_seeds[j] + lambda_matrix[j, i] * growth_directions[j]

                # Check if both are in the melt pool foot print.
                in_domain_i = _in_meltpool(pi[1], pi[2], thermal_layers, layer_index, bd_increment, simulation_width, z_ultimate)
                in_domain_j = _in_meltpool(pj[1], pj[2], thermal_layers, layer_index, bd_increment, simulation_width, z_ultimate)

                if in_domain_i and in_domain_j:

                    # Compute local thermal gradient at the conflict coordinates
                    Gi = compute_normalized_gradient_at_point(pi[1], pi[2], thermal_layers[layer_index], meltpool_z)[0]
                    Gj = compute_normalized_gradient_at_point(pj[1], pj[2], thermal_layers[layer_index], meltpool_z)[0]

                    # Most favorable dendrite oriented (Walton and Chalmers theory)
                    si = np.dot(growth_directions[i], Gi)
                    sj = np.dot(growth_directions[j], Gj)

                    if si >= sj:  # Index 'i' better aligned with thermal gradient: 'j' stops growing
                        alive[j] = False
                        dendrite_lengths[j] = lamb_matrix_copy[j, i]
                        lambda_matrix[j, :] = np.nan
                        lambda_matrix[:, j] = np.nan
                    else:
                        alive[i] = False
                        dendrite_lengths[i] = lamb_matrix_copy[i, j]
                        lambda_matrix[i, :] = np.nan
                        lambda_matrix[:, i] = np.nan

                    # data tracking 
                    # id_grain[i] make te correspondence between local id and global id of the whole simulation
                    grain_conflict_history[global_id[i]]["opponent"].append(int(global_id[j]))
                    grain_conflict_history[global_id[i]]["coordinates conflict opponent"].append((pj))
                    grain_conflict_history[global_id[i]]["direction conflict opponent"].append(growth_directions[j].tolist())
                    grain_conflict_history[global_id[i]]["gradient conflict opponent"].append(Gj.tolist())
                    grain_conflict_history[global_id[i]]["coordinates conflict"].append((pi))
                    grain_conflict_history[global_id[i]]["gradient conflict"].append(Gi.tolist())
                    #
                    grain_conflict_history[global_id[j]]["opponent"].append(int(global_id[i]))
                    grain_conflict_history[global_id[j]]["coordinates conflict opponent"].append((pi))
                    grain_conflict_history[global_id[j]]["direction conflict opponent"].append(growth_directions[i].tolist())
                    grain_conflict_history[global_id[j]]["gradient conflict opponent"].append(Gi.tolist())
                    grain_conflict_history[global_id[j]]["coordinates conflict"].append((pj))
                    grain_conflict_history[global_id[j]]["gradient conflict"].append(Gj.tolist())

                elif in_domain_i and not in_domain_j:
                    alive[j] = False
                    dendrite_lengths[j] = lamb_matrix_copy[j, i]
                    lambda_matrix[j, :] = np.nan
                    lambda_matrix[:, j] = np.nan

                elif in_domain_j and not in_domain_i:
                    alive[i] = False
                    dendrite_lengths[i] = lamb_matrix_copy[i, j]
                    lambda_matrix[i, :] = np.nan
                    lambda_matrix[:, i] = np.nan

                else:
                    alive[i] = False
                    alive[j] = False
                    dendrite_lengths[i] = lamb_matrix_copy[i, j]
                    dendrite_lengths[j] = lamb_matrix_copy[j, i]
                    lambda_matrix[i, :] = np.nan
                    lambda_matrix[:, i] = np.nan
                    lambda_matrix[j, :] = np.nan
                    lambda_matrix[:, j] = np.nan

    return dendrite_lengths, alive, grain_conflict_history


def _in_meltpool(y, z, thermal_layers: dict, layer_index: int, bd_increment:float, simulation_width:float, z_ultimate:float):
    """
    Helper to check if (y, z) lies within current meltpool domain.

    Parameters
    ----------
        y, z (float): Coordinates to evaluate
        thermal_layers (dict): Full thermal history.
        layer_index (int): Index of the current layer.
        bd_increment (float): building direction increment (constant all along the wall)
        simulation_width (float): width of the simulation, periodic-like aspect along x-axis.
        z_ultimate (float): Top limit of the simulation.

    Returns
    -------
        (bool): Inside or not the melt pool footprint.
    """

    h_inf = thermal_layers[layer_index]['height']
    w_inf = thermal_layers[layer_index]['width']
    z0_inf = layer_index * bd_increment + thermal_layers[0]['height']

    z_check_inf = z0_inf - h_inf*np.sqrt(1 - ((2*y)/w_inf)**2) < z
    y_check = (simulation_width/2 > y) & (-simulation_width/2 < y)

    if layer_index < len(thermal_layers)-1:
        z0_sup = (layer_index+1)*bd_increment + thermal_layers[0]['height']
        h_sup=thermal_layers[layer_index+1]['height']
        w_sup=thermal_layers[layer_index+1]['width']
        z_check_sup = (z0_sup - h_sup*np.sqrt(1 - ((2*y)/w_sup)**2)) > z

    else:
        z_check_sup = z < z_ultimate

    return y_check & z_check_inf & z_check_sup


def resolve_batch_overlaps(data, batch_size, coords_all, ids_all, simulation_length):
    """
    Merge overlapping batch data by keeping only the local-domain value.

    Returns
    -------
        lengths : ndarray
        growth_dirs : ndarray
    """
    result_lengths = np.full(ids_all.shape, np.nan)
    result_dirs = np.zeros((ids_all.shape[0], 3))
    x_all = coords_all[:, 0]
    boundaries = np.arange(0, simulation_length + batch_size, batch_size)

    for i, (x, gid) in enumerate(zip(x_all, ids_all)):
        area = np.argmax(boundaries - x > 0) - 1
        idx_local = np.where(data[area]['id'] == gid)[0]
        if idx_local.size == 0:
            continue
        result_lengths[i] = data[area]['length dendrites'][idx_local][0]
        result_dirs[i] = data[area]['growth direction'][idx_local][0]

    return result_lengths, result_dirs
