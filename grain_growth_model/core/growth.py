import numpy as np
from tqdm import tqdm
from os.path import join

from grain_growth_model.core.meltpool import (
    compute_normalized_gradient_at_point,
    find_meltpool_exit_point,
    find_simulation_top_exit_point,
    is_point_inside_meltpool
)


def compute_easy_growth_directions(euler_angles, thermal_gradient, crystal_directions, y_positions, simulation_width):
    """
    Compute the most favorable growth direction (e.g. <001>) based on the local thermal gradient. 
    If easy growth direction is pointing out outward, the second or third better orientation is set.

    Parameters
    ----------
        euler_angles : ndarray (N, 3), Euler-Bunge angles for each grain.
        thermal_gradient : ndarray (N, 3), Local normalized thermal gradients.
        crystal_directions : ndarray (M, 3), Possible crystal growth directions (e.g., <001> family).
        y_positions : ndarray (N,) Y-positions to handle wall constraints.
        simulation_width : float, width of the simulation (y=0 is the median plane)

    Returns
    -------
        ndarray, Best growth direction (OOI) for each grain, shape (N, 3).
    """
    N = euler_angles.shape[0]
    M = crystal_directions.shape[0]

    # alpha, beta, gamma = euler_angles.T
    # ca: cos(alpha)/ sb: sin(beta)/ ...
    ca, cb, cg = np.cos(euler_angles).T
    sa, sb, sg = np.sin(euler_angles).T

    # rotation matrix according to Euler-Bunges convention
    R = np.zeros((N, 3, 3), dtype=np.float32)
    R[:, 0, 0] = ca * cg - sa * cb * sg
    R[:, 0, 1] = -ca * sg - sa * cb * cg
    R[:, 0, 2] = sa * sb

    R[:, 1, 0] = sa * cg + ca * cb * sg
    R[:, 1, 1] = -sa * sg + ca * cb * cg
    R[:, 1, 2] = -ca * sb

    R[:, 2, 0] = sb * sg
    R[:, 2, 1] = sb * cg
    R[:, 2, 2] = cb


    rotated_directions = np.einsum('nij,mj->nmi', R, crystal_directions)
    scalar_products = np.einsum('ni,nmi->nm', thermal_gradient, rotated_directions)
    sorted_indices = np.argsort(scalar_products, axis=1)[:, ::-1]

    OOI = np.array([
        rotated_directions[i, sorted_indices[i, 0]] for i in range(N)
    ])

    # Checking grains on simulation walls to force growth inwards
    outwards_growth = (
        ((y_positions == +0.5 * simulation_width) & (OOI[:, 1] > 0)) |
        ((y_positions == -0.5 * simulation_width) & (OOI[:, 1] < 0))
    )
    if np.any(outwards_growth):
        second = np.array([
            rotated_directions[i, sorted_indices[i, 1]] for i in range(N)
        ])
        third = np.array([
            rotated_directions[i, sorted_indices[i, 2]] for i in range(N)
        ])
        OOI[outwards_growth] = second[outwards_growth]
        still_bad = (
            ((y_positions == +0.5 * simulation_width) & (OOI[:, 1] > 0)) |
            ((y_positions == -0.5 * simulation_width) & (OOI[:, 1] < 0))
        )
        OOI[still_bad] = third[still_bad]

    return OOI


def simulate_grain_growth(
    coordinates_seeds: np.ndarray,
    growth_directions: np.ndarray,
    grain_ids: np.ndarray,
    gradients: np.ndarray,
    euler_angles: np.ndarray,
    dendrite_lengths: np.ndarray,
    thermal_layers: dict,
    simulation_length:float,
    current_z: float,
    z_ultimate:float,
    inter_seeds_distance:float,
    layer_index: int,
    noise_neper:float,
    output_folder: str
):
    """
    Simulate the growth of all grains in the current layer using their preferred direction.

    Parameters
    ----------
        coordinates_seeds (ndarray) (N, 3): Starting positions of grains.
        growth_directions (ndarray) (N, 3): Easy growth direction vectors (OOI).
        grain_ids (ndarray) (N,): Unique grain identifiers.
        gradients (ndarray) (N, 3): Initial thermal gradient vectors at seed points.
        euler_angles (ndarray) (N, 3): Crystal orientations.
        dendrite_lengths (ndarray) (N,): Maximum growth length per grain.
        thermal_layers (dict): Full thermal history.
        simulation_length (float): length of the simulation, periodic-like aspect along x-axis.
        current_z (float): Current melt pool height.
        z_ultimate (float): Top limit of the simulation.
        inter_seeds_distance (float): distance between two consecutive seeds
        layer_index (int): Index of the current layer.
        noise_neper (float): Noise for tesselation construction (Neper)
        output_folder (str): Output path for coordinate, orientation, and ID saving.

    Returns
    -------
        (int): Number of grains that reached the top of the simulation domain.
    """

    reach_top = 0
    is_last = layer_index == (len(thermal_layers) - 1)

    path_pos = join(output_folder, f"coordinates/coordinates_{layer_index}.txt")
    path_ori = join(output_folder, f"orientations/orientations_{layer_index}.txt")
    path_idx = join(output_folder, f"index/index_{layer_index}.txt")

    path_top_pos = join(output_folder, f"coordinates/coordinates_{layer_index}_{layer_index+1}.txt")
    path_top_ori = join(output_folder, f"orientations/orientations_{layer_index}_{layer_index+1}.txt")
    path_top_idx = join(output_folder, f"index/index_{layer_index}_{layer_index+1}.txt")

    threshold = 1e-3
    next_layer_MP_description = thermal_layers[layer_index + 1] if not is_last else None

    with open(path_pos, 'w') as pos_f, open(path_ori, 'w') as ori_f, open(path_idx, 'w') as id_f, \
         open(path_top_pos, 'w') as top_pos_f, open(path_top_ori, 'w') as top_ori_f, open(path_top_idx, 'w') as top_idx_f:

        for i in tqdm(range(grain_ids.size), desc="Grain growth"):

            xg, yg, zg = coordinates_seeds[i]  # Initial seed coordinates.
            grad = gradients[i]  # Initial unary thermal gradient direction for epitaxial growth.
            ooi = growth_directions[i]  # Easy growth direction of the dendrite (Computed at the beginning and still the same until the end).
            g_id = grain_ids[i]  # Index of the grain, every seeds sharing the same index belong to the same grain.
            o1, o2, o3 = euler_angles[i]  # Euler-Bunges angles of the crystal orientation of the grain.
            length_max = dendrite_lengths[i]  # Length of the dendrite computed during the competitive stage of the model. The length of the grain is related to this datum.
            growth = 0  # indicator of grain growth.

            while growth < length_max:

                xn, yn, zn = (xg, yg, zg) + np.random.uniform(0, noise_neper, 3)
                pos_f.write(f"{xn} {yn} {zn}\n")
                ori_f.write(f"{o1} {o2} {o3}\n")
                id_f.write(f"{g_id}\n")

                cos_theta = np.clip(abs(grad[0]), threshold, None)  # Low threshold to help worst oriented grains starting growing
                d_push = inter_seeds_distance * cos_theta
                xg, yg, zg = (xg, yg, zg) + d_push * grad

                # Keep coordinates in the simulation windows. Periodic-like aspect along x-axis
                if xg < 0:
                    xg += simulation_length
                elif xg > simulation_length:
                    xg -= simulation_length

                # First stop condition: last layer of the growth --> flat upper limits reached of the domain
                if is_last and zg > z_ultimate:
                    stop_point = find_simulation_top_exit_point(
                        point=np.array([xg, yg, zg]),
                        direction=grad,
                        z_top=z_ultimate
                    )
                    top_pos_f.write(f"{' '.join(map(str, stop_point))}\n")
                    top_ori_f.write(f"{o1} {o2} {o3}\n")
                    top_idx_f.write(f"{g_id}\n")
                    break

                if not is_last:
                    # Check if the actual position of the seed won't be remelted by the next layer.
                    if is_point_inside_meltpool(
                        y=yg,
                        z=zg,
                        half_width=next_layer_MP_description['width'] / 2,
                        depth=next_layer_MP_description['height'],
                        yc=0,
                        zc=current_z + thermal_layers[0]['height']
                    ):
                        # If so, then computation of the intersection point between the growth direction \
                        # and the foot print of the pass of the next melt pool
                        x_stop, y_stop, z_stop = find_meltpool_exit_point(
                            point=np.array([xg, yg, zg]),
                            direction=grad,
                            meltpool=next_layer_MP_description,
                            z0=current_z + thermal_layers[0]['height']
                        )
                        if x_stop > simulation_length:
                            x_stop-=simulation_length
                        elif x_stop < 0:
                            x_stop+=simulation_length

                        top_pos_f.write(f"{x_stop} {y_stop} {z_stop}\n")
                        top_ori_f.write(f"{o1} {o2} {o3}\n")
                        top_idx_f.write(f"{g_id}\n")

                        reach_top += 1
                        break

                # Update growth
                step = d_push / max(np.dot(grad, ooi), 1e-8)
                growth += step

                # Update gradient
                grad = compute_normalized_gradient_at_point(
                    y=yg, z=zg,
                    meltpool=thermal_layers[layer_index],
                    z_center=current_z
                )[0]

    return reach_top
