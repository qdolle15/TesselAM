import numpy as np

from grain_growth_model.core.meltpool import ellipse_arc_length, is_point_inside_meltpool


def generate_initial_substrate_grains(
        thermal_layers: dict, 
        simulation_length: float,
        simulation_width: float,
        mean_grain_diameter:float
        ):
    """
    Generate randomly distributed grains within the initial substrate volume.
    Any grains located in the remelted zone (first meltpool) are removed.
    Purely esthetic.

    Parameters
    ----------
        thermal_layers (dict): Full thermal history.
        simulation_length (float): length of the simulation, periodic-like aspect along x-axis.
        simulation_width (float): width of the simulation.
        mean_grain_diameter (float): Mean grain size for substrate initialization.

    Returns
    -------
        tuple of (ndarray, ndarray, ndarray) Coordinates, orientations and indexes of non-remelted substrate grains.
    """
    depth = thermal_layers['height']
    width = thermal_layers['width']

    XMIN, XMAX = 0, simulation_length
    YMIN, YMAX = -simulation_width / 2, simulation_width / 2
    ZMIN, ZMAX = 0, depth

    num_grains = int(((XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN)) / mean_grain_diameter**3)
    X = np.random.uniform(XMIN, XMAX, num_grains)
    Y = np.random.uniform(YMIN, YMAX, num_grains)
    Z = np.random.uniform(ZMIN, ZMAX, num_grains)

    orientations = np.random.uniform(-np.pi, np.pi, size=(num_grains, 3))

    # Filter out remelted grains
    is_inside = is_point_inside_meltpool(
        y=Y, z=Z,
        half_width=width / 2,
        depth=depth,
        yc=0,
        zc=ZMAX
    )

    return np.vstack((X[~is_inside], Y[~is_inside], Z[~is_inside])).T, orientations[~is_inside, :], -1 * np.ones(len(X[~is_inside]), dtype=int)


def initialize_interface_grains(
        thermal_layers: dict, 
        layer_index: int, 
        known_coords: np.ndarray,
        known_orientations: np.ndarray,
        known_ids: np.ndarray,
        simulation_length: float,
        simulation_width: float,
        mean_grain_diameter:float,
        bd_increment: float,
        z_ultimate: float,
        No_grains_created: int

        ):
    """
    Initialize grain seeds on the meltpool interface depending on whether epitaxial growth applies.
    Add randomly distributed nucleated grains on the lateral side of the bead (call wall).

    Parameters
    ----------
        thermal_layers (dict): Dictionary containing meltpool descriptions at each layer.
        layer_index (int): Index of the current layer.
        known_coords (ndarray or None): Interface coordinates from previous layer.
        known_orientations (ndarray or None): Interface orientations from previous layer.
        known_ids (ndarray or None): Interface indexes from previous layer.
        simulation_length (float): length of the simulation.
        mean_grain_diameter (float): Mean grain size for substrate initialization.
        bd_increment (float): building direction increment (constant all along the wall)
        simulation_width (float): width of the simulation.
        z_ultimate (float): Top limit of the simulation.
        No_grains_created (int): Number of grains already created along the simulation.

    Returns
    -------
        (tuple): Coordinates, orientations, and number of new grains created.
    """
    y0 = 0
    z0 = layer_index * bd_increment + thermal_layers[0]['height']
    k_inside = 0.999

    current_depth = thermal_layers[layer_index]['height']
    current_width = thermal_layers[layer_index]['width']
    theta_limit = np.arcsin(simulation_width / current_width)

    if known_coords is None:
        arc_length = ellipse_arc_length(current_width / 2, current_depth, -theta_limit, theta_limit)
        interface_area = arc_length * simulation_length
        num_interface_grains = int(interface_area / mean_grain_diameter**2)

        theta = np.random.uniform(-theta_limit, theta_limit, size=num_interface_grains)
        x_coords = np.random.uniform(0, simulation_length, num_interface_grains)
        y_coords = y0 + (current_width / 2) * k_inside * np.sin(theta)
        z_coords = z0 - current_depth * k_inside * np.cos(theta)
        orientations = np.random.uniform(-np.pi, np.pi, size=(num_interface_grains, 3))

        z_low = z0 - current_depth * k_inside * np.cos(theta_limit)
        if layer_index + 1 < len(thermal_layers):
            next_width = thermal_layers[layer_index + 1]['width']
            next_depth = thermal_layers[layer_index + 1]['height']
            next_theta_limit = np.arcsin(simulation_width / next_width)
            z_high = z0 + bd_increment - next_depth * k_inside * np.cos(next_theta_limit)
        else:
            z_high = z_ultimate

        num_wall_grains = int((z_high - z_low) * simulation_length / mean_grain_diameter**2)

        x1 = np.random.uniform(0, simulation_length, num_wall_grains)
        y1 = np.ones(num_wall_grains) * (simulation_width / 2)
        z1 = np.random.uniform(z_low, z_high, num_wall_grains)

        x2 = np.random.uniform(0, simulation_length, num_wall_grains)
        y2 = -np.ones(num_wall_grains) * (simulation_width / 2)
        z2 = np.random.uniform(z_low, z_high, num_wall_grains)

        wall_orientations = np.random.uniform(-np.pi, np.pi, size=(2 * num_wall_grains, 3))

        all_coords = np.hstack((
            np.vstack((x_coords, y_coords, z_coords)),
            np.vstack((x1, y1, z1)),
            np.vstack((x2, y2, z2))
        )).T
        all_orientations = np.vstack((orientations, wall_orientations))
        new_grains = num_interface_grains + 2 * num_wall_grains

    else:
        if layer_index + 1 < len(thermal_layers):
            next_width = thermal_layers[layer_index + 1]['width']
            next_depth = thermal_layers[layer_index + 1]['height']
            next_theta_limit = np.arcsin(simulation_width / next_width)
            z_high = z0 + bd_increment - next_depth * k_inside * np.cos(next_theta_limit)
        else:
            z_high = z_ultimate

        z_low = z0 - current_depth * k_inside * np.cos(theta_limit)
        num_wall_grains = int((z_high - z_low) * simulation_length / mean_grain_diameter**2)

        x1 = np.random.uniform(0, simulation_length, num_wall_grains)
        y1 = np.ones(num_wall_grains) * (simulation_width / 2)
        z1 = np.random.uniform(z_low, z_high, num_wall_grains)

        x2 = np.random.uniform(0, simulation_length, num_wall_grains)
        y2 = -np.ones(num_wall_grains) * (simulation_width / 2)
        z2 = np.random.uniform(z_low, z_high, num_wall_grains)

        wall_orientations = np.random.uniform(-np.pi, np.pi, size=(2 * num_wall_grains, 3))

        all_coords = np.hstack((
            known_coords.T,
            np.vstack((x1, y1, z1)),
            np.vstack((x2, y2, z2))
        )).T
        all_orientations = np.vstack((known_orientations, wall_orientations))
        new_grains = 2 * num_wall_grains

    new_ids = np.arange(No_grains_created, No_grains_created + new_grains, dtype=int)
    all_ids = np.concatenate([known_ids, new_ids])

    return all_coords, all_orientations, all_ids, new_grains
