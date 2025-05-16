import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os


def plot_thermal_history_cross_section(
        thermal_layers:dict, 
        bd_increment: float,
        simulation_width: float,
        z_ultimate:float,
        save_path:str, 
        save:bool
        ):
    
    """
    Plot a cross-sectional view of the thermal history across the meltpool height.

    Parameters
    ----------
        thermal_layers (dict): Dictionary of meltpool parameters.
        bd_increment (float): building direction increment (constant all along the wall)
        simulation_width (float): width of the simulation.
        z_ultimate (float): Virtual top limit of the simulation to avoid singularities from the model.
        save_path (str): Directory to save image if save=True.
        save (bool): If True, saves the plot instead of displaying it.
    """

    # Setup colors and gradients
    base_colors = list(mcolors.BASE_COLORS.values())
    n_lines = 10
    log_space = np.logspace(0, -2, num=n_lines)
    alpha_space = 1 - log_space  # transparency control

    # Setup Y axis (track width) and substrate shape
    y_pos = np.linspace(-simulation_width / 2, simulation_width / 2, 5000)
    z_laser = thermal_layers[0]['height']
    z_top_substrate_limit = z_laser*(1-np.sqrt(1 - (2*y_pos/thermal_layers[0]['width'])**2))

    # Compute vertical limits
    n_layers = len(thermal_layers)
    y_max = (n_layers - 1) * bd_increment + thermal_layers[0]['height']
    y_limits = [i * bd_increment + thermal_layers[0]['height'] for i in range(n_layers)]
    y_limits.insert(0, 0)

    fig, ax = plt.subplots()
    ax.fill_between(y_pos, 0, z_top_substrate_limit, color='black', alpha=0.5)
    for idx, layer_thermal_description in thermal_layers.items():
        color = base_colors[idx % len(base_colors)]

        l = layer_thermal_description['length']
        w = layer_thermal_description['width']
        d = layer_thermal_description['height']

        # Heat source impact (dot marker at top of meltpool)
        ax.scatter(0, idx*bd_increment + thermal_layers[0]['height'], color=color)

        # Plot isothermal contours within the meltpool ellipse
        for t in alpha_space[:-1] * l:
            with np.errstate(invalid='ignore'):  # suppress sqrt domain warnings
                z_curve = z_laser - d * np.sqrt(1 - (t / l) ** 2 - (2 * y_pos / w) ** 2)
                ax.plot(y_pos, z_curve, color=color, alpha=1 - t / l)

        z_laser += bd_increment


    # Guide lines for boundaries
    ax.vlines(x=[-simulation_width / 2, simulation_width / 2], ymin=0, ymax=y_max, color='k', linewidth=0.8)
    ax.hlines(y=y_limits, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='gray', linewidth=0.6)
    ax.hlines(y=z_ultimate, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='k', ls='--', label='Z ultimate')

    # Plot formatting
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("Y axis (track width)", fontsize=8)
    ax.set_ylabel("Z axis (height)", fontsize=8)

    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    if save:
        plt.savefig(f"{save_path}/meltpool_profiles.png", bbox_inches="tight", dpi=300)
    else:
        plt.show()

    plt.close()


def show_domains(thermal_history, working_dir, save, domain, cut_view):
    """
    Show visual cuts (e.g., XZ, YZ) of meltpool and domain of interest.

    Parameters
    ----------
    thermal_history : dict
    working_dir : str
    save : bool
    domain : dict
        Keys: x_min, x_max, y_min, y_max, z_min, z_max
    cut_view : list of (str, float)
    """
    fig, ax = plt.subplots(figsize=(6, 4))

    for view, value in cut_view:
        for idx, layer in thermal_history.items():
            zc = idx * thermal_history[0]["height"]
            height = layer["height"]
            width = layer["width"]
            length = layer["length"]

            theta = np.linspace(0, 2 * np.pi, 200)

            if view.upper() == "XZ":
                x = (length / 2) * np.cos(theta)
                z = zc - height * np.sin(theta)
                ax.plot(x, z, label=f"Layer {idx}")
                ax.set_xlabel("x")
                ax.set_ylabel("z")

            elif view.upper() == "YZ":
                y = (width / 2) * np.cos(theta)
                z = zc - height * np.sin(theta)
                ax.plot(y, z, label=f"Layer {idx}")
                ax.set_xlabel("y")
                ax.set_ylabel("z")

            elif view.upper() == "XY":
                x = (length / 2) * np.cos(theta)
                y = (width / 2) * np.sin(theta)
                ax.plot(x, y, label=f"Layer {idx}")
                ax.set_xlabel("x")
                ax.set_ylabel("y")

    ax.set_title(f"{view} Section — Domain")
    ax.axis("equal")
    ax.legend()
    if save:
        os.makedirs(working_dir, exist_ok=True)
        plt.savefig(f"{working_dir}/domain_{view}.png", dpi=300)
    else:
        plt.show()
    plt.close()


def visualize_growth_vectors_2D(positions, directions, title="Growth directions", save_path=None, scale=0.03):
    """
    Display 2D arrows showing grain growth vectors.

    Parameters
    ----------
    positions : ndarray (N, 3)
        Coordinates of grains (x, y, z).
    directions : ndarray (N, 3)
        Direction vectors for each grain.
    title : str
        Title of the plot.
    save_path : str or None
        If given, path to save the image.
    scale : float
        Arrow scale.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.quiver(
        positions[:, 0], positions[:, 2],
        directions[:, 0], directions[:, 2],
        angles='xy', scale_units='xy', scale=1/scale, color="steelblue", width=0.002
    )
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(title)
    ax.axis("equal")
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()
    plt.close()


def visualize_seeds_at_interface(
        thermal_layers:dict, 
        coo_substrate:np.ndarray,
        coo_interface:np.ndarray,
        width_simulation:float,
        length_simulation:float,
        save:bool,
        save_path:str
        ):
    
    projections = [('X', 'Y', 0, 1), ('X', 'Z', 0, 2), ('Y', 'Z', 1, 2)]
    zc = thermal_layers['height']
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    for ax, (label_i, label_j, i, j) in zip(axs, projections):


        # Domain bounding box
        if label_i  == 'X' and label_j == 'Y':
            ax.scatter(coo_interface[:, i], coo_interface[:, j], color='blue', label='Interface', s=10, zorder=10)
            ax.scatter(coo_substrate[:, i], coo_substrate[:, j], color='red', label='Substrate', s=10)
            ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
            ax.axvline(length_simulation, color='k', linestyle='--', linewidth=0.5)
            
        if label_i  == 'X' and label_j == 'Z':
            ax.scatter(coo_substrate[:, i], coo_substrate[:, j], color='red', label='Substrate', s=10, zorder=10)
            ax.scatter(coo_interface[:, i], coo_interface[:, j], color='blue', label='Interface', s=10)
            ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
            ax.axvline(-width_simulation / 2, color='k', linestyle='--', linewidth=0.5)


        # Meltpool profiles
        if label_i == 'Y' and label_j == 'Z':
            ax.scatter(coo_substrate[:, i], coo_substrate[:, j], color='red', label='Substrate', s=10)
            ax.scatter(coo_interface[:, i], coo_interface[:, j], color='blue', label='Interface', s=10)
            y = np.linspace(-width_simulation / 2, width_simulation / 2, 500)
            z = zc - thermal_layers["height"] * np.sqrt(1 - (2 * y / thermal_layers["width"]) ** 2)
            z[np.isnan(z)] = np.nan
            ax.plot(y, z, linewidth=1, linestyle='--')

        ax.set_xlabel(f"{label_i}")
        ax.set_ylabel(f"{label_j}")
        ax.set_title(f"{label_i}{label_j} projection")
        ax.set_aspect('equal')
        ax.grid(True)

    plt.tight_layout()
    if save:
        plt.savefig(f"{save_path}/meltpool_profiles.png", bbox_inches="tight", dpi=300)
    else:
        plt.show()


def visualize_3d_growth_directions(
    display_ratio: float,
    substrate_coords: np.ndarray = None,
    interface_coords: np.ndarray = None,
    growth_time: np.ndarray = None,
    gradient_vectors: np.ndarray = None,
    growth_directions: np.ndarray = None,
    arrow_length: float = 0.1
):
    """
    Visualize a 3D scatter plot of the interface and substrate grains,
    along with growth directions and thermal gradients.

    Parameters
    ----------
    display_ratio : float
        Fraction of interface points to display (e.g. 0.05 for 5%).
    substrate_coords : np.ndarray, optional
        (N, 3) array of grain seed coordinates in the substrate.
    interface_coords : np.ndarray
        (N, 3) array of grain seed coordinates at the interface.
    growth_time : np.ndarray, optional
        (N,) array of growth times for the interface grains.
    gradient_vectors : np.ndarray, optional
        (N, 3) array of thermal gradient vectors at each interface grain.
    growth_directions : np.ndarray, optional
        (N, 3) array of preferred growth directions (OOI vectors).
    arrow_length : float
        Length of quiver arrows for gradients and directions.
    """
    if interface_coords is None:
        raise ValueError("interface_coords must be provided for visualization.")

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Sample interface points
    mask_int = np.arange(0, len(interface_coords), int(1/display_ratio))
    x_int, y_int, z_int = interface_coords[mask_int].T

    # Plot interface points (colored if growth_time is provided)
    if growth_time is not None:
        ax.scatter(x_int, y_int, z_int, c=growth_time[mask_int], cmap='viridis', marker='o')
    else:
        ax.scatter(x_int, y_int, z_int, color='red', marker='o', label='Interface seeds')

    # Plot substrate points
    if substrate_coords is not None:
        n_sub = int(display_ratio * len(substrate_coords))
        mask_sub = np.random.choice(len(substrate_coords), size=n_sub, replace=False)
        x_sub, y_sub, z_sub = substrate_coords[mask_sub].T
        ax.scatter(x_sub, y_sub, z_sub, color='blue', marker='o', label='Substrate seeds')

    # Plot growth directions (OOI)
    if growth_directions is not None:
        u, v, w = growth_directions[mask_int].T
        ax.quiver(x_int, y_int, z_int, u, v, w, length=arrow_length, color='black', normalize=True, label='Growth direction')

    # Plot thermal gradients
    if gradient_vectors is not None:
        gx, gy, gz = gradient_vectors[mask_int].T
        ax.quiver(x_int, y_int, z_int, gx, gy, gz, length=arrow_length, color='green', normalize=True, label='Thermal gradient')

    # Colorbar
    if growth_time is not None:
        mappable = plt.cm.ScalarMappable(cmap='viridis')
        mappable.set_array(growth_time[mask_int])
        fig.colorbar(mappable, ax=ax, shrink=0.75, label='Growth time')

    # Axis settings
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("3D Grain Seeds and Growth Vectors")
    ax.set_box_aspect([1, 1, 1])
    ax.legend()
    plt.tight_layout()
    plt.show()
