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


def show_domains(thermal_history, simulation_length, simulation_width, bd_increment, z_ultimate, domain_path_dir, save, domain, cut_view):
    """
    Show visual cuts (e.g., XZ, YZ) of meltpool and domain of interest.

    Parameters
    ----------
    thermal_history : dict
    domain_path_dir : str
    save : bool
    domain : dict
        Keys: x_min, x_max, y_min, y_max, z_min, z_max
    cut_view : list of (str, float)
    """

    YMAX = (len(thermal_history) - 1) * bd_increment + thermal_history[0]['height']
    YLIM = [i * bd_increment + thermal_history[0]['height'] for i in range(len(thermal_history))]
    YLIM.insert(0, 0)
    colors = mcolors.BASE_COLORS
    col_keys = list(colors.keys())
    z_source = thermal_history[0]['height']
    ypos = np.linspace(-simulation_width / 2, simulation_width / 2, 5000, endpoint=True)
    zsubstrat = z_source - thermal_history[0]['height'] * np.sqrt(1 - (2 * ypos / thermal_history[0]['width']) ** 2)

    # Création des subplots pour la vue de face et la vue de côté
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    plt.subplots_adjust(wspace=0)

    # Substrat
    ax1.fill_between(ypos, 0, zsubstrat, color='black', alpha=0.5)  
    
    for layer in range(len(thermal_history)):
        col = colors[col_keys[layer % len(thermal_history)]]
        w = thermal_history[layer]['width']
        d = thermal_history[layer]['height']
        zpos = z_source - d * np.sqrt(1 - (2 * ypos / w) ** 2)
        z_source += bd_increment

        ax1.scatter(0, layer * bd_increment + thermal_history[0]['height'], color=col)
        ax1.plot(ypos, zpos, c=col)

        ax2.hlines(y=layer * bd_increment + thermal_history[0]['height'], xmin=0, xmax=simulation_length, color='k')
        # ax2.scatter(LENGTH_SIMULATION/2, layer * BD_INCREMENTS + thermal_history[0]['height'], color=col)
        ax2.hlines(y=np.min(zpos), xmin=0, xmax=simulation_length, color=col)
        ax2.hlines(y=np.max(zpos), xmin=0, xmax=simulation_length, color=col)
        ax2.fill_between([0, simulation_length], np.min(zpos), np.max(zpos), color=col, alpha=0.2)
        
    # Bordures largeur 
    ax1.vlines(x=[-simulation_width / 2, simulation_width / 2], ymin=0, ymax=YMAX, color='k')    
    ax1.hlines(y=YLIM, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='k')
    ax1.hlines(y=z_ultimate, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='k', ls='--')
    # Plan de coupe
    ax1.vlines(x=0, ymin=-0.05, ymax=YMAX+0.05, color='silver', ls='dashdot')

    ax1.set_ylim(-0.05, YMAX + 0.05)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlabel('Y axis')
    ax1.set_ylabel('Z axis')
    ax1.set_title('Front View')
    ax1.set_aspect('equal', adjustable='box')

    ax2.set_ylim(-0.05, YMAX + 0.05)
    ax2.vlines(x=[0, simulation_length], ymin=0, ymax=YMAX, color='k')
    ax2.hlines(y=z_ultimate, xmin=0, xmax=simulation_length, color='k', ls='--')
    # ax2.vlines(x=LENGTH_SIMULATION/2, ymin=-0.05, ymax=YMAX+0.05, color='k', ls='dotted')
    ax2.set_xlabel('X axis')
    ax2.set_ylabel('Z axis')
    ax2.set_title('Side View')

    # Domaine
    xmax = domain['x_max']
    xmin = domain['x_min']
    ymax = domain['y_max']
    ymin = domain['y_min']
    zmax = domain['z_max']
    zmin = domain['z_min']
    color_domain='chartreuse'
    ax1.hlines(y=[zmin, zmax], xmin=ymin, xmax=ymax, color=color_domain)
    ax1.vlines(x=[ymin, ymax], ymin=zmin, ymax=zmax, color=color_domain)
    ax2.hlines(y=[zmin, zmax], xmin=xmin, xmax=xmax, color=color_domain)
    ax2.vlines(x=[xmin, xmax], ymin=zmin, ymax=zmax, color=color_domain)

    # Plans de coupe
    color_cut='chartreuse'
    for plane, position in cut_view:
        if plane == 'XZ':
            y_cut = ymin + position * (ymax - ymin)
            ax1.vlines(x=y_cut, ymin=zmin, ymax=zmax, color=color_cut, ls='dashdot', lw=1.5)
        elif plane == 'YZ':
            x_cut = xmin + position * (xmax - xmin)
            ax2.vlines(x=x_cut, ymin=zmin, ymax=zmax, color=color_cut, ls='dashdot', lw=1.5)
        elif plane == 'XY':
            z_cut = zmin + position * (zmax - zmin)
            ax1.hlines(y=z_cut, xmin=ymin, xmax=ymax, color=color_cut, ls='dashdot', lw=1.5)
            ax2.hlines(y=z_cut, xmin=xmin, xmax=xmax, color=color_cut, ls='dashdot', lw=1.5)

    # Sauvegarde ou affichage des graphiques
    plt.tight_layout()
    if save:
        plt.savefig(f'{domain_path_dir}/domain_selection.png')
    else:
        plt.show()
    plt.close(fig)


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
