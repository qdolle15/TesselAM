import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Dict
import os

## Check thermal history
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
    base_colors = list(mcolors.BASE_COLORS.values())[:-1]
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


## Visualization domain segmentation 
def generate_subdomains(global_domain, axis, num_subdomains, overlap_ratio=0.1):
    """
    Generate subdomains of equal size with overlap for stitching.

    Parameters:
        global_domain (dict): Global domain (e.g., {"x_min": 0, "x_max": 2, ...}).
        axis (str): Direction of splitting ("x", "y", or "z").
        num_subdomains (int): Number of subdomains to create.
        overlap_ratio (float): Overlap ratio between subdomains (e.g., 0.1 for 10%).

    Returns:
        list: List of subdomain dictionaries.
    """
    min_key = f"{axis}_min"
    max_key = f"{axis}_max"
    min_val = global_domain[min_key]
    max_val = global_domain[max_key]
    total_length = max_val - min_val

    # Calculate the length of each subdomain (including overlap)
    subdomain_length = total_length / num_subdomains
    overlap = subdomain_length * overlap_ratio
    step = subdomain_length - overlap

    subdomains = []
    for i in range(num_subdomains):
        start = min_val + i * step
        end = start + subdomain_length

        # Create the subdomain
        subdomain = global_domain.copy()
        subdomain[min_key] = start
        subdomain[max_key] = end

        subdomains.append(subdomain)

    return subdomains

def generate_segmented_cut_views(
    global_domain: dict,
    axis: str,
    num_subdomains: int,
    overlap_ratio: float = 0.1,
    planes: list = ["XZ"],  # Par défaut, on utilise des plans XZ
    rel_positions: list = [0.25, 0.5, 0.75]
) -> list:
    """
    Generate CUT_VIEWS by segmenting a global domain into subdomains with overlap.

    Args:
        global_domain: Global domain (e.g., {"x_min": 0, "x_max": 2, ...}).
        axis: Axis along which to split the domain ("x", "y", or "z").
        num_subdomains: Number of subdomains to create.
        overlap_ratio: Overlap ratio between subdomains (e.g., 0.1 for 10%).
        planes: List of planes for cuts (e.g., ["XZ"]).
        rel_positions: Relative positions for cuts (e.g., [0.25, 0.5, 0.75]).

    Returns:
        List of CUT_VIEWS dictionaries.
    """
    # Générer les sous-domaines
    subdomains = generate_subdomains(global_domain, axis, num_subdomains, overlap_ratio)

    # Générer les CUT_VIEWS pour chaque sous-domaine
    cut_views = []
    for subdomain in subdomains:
        cut_views.append({
            "domain": subdomain,
            "plans": [(plane, pos) for plane in planes for pos in rel_positions]
        })

    return cut_views

def show_domains(
    thermal_history, simulation_length, simulation_width, bd_increment, z_ultimate,
    domain_path_dir, save, domain, cut_view, global_domain=None
):
    """
    Show visual cuts (e.g., XZ, YZ) of meltpool, global domain, and subdomains.

    Parameters:
        thermal_history (dict): Thermal history data.
        simulation_length (float): Length of the simulation.
        simulation_width (float): Width of the simulation.
        bd_increment (float): Build increment.
        z_ultimate (float): Ultimate Z height.
        domain_path_dir (str): Path to save the plot.
        save (bool): Whether to save the plot.
        domain (dict): Subdomain to visualize (keys: x_min, x_max, y_min, y_max, z_min, z_max).
        cut_view (list): List of cut planes (e.g., [("XZ", 0.5), ("YZ", 0.5)]).
        global_domain (dict, optional): Global domain to visualize.
    """
    # Calculate Y limits
    YMAX = (len(thermal_history) - 1) * bd_increment + thermal_history[0]['height']
    YLIM = [i * bd_increment + thermal_history[0]['height'] for i in range(len(thermal_history))]
    YLIM.insert(0, 0)

    # Initialize plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    plt.subplots_adjust(wspace=0)

    # --- Substrate ---
    ypos = np.linspace(-simulation_width / 2, simulation_width / 2, 5000, endpoint=True)
    zsubstrat = thermal_history[0]['height'] - thermal_history[0]['height'] * np.sqrt(1 - (2 * ypos / thermal_history[0]['width']) ** 2)
    ax1.fill_between(ypos, 0, zsubstrat, color='black', alpha=0.5)

    # --- Thermal history layers ---
    colors = list(mcolors.BASE_COLORS.values())
    z_source = thermal_history[0]['height']
    for layer in range(len(thermal_history)):
        col = colors[layer % len(colors)]
        w = thermal_history[layer]['width']
        d = thermal_history[layer]['height']
        zpos = z_source - d * np.sqrt(1 - (2 * ypos / w) ** 2)
        z_source += bd_increment

        ax1.scatter(0, layer * bd_increment + thermal_history[0]['height'], color=col)
        ax1.plot(ypos, zpos, c=col)
        ax2.hlines(y=layer * bd_increment + thermal_history[0]['height'], xmin=0, xmax=simulation_length, color='k')
        ax2.hlines(y=np.min(zpos), xmin=0, xmax=simulation_length, color=col)
        ax2.hlines(y=np.max(zpos), xmin=0, xmax=simulation_length, color=col)
        ax2.fill_between([0, simulation_length], np.min(zpos), np.max(zpos), color=col, alpha=0.2)

    # --- Simulation domain ---
    ax1.vlines(x=[-simulation_width / 2, simulation_width / 2], ymin=0, ymax=YMAX, color='k')
    ax1.hlines(y=YLIM, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='k')
    ax1.hlines(y=z_ultimate, xmin=-simulation_width / 2, xmax=simulation_width / 2, color='k', ls='--')
    ax2.vlines(x=[0, simulation_length], ymin=0, ymax=YMAX, color='k')
    ax2.hlines(y=z_ultimate, xmin=0, xmax=simulation_length, color='k', ls='--')

    # --- Visualisation domain ---
    if global_domain:
        global_color = 'chartreuse'
        ax1.hlines(y=[global_domain['z_min'], global_domain['z_max']],
                  xmin=global_domain['y_min'], xmax=global_domain['y_max'], color=global_color, ls='--')
        ax1.vlines(x=[global_domain['y_min'], global_domain['y_max']],
                  ymin=global_domain['z_min'], ymax=global_domain['z_max'], color=global_color, ls='--')
        ax2.hlines(y=[global_domain['z_min'], global_domain['z_max']],
                  xmin=global_domain['x_min'], xmax=global_domain['x_max'], color=global_color, ls='--')
        ax2.vlines(x=[global_domain['x_min'], global_domain['x_max']],
                  ymin=global_domain['z_min'], ymax=global_domain['z_max'], color=global_color, ls='--')

    # --- Subdomain ---
    subdomain_color = 'chartreuse'
    ax1.hlines(y=[domain['z_min'], domain['z_max']],
              xmin=domain['y_min'], xmax=domain['y_max'], color=subdomain_color)
    ax1.vlines(x=[domain['y_min'], domain['y_max']],
              ymin=domain['z_min'], ymax=domain['z_max'], color=subdomain_color)
    ax2.hlines(y=[domain['z_min'], domain['z_max']],
              xmin=domain['x_min'], xmax=domain['x_max'], color=subdomain_color)
    ax2.vlines(x=[domain['x_min'], domain['x_max']],
              ymin=domain['z_min'], ymax=domain['z_max'], color=subdomain_color)

    # --- Cut planes ---
    cut_color = 'red'
    for plane, position in cut_view:
        if plane == 'XZ':
            y_cut = domain['y_min'] + position * (domain['y_max'] - domain['y_min'])
            ax1.vlines(x=y_cut, ymin=domain['z_min'], ymax=domain['z_max'], color=cut_color, ls='dashdot', lw=1.5)
        elif plane == 'YZ':
            x_cut = domain['x_min'] + position * (domain['x_max'] - domain['x_min'])
            ax2.vlines(x=x_cut, ymin=domain['z_min'], ymax=domain['z_max'], color=cut_color, ls='dashdot', lw=1.5)
        elif plane == 'XY':
            z_cut = domain['z_min'] + position * (domain['z_max'] - domain['z_min'])
            ax1.hlines(y=z_cut, xmin=domain['y_min'], xmax=domain['y_max'], color=cut_color, ls='dashdot', lw=1.5)
            ax2.hlines(y=z_cut, xmin=domain['x_min'], xmax=domain['x_max'], color=cut_color, ls='dashdot', lw=1.5)

    # --- Plot formatting ---
    ax1.set_ylim(-0.05, YMAX + 0.05)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlabel('Y axis')
    ax1.set_ylabel('Z axis')
    ax1.set_title('Front View')

    ax2.set_ylim(-0.05, YMAX + 0.05)
    ax2.set_xlabel('X axis')
    ax2.set_ylabel('Z axis')
    ax2.set_title('Side View')

    # --- Save or show ---
    plt.tight_layout()
    if save:
        plt.savefig(f'{domain_path_dir}')
    else:
        plt.show()
    plt.close(fig)

## Debugging zone
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
