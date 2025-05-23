import numpy as np

from grain_growth_model.neper.neper_tools import create_sub_selection, read_data, transform_data
from grain_growth_model.neper.neper_command import create_tessellation
from grain_growth_model.neper.neper_visualization import (
    extract_and_plot_slice, color_IPF,
    image_cross_section, image_IPF_triangle, pole_figure
)

def EBSD_like(No_layer:int, domain: dict, cut_view: list, save: bool, data_path_dir: str, domain_path_dir: str) -> dict:
    """
    Generate EBSD-like maps and pole figures using Neper tessellation.

    Parameters
    ----------
        No_layer (int): Number of layers to consider in the simulation (including substrate).
        domain (dict): Spatial domain to crop, with keys: x_min, x_max, y_min, y_max, z_min, z_max.
        cut_view (list of tuple): List of planes and relative positions (e.g., [("YZ", 0.5), ("XZ", 0.25)]).
        save (bool) Whether to save images (True) or display them (False).
        data_path_dir (str) Path to full simulation output files (for sub-selection).
        domain_path_dir (str): Output path for saving EBSD-like visualizations.

    Returns
    -------
        (dict) Dictionary containing tessellation information and statistics.
    """

    # Crop seeds in domain
    No_seeds, grains_per_layers = create_sub_selection(
        No_layer=No_layer,
        domain=domain,
        data_path=data_path_dir,
        domain_path_dir=domain_path_dir
    )
    # Create Neper tessellation of this sub-domain
    cmd, voxels, time_tess = create_tessellation(
        nbr_seeds=No_seeds,
        domain=domain,
        domain_path_dir=domain_path_dir
    )

    # Load crystal orientations
    angles = np.loadtxt(f'{domain_path_dir}/sub_ori.txt', delimiter=' ')

    # Load 3D tessellation
    tesr_raw = read_data(f'{domain_path_dir}/tessellation_3d.tesr')
    tesr_flatten = transform_data(tesr_raw)
    tesr_3d = np.transpose(
        tesr_flatten.reshape((voxels[2], voxels[1], voxels[0])),
        (0, 1, 2)
    )
    # The grain is made up of voxels, each voxel with the identifier of the grain to which it belongs. \
    # This identifier corresponds to the line in the 'sub_ori.txt' file.

    # Loop over requested 2D views
    for plane, position in cut_view:
        direction = "z"
        suf = f"{int(position * 100):03d}"

        # Cut the 3D tessellation
        tesr_cut_2d = extract_and_plot_slice(
            arr=tesr_3d,
            plane=plane,
            position=position
        )
        cut_angles_flatten = angles[tesr_cut_2d.ravel() - 1]

        # Color the section with IPF coloring
        array_EBSD_like_colored = color_IPF(
            flatten_ori=cut_angles_flatten,
            direction=direction,
            plane_dimension=tesr_cut_2d.shape
        )
        
        # Save or show results
        image_cross_section(
            arr=array_EBSD_like_colored,
            save=save,
            path_save=f'{domain_path_dir}/EBSD_{direction}_{plane}_{suf}.png'
        )
        image_IPF_triangle(
            flatten_ori=cut_angles_flatten,
            direction=direction,
            save=save,
            path_save=f'{domain_path_dir}/IPF_{direction}_triangle_{plane}_{suf}.png'
        )
        pole_figure(
            flatten_ori=cut_angles_flatten,
            save=save,
            path_save=f'{domain_path_dir}/PF_{direction}_{plane}_{suf}.png'
        )

    return {
        'voxels': voxels,
        'N': No_seeds,
        'cmd': cmd,
        'work_path': domain_path_dir,
        'time_tess': time_tess,
        'grains_per_layers': grains_per_layers
    }

