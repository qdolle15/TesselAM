import numpy as np

from grain_growth_model.neper.tools import create_sub_selection, read_data, transform_data
from grain_growth_model.neper.cmd import create_tessellation
from grain_growth_model.neper.neper_visualization import (
    extract_and_plot_slice, color_IPF,
    image_cross_section, image_IPF_triangle, pole_figure
)

def EBSD_like(domain: dict, cut_view: dict, save: bool, working_dir: str, global_dir: str):
    """
    Generate EBSD-like maps and pole figures using Neper tessellation.

    Parameters
    ----------
    domain : dict
        Domain bounding box with keys x_min, x_max, y_min, y_max, z_min, z_max.
    cut_view : dict
        Dictionary of plane:position (e.g. {'YZ': 0.5}) to cut for EBSD-like maps.
    save : bool
        If True, save all generated images.
    working_dir : str
        Directory to save intermediate and final results.
    global_dir : str
        Directory where global simulation data is stored.

    Returns
    -------
    dict
        Information about tessellation (voxels, command, seeds, time, etc.).
    """
    N, grains_per_layers = create_sub_selection(domain, working_dir, global_dir)
    cmd, voxels, time_tess = create_tessellation(N, domain, working_dir)

    angles = np.loadtxt(f'{working_dir}/sub_ori.txt', delimiter=' ')
    tesr_raw = read_data(f'{working_dir}/tessellation_3d.tesr')
    tesr_flatten = transform_data(tesr_raw)
    tesr_3d = np.transpose(tesr_flatten.reshape((voxels[2], voxels[1], voxels[0])), (0, 1, 2))

    for plane, position in cut_view.items():
        direction = "z"
        suf = str(position).replace('.', 'o')
        tesr_cut_2d = extract_and_plot_slice(tesr_3d, plane=plane, position=position)
        cut_angles_flatten = angles[tesr_cut_2d.ravel() - 1]

        colors = color_IPF(flatten_ori=cut_angles_flatten, direction=direction, plane_dimension=tesr_cut_2d.shape)

        image_cross_section(colors, save, f'{working_dir}/EBSD_{direction}_{plane}_{suf}.pdf')
        image_IPF_triangle(cut_angles_flatten, save, direction, f'{working_dir}/IPF_{direction}_triangle_{plane}_{suf}.pdf')
        pole_figure(cut_angles_flatten, f'{working_dir}/PF_{direction}_{plane}_{suf}.pdf', save)

    return {
        'voxels': voxels,
        'N': N,
        'cmd': cmd,
        'work_path': working_dir,
        'time_tess': time_tess,
        'grains_per_layers': grains_per_layers
    }

