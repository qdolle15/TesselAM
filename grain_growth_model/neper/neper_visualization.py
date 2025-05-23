import numpy as np
import matplotlib.pyplot as plt
from orix.quaternion import Orientation, symmetry
from orix.vector import Vector3d
from orix.plot import IPFColorKeyTSL
from orix.crystal_map import Phase
from orix.vector import Miller


def extract_and_plot_slice(arr: np.ndarray, plane: str, position: float) -> np.ndarray:
    """
    Extract a 2D slice from a 3D volume along a specified plane at a relative position.

    Parameters
    ----------
    arr : np.ndarray
        3D volume array of shape (Z, Y, X).
    plane : str
        Slice plane: 'XY', 'XZ', or 'YZ'.
    position : float
        Relative position in [0, 1] within the volume.

    Returns
    -------
    slice_2d : np.ndarray
        2D array slice.
    """
    assert 0.0 <= position <= 1.0, "Position must be between 0 and 1."

    shape = arr.shape
    planes = {
        'XY': arr[int(position * shape[0]), :, :],
        'XZ': arr[:, int(position * shape[1]), :],
        'YZ': arr[:, :, int(position * shape[2])]
    }

    if plane not in planes:
        raise ValueError("Invalid plane. Choose from 'XY', 'XZ', or 'YZ'.")

    return planes[plane]


def color_IPF(flatten_ori: np.ndarray, direction: str, plane_dimension: tuple) -> np.ndarray:
    """
    Convert orientations to IPF-colored RGB image for a given direction.

    Parameters
    ----------
    flatten_ori : np.ndarray
        Array of Euler angles in degrees (N, 3).
    direction : str
        Projection direction: 'x', 'y', or 'z'.
    plane_dimension : tuple
        Dimensions of the 2D plane (height, width).

    Returns
    -------
    np.ndarray
        RGB array of shape (H, W, 3).
    """
    assert direction in ("x", "y", "z"), "Direction must be 'x', 'y', or 'z'."
    v = {"x": Vector3d([1, 0, 0]), "y": Vector3d([0, 1, 0]), "z": Vector3d([0, 0, 1])}[direction]

    ori = Orientation.from_euler(flatten_ori, degrees=False)
    ipfkey = IPFColorKeyTSL(symmetry.Oh, direction=v)
    ori.symmetry = ipfkey.symmetry
    rgb = ipfkey.orientation2color(ori)
    return rgb.reshape((plane_dimension[0], plane_dimension[1], 3))


def image_cross_section(arr: np.ndarray, save: bool, path_save: str = None):
    """
    Display or save an RGB image of a cross-section.

    Parameters
    ----------
    arr : np.ndarray
        RGB array (H, W, 3).
    save : bool
        If True, saves the image to `path_save`.
    path_save : str, optional
        Path where to save the image if `save` is True.
    """
    fig, ax = plt.subplots()
    ax.imshow(arr, origin='lower')
    ax.axis('off')

    if save and path_save:
        plt.savefig(path_save, bbox_inches='tight', dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)


def image_IPF_triangle(flatten_ori: np.ndarray, direction: str, save: bool, path_save: str = None):
    """
    Plot an IPF triangle from orientation data.

    Parameters
    ----------
    flatten_ori : np.ndarray
        Array of Euler angles in degrees (N, 3).
    direction : str
        Projection direction: 'x', 'y', or 'z'.
    save : bool
        If True, saves the plot.
    path_save : str, optional
        Output path for saving.
    """
    assert direction in ("x", "y", "z"), "Direction must be 'x', 'y', or 'z'."
    v = {"x": Vector3d([1, 0, 0]), "y": Vector3d([0, 1, 0]), "z": Vector3d([0, 0, 1])}[direction]

    ori = Orientation.from_euler(flatten_ori, degrees=False)
    ipfkey = IPFColorKeyTSL(symmetry.Oh, direction=v)
    ori.symmetry = ipfkey.symmetry
    rgb = ipfkey.orientation2color(ori)

    ori.scatter("ipf", c=rgb, direction=v)
    if save and path_save:
        plt.savefig(path_save, bbox_inches='tight', dpi=300)
        plt.close()
    else:
        plt.show()


def pole_figure(flatten_ori: np.ndarray, path_save: str, save: bool):
    """
    Plot pole figures for <100>, <101>, <111> directions.

    Parameters
    ----------
    flatten_ori : np.ndarray
        Euler angles in degrees.
    path_save : str
        Path where to save the figure.
    save : bool
        If True, saves the figure.
    """
    ori = Orientation.from_euler(flatten_ori, degrees=False, symmetry=symmetry.Oh)

    phase = Phase(point_group="m-3m")
    t_100 = Miller(xyz=[1, 0, 0], phase=phase).symmetrise(unique=True)
    t_101 = Miller(xyz=[1, 0, 1], phase=phase).symmetrise(unique=True)
    t_111 = Miller(xyz=[1, 1, 1], phase=phase).symmetrise(unique=True)

    v100 = ori.inv().outer(t_100)
    v101 = ori.inv().outer(t_101)
    v111 = ori.inv().outer(t_111)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": "stereographic"})
    axs[0].pole_density_function(v100, cmap='jet', sigma=5)
    axs[0].set(title='<100>')

    axs[1].pole_density_function(v101, cmap='jet', sigma=5)
    axs[1].set(title='<101>')

    axs[2].pole_density_function(v111, cmap='jet', sigma=5)
    axs[2].set(title='<111>')

    for ax in axs:
        ax.set_labels("PD", "ND", None)

    plt.tight_layout()
    if save:
        plt.savefig(path_save, dpi=300)
        plt.close()
    else:
        plt.show()
