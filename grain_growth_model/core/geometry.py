import numpy as np


def compute_min_distances_between_lines(
        coordinates_seeds: np.ndarray, 
        growth_directions: np.ndarray, 
        mean_grain_diameter:float
        ):
    """
    Computes shortest distances and parameters between all pairs of growth lines.

    Parameters
    ----------
        coordinates_seeds (ndarray) (N, 3): Starting positions of grains.
        growth_directions (ndarray) (N, 3): Easy growth direction vectors (OOI).
        mean_grain_diameter (float): Mean grain size for substrate initialization.

    Returns
    -------
        distances_normalized (ndarray) (N, N): Normalized minimal distance between lines.
        coefficients : (ndarray) (N, N): Parameters λ of intersection along each direction.
    """
    N = coordinates_seeds.shape[0]
    d1 = growth_directions[:, None, :]  # shape (N,1,3)
    d2 = growth_directions[None, :, :]  # shape (1,N,3)

    n = np.cross(d1, d2, axis=2)
    norm_n = np.linalg.norm(n, axis=2)

    diff = coordinates_seeds[:, None, :] - coordinates_seeds[None, :, :]  # shape (N,N,3)
    n_dot_diff = np.einsum('ijk,ijk->ij', n, diff)

    with np.errstate(divide='ignore', invalid='ignore'):
        distances = np.abs(n_dot_diff) / norm_n
        distances[norm_n == 0] = 0.0
        distances_normalized = distances / mean_grain_diameter

        n1 = np.cross(n, d2, axis=2)
        n1_dot_diff = np.einsum('ijk,ijk->ij', n1, diff)
        n1_dot_d1 = np.einsum('ijk,ijk->ij', n1, d1)

        coefficients = n1_dot_diff / n1_dot_d1
        coefficients[n1_dot_d1 == 0] = np.nan

    return distances_normalized, coefficients