import numpy as np
from scipy.integrate import quad

def ellipse_arc_length(x_axis, y_axis, theta_min, theta_max):
    """
    Computes the arc length of an ellipse between two angles.

    Parameters
    ----------
        x_axis (float) Semi-major axis.
        y_axis (float): Semi-minor axis.
        theta_min, theta_max  (float) Angular bounds in radians.

    Returns
    -------
        (float): Elliptical arc length.
    """
    def integrand(theta, x_axis, y_axis):
        return np.sqrt((x_axis * np.sin(theta))**2 + (y_axis * np.cos(theta))**2)
    
    length, _ = quad(integrand, theta_min, theta_max, args=(x_axis, y_axis))
    return length


def find_x_on_ellipsoid(y, z, xc, yc, zc, x_axis, y_axis, z_axis):
    """
    Given y and z, solve for x on an ellipsoid centered at (x_axis, y_axis, z_axis).

    Parameters
    ----------
        y, z (float or array-like): Coordinates of evaluation.
        xc, yc, zc (float): Center of the ellipse coordinates.
        x_axis, y_axis, z_axis (float): Half axis of the ellipse.

    Returns
    -------
        x1, x2 (tuple): Two possible x solutions on the ellipsoid.
    """
    term = 1 - ((y - yc)**2 / y_axis**2) - ((z - zc)**2 / z_axis**2)
    if np.any(term < 0):
        raise ValueError("Point lies outside the ellipsoid domain.")
    sqrt_term = np.sqrt(term)
    return xc + x_axis * sqrt_term, xc - x_axis * sqrt_term


def compute_normalized_gradient_at_point(y, z, meltpool, z_center, x_center=0.0, y_center=0.0):
    """
    Compute unit gradient vector of the meltpool surface at point (y,z).

    Parameters
    ----------
        y, z (float or array-like): Coordinates of evaluation.
        meltpool (dict): Meltpool description with keys: 'length', 'width', 'height', 'PD'.
        z_center (float): Center position in z.
        x_center, y_center (float): Center positions in x and y.

    Returns
    -------
        (ndarray): Normalized gradient vectors, shape (N, 3).
    """
    strategy = meltpool['PD']
    x_half_axis = meltpool['length']
    y_half_axis = meltpool['width'] / 2
    z_half_axis = meltpool['height']

    x1, x2 = find_x_on_ellipsoid(y, z, x_center, y_center, z_center, x_half_axis, y_half_axis, z_half_axis)
    x = x2 if strategy > 0 else x1

    grad_x = 2 * (x - x_center) / x_half_axis**2
    grad_y = 2 * (y - y_center) / y_half_axis**2
    grad_z = 2 * (z - z_center) / z_half_axis**2

    norm = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    grad = -np.vstack((grad_x, grad_y, grad_z)) / norm
    return grad.T


def compute_parametric_time(x, y, z, meltpool, zc, simulation_length):
    """
    Compute parametric time t used in a spherical meltpool description.

    Parameters
    ----------
        x, y, z (float or array-like): Coordinates of evaluation.
        meltpool (dict): Meltpool description with keys: 'length', 'width', 'height', 'PD'.
        zc (float): z-coordinates of the ellipse center.
        simulation_length (float): Length of the simulation

    Returns
    -------
        (float or ndarray): Computed value of 't'.
    """
    x_reference = 0.0
    strategy = meltpool['PD']
    x_half_axis = meltpool['length']
    y_half_axis = meltpool['width'] / 2
    z_half_axis = meltpool['height']

    theta = -np.arcsin((z - zc) / z_half_axis)
    cos_theta = np.cos(theta)
    arg = y / (y_half_axis * cos_theta)
    arg = np.clip(arg, -1.0, 1.0)
    phi = -np.arccos(arg)

    if strategy < 0:
        phi = abs(phi)
        x_reference = simulation_length

    t = x - x_reference - x_half_axis * np.cos(theta) * np.sin(phi)
    return abs(t)


def is_point_inside_meltpool(y, z, half_width, depth, yc=0.0, zc=0.0):
    """
    Checks if a (y,z) point is inside a 2D elliptical meltpool cross-section.

    Parameters
    ----------
        y, z (float or ndarray): Coordinates to check.
        half_width (float): Half-width of the meltpool (semi-axis).
        depth (float): Depth of the meltpool (semi-axis).
        yc, zc (float): Center coordinates of the ellipse.

    Returns
    -------
        (bool or ndarray): True if the point lies within the meltpool.
    """
    return ((y - yc) / half_width)**2 + ((z - zc) / depth)**2 < 1


def find_meltpool_exit_point(point, direction, meltpool, z0):
    """
    Find the intersection point between a growth direction and the meltpool boundary.

    Parameters
    ----------
        point (ndarray) (3,): Starting coordinates.
        direction (ndarray) (3,): Growth direction vector.
        meltpool (dict): Geometry of the meltpool.
        z0 (float): z-position of the meltpool center.

    Returns
    -------
        (tuple) (x, y, z): exit point.
    """
    _, dy, dz = direction
    xa, ya, za = point
    w, d = meltpool['width'], meltpool['height']

    A = (2 * dy / w)**2 + (dz / d)**2
    B = 2 * ((4 * ya * dy / w**2) + (2 * (za - z0) * dz / d**2))
    C = (2 * ya / w)**2 + ((za - z0) / d)**2 - 1

    delta = B**2 - 4 * A * C
    if delta <= 0:
        return tuple(point)

    t1 = (-B - np.sqrt(delta)) / (2 * A)
    t2 = (-B + np.sqrt(delta)) / (2 * A)
    t = min(t1, t2)
    return tuple(np.asarray(point) + t * np.asarray(direction))


def find_simulation_top_exit_point(point, direction, z_top):
    """
    Computes the point where a vector intersects the simulation top surface (plane).

    Parameters
    ----------
        point (ndarray) (3, ): Starting coordinates.
        direction (ndarray) (3, ): Direction vector.
        z_top (float): Top surface Z coordinate.

    Returns
    -------
        tuple (x, y, z): intersection point on top plane.
    """
    dz = direction[2]
    if dz == 0:
        return tuple(point)
    t = (z_top - point[2]) / dz
    return tuple(np.asarray(point) + t * np.asarray(direction))
