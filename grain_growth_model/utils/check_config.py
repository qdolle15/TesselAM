from typing import Protocol
import numpy as np


class ConfigProtocol(Protocol):
    START_SIMULATION: bool
    LENGTH_SIMULATION: float
    WIDTH_SIMULATION: float
    D0: float
    D_seeds: float
    THRESHOLD: float
    NOISE_NEPER: float
    BD_INCREMENTS: float
    Z_ULTIMATE: float
    THERMAL_HISTORY: dict
    EGD_FAMILY: np.ndarray
    CUT_VIEWS: list
    START: bool
    RANDOM_SEED: int


def validate_config_module(config_module, required_fields=None):
    """
    Validate that a configuration module contains all required attributes.

    Parameters
    ----------
        config_module (module) Python module object loaded dynamically (e.g. via load_config_module()).
        required_fields (list[str] or None) List of required attribute names. If None, uses default FAST-MMAM config fields.

    Raises
    ------
        (AttributeError) If any required field is missing from the config module.
    """
    if required_fields is None:
        required_fields = [
            "START_SIMULATION",
            "LENGTH_SIMULATION",
            "WIDTH_SIMULATION",
            "D0",
            "D_seeds",
            "THRESHOLD",
            "NOISE_NEPER",
            "BD_INCREMENTS",
            "Z_ULTIMATE",
            "THERMAL_HISTORY",
            "EGD_FAMILY",
            "CUT_VIEWS",
            "START",
            "RANDOM_SEED"
        ]

    missing = [field for field in required_fields if not hasattr(config_module, field)]

    if missing:
        raise AttributeError(f"Config file is missing required fields: {missing}")
    else:
        print("Configuration file validated successfully.")


def max_array_size_from_bytes(max_bytes: int, dtype=np.float32) -> int:
    """
    Compute the maximum dimension N for an N x N array of a given dtype
    that fits within a byte size limit.

    Parameters
    ----------
    max_bytes : int
        Maximum allowed memory size in bytes.
    dtype : np.dtype
        NumPy data type (default: np.float32).

    Returns
    -------
    int
        Maximum allowed value of N for an (N x N) array.
    """
    bytes_per_element = np.dtype(dtype).itemsize
    max_elements = max_bytes // bytes_per_element
    max_side_length = int(np.floor(np.sqrt(max_elements)))
    return max_side_length
