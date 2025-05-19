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
    SUBSTRATE_BBOX_EXCESS: float
    BD_INCREMENTS: float
    Z_ULTIMATE: float
    THERMAL_HISTORY: dict
    EGD_FAMILY: np.ndarray
    DOMAINS: list
    CUT_VIEWS: list
    START: bool

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
            "SUBSTRATE_BBOX_EXCESS",
            "BD_INCREMENTS",
            "Z_ULTIMATE",
            "THERMAL_HISTORY",
            "EGD_FAMILY",
            "DOMAINS",
            "CUT_VIEWS",
            "START",
        ]

    missing = [field for field in required_fields if not hasattr(config_module, field)]

    if missing:
        raise AttributeError(f"Config file is missing required fields: {missing}")
    else:
        print("Configuration file validated successfully.")
