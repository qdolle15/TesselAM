import os
import numpy as np
from datetime import datetime


def initialize_output_directory(label=None, base_folder="outputs"):
    """
    Create and return a fresh output directory for the current simulation.

    Parameters
    ----------
        label (str or None) Optional subfolder name. If None, a timestamped folder is created.
        base_folder (str): Root permanent folder for all simulations.

    Returns
    -------
        (dict): Dictionary with all created paths (root, coordinates/, etc.)
    """
    os.makedirs(base_folder, exist_ok=True)

    if label is None:
        label = "simulation"

    # Base timestamp
    base_name = f"{label}_0001"
    output_root = os.path.join(base_folder, base_name)

    # Check for existing folders with similar names and increment if needed
    counter = 2
    while os.path.exists(output_root):
        base_name = f"{label}_{counter:04d}"
        output_root = os.path.join(base_folder, base_name)
        counter += 1

    output_root = os.path.join(base_folder, base_name)
    os.makedirs(output_root)

    subfolders = ["data_layers", "preliminary_results", 'results']
    paths = {"root": output_root}

    for sub in subfolders:
        path = os.path.join(output_root, sub)
        os.makedirs(path, exist_ok=True)
        paths[sub] = path

    return paths


def generate_report(data: dict, config_module, output_file: str):
    """
    Generate a plain text report including global configuration and simulation statistics.

    Parameters
    ----------
    data : dict
        Nested dictionary with numerical tracking information (global and per-layer).
    config_module : module
        Configuration module used for the simulation (from load_config_module).
    output_file : str
        Path to the output text file.
    """
    with open(output_file, 'w') as f:
        f.write("# FAST-MMAM Simulation Report\n\n")

        # Configuration summary
        f.write("## Configuration Parameters\n")
        config_keys = [
            "LENGTH_SIMULATION", "WIDTH_SIMULATION", "D0", "D_seeds",
            "BD_INCREMENTS", "Z_ULTIMATE", "START_SIMULATION",
        ]
        for key in config_keys:
            value = getattr(config_module, key, "N/A")
            f.write(f"{key:<25}: {value}\n")

        # Extra fields (thermal layers, EGD family)
        f.write(f"\nTHERMAL_HISTORY (n layers): {len(config_module.THERMAL_HISTORY)}\n")
        f.write(f"EGD_FAMILY shape            : {getattr(config_module, 'EGD_FAMILY', np.array([])).shape}\n")
        f.write(f"DOMAINS defined             : {len(getattr(config_module, 'DOMAINS', []))}\n")
        f.write(f"CUT_VIEWS                   : {getattr(config_module, 'CUT_VIEWS', [])}\n")

        f.write("\n## Simulation Summary\n")
        for key, val in data.items():
            if isinstance(val, dict):
                f.write(f"\n--- Layer {key} ---\n")
                for subkey, subval in val.items():
                    f.write(f"{subkey:<30}: {subval}\n")
            else:
                f.write(f"{key:<30}: {val}\n")


def generate_visualization_report(output_file, data, domains, plans):
    """
    Save a summary of the visualization settings and selections.

    Parameters
    ----------
    output_file : str
        Output file path.
    data : dict
        Dictionary describing visualized domains.
    domains : list of dict
        Spatial selection areas.
    plans : list of tuple
        Cut views requested (e.g., 'XZ', 0.5).
    """
    with open(output_file, 'w') as f:
        f.write("# FAST-MMAM Visualization Summary\n\n")
        for idx, domain in data.items():
            f.write(f"--- Domain {idx + 1} ---\n")
            for key, val in domains[idx].items():
                f.write(f"{key:<20}: {val}\n")
            f.write("Cut views:\n")
            for plane, value in plans:
                f.write(f"  {plane} at {value}\n")
            f.write("\n")
