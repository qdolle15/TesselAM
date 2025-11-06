import os
from typing import Dict
from datetime import datetime
import importlib.util
import argparse
import re
from typing import Tuple, List

# Inputs
def load_config_module(config_path: str):
    """Safely load a Python configuration file."""
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    spec = importlib.util.spec_from_file_location("custom_config", config_path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config

def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Visualize the results.")
    parser.add_argument("config_path", help="Path to the configuration file.")
    parser.add_argument("--root", help="Root directory for the simulation.")
    parser.add_argument("--data", help="Directory for intermediate data.")
    parser.add_argument("--results", help="Directory for visualizations.")
    parser.add_argument("--reports", help="Directory for reports.")
    return parser.parse_args()


# 
def initialize_output_directory(label: str, base_folder: str = "outputs", force_new: bool = False) -> Dict[str, str]:
    """
    Create and return a directory for the current simulation.

    If force_new is False and the directory does not exist, a FileNotFoundError is raised.
    If force_new is True, a new directory is created, possibly with an incremented suffix.

    Parameters:
        label (str): Name of the simulation (e.g., config file name).
                     If None, a timestamped folder is created.
        base_folder (str): Root folder for all simulations (default: "outputs").
        force_new (bool): If True, create a new directory even if the label already exists.
                          If False, reuse the existing directory or raise an error if it does not exist.

    Returns:
        dict: Dictionary with paths to the created directories:
              - "root": Root directory for the simulation.
              - "data": Directory for intermediate data (e.g., `.npz` files).
              - "results": Directory for visualizations (e.g., images).
              - "reports": Directory for reports and metrics.

    Raises:
        FileNotFoundError: If force_new is False and the directory does not exist.
    """
    # Create the base folder if it doesn't exist
    os.makedirs(base_folder, exist_ok=True)

    # Initialize the output root path
    output_root = os.path.join(base_folder, label)

    if force_new:
        # If force_new is True, find a new directory name
        counter = 1
        while os.path.exists(output_root):
            output_root = os.path.join(base_folder, f"{label}_{counter:04d}")
            counter += 1
    else:
        # If force_new is False, check if the directory exists
        if not os.path.exists(output_root):
            raise FileNotFoundError(f"The directory {output_root} does not exist.")

    # Create the root directory if it doesn't exist (only if force_new is True or directory already exists)
    os.makedirs(output_root, exist_ok=True)

    # Create subdirectories
    subfolders = ["data", "results"]
    paths = {"root": output_root}
    for subfolder in subfolders:
        path = os.path.join(output_root, subfolder)
        os.makedirs(path, exist_ok=True)
        paths[subfolder] = path

    return paths

def get_domain_output_directory(root_results_dir: str, domain_index: int, plane: str) -> str:
    """
    Create a subdirectory for a specific domain's results.

    Args:
        root_results_dir: Path to the root results directory.
        domain_index: Index of the domain (e.g., 1, 2, 3).
        plane: Cut plane (e.g., "XZ", "YZ").

    Returns:
        str: Path to the domain-specific subdirectory.
    """
    domain_dir = os.path.join(root_results_dir, f"domain_{domain_index}__{plane}")
    os.makedirs(domain_dir, exist_ok=True)
    return domain_dir


# Reports
def generate_simulation_report(data: dict, config_module, output_file: str):
    """
    Generate a well-formatted plain text report for FAST-MMAM, including configuration,
    per-layer statistics, and competitive growth metrics.

    Parameters
    ----------
    data : dict
        Nested dictionary with per-layer tracking information.
    config_module : module
        Configuration module used for the simulation.
    output_file : str
        Path to the output .txt report.
    """
    def section(title):
        return f"\n\t--------\n\t{title}\n\t--------\n"

    def subsection(title):
        return f"{title}:\n" + "-" * (len(title)+1) + "\n"

    with open(output_file, "w") as f:
        f.write("# TesselAM Simulation Report\n")
        f.write("=" * 29 + "\n\n")
        f.write(f"{'Date simulation':<37}: {data['date start simulation']}\n")
        f.write(f"{'Elapsed time:':<37}: {data['tot. simulation time']} sec.\n")
        f.write(f"{'Number of layers':<37}: {len(config_module.THERMAL_HISTORY)}\n")
        f.write(f"{'Simulation length':<37}: {config_module.LENGTH_SIMULATION} mm\n")
        f.write(f"{'Simulation width':<37}: {config_module.WIDTH_SIMULATION} mm\n")
        f.write(f"{'Building direction increment':<37}: {config_module.BD_INCREMENTS} mm\n")
        f.write(f"{'Mean grain size (D0)':<37}: {config_module.D0} mm\n")
        f.write(f"{'Threshold for conflict detection':<37}: {config_module.THRESHOLD}\n")
        f.write(f"{'Randomness number':<37}: {config_module.RANDOM_SEED}\n")
        f.write("-" * 29 + "\n")
        f.write(f"{'Number of grains into the substrate':<37}: {data['No substrate grains generated']:,}".replace(",", " ") + "\n")
        f.write(f"{'Number of grains created':<37}: {data["No of grains generated"]:,}".replace(",", " ") + " (Except substrate)\n")
        f.write(f"{'Number of seeds created':<37}: {data["No of seeds generated"]:,}".replace(",", " ") + "\n")           
        f.write("=" * 29 + "\n\n")

        # Per layers
        for layer_idx in sorted(k for k in data if isinstance(k, int)):
            layer_data = data[layer_idx]

            f.write(section(f"Layer {layer_idx + 1}"))

            # --- Thermal activity
            thermal = config_module.THERMAL_HISTORY[layer_idx]
            f.write(subsection("Meltpool description"))
            f.write(f"{'length':<20}: {thermal['length']:>6} mm\n")
            f.write(f"{'width':<20}: {thermal['width']:>6} mm\n")
            f.write(f"{'depth':<20}: {thermal['height']:>6} mm\n")
            f.write(f"{'epitaxial growth':<20}:    {thermal['epitaxy']}\n")
            f.write(f"{'direction':<20}:    {'+x' if thermal['PD'] == 1 else '-x'}\n\n")

            # --- Grains
            f.write(subsection("Grains"))
            f.write(f"{'grains from the previous layer':<35}: {layer_data['grains for epitaxy']:,}".replace(",", " ") + "\n")
            f.write(f"{'grains nucleated':<35}: {layer_data['grains nucleated']:,}".replace(",", " ") + "\n")
            f.write(f"{'grains for growing':<35}: {layer_data['total grains for growth']:,}".replace(",", " ") + "\n")
            f.write(f"{'grains reaching top':<35}: {layer_data['reached_top']:,}".replace(",", " ") + "\n\n")

            # --- Conflicts
            f.write(subsection("Competitive growth"))
            f.write(f"batch size : {layer_data['batch size']:,}".replace(",", " ") + "\n")
            f.write(f"overlap ratio between batches : {100*layer_data['overlap ratio']:.0f} %\n")
            f.write(f"length of batch : {layer_data['length of one batch']} mm\n")
            f.write(f"number of batches: {layer_data['number of batches']}\n\n")
            # f.write(f"Nombre longueurs dendrites == 0 : {layer_data.get('null_lengths', 'N/A')}\n")
            # f.write(f"Nombre longueurs dendrites inchangés : {layer_data.get('unchanged_lengths', 'N/A')}\n")

            
            # --- Statistics
            f.write(subsection("Statistics"))
            f.write(f"Elapsed time for layer : {layer_data['layer simu. duration']} sec.\n")
            # if "statistics" in layer_data:
            #     f.write("+--------------------------------+------------------+\n")
            #     for key, val in layer_data["statistics"].items():
            #         f.write(f"| {key:<30} | {val:<16} |\n")
            #     f.write("+--------------------------------+------------------+\n\n")

def generate_visualization_report(output_file: str, data: dict, time_information: dict, cut_views: list[dict]):
    """
    Save a structured summary of the visualization settings and domains.

    Parameters
    ----------
    output_file : str
        Path to the output text file.
    data : dict
        Dictionary containing visualization metadata or results per domain.
    time_information : dict
        Dictionnary containing the time elapsed for all the visualization steps and for each domain
    cut_views : list of dict
        List of domain + plane/position selections. Each item must have:
            - "domain": dict with keys x_min, x_max, y_min, y_max, z_min, z_max
            - "plans": list of tuples (plane, relative position)
    """

    def box_title(text: str) -> str:
        width = 16
        top = "+" + "-" * (width - 2) + "+"
        mid = "|" + text.center(width - 2) + "|"
        return f"{top}\n{mid}\n{top}"


    with open(output_file, 'w') as f:
        f.write("# TesselAM Visualization Report\n\n")
        f.write(f"Total domains visualized: {len(cut_views)}\n")
        f.write(f"Total time elapsed: {time_information['total']} s\n\n")

        for idx, cut_view in enumerate(cut_views):
            domain = cut_view["domain"]
            planes = cut_view["plans"]
            f.write(box_title(f"Domain {idx + 1}"))

            f.write(f"\nTime elapsed for the domain: {time_information[idx]} s\n")

            f.write("\nSpatial bounds:\n")
            for key in ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]:
                f.write(f"  {key:<10}: {domain.get(key, 'N/A')}\n")

            f.write("\nRequested cut views:\n")
            for plane, value in planes:
                f.write(f"  - {plane} at relative position {value:.2f}\n")

            # Optional: write info from EBSD_like or Neper results
            if idx in data:
                info = data[idx]
                f.write("\nTessellation Info:\n")
                f.write(f"  Status               : {'Aborted due to size file' if not info['status'] else 'Clear'}\n")
                f.write(f"  Number of voxels     : {info.get('voxels', 'N/A')}\n")
                f.write(f"  Number of seeds      : {info.get('N', 'N/A')}\n")
                f.write(f"  Number of grains/layer: {info.get('grains_per_layers', 'N/A')}\n")
                f.write(f"  Elapsed tessellation time: {round(info.get('time_tess', 0), 3)} s\n")
                f.write(f"  Command used         :\n    {info.get('cmd', '').strip()}\n")

            f.write("\n")


# Compare config and reports
def compare_simulation_data(report_path: str, config_path: str) -> Tuple[bool, List[str]]:
    """
    Compare simulation data between a report text file and a configuration Python file.
    Returns a tuple with a boolean indicating if all data match and a list of error messages.

    Args:
        report_path (str): Path to the report text file.
        config_path (str): Path to the configuration Python file.

    Returns:
        Tuple[bool, List[str]]: A tuple with a boolean indicating if all data match and a list of error messages.
    """
    if not os.path.exists(report_path):
        return False, [f"Report file not found: {report_path}"]

    # Load the configuration file
    spec = importlib.util.spec_from_file_location("config", config_path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)

    # Read the report file
    with open(report_path, 'r') as file:
        report_content = file.read()

    errors = []

    try:
        # Extract relevant data from the report
        report_data = {}
        try:
            report_data['Simulation length'] = extract_value(report_content, r'Simulation length\s*:\s*([\d.]+)\s*mm')
        except ValueError as e:
            errors.append(str(e))

        try:
            report_data['Simulation width'] = extract_value(report_content, r'Simulation width\s*:\s*([\d.]+)\s*mm')
        except ValueError as e:
            errors.append(str(e))

        try:
            report_data['Building direction increment'] = extract_value(report_content, r'Building direction increment\s*:\s*([\d.]+)\s*mm')
        except ValueError as e:
            errors.append(str(e))

        try:
            report_data['Mean grain size (D0)'] = extract_value(report_content, r'Mean grain size \(D0\)\s*:\s*([\d.]+)\s*mm')
        except ValueError as e:
            errors.append(str(e))

        try:
            report_data['Threshold for conflict detection'] = extract_value(report_content, r'Threshold for conflict detection\s*:\s*([\d.]+)')
        except ValueError as e:
            errors.append(str(e))

        try:
            report_data['Randomness number'] = extract_value(report_content, r'Randomness number\s*:\s*([\d.]+)')
        except ValueError as e:
            errors.append(str(e))

        if not errors:
            # Compare relevant data
            if abs(report_data['Simulation length'] - config.LENGTH_SIMULATION) > 1e-6:
                errors.append(f"Simulation length mismatch: Report={report_data['Simulation length']}, Config={config.LENGTH_SIMULATION}")

            if abs(report_data['Simulation width'] - config.WIDTH_SIMULATION) > 1e-6:
                errors.append(f"Simulation width mismatch: Report={report_data['Simulation width']}, Config={config.WIDTH_SIMULATION}")

            if abs(report_data['Building direction increment'] - config.BD_INCREMENTS) > 1e-6:
                errors.append(f"Building direction increment mismatch: Report={report_data['Building direction increment']}, Config={config.BD_INCREMENTS}")

            if abs(report_data['Mean grain size (D0)'] - config.D0) > 1e-6:
                errors.append(f"Mean grain size (D0) mismatch: Report={report_data['Mean grain size (D0)']}, Config={config.D0}")

            if abs(report_data['Threshold for conflict detection'] - config.THRESHOLD) > 1e-6:
                errors.append(f"Threshold for conflict detection mismatch: Report={report_data['Threshold for conflict detection']}, Config={config.THRESHOLD}")

            if abs(report_data['Randomness number'] - config.RANDOM_SEED) > 1e-6:
                errors.append(f"Randomness number mismatch: Report={report_data['Randomness number']}, Config={config.RANDOM_SEED}")

    except Exception as e:
        errors.append(f"Unexpected error during comparison: {e}")

    return (len(errors) == 0, errors)

def extract_value(content: str, pattern: str) -> float:
    """
    Extract a value from text content using a regex pattern.
    Handles values with units (e.g., "3 mm" -> 3).

    Args:
        content (str): Text content to search.
        pattern (str): Regex pattern to find the value.

    Returns:
        float: Extracted value.

    Raises:
        ValueError: If the value cannot be found.
    """
    match = re.search(pattern, content)
    if not match:
        raise ValueError(f"Value not found for pattern: {pattern}")

    # Extract the numeric part, ignoring any units
    value_str = match.group(1)
    numeric_value = re.sub(r'[^\d.]', '', value_str)  # Remove non-numeric characters
    return float(numeric_value)