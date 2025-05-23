import os


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

    subfolders = ["data", "results", 'reports']
    paths = {"root": output_root}

    for sub in subfolders:
        path = os.path.join(output_root, sub)
        os.makedirs(path, exist_ok=True)
        paths[sub] = path

    return paths


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
        f.write("# FAST-MMAM Simulation Report\n")
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


def generate_visualization_report(output_file: str, data: dict, cut_views: list[dict]):
    """
    Save a structured summary of the visualization settings and domains.

    Parameters
    ----------
    output_file : str
        Path to the output text file.
    data : dict
        Dictionary containing visualization metadata or results per domain.
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
        f.write("# FAST-MMAM Visualization Report\n\n")
        f.write(f"Total domains visualized: {len(cut_views)}\n\n")

        for idx, cut_view in enumerate(cut_views):
            domain = cut_view["domain"]
            planes = cut_view["plans"]
            f.write(box_title(f"Domain {idx + 1}"))

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
                f.write(f"  Number of voxels     : {info.get('voxels', 'N/A')}\n")
                f.write(f"  Number of seeds      : {info.get('N', 'N/A')}\n")
                f.write(f"  Number of grains/layer: {info.get('grains_per_layers', 'N/A')}\n")
                f.write(f"  Elapsed tessellation time: {round(info.get('time_tess', 0), 3)} s\n")
                f.write(f"  Command used         :\n    {info.get('cmd', '').strip()}\n")

            f.write("\n")