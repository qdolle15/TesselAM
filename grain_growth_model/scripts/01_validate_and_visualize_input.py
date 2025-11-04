import os
import numpy as np

from grain_growth_model.utils.check_config import ConfigProtocol, validate_config_module
from grain_growth_model.utils.io import load_config_module, parse_arguments
from grain_growth_model.utils.visualization import plot_thermal_history_cross_section, show_domains

def validate_thermal_history(config: ConfigProtocol) -> bool:
    """
    Validate the thermal history configuration.

    Args:
        config: Configuration object.

    Returns:
        bool: True if validation passes, False otherwise.
    """
    START = True

    # Check if meltpool width is larger than simulation width
    try:
        assert np.all(
            np.asarray([config.THERMAL_HISTORY[i]['width'] > config.WIDTH_SIMULATION for i in config.THERMAL_HISTORY.keys()])
        )
    except AssertionError:
        print("Error: The width of the meltpool is larger than the simulation width.")
        START = False

    # Check if corners are well remelted
    try:
        assert np.all(
            np.asarray([
                config.BD_INCREMENTS < np.sqrt(1 - (config.WIDTH_SIMULATION / config.THERMAL_HISTORY[i]['width']) ** 2) * config.THERMAL_HISTORY[i]['height']
                for i in range(1, len(config.THERMAL_HISTORY))
            ])
        )
    except AssertionError:
        print("Error: Areas on the corners aren't well remelted.")
        START = False

    return START

def main():
    # Load and validate config
    args = parse_arguments()
    config: ConfigProtocol = load_config_module(args.config_path)
    validate_config_module(config)

    # Set working directory
    working_paths = {
        "root": args.root,
        "data": args.data,
        "results": args.results,
        "reports": args.reports
    }

    checking_results_dir = os.path.join(working_paths["results"], "checking")
    os.makedirs(checking_results_dir, exist_ok=True)

    # Plot thermal history
    plot_thermal_history_cross_section(
        thermal_layers=config.THERMAL_HISTORY,
        bd_increment=config.BD_INCREMENTS,
        simulation_width=config.WIDTH_SIMULATION,
        z_ultimate=config.Z_ULTIMATE,
        save_path=checking_results_dir,
        save=True
    )

    # Check if thermal history is consistent
    consistent_thermal_history = validate_thermal_history(config=config)
    if not consistent_thermal_history:
        raise ValueError("Thermal history is not consistent.")
    
    # Show domains for each cut view
    for i, cut_view in enumerate(config.CUT_VIEWS):
        domain = cut_view["domain"]
        plans = cut_view["plans"]

        domain_path_dir = os.path.join(
            checking_results_dir, 
            f"domain_{i+1}__{plans[0][0]}.png"
            )

        show_domains(
            thermal_history=config.THERMAL_HISTORY,
            simulation_length=config.LENGTH_SIMULATION,
            simulation_width=config.WIDTH_SIMULATION,
            bd_increment=config.BD_INCREMENTS,
            z_ultimate=config.Z_ULTIMATE,
            domain_path_dir=domain_path_dir,
            save=True,
            domain=domain,
            cut_view=plans,
            global_domain=config.GLOBAL_DOMAIN
        )

    print(f"\nInput validation and visualization complete. Results saved in:\n{checking_results_dir}")


if __name__ == "__main__":
    main()
