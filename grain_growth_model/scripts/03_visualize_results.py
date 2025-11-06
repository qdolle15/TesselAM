import os
import sys
import time
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

from grain_growth_model.neper.runner import EBSD_like
from grain_growth_model.utils.visualization import show_domains
from grain_growth_model.utils.check_config import ConfigProtocol, validate_config_module
from grain_growth_model.neper.neper_visualization import load_images_from_folders, manual_stitch
from grain_growth_model.utils.io import load_config_module, parse_arguments, compare_simulation_data, generate_visualization_report

def main():

    try:
        # Load and validate config
        args = parse_arguments()
        config: ConfigProtocol = load_config_module(args.config_path)
        validate_config_module(config)

        # Set working directory
        working_paths = {
            "root": args.root,
            "data": args.data,
            "results": args.results,
        }

        # Check if the report file exists
        report_path = os.path.join(working_paths["data"], "report_simulation.txt")
        if not os.path.exists(report_path):
            raise FileNotFoundError(f"Report file not found: {report_path}")

        # Compare simulation data
        success, errors = compare_simulation_data(report_path=report_path, config_path=args.config_path)
        if not success:
            print("\nDifferences between the config and the report:")
            for error in errors:
                print(f"  - {error}")
            raise ValueError("\nSimulation results may not match with what is asked.\n")

        # Create sub repository
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = os.path.join(working_paths["results"], f"results_{timestamp}")
        os.makedirs(results_dir, exist_ok=True)

        # Visualization for each domain
        time_start_visualization_module = time.time()
        visualization_results = {}
        time_information = {}

        for i, cut_view in enumerate(config.CUT_VIEWS):

            time_start_layer_visualization = time.time()

            domain = cut_view["domain"]
            plans = cut_view["plans"]

            # Create a subdirectory for this domain
            domain_path_dir = os.path.join(results_dir, f"domain_{i+1}__{plans[0][0]}")
            os.makedirs(domain_path_dir, exist_ok=True)

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

            # Execute EBSD_like and save the results
            information = EBSD_like(
                No_layer=len(config.THERMAL_HISTORY),
                domain=domain,
                cut_view=plans,
                save=True,
                data_path_dir=working_paths["data"],
                domain_path_dir=domain_path_dir
            )
            visualization_results[i] = information
            time_information[i] = round(time.time() - time_start_layer_visualization, 2)

        time_information['total'] = round(time.time() - time_start_visualization_module, 2)
        print(f"\nVisualization complete. Results saved in: {results_dir}")

        # Generate reports
        report_file = os.path.join(results_dir, "visualization.txt")
        generate_visualization_report(
            data=visualization_results,
            cut_views=config.CUT_VIEWS,
            output_file=report_file,
            time_information=time_information
        )

        # Stiching Domains
        print(f"\nStitchin images from {config.N_SUBDOMAINS} domains")
        for cut_view in ['025', '050', '075']:
            images = load_images_from_folders(
                base_path=results_dir,
                num_subdomains=config.N_SUBDOMAINS,
                filename=f"EBSD_z_XZ_{cut_view}.npy",
                )
            stitched_images = manual_stitch(images=images)

            np.save(os.path.join(results_dir, f"stitched_image_{cut_view}.npy"), stitched_images)
            plt.imsave(os.path.join(results_dir, f"stitched_image_{cut_view}.png"), stitched_images)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
