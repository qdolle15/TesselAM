import os
import sys
import subprocess
import argparse
from datetime import datetime
from typing import Dict

from grain_growth_model.utils.io import initialize_output_directory

def run_script(script_name: str, config_path: str, working_paths: Dict[str, str]) -> None:
    """
    Execute a Python script with the working paths.

    Args:
        script_name (str): Name of the script to run.
        config_path (str): Path to the configuration file.
        working_paths (Dict[str, str]): Dictionary containing paths to the working directories.
    """
    project_root = os.getcwd()
    script_path = os.path.join(project_root, "grain_growth_model", "scripts", script_name)
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"Script not found: {script_path}")

    # Pass the working paths as arguments
    cmd = [sys.executable, script_path, config_path]
    for key, path in working_paths.items():
        cmd.extend([f"--{key}", path])

    subprocess.run(cmd, check=True)

def main():
    """
    Main function to execute the simulation workflow.
    """
    parser = argparse.ArgumentParser(description="Run simulation workflow.")
    parser.add_argument("-m", "--mode", required=True,
                        help="Mode: C (check), S (simulation), V (visualization), CSV (all)")
    parser.add_argument("-f", "--file", required=True, help="Path to the configuration file")
    parser.add_argument("-o", "--output", help="Name for the output directory")
    args = parser.parse_args()

    # Determine the output directory name
    if args.output is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_label = f"simulation_{timestamp}"
    else:
        output_label = args.output

    try:
        # Initialize the output directory based on the mode
        if 'C' in args.mode:
            paths_working_dir = initialize_output_directory(label=output_label, force_new=True)
        else:
            paths_working_dir = initialize_output_directory(label=output_label, force_new=False)

        # Run the scripts according to the mode
        if 'C' in args.mode:
            run_script("01_validate_and_visualize_input.py", args.file, paths_working_dir)
        if 'S' in args.mode:
            run_script("02_run_simulation.py", args.file, paths_working_dir)
        if 'V' in args.mode:
            run_script("03_visualize_results.py", args.file, paths_working_dir)

    except subprocess.CalledProcessError as e:
        print(f"Error: Script execution failed with error: {e}")
        sys.exit(1)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()


# Exemples:
# python3 main.py -m C -f configs/config_article.py
# python3 main.py -m CS -f configs/config_article.py
# python3 main.py -m V -f configs/config_article.py
# python3 main.py -m CSV -f configs/config_article.py
# python3 main.py -m CSV -f configs/config_article.py -o ma_simulation