# source .venv/bin/activate

import os
import sys
import time
import numpy as np
import warnings
import importlib.util
from tqdm import tqdm

from grain_growth_model.utils.check_config import ConfigProtocol, validate_config_module
from grain_growth_model.utils.io import initialize_output_directory, generate_simulation_report, generate_visualization_report
from grain_growth_model.utils.visualization import plot_thermal_history_cross_section, show_domains, visualize_seeds_at_interface, visualize_3d_growth_directions
from grain_growth_model.core.grains import generate_initial_substrate_grains, initialize_interface_grains
from grain_growth_model.core.meltpool import compute_normalized_gradient_at_point
from grain_growth_model.core.growth import compute_easy_growth_directions, simulate_grain_growth
from grain_growth_model.core.competition import split_seeds_into_batches, resolve_growth_conflicts, resolve_batch_overlaps
from grain_growth_model.neper.runner import EBSD_like

np.seterr(divide='ignore', invalid='ignore')
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, message='All-NaN slice encountered')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in sqrt')


def load_config_module(config_path:str):
    """Safely load a Python configuration file."""
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    spec = importlib.util.spec_from_file_location("custom_config", config_path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def main(config_path):

    # Load user configuration
    config: ConfigProtocol = load_config_module(config_path)
    validate_config_module(config)

    if not config.START:
        print("Simulation can't start due to config file")
        return
    

    # Extract label from config filename
    label = os.path.splitext(os.path.basename(config_path))[0]
    working_paths = initialize_output_directory(label=label, base_folder="outputs")

    time_start_simulation = time.time()
    debug = False
    tracking = {}
    tracking['date start simulation'] = time.ctime(time.time())
    total_grains_created = 0
    total_seeds_created = 0
    laser_z = config.THERMAL_HISTORY[0]["height"]

    plot_thermal_history_cross_section(
                thermal_layers=config.THERMAL_HISTORY, 
                bd_increment=config.BD_INCREMENTS,
                simulation_width=config.WIDTH_SIMULATION,
                z_ultimate=config.Z_ULTIMATE,
                save_path=working_paths["results"],
                save=True
    )

    # Substrate management
    substrate_coords, substrate_oris, substrate_ids = generate_initial_substrate_grains(
        thermal_layers=config.THERMAL_HISTORY[0],
        simulation_length=config.LENGTH_SIMULATION,
        simulation_width=config.WIDTH_SIMULATION,
        mean_grain_diameter=config.D0
        )
    np.savez(
        file=os.path.join(working_paths["data"], "substrate.npz"), 
        coords=substrate_coords,
        oris=substrate_oris,
        indexes=substrate_ids
        )
    tracking['No substrate grains generated'] = len(substrate_ids)

    # Beads management
    for layer_idx in range(len(config.THERMAL_HISTORY)):
        print(f"\n=== Processing Layer {layer_idx} ===")
        tracking[layer_idx] = {}
        time_start_layer_simulation = time.time()

        layer_meltpool_description: dict = config.THERMAL_HISTORY[layer_idx]
        is_epitaxy = layer_meltpool_description.get("epitaxy", False)

        if is_epitaxy and layer_idx > 0:
            interface_data = np.load(os.path.join(working_paths["data"], f"interface_{layer_idx-1}_{layer_idx}.npz"))
            interface_coords = interface_data['coords']
            interface_oris = interface_data['oris']
            interface_ids = interface_data['indexes']
        else:
            interface_coords, interface_oris, interface_ids = None, None, np.array([], dtype=int)

        extended_interface_coords, extended_interface_oris, extended_interface_ids, \
            grains_generated = initialize_interface_grains(
                thermal_layers=config.THERMAL_HISTORY,
                layer_index=layer_idx,
                known_coords=interface_coords,
                known_orientations=interface_oris,
                known_ids=interface_ids,
                simulation_length=config.LENGTH_SIMULATION,
                simulation_width=config.WIDTH_SIMULATION,
                mean_grain_diameter=config.D0,
                bd_increment=config.BD_INCREMENTS,
                z_ultimate=config.Z_ULTIMATE,
                No_grains_created=total_grains_created
        )

        if debug:
            visualize_seeds_at_interface(
                thermal_layers=config.THERMAL_HISTORY[0], 
                coo_substrate=substrate_coords,
                coo_interface=extended_interface_coords,
                width_simulation=config.WIDTH_SIMULATION,
                length_simulation=config.LENGTH_SIMULATION,
                save=True,
                save_path=working_paths["results"]
            )

        # Creating indexes that have never been attributed before for following \
        # potential grain growth evolution across layers. 
        total_grains_created += grains_generated
        tracking[layer_idx]['grains for epitaxy'] = len(interface_ids)
        tracking[layer_idx]['grains nucleated'] = grains_generated
        tracking[layer_idx]['total grains for growth'] = len(extended_interface_ids)

        extended_gradients_at_interface = compute_normalized_gradient_at_point(
            y=extended_interface_coords[:, 1],
            z=extended_interface_coords[:, 2],
            meltpool=layer_meltpool_description,
            z_center=laser_z
        )

        batch_size=6000
        overlap_ratio=0.3
        batch_len, batch_data = split_seeds_into_batches(
            size_limit=batch_size,
            overlap_ratio=overlap_ratio,
            simulation_length=config.LENGTH_SIMULATION,
            coordinates_seeds=extended_interface_coords,
            euler_angles=extended_interface_oris,
            grain_ids=extended_interface_ids,
            gradients=extended_gradients_at_interface,
        )

        tracking[layer_idx]['batch size'] = batch_size
        tracking[layer_idx]['overlap ratio'] = overlap_ratio
        tracking[layer_idx]['length of one batch'] = batch_len
        tracking[layer_idx]['number of batches'] = len(batch_data)

        for batch_id, data in tqdm(batch_data.items(), desc="Batches"):

            sub_extended_interface_coords = data["coordinates"]
            sub_extended_interface_oris = data["orientations"]
            sub_extended_interface_thermal_gradient = data["gradient interface"]
            sub_extended_interface_id = data["id"]

            sub_dendrite_growth_dirs = compute_easy_growth_directions(
                euler_angles=sub_extended_interface_oris,
                thermal_gradient=sub_extended_interface_thermal_gradient,
                crystal_directions=config.EGD_FAMILY,
                y_positions=sub_extended_interface_coords[:, 1],
                simulation_width=config.WIDTH_SIMULATION
            )
            data["growth direction"] = sub_dendrite_growth_dirs

            if debug:
                visualize_3d_growth_directions(
                    display_ratio=0.05,
                    substrate_coords=None,
                    interface_coords=sub_extended_interface_coords,
                    growth_time=None,
                    gradient_vectors=sub_extended_interface_thermal_gradient,
                    growth_directions=sub_dendrite_growth_dirs,
                    arrow_length=0.1
                )

            lengths, mask, conflict_history = resolve_growth_conflicts(
                growth_directions=sub_dendrite_growth_dirs,
                coordinates_seeds=sub_extended_interface_coords,
                global_id=sub_extended_interface_id,
                thermal_layers=config.THERMAL_HISTORY,
                layer_index=layer_idx,
                meltpool_z=laser_z,
                threshold=config.THRESHOLD,
                mean_grain_diameter=config.D0,
                bd_increment=config.BD_INCREMENTS,
                simulation_width=config.WIDTH_SIMULATION,
                z_ultimate=config.Z_ULTIMATE,
            )
            data["length dendrites"] = lengths

        dendrite_lengths, extended_dendrite_growth_dirs = resolve_batch_overlaps(
            data=batch_data,
            batch_size=batch_len,
            coords_all=extended_interface_coords,
            ids_all=extended_interface_ids,
            simulation_length=config.LENGTH_SIMULATION
        )

        reached_top, new_seeds = simulate_grain_growth(
            coordinates_seeds=extended_interface_coords,
            growth_directions=extended_dendrite_growth_dirs,
            grain_ids=extended_interface_ids,
            gradients=extended_gradients_at_interface,
            euler_angles=extended_interface_oris,
            dendrite_lengths=dendrite_lengths,
            thermal_layers=config.THERMAL_HISTORY,
            current_z=laser_z,
            simulation_length=config.LENGTH_SIMULATION,
            z_ultimate=config.Z_ULTIMATE,
            inter_seeds_distance=config.D_seeds,
            layer_index=layer_idx,
            noise_neper=config.NOISE_NEPER,
            saving_path=working_paths["data"]
        )

        total_seeds_created += new_seeds
        tracking[layer_idx]["reached_top"] = reached_top
        tracking[layer_idx]["layer simu. duration"] = round(time.time() - time_start_layer_simulation, 2)
        laser_z += config.BD_INCREMENTS

    tracking["No of grains generated"] = total_grains_created
    tracking["No of seeds generated"] = total_seeds_created

    # Results visualization
    visualization_results = {}
    for i, cut_view in enumerate(config.CUT_VIEWS):

        domain = cut_view["domain"]
        plans = cut_view["plans"]

        domain_path_dir = os.path.join(working_paths["results"], f"domain_{i+1}__{plans[0][0]}")
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
            cut_view=plans
            )
        information = EBSD_like(
            No_layer=len(config.THERMAL_HISTORY),
            domain=domain, 
            cut_view=plans,
            save=True,
            data_path_dir=working_paths['data'], 
            domain_path_dir=domain_path_dir
            )

        visualization_results[i] = information

    simulation_duration = round(time.time() - time_start_simulation, 2)
    tracking['tot. simulation time'] = simulation_duration

    # Writing reports
    generate_simulation_report(
        data=tracking, 
        config_module=config, 
        output_file=os.path.join(working_paths["reports"], "report_simulation.txt"),
        )

    generate_visualization_report(
        data=visualization_results,
        cut_views=config.CUT_VIEWS,
        output_file=os.path.join(working_paths["reports"], "report_visualization.txt"),
    )

    print(f"\nSimulation complete. Results saved in: {working_paths["root"]}")
    print(f"\nElapsed time: {simulation_duration} seconds")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please provide the path to a config file.")
        print("Usage: python main.py configs/config_development.py")
        sys.exit(1)

    main(config_path=sys.argv[1])
    # main(config_path="./configs/config_development.py")
    