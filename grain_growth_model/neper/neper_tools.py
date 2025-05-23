import os
import re
import numpy as np

def create_sub_selection(No_layer:int, domain: dict, data_path: str, domain_path_dir: str):
    """
    Filter all seeds to keep only those inside the analysis domain.

    Returns
    -------
    int : Total number of filtered seeds
    np.ndarray : Grain count per layer
    """
    coo_filename = f"{domain_path_dir}/sub_coo.txt"
    ori_filename = f"{domain_path_dir}/sub_ori.txt"
    id_filename = f"{domain_path_dir}/sub_id.txt"

    grains_per_layers = []
    N_seeds = 0

    with open(coo_filename, 'w') as final_coo, open(ori_filename, 'w') as final_ori, open(id_filename, 'w') as final_id:
        for lay in np.arange(-1, No_layer):
            if lay == -1:
                substrate_data = np.load(os.path.join(data_path, f"substrate.npz"))
                coo = substrate_data['coords']
                ori = substrate_data['oris']
                ids = substrate_data['indexes']
            else:
                layer_data = np.load(os.path.join(data_path, f"layer_{lay}.npz"))
                coo = layer_data['coords']
                ori = layer_data['oris']
                ids = layer_data['indexes']

            selection = (
                (coo[:, 0] > domain['x_min']) & (coo[:, 0] < domain['x_max']) &
                (coo[:, 1] > domain['y_min']) & (coo[:, 1] < domain['y_max']) &
                (coo[:, 2] > domain['z_min']) & (coo[:, 2] < domain['z_max'])
            )

            sub_coo = np.vstack((coo[selection, 0] - domain['x_min'],
                                 coo[selection, 1] - domain['y_min'],
                                 coo[selection, 2] - domain['z_min'])).T
            np.savetxt(final_coo, sub_coo)
            np.savetxt(final_ori, ori[selection])
            np.savetxt(final_id, ids[selection])
            
            grains_per_layers.append(np.unique(ids[selection]).size)
            N_seeds += len(sub_coo)

    return N_seeds, np.asarray(grains_per_layers)


def read_data(file_path, start_marker='**data', end_marker='***end'):
    with open(file_path, 'r') as f:
        content = f.read()
    pattern = re.compile(f'{re.escape(start_marker)}(.*?){re.escape(end_marker)}', re.DOTALL)
    match = pattern.search(content)
    if not match:
        return []
    return match.group(1).strip().split('\n')


def transform_data(data):
    data_flat = []
    for line in data:
        if line != 'ascii':
            data_flat.extend([int(x) for x in line.split()])
    return np.array(data_flat, dtype=np.uint32)
