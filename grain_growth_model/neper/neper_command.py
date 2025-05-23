import os
import time

def create_tessellation(nbr_seeds: int, domain: dict, domain_path_dir: str):
    """
    Run Neper to create a raster tessellation from filtered seeds.
    """
    x_dom = domain['x_max'] - domain['x_min']
    y_dom = domain['y_max'] - domain['y_min']
    z_dom = domain['z_max'] - domain['z_min']
    cube_size = 1e-3

    tesr_X = int(round(x_dom / cube_size))
    tesr_Y = int(round(y_dom / cube_size))
    tesr_Z = int(round(z_dom / cube_size))

    command = (
        f"neper -T -dim 3 -n {nbr_seeds} "
        f"-domain 'cube({x_dom},{y_dom},{z_dom})' "
        f"-morpho 'voronoi' "
        f"-morphooptiini 'coo:file(./{domain_path_dir}/sub_coo.txt)' "
        f"-ori 'file(./{domain_path_dir}/sub_ori.txt)' "
        f"-format 'tesr' -tesrformat 'ascii' "
        f"-tesrsize '{tesr_X}:{tesr_Y}:{tesr_Z}' "
        f"-o tessellation_3d"
    )

    start = time.time()
    os.system(command)
    elapsed = time.time() - start
    os.system(f"mv tessellation_3d.* {domain_path_dir}")
    
    return command, (tesr_X, tesr_Y, tesr_Z), elapsed

