# Grain Growth Model

**TesselAM** stands for **Tessellation-inspired Simulation for Additive Manufacturing**. 

It is a modular, extensible Python framework for simulating **competitive grain growth** during **directional solidification** such as in **additive manufacturing (AM)** or **welding**. The objective is to deliver **physically motivated**, **statistically meaningful**, and **visualization-ready** predictions of microstructural evolution, while maintaining **computational efficiency**.

This framework is not intended to replace existing models such as **Phase-Field** or **Cellular Automaton** approaches. Instead, **TesselAM** offers a complementary perspective by leveraging **upscaled physical observations** to enable **faster** and **less resource-intensive simulations**, particularly suitable for **in situ microstructure monitoring** and rapid design exploration.

---

## Features

- **3D melt pool modeling** with quarter of ellipsoid, layer by layer.

- **Microstructure as tessellation**: the domain is filled with small elements called *seeds*, each representing a portion of a grain.
  - Each *seed* is defined by **7 degrees of freedom**: 3 coordinates (position), 3 orientation angles (crystallographic), and 1 grain index.
  - *Seeds* belonging to the same grain share the same orientation and index.

- **Grain construction**:
  - *Seeds* are iteratively positioned along the **thermal gradient direction**, starting from an initial position.
  - Growth proceeds as long as the **life expectancy** of the grain permits it (defined by competitive interactions).

- **Two directions of growth** are considered:
  - **Dendrite growth**: along the dendrite’s *easy growth direction* (EGD), intrinsic to its crystallographic orientation.
  - **Grain growth**: driven by the local thermal gradient field.

- **Two growth stages** are modeled:
  1. **Competitive growth stage**:
     - All dendrites grow simultaneously along their EGD.
     - Potential conflicts are detected as **minimal distances between dendrite trajectories** below a user-defined threshold.
     - The **Walton & Chalmers criterion** is used: the dendrite most aligned with the thermal gradient wins the competition and continues to grow.
     - The segment up to the lost conflict defines the **maximum extent (life)** of each grain.
  2. **Grain growth stage**:
     - Grains grow iteratively along the thermal gradient until their life expectancy is reached.

- **Epitaxial growth & CET (Columnar-to-Equiaxed Transition)**:
  - Controlled via configuration (e.g. thermal profile & Hunt criterion).
  - If CET is True: a new interface is initialized with randomly oriented grains.
  - If CET is False: grains that reach the meltpool top can regrow in the next layer.

- **Batch-based conflict resolution**:
  - To handle large NxN combinations in the competitive stage, the domain is sliced along its length and treated **batch by batch** to reduce memory usage.

- **Post-processing**:
  - Automatic export of **coordinates**, **orientations**, and **grain indices** for each layer and interface.

- **EBSD-like visualization**:
  - Generation of 2D cross-sections colored by crystallographic orientation using **Neper** and **ORIX**.

---

## Project Structure

```text
TesselAM/
|
├── configs/                             # Configuration scripts for multiple simulations (Respect name and format of this file)
│   ├── config_1.py
│   ├── config_2.py
│   └── ...                       
|
├── grain_growth_model/
│   ├── __init__.py
│   ├── core/                            # Main algorithms and meltpool geometry
│   ├── analysis/                        # Statistical post-processing
│   ├── neper/                           # Tessellation-based EBSD visualization
│   └── utils/                           # I/O, visualization tools, configuration checks
|
├── outputs/                             # <--- permanent output folder
│   └── name_configuration_file_XXXXX/
│       ├── data/
│       ├── report/
│       └── results/
|
├── main.py                               # Main simulation driver script (root level)
├── requirements.txt                      # List of dependencies
├── setup.py                              # Installable package
└── README.md

```

## Quick Start

### 1. Install dependencies
You can use `pip` with a virtual environment:
```bash
pip install -r requirements.txt
```

### 2. Configure your simulation
Edit or duplicate any file in configs/. Each config defines:
* Simulation domain size
* Meltpool thermal profiles (length, width, depth)
* Growth thresholds
* Crystallographic directions


### 3. Run the model
```bash
python main.py configs/{name_of_the_python_configuration_file}.py
```

### 4. Visualize your results
results saved in 'outputs/name_configuration_file_XXXXX/' containing:
* Report of the simulation
* Report of the visualization part
* sub repositories for each domain with the EBSD like images / IPF triangles / PF

## Dependencies
| Library     | Use                         |
|-------------|-----------------------------|
| `numpy`     | Arrays, math                |
| `scipy`     | Integration, geometry       |
| `matplotlib`| Plotting                    |
| `tqdm`      | Progress bars               |
| `orix`      | Crystal orientation handling (IPF, PF) |
| `neper`     | 3D tessellation visualization (external, see below) |


## Install / Uninstall
* Execute the following command line at the level of setup.py
```bash
pip install -e .
```
* Uninstall through the following command line:
```bash
pip uninstall TesselAM
```
* Check the presence using:
```bash
pip list | TesselAM
```

## EBSD-like Visualization
**Neper** is required to generate grain tessellations.

Official install guide for ubuntu:  
https://neper.info/doc/tutorials/install_ubuntu22.html#installation-ubuntu-22

Once installed, TesselAM will:
- extract seeds in a subdomain of the simulation domain
- run a 3D raster tessellation with Neper
- Reconstruct a 3D arrays in python with id of the voxels corresponding to grain id
- Get the IPF color associated to the Euler-Bunges angles thanks to 'Orix'
- extract 2D planes
- visualize orientations using matplotlib

## License / Citation
Distributed for academic use. Please cite the author or related publication if used in a research project.
DOI of the associated publication: https://doi.org/10.1016/j.commatsci.2024.113112

Citation for the python framework: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17482024.svg)](https://doi.org/10.5281/zenodo.17482024)

## Author
Developed by Quentin Dollé. For questions or contributions, open an issue or contact me directly.
mel: quentin.dolle@polytechnique.edu  
