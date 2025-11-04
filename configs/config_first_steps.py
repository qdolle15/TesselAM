import numpy as np
from grain_growth_model.utils.visualization import generate_segmented_cut_views

# ========== Simulation Parameters ==========

START_SIMULATION = True
RANDOM_SEED = 15
np.random.seed(RANDOM_SEED)

# Dimensions of the simulation domain
LENGTH_SIMULATION = 3.1     # x-axis (mm)
WIDTH_SIMULATION = 0.3        # y-axis (mm)
BD_INCREMENTS = 0.15         # z-step between meltpool layers

# Grain spacing and diameter
D0 = 0.02
D_seeds = 0.3

# Grain competition threshold
THRESHOLD = 1.2

# Add random noise to seeds to avoid singularities
NOISE_NEPER = 0


# Crystallographic directions (e.g. <100>)
EGD_FAMILY = np.array([
    [1,  0,  0],
    [-1, 0,  0],
    [0,  1,  0],
    [0, -1,  0],
    [0,  0,  1],
    [0,  0, -1]
])

# ========== Thermal History ==========

THERMAL_HISTORY = {
    0:{
        'length':21,
        'width':0.35,
        'height':0.2,
        'PD':+1,
        'epitaxy':False
    },
    1:{
        'length':23,
        'width':0.45,
        'height':0.222,
        'PD':-1,
        'epitaxy':False
    },
}

# --------------------------------------------

# Depth limit (for top plane)
MARGIN = 0.3
Z_ULTIMATE = (
    THERMAL_HISTORY[0]["height"] +
    (len(THERMAL_HISTORY) - 1 - MARGIN) * BD_INCREMENTS
)


# ========== Global Domain ==========
GLOBAL_DOMAIN = {
    "x_min": 0.5, "x_max": 2.5,
    "y_min": -.05, "y_max": .05,
    "z_min": 0., "z_max": .32
}

# ========== Generate CUT_VIEWS automatically ==========
CUT_VIEWS = generate_segmented_cut_views(
    global_domain=GLOBAL_DOMAIN,
    axis="x",                 # Split along the x-axis
    num_subdomains=4,         # Create 3 subdomains
    overlap_ratio=0.1,        # 10% overlap between subdomains
    rel_positions=[0.25, 0.5, 0.75]
)
