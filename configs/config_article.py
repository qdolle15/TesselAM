import numpy as np
from grain_growth_model.utils.visualization import generate_segmented_cut_views

# ========== Simulation Parameters ==========

START_SIMULATION = True
RANDOM_SEED = 15
np.random.seed(RANDOM_SEED)

# Dimensions of the simulation domain
LENGTH_SIMULATION = 1.5       # x-axis (mm)
WIDTH_SIMULATION = 3        # y-axis (mm)
BD_INCREMENTS = 0.81         # z-step between meltpool layers

# Grain spacing and diameter
D0 = 0.01
D_seeds = 0.15

# Grain competition threshold
THRESHOLD = 1.

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
        'length':24.38,
        'width':10,
        'height':3.06,
        'PD':+1,
        'epitaxy':True
    },
    1:{
        'length':36.67,
        'width':10,
        'height':2.45,
        'PD':-1,
        'epitaxy':True
    },
    2:{
        'length':24.96,
        'width':10,
        'height':2.91,
        'PD':+1,
        'epitaxy':True
    },
    3:{
        'length':28.28,
        'width':10,
        'height':1.92+0.1,
        'PD':-1,
        'epitaxy':True
    },
    4:{
        'length':26.09,
        'width':10,
        'height':2.80,
        'PD':+1,
        'epitaxy':True
    },
    5:{
        'length':22.44,
        'width':10,
        'height':1.84+0.1,
        'PD':-1,
        'epitaxy':True
    },
    6:{
        'length':24.46,
        'width':10,
        'height':2.73,
        'PD':+1,
        'epitaxy':True
    },
    7:{
        'length':24.20,
        'width':10,
        'height':1.98,
        'PD':-1,
        'epitaxy':True
    },
    8:{
        'length':21.54,
        'width':10,
        'height':2.71,
        'PD':+1,
        'epitaxy':True
    },
    9:{
        'length':22.23,
        'width':10,
        'height':1.99+0.1,
        'PD':-1,
        'epitaxy':True
    },
    10:{
        'length':23.20,
        'width':10,
        'height':2.76,
        'PD':+1,
        'epitaxy':True
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
    "x_min": 0.0, "x_max": 1.5,
    "y_min": -.5, "y_max": .5,
    "z_min": 1.6, "z_max": 3.6
}

# ========== Generate CUT_VIEWS automatically ==========
N_SUBDOMAINS = 7
OVERLAPING_DOMAINS = 0.2

CUT_VIEWS = generate_segmented_cut_views(
    global_domain=GLOBAL_DOMAIN,
    axis="x",                          # Split along the x-axis
    num_subdomains=N_SUBDOMAINS,       # Create 5 subdomains
    overlap_ratio=OVERLAPING_DOMAINS,    # 20% overlap between subdomains
    rel_positions=[0.25, 0.5, 0.75]
)
