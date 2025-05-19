import numpy as np

# ========== Simulation Parameters ==========

START_SIMULATION = True
RANDOM_SEED = 15
np.random.seed(RANDOM_SEED)

# Dimensions of the simulation domain
LENGTH_SIMULATION = 1.0       # x-axis (mm)
WIDTH_SIMULATION = 0.3        # y-axis (mm)
BD_INCREMENTS = 0.35          # z-step between meltpool layers

# Grain spacing and diameter
D0 = 0.015
D_seeds = 0.01

# Grain competition threshold
THRESHOLD = 1

# Add random noise to growth
NOISE_NEPER = 0

# Substrate thickness control
SUBSTRATE_BBOX_EXCESS = 1.05

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
        'length':16,
        'width':6,
        'height':1.2,
        'PD':1,
        'epitaxy':False
    },
    1:{
        'length':17,
        'width':6,
        'height':1.3,
        'PD':1,
        'epitaxy':True
    },
    2:{
        'length':18,
        'width':6,
        'height':1.4,
        'PD':1,
        'epitaxy':True
    },
}

# --------------------------------------------
START = True
try:
    assert np.all(
        np.asarray([THERMAL_HISTORY[i]['width'] > WIDTH_SIMULATION for i in THERMAL_HISTORY.keys()])
        )
except AssertionError as e:
    print("The width of the meltpool is larger than the simulation width")
    START=False

try:
    assert np.all(
            np.asarray(
                [BD_INCREMENTS < np.sqrt(1-(WIDTH_SIMULATION/THERMAL_HISTORY[i]['width'])**2)*THERMAL_HISTORY[i]['height'] for i in range(1, len(THERMAL_HISTORY))]
                )
        )
except AssertionError as e:
    print("Area on the corner aren't well remelted")
    START=False
# --------------------------------------------

# Depth limit (for top plane)
MARGIN = 0.3
Z_ULTIMATE = (
    THERMAL_HISTORY[0]["height"] +
    (len(THERMAL_HISTORY) - 1 - MARGIN) * BD_INCREMENTS
)

# ========== Visualization Domains ==========

DOMAINS = [
    {"x_min": 0.1, "x_max": 0.9, "y_min": -0.05, "y_max": 0.05, "z_min": 0, "z_max": 0.3}
]

CUT_VIEWS = [
    ("XZ", 0.5),
]
