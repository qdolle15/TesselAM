import numpy as np

# ========== Simulation Parameters ==========

START_SIMULATION = True
RANDOM_SEED = 108
np.random.seed(RANDOM_SEED)

# Dimensions of the simulation domain
LENGTH_SIMULATION = 1.5       # x-axis (mm)
WIDTH_SIMULATION = 0.3        # y-axis (mm)
BD_INCREMENTS = 0.15          # z-step between meltpool layers

# Grain spacing and diameter
D0 = 0.01
D_seeds = 0.1

# Grain competition threshold
THRESHOLD = 1

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
        'length':12,
        'width':0.35,
        'height':0.2,
        'PD':+1,
        'epitaxy':False
    },
    1:{
        'length':14,
        'width':0.45,
        'height':0.22,
        'PD':-1,
        'epitaxy':True
    },
    2:{
        'length':16,
        'width':0.55,
        'height':0.25,
        'PD':+1,
        'epitaxy':True
    },
    3:{
        'length':18,
        'width':0.65,
        'height':0.29,
        'PD':-1,
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
CUT_VIEWS = [
    {
        "domain" : {"x_min": 0.4, "x_max": 1.1, "y_min": -0.02, "y_max": 0.02, "z_min": 0, "z_max": 0.3},
        "plans" : [
            ("XZ", 0.25), 
            ("XZ", 0.50),
            ("XZ", 0.75),
        ]
    },
    {
        "domain" : {"x_min": 0.7, "x_max": 0.8, "y_min": -0.15, "y_max": 0.15, "z_min": 0, "z_max": 0.3},
        "plans" : [
            ("YZ", 0.25), 
            ("YZ", 0.50),
            ("YZ", 0.75),
        ]
    },  
]
