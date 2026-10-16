"""Frozen sample definitions (cutout origins, snapshots, channel blocks) and
file naming. Values are those fixed in the analysis plans (see prereg/)."""
import os
import numpy as np
import h5py

CUTOUT = 128

# Isotropic windows: JHTDB dataset, Reynolds number, Kolmogorov length, snapshot
# index passed to getCutout (1-based), and cutout origins (x, y, z; 1-based).
WINDOWS = {
    "stage1_611": dict(dataset="isotropic4096", Re=611, N=4096, eta=1.3844e-3, nu=1.732e-4, snapshot=[1] * 12, origins=[
        (3776, 3851, 3102), (3074, 1355, 3727), (3786, 525, 479), (663, 1010, 1237), (622, 3449, 3340), (3236, 3936, 1127),
        (2389, 3127, 46), (3329, 1716, 217), (592, 3440, 1150), (318, 2403, 545), (1734, 2478, 105), (2171, 3869, 3121)]),
    "stage1_1300": dict(dataset="isotropic8192", Re=1300, N=8192, eta=5.00e-4, nu=4.385e-5, snapshot=[1] * 12, origins=[
        (2322, 3143, 867), (4932, 4272, 6794), (6774, 466, 7833), (5017, 2465, 6944), (355, 4403, 3939), (4432, 6810, 3255),
        (2465, 4028, 4256), (7036, 3870, 1643), (2080, 1810, 6660), (6587, 7647, 5460), (3893, 3784, 1030), (5106, 3254, 7070)]),
    "holdout_611": dict(dataset="isotropic4096", Re=611, N=4096, eta=1.3844e-3, nu=1.732e-4, snapshot=[1] * 12, origins=[
        (3275, 3767, 1618), (1764, 3092, 2503), (473, 2180, 1145), (2361, 3506, 2849), (2293, 3865, 2074), (954, 30, 714),
        (3771, 2531, 837), (2519, 43, 2483), (1305, 1707, 262), (1646, 319, 663), (1251, 3026, 2596), (1860, 1710, 1899)]),
    "holdout_1300": dict(dataset="isotropic8192", Re=1300, N=8192, eta=5.00e-4, nu=4.385e-5, snapshot=[4] * 12, origins=[
        (4449, 4284, 2942), (2478, 3851, 1593), (4982, 674, 2432), (3613, 1944, 1588), (6615, 4099, 1150), (2958, 767, 6665),
        (6447, 6200, 6973), (4622, 3833, 4708), (861, 4707, 6583), (6176, 7580, 5400), (525, 5502, 5098), (4543, 357, 5963)]),
    "extension_433": dict(dataset="isotropic1024coarse", Re=433, N=1024, eta=2.87e-3, nu=1.85e-4,
                          snapshot=[1] * 4 + [2001] * 4 + [4001] * 4, origins=[
        (698, 679, 259), (64, 683, 391), (859, 439, 35), (246, 343, 324), (657, 813, 293), (349, 585, 27),
        (780, 213, 230), (881, 504, 406), (90, 536, 52), (457, 127, 614), (603, 161, 870), (864, 707, 591)]),
}

# Order in which gate-passing cutouts are processed; it fixes the random seeds
# of the reference ensembles (see dense.py).
PROCESSING_ORDER = ["stage1_611", "stage1_1300", "holdout_611", "holdout_1300", "extension_433"]

# Channel flow (JHTDB "channel", Re_tau ~ 1000): wall-normal bands and blocks
# (time, x centre, z centre). Stage 1 and stage 2 use different blocks.
RE_TAU = 1000.0
YPLUS = [15, 30, 100, 300, 1000]
CHANNEL_STAGES = {
    "stage1": dict(seed=20261012, blocks=[
        (1.0, 7.0611, 3.853), (3.0, 14.6972, 7.0311), (5.0, 3.4486, 2.5809), (7.0, 1.2899, 3.0021),
        (9.0, 13.9097, 5.8508), (11.0, 23.0263, 1.5872), (13.0, 6.4042, 6.7532), (15.0, 8.0422, 4.3598),
        (17.0, 17.5063, 3.1113), (19.0, 13.3004, 2.8835), (21.0, 14.0915, 9.0557), (23.0, 14.6861, 0.8295)]),
    "stage2": dict(seed=20261013, blocks=None),   # generated from the seed, see channel_blocks()
}


def channel_blocks(stage):
    cfg = CHANNEL_STAGES[stage]
    if cfg["blocks"] is not None:
        return cfg["blocks"]
    rng = np.random.default_rng(cfg["seed"] + 1000)
    return [(2.0 + 2.0 * k, float(rng.uniform(0.3, 8 * np.pi - 0.3)), float(rng.uniform(0.2, 3 * np.pi - 0.2)))
            for k in range(12)]


def channel_points(stage, band, block):
    """Time and 2000 sample points (x, y, z) for one channel block."""
    cfg = CHANNEL_STAGES[stage]
    yp = YPLUS[band]
    rng = np.random.default_rng(cfg["seed"] + 2000 + 100 * band + block)
    t, xc, zc = channel_blocks(stage)[block]
    y = rng.uniform(-0.1, 0.1, 2000) if yp >= 1000 else -1 + rng.uniform(0.9 * yp, 1.1 * yp, 2000) / RE_TAU
    x = xc + rng.uniform(-0.2, 0.2, 2000)
    z = zc + rng.uniform(-0.1, 0.1, 2000)
    return t, np.stack([x, y, z], 1)


def cutout_path(data_dir, window, k):
    return os.path.join(data_dir, "%s_c%02d.h5" % (window, k))


def channel_path(data_dir, stage):
    return os.path.join(data_dir, "channel_%s.h5" % stage)


def load_cutout(path):
    """Velocity (3, n, n, n) with component i along array axis i, and grid spacing.
    Files store u, v, w (x, y, z components) on arrays ordered (z, y, x)."""
    with h5py.File(path, "r") as f:
        u = np.stack([f["u"][:], f["v"][:], f["w"][:]]).astype(np.float64)
        attrs = dict(f.attrs)
    return u[[2, 1, 0]], float(attrs["dx"] if "dx" in attrs else attrs["dx_phys"]), attrs
