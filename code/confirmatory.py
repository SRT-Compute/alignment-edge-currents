"""Confirmatory sign and organisation test (2000 preregistered points per cutout).

For each window: mean edge currents over the sampling points of all gate-passing
cutouts, organisation O = |<J>| / <|J>|, and the verdict. A stage is REPLICATED
if all three predicted signs and all three O >= 0.05 hold in both of its windows.
"""
import numpy as np

import instrument as I

PREDICTED_SIGNS = np.array([-1, +1, -1])
O_MIN = 0.05
G_MIN = 0.01


def sampling_points():
    return np.random.default_rng(20260811).integers(32, 128 - 32, size=(4000, 3))[:2000]


def cutout_currents(u, dx):
    pts = sampling_points()
    A = I.gradient(u, dx)
    S, w = I.strain_vorticity(A)
    Sp = S[:, :, pts[:, 0], pts[:, 1], pts[:, 2]].transpose(2, 0, 1)
    wv = w[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
    J, g, st = I.edge_currents(Sp, wv)
    keep = g > G_MIN
    return J[keep], st["c"][keep], I.chi(st)[keep]


def window_result(blocks, chis):
    flat = np.concatenate(blocks)
    mean = flat.mean(0)
    O = np.abs(mean) / np.abs(flat).mean(0)
    rng = np.random.default_rng(11)
    n = len(blocks)
    boot = np.empty((2000, 3))
    for b in range(2000):
        idx = rng.integers(0, n, n)
        boot[b] = np.concatenate([blocks[i] for i in idx]).sum(0) / sum(len(blocks[i]) for i in idx)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    checks = [{"edge": I.EDGE_LABELS[i], "sign_ok": bool(np.sign(mean[i]) == PREDICTED_SIGNS[i]),
               "O_ok": bool(O[i] >= O_MIN)} for i in range(3)]
    return {"n_cutouts": n, "n_points": int(len(flat)), "mean": mean.tolist(), "O": O.tolist(),
            "median": np.median(flat, 0).tolist(), "ci95_lo": lo.tolist(), "ci95_hi": hi.tolist(),
            "chi_median": np.median(np.concatenate(chis), 0).tolist(), "checks": checks,
            "window_pass": all(c["sign_ok"] and c["O_ok"] for c in checks)}
