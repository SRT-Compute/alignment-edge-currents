"""Dense per-cutout analysis: all interior points of one 128^3 cutout.

Each function returns block statistics {n, s = sum J, a = sum |J|}; window
means are pooled ratios over blocks (stats.py). The block seed fixes the random
draws of the reference ensembles and the subsamples:

    seed = 20261012 + position + 100 * WINDOW_CODE[window],

where `position` is the index of the cutout in the processing order of the
gate-passing cutouts (datasets.PROCESSING_ORDER, cutouts in index order).
"""
import math

import numpy as np
from scipy import ndimage

import instrument as I

SEED = 20261012
WINDOW_CODE = {"stage1_611": 182, "stage1_1300": 502, "holdout_611": 523, "holdout_1300": 755, "extension_433": 548}
MARGIN = 8
G_MIN = 0.01
G_MIN_LEVELS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.2]
FILTER_LEVELS = [0, 2, 4, 8, 16, 32, 64]                  # Delta / eta; 0 = unfiltered
ENSTROPHY_BINS = [0.0, 0.5, 1.0, 2.0, 4.0, np.inf]       # |omega|^2 / <|omega|^2>_block
N_SUBSAMPLE = 50000
N_GAUSSIAN = 200000
CORR_LAGS = 30


def block_seed(window, position):
    return SEED + position + 100 * WINDOW_CODE[window]


def bstat(J):
    return {"n": int(len(J)), "s": J.sum(0), "a": np.abs(J).sum(0)}


def gaussian_reference(n, ss, w2, rng):
    """Isotropic Gaussian gradient states: traceless symmetric S and independent
    omega, with <S:S> = ss and <|omega|^2> = w2."""
    X = rng.standard_normal((n, 3, 3))
    S = 0.5 * (X + X.transpose(0, 2, 1))
    S -= (np.trace(S, axis1=1, axis2=2) / 3)[:, None, None] * np.eye(3)
    S *= math.sqrt(ss / (S ** 2).sum((1, 2)).mean())
    w = rng.standard_normal((n, 3)) * math.sqrt(w2 / 3)
    return S, w


def invariance_check(S, w, rng):
    """Largest change of J under a random proper rotation and a reflection,
    relative to <|J|>, on up to 10 000 points with g > G_MIN."""
    J0, g0, _ = I.edge_currents(S[:40000], w[:40000])
    keep = np.where(g0 > G_MIN)[0][:10000]
    S, w, J0 = S[keep], w[keep], J0[keep]
    worst = 0.0
    for det in (+1, -1):
        Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
        if np.linalg.det(Q) * det < 0:
            Q[:, 0] *= -1
        S2 = np.einsum("ik,nkl,jl->nij", Q, S, Q)
        w2 = np.linalg.det(Q) * (w @ Q.T)                  # pseudovector
        J1, _, _ = I.edge_currents(S2, w2)
        worst = max(worst, float((np.abs(J1 - J0).max(0) / np.abs(J0).mean(0)).max()))
    return worst


def analyse_cutout(u, dx, eta, seed):
    rng = np.random.default_rng(seed)
    out = {}
    S, w = I.interior_points(u, dx, MARGIN)
    J, g, st = I.edge_currents(S, w)
    om2 = st["om2"]
    out["all"] = {"n": len(S), "ss": float((S ** 2).sum()), "w2": float(om2.sum())}
    if (seed - SEED) % 7 == 0:
        out["invariance"] = invariance_check(S, w, rng)

    out["gmin"] = {str(gm): bstat(J[g > gm]) for gm in G_MIN_LEVELS}
    sel = g > G_MIN
    Jb = J[sel]
    out["base"] = bstat(Jb)
    idx = rng.choice(len(Jb), size=min(N_SUBSAMPLE, len(Jb)), replace=False)
    out["sub"] = Jb[idx].astype(np.float32)

    r = om2[sel] / om2.mean()
    out["enstrophy_bins"] = [bstat(Jb[(r >= ENSTROPHY_BINS[i]) & (r < ENSTROPHY_BINS[i + 1])])
                             for i in range(len(ENSTROPHY_BINS) - 1)]

    S_, w_ = S[sel], w[sel]
    Om = np.zeros_like(S_)
    Om[:, 0, 1], Om[:, 0, 2], Om[:, 1, 2] = -0.5 * w_[:, 2], 0.5 * w_[:, 1], -0.5 * w_[:, 0]
    Om[:, 1, 0], Om[:, 2, 0], Om[:, 2, 1] = -Om[:, 0, 1], -Om[:, 0, 2], -Om[:, 1, 2]
    Q = 0.5 * (0.5 * om2[sel] - (S_ ** 2).sum((1, 2)))
    R = -np.linalg.det(S_ + Om)
    quadrants = {"Q+R-": (Q > 0) & (R <= 0), "Q+R+": (Q > 0) & (R > 0),
                 "Q-R+": (Q <= 0) & (R > 0), "Q-R-": (Q <= 0) & (R <= 0)}
    out["qr"] = {k: bstat(Jb[m]) for k, m in quadrants.items()}
    del Om, Q, R, S_, w_

    if (seed - SEED) % 5 == 0:
        L = 128 - 2 * MARGIN
        F = np.where(sel, J[:, 0], 0.0).reshape(L, L, L)
        F = F - F.mean()
        var = (F ** 2).mean()
        out["autocorr_12"] = [1.0] + [float((F[:-k] * F[k:]).mean() / var) for k in range(1, CORR_LAGS + 1)]
        del F

    SG, wG = gaussian_reference(N_GAUSSIAN, (S ** 2).sum((1, 2)).mean(), om2.mean(), rng)
    JG, gG, _ = I.edge_currents(SG, wG)
    out["ref_gaussian"] = bstat(JG[gG > G_MIN])
    del SG, wG, JG

    v = rng.standard_normal((len(S), 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    Jx, gx, _ = I.edge_currents(S, v * np.sqrt(om2)[:, None])
    out["ref_random_direction"] = bstat(Jx[gx > G_MIN])
    perm = rng.permutation(len(S))
    Jo, go, _ = I.edge_currents(S, st["xi"] * np.sqrt(om2[perm])[:, None])
    out["ref_permuted_enstrophy"] = bstat(Jo[go > G_MIN])
    del Jx, Jo, v, st, J, S, w

    for order in (4, 8):
        S2, w2 = I.interior_points(u, dx, MARGIN, order)
        J2, g2, _ = I.edge_currents(S2, w2)
        out["FD%d" % order] = bstat(J2[g2 > G_MIN])
    S2, w2 = I.interior_points(np.ascontiguousarray(u[:, ::2, ::2, ::2]), 2 * dx, 4)
    J2, g2, _ = I.edge_currents(S2, w2)
    out["DX2"] = bstat(J2[g2 > G_MIN])
    del S2, w2, J2

    out["filter"] = {}
    for D in FILTER_LEVELS:
        if D == 0:
            out["filter"]["0"] = out["base"]
            continue
        sigma = (D * eta / math.sqrt(12)) / dx
        m = int(math.ceil(3 * sigma)) + 4
        if 128 - 2 * m < 16:
            out["filter"][str(D)] = None
            continue
        uf = np.stack([ndimage.gaussian_filter(u[i], sigma, mode="nearest", truncate=3.0) for i in range(3)])
        S2, w2 = I.interior_points(uf, dx, m)
        J2, g2, _ = I.edge_currents(S2, w2)
        out["filter"][str(D)] = bstat(J2[g2 > G_MIN])
    return out
