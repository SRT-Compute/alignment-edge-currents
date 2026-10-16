"""Quantities derived from the dense estimator: normalised currents, ratios to
the Gaussian reference, edge asymmetry, net drifts N_i, cycle component K, and
the alignment-induced rate of change of the stretching efficiency.

    N_i = sum_j <J[i<-j]>          (net drift into weight i)
    K   = (-<J12> - <J23> + <J13>) / 3   (cycle component of the oriented edge currents)
    dsigma/dt|_RE = sum_{i<j} (lambda_i - lambda_j) J[i<-j],  sigma = sum_i lambda_i c_i
"""
import math

import numpy as np

import dense as DN
import instrument as I

N_BOOT = 10000


def _arrays(bl, key):
    return (np.array([b[key]["s"] for b in bl]), np.array([b[key]["a"] for b in bl]),
            np.array([b[key]["n"] for b in bl], float))


def _quantities(bl, idx):
    sD, aD, nD = _arrays(bl, "base")
    sG, _, nG = _arrays(bl, "ref_gaussian")
    ss = np.array([b["all"]["ss"] for b in bl])
    na = np.array([b["all"]["n"] for b in bl], float)
    mD = sD[idx].sum(0) / nD[idx].sum()
    mG = sG[idx].sum(0) / nG[idx].sum()
    srms = math.sqrt(ss[idx].sum() / na[idx].sum())
    out = {"mean": mD, "O": np.abs(mD) / (aD[idx].sum(0) / nD[idx].sum()), "norm": mD / srms, "norm_gaussian": mG / srms,
           "ratio_to_gaussian": mD / mG, "asymmetry": np.array([mD[0] / mD[2]]), "asymmetry_gaussian": np.array([mG[0] / mG[2]])}
    for tag, m in (("", mD / srms), ("_gaussian", mG / srms)):
        oriented = np.array([-m[0], -m[2], m[1]])                       # 1->2, 2->3, 3->1
        out["K" + tag] = np.array([oriented.mean()])
        out["N" + tag] = np.array([m[0] + m[1], -m[0] + m[2], -m[1] - m[2]])
    out["K_ratio"] = out["K"] / out["K_gaussian"]
    return out


def edge_balance(W):
    """W: dict window -> list of per-cutout results. Point estimates and 95 %
    block-bootstrap intervals; windows are processed in a fixed order with one
    shared random generator."""
    rng = np.random.default_rng(20261012)
    res = {}
    for w in ("holdout_611", "holdout_1300", "pooled_611", "pooled_1300"):
        bl = W[w]
        n = len(bl)
        full = _quantities(bl, np.arange(n))
        boot = [_quantities(bl, rng.integers(0, n, n)) for _ in range(N_BOOT)]
        res[w] = {"n_cutouts": n}
        for key, v in full.items():
            arr = np.array([x[key] for x in boot])
            lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
            res[w][key] = {"est": v.tolist(), "lo": lo.tolist(), "hi": hi.tolist()}
    return res


def stretching_efficiency(files_by_window, load):
    """files_by_window: ordered dict window -> list of cutout paths (gate-passing).
    Normalised by S_rms^2 of the same cutouts."""
    rng = np.random.default_rng(20261012)
    res = {}
    for w, files in files_by_window.items():
        rows = []
        for p in files:
            u, dx, _ = load(p)
            S, wv = I.interior_points(u, dx, DN.MARGIN)
            J, g, st = I.edge_currents(S, wv)
            sel = g > DN.G_MIN
            lam = st["lam"]
            gaps = np.stack([lam[:, 0] - lam[:, 1], lam[:, 0] - lam[:, 2], lam[:, 1] - lam[:, 2]], 1)
            rate = (gaps * J).sum(1)[sel]
            sigma = (lam * st["c"]).sum(1)[sel]
            ss = (S ** 2).sum((1, 2))
            SG, wG = DN.gaussian_reference(200000, ss.mean(), st["om2"].mean(), rng)
            JG, gG, stG = I.edge_currents(SG, wG)
            selG = gG > DN.G_MIN
            lG = stG["lam"]
            gapsG = np.stack([lG[:, 0] - lG[:, 1], lG[:, 0] - lG[:, 2], lG[:, 1] - lG[:, 2]], 1)
            rateG = (gapsG * JG).sum(1)[selG]
            rows.append((rate.sum(), sel.sum(), sigma.sum(), rateG.sum(), selG.sum(), ss.sum(), len(ss),
                         (gaps[sel] * J[sel]).sum(0)))
        R = np.array([r[:7] for r in rows], float)
        E = np.array([r[7] for r in rows])

        def est(idx):
            s2 = R[idx, 5].sum() / R[idx, 6].sum()
            return np.array([R[idx, 0].sum() / R[idx, 1].sum() / s2, R[idx, 3].sum() / R[idx, 4].sum() / s2,
                             R[idx, 2].sum() / R[idx, 1].sum() / math.sqrt(s2)] + list(E[idx].sum(0) / R[idx, 1].sum() / s2))
        full = est(np.arange(len(rows)))
        boot = np.array([est(rng.integers(0, len(rows), len(rows))) for _ in range(N_BOOT)])
        lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
        names = ["rate", "rate_gaussian", "sigma", "rate_edge_12", "rate_edge_13", "rate_edge_23"]
        res[w] = {nm: [float(full[i]), float(lo[i]), float(hi[i])] for i, nm in enumerate(names)}
    return res
