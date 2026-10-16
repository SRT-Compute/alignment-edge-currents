"""Channel-flow analysis.

Stage 1: signs and organisation of the edge currents per wall-normal band, for
the full gradient and with the band-mean gradient removed.
Stage 2 (new blocks): currents normalised by the local S_rms, compared with the
isotropic reference values r = (-0.1742, +0.0474, -0.0154):
  K1  for y+ = 100, 300, 1000 and edges 1<-2, 1<-3 the 95 % interval overlaps [0.75 r, 1.25 r];
  K2  at y+ = 15 the upper bound of |J|/S_rms is below 0.5 |r| for edges 1<-2 and 1<-3.
"""
import math

import h5py
import numpy as np

import datasets as D
import gates as G
import instrument as I
import stats as SX

ISOTROPIC_REFERENCE = np.array([-0.1742, 0.0474, -0.0154])
G_MIN = 0.01


def _block(A):
    S = 0.5 * (A + A.transpose(0, 2, 1))
    w = np.stack([A[:, 2, 1] - A[:, 1, 2], A[:, 0, 2] - A[:, 2, 0], A[:, 1, 0] - A[:, 0, 1]], 1)
    return S, w


def stage1(path):
    out = {}
    with h5py.File(path, "r") as f:
        for b, yp in enumerate(D.YPLUS):
            gates, kept = [], []
            for k in range(12):
                A = f["b%d_k%02d" % (b, k)][:].reshape(-1, 3, 3)
                g = G.channel(A)
                gates.append(dict(g, block=k))
                if g["pass"]:
                    kept.append(A)
            nfail = sum(not g["pass"] for g in gates)
            res = {"yplus": yp, "gates": gates, "n_excluded": nfail}
            if nfail > G.MAX_FAILED or not kept:
                res["verdict"] = "GATE_FAILURE"
                out[str(yp)] = res
                continue
            A_mean = np.concatenate(kept).mean(0)
            for variant in ("full", "fluctuating"):
                blocks = []
                for A in kept:
                    S, w = _block(A - A_mean if variant == "fluctuating" else A)
                    J, g, _ = I.edge_currents(S, w)
                    blocks.append({"n": int((g > G_MIN).sum()), "s": J[g > G_MIN].sum(0), "a": np.abs(J[g > G_MIN]).sum(0)})
                res[variant] = SX.cell(blocks)
            S_all, _ = _block(np.concatenate(kept))
            res["S_rms_local"] = math.sqrt((S_all ** 2).sum((1, 2)).mean())
            res["mean_gradient"] = A_mean
            res["verdict"] = "HOLDS" if res["full"] is not None and res["full"]["pass"].all() else "FAILS"
            out[str(yp)] = res
    return out


def stage2(path):
    rng = np.random.default_rng(7)
    bands = {}
    with h5py.File(path, "r") as f:
        for b, yp in enumerate(D.YPLUS):
            rows, gvals = [], []
            for k in range(12):
                A = f["b%d_k%02d" % (b, k)][:].reshape(-1, 3, 3)
                g = G.channel(A)
                gvals.append(g["G2"])
                if not g["pass"]:
                    continue
                S, w = _block(A)
                J, gg, _ = I.edge_currents(S, w)
                sel = gg > G_MIN
                rows.append((J[sel].sum(0), sel.sum(), (S ** 2).sum(), len(S)))
            if sum(v > G.G2_MAX_CHANNEL for v in gvals) > G.MAX_FAILED:
                bands[yp] = {"verdict": "GATE_FAILURE", "G2": gvals}
                continue
            s = np.array([r[0] for r in rows]); n = np.array([r[1] for r in rows], float)
            ss = np.array([r[2] for r in rows]); na = np.array([r[3] for r in rows], float)
            est = (s.sum(0) / n.sum()) / math.sqrt(ss.sum() / na.sum())
            idx = rng.integers(0, len(rows), (10000, len(rows)))
            boot = (s[idx].sum(1) / n[idx].sum(1)[:, None]) / np.sqrt(ss[idx].sum(1) / na[idx].sum(1))[:, None]
            lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
            bands[yp] = {"J_norm": est, "lo": lo, "hi": hi, "n_blocks": len(rows), "G2": gvals,
                         "ratio_to_isotropic": est / ISOTROPIC_REFERENCE}
    r = ISOTROPIC_REFERENCE
    K1 = []
    for yp in (100, 300, 1000):
        for i in (0, 1):
            c = bands.get(yp, {})
            if "J_norm" not in c:
                K1.append(False); continue
            a, b = sorted([0.75 * r[i], 1.25 * r[i]])
            K1.append(bool(c["hi"][i] >= a and c["lo"][i] <= b))
    c15 = bands.get(15, {})
    K2 = [bool("J_norm" in c15 and max(abs(c15["lo"][i]), abs(c15["hi"][i])) < 0.5 * abs(r[i])) for i in (0, 1)]
    verdict = "CONFIRMED" if all(K1) and all(K2) else ("PARTIAL" if all(K1) or all(K2) else "NOT_CONFIRMED")
    return {"reference": r, "bands": bands, "K1": K1, "K2": K2, "verdict": verdict}
