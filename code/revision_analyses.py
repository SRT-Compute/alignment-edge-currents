"""Analyses added in revision (analysis plan 5, prereg/PREREG_5_revision_analyses.md).

A  (Reviewer 4, comment 3) spatial dependence between cutouts: geometry, distance-correlation with permutation
   test, grouping (single linkage on edge gap < d*) and greedy thinning, recomputed uncertainties.
B1 (Reviewer 4, comment 4) N_2 = -<J12> + <J23> under trimming (0.1 % and 1 %), DNS against the Gaussian reference.
B2 (Reviewer 4, comment 4) characterisation of the largest |J23| points: enstrophy, eigenvalue gaps, alignment, (Q, R).

Inputs: the per-cutout cache and result files written by run_all.py (--results) and the cutout files (--data-dir,
used only in B2). The repository modules are used unchanged. Run after run_all.py:

    cd code
    python revision_analyses.py --results ../results --data-dir ../data --out ../results
"""
import argparse
import itertools
import json
import math
import os
import pickle
import time

import numpy as np
from scipy import stats as sst

import datasets as D
import dense as DN
import derived as DV
import instrument as I
import stats as SX

RESULTS = os.path.join("..", "results")
CACHE = os.path.join(RESULTS, "cache")
DATA = os.path.join("..", "data")

D_STAR = 1.3             # primary limit for the edge gap between cutouts (code units, ~ L)
D_STAR_2 = 2.6           # secondary limit (thinning only; grouping merges everything into one group)
L_INT = {"611": 1.2983, "1300_s0": 1.2438, "1300_s3": 1.3192}   # JHTDB L1 (see analysis plan 5)
N_PERM = 10000
N_BOOT = 10000
N_BOOT_TRIM = 2000
TRIMS = [0.001, 0.01]            # primary trimming levels
TRIMS_DESC = [0.05, 0.10]        # descriptive
N_SUB_G = 50000                  # Gaussian states per block for trimming (= DNS subsample)
SEED_G = 20261103
SEED_BOOT = 20261104
SEED_PERM = 20261105
LOG_BINS = np.linspace(-30.0, 15.0, 4501)

SETS = {"611": ("stage1_611", "holdout_611"), "1300_s0": ("stage1_1300",), "1300_s3": ("holdout_1300",)}
WINDOWS = {"pooled_611": ("stage1_611", "holdout_611"), "pooled_1300": ("stage1_1300", "holdout_1300"),
           "holdout_611": ("holdout_611",), "holdout_1300": ("holdout_1300",)}
SET_OF = {"stage1_611": "611", "holdout_611": "611", "stage1_1300": "1300_s0", "holdout_1300": "1300_s3"}


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, (np.integer, int)) and not isinstance(o, bool):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def load_blocks():
    """Gate-passing cutouts in the processing order of run_all.py, with their processing position."""
    gates = json.load(open(os.path.join(RESULTS, "gates.json"), encoding="utf-8"))
    out, pos = [], 0
    for window in D.PROCESSING_ORDER:
        g = gates.get(window)
        if g is None or g["verdict"] != "GATE_PASS":
            continue
        for c in g["cutouts"]:
            if not c["pass"]:
                continue
            k = c["cutout"]
            r = pickle.load(open(os.path.join(CACHE, "%s_c%02d.pkl" % (window, k)), "rb"))
            w = D.WINDOWS[window]
            out.append({"window": window, "k": k, "pos": pos, "origin": w["origins"][k], "N": w["N"], "r": r})
            pos += 1
    return out


def check_cache(blocks):
    """Largest deviation between the cache's normalised mean currents and the published ones."""
    pub = json.load(open(os.path.join(RESULTS, "dense_results.json"), encoding="utf-8"))
    pub = pub["significance"]["normalised"]
    worst = 0.0
    for key, wins in WINDOWS.items():
        bl = [b["r"] for b in blocks if b["window"] in wins]
        est = sum(b["base"]["s"] for b in bl) / sum(b["base"]["n"] for b in bl)
        srms = math.sqrt(sum(b["all"]["ss"] for b in bl) / sum(b["all"]["n"] for b in bl))
        worst = max(worst, float(np.max(np.abs(est / srms - np.array(pub[key]["mean_norm"])))))
    return worst


# ---------------------------------------------------------------- A: geometry and dependence

def _d(a, b, N):
    return [min(abs(a[i] - b[i]), N - abs(a[i] - b[i])) for i in range(3)]


def edge_gap(a, b, N):
    """Periodic distance between the boundaries of two cutouts (code units; 0 if they touch or overlap)."""
    return math.sqrt(sum(max(0, d - D.CUTOUT) ** 2 for d in _d(a, b, N))) * 2 * math.pi / N


def centre_dist(a, b, N):
    return math.sqrt(sum(d ** 2 for d in _d(a, b, N))) * 2 * math.pi / N


def geometry(bl, L):
    n = len(bl)
    gaps = np.full((n, n), np.inf)
    cents = np.full((n, n), np.inf)
    for i, j in itertools.combinations(range(n), 2):
        gaps[i, j] = gaps[j, i] = edge_gap(bl[i]["origin"], bl[j]["origin"], bl[i]["N"])
        cents[i, j] = cents[j, i] = centre_dist(bl[i]["origin"], bl[j]["origin"], bl[i]["N"])
    iu = np.triu_indices(n, 1)
    nn = gaps.min(1)
    return {"n": n, "L": L, "cutout_side": D.CUTOUT * 2 * math.pi / bl[0]["N"],
            "edge_gap_min": float(gaps[iu].min()), "edge_gap_median": float(np.median(gaps[iu])),
            "centre_min": float(cents[iu].min()), "centre_median": float(np.median(cents[iu])),
            "nn_edge_gap_min": float(nn.min()), "nn_edge_gap_median": float(np.median(nn)),
            "nn_edge_gap_max": float(nn.max()), "n_pairs": int(len(iu[0])),
            "n_pairs_gap_below_dstar": int((gaps[iu] < D_STAR).sum()),
            "in_L": {"edge_gap_min": float(gaps[iu].min() / L), "nn_edge_gap_median": float(np.median(nn) / L),
                     "centre_median": float(np.median(cents[iu]) / L)},
            "pairs": [{"a": bl[i]["window"] + "_c%02d" % bl[i]["k"], "b": bl[j]["window"] + "_c%02d" % bl[j]["k"],
                       "edge_gap": float(gaps[i, j]), "centre": float(cents[i, j])} for i, j in zip(*iu)]}, gaps


def dependence(bl, gaps, rng):
    """Correlation between block means as a function of distance. z = standardised block mean per edge;
    rho_near = mean of z_i z_j over pairs with edge gap < D_STAR; Spearman(edge gap, z_i z_j)."""
    n = len(bl)
    m = np.array([b["r"]["base"]["s"] / b["r"]["base"]["n"] for b in bl])
    srms = np.array([math.sqrt(b["r"]["all"]["ss"] / b["r"]["all"]["n"]) for b in bl])
    iu = np.triu_indices(n, 1)
    g = gaps[iu]
    near = g < D_STAR
    rank_g = sst.rankdata(g)
    perms = [rng.permutation(n) for _ in range(N_PERM)]

    def spear(prod):
        rp = np.apply_along_axis(sst.rankdata, 0, prod)
        a = rank_g - rank_g.mean()
        b = rp - rp.mean(0)
        return (a[:, None] * b).sum(0) / np.sqrt((a ** 2).sum() * (b ** 2).sum(0))

    out = {"n_blocks": n, "n_pairs": int(len(g)), "n_near_pairs": int(near.sum())}
    for name, x in (("raw", m), ("normalised", m / srms[:, None])):
        z = (x - x.mean(0)) / x.std(0, ddof=1)

        def stat(zz):
            prod = zz[iu[0]] * zz[iu[1]]
            rn = prod[near].mean(0) if near.any() else np.full(3, np.nan)
            return rn, float(np.mean(rn)), spear(prod)

        r0, c0, s0 = stat(z)
        nr, nc, ns = [], [], []
        for p in perms:
            a, b, c = stat(z[p])
            nr.append(a); nc.append(b); ns.append(c)
        nr, nc, ns = np.array(nr), np.array(nc), np.array(ns)
        out[name] = {"rho_near": r0, "rho_near_combined": c0, "spearman_gap_vs_product": s0,
                     "p_rho_near": (1 + (nr >= r0 - 1e-12).sum(0)) / (1 + N_PERM),
                     "p_rho_near_combined": float((1 + (nc >= c0 - 1e-12).sum()) / (1 + N_PERM)),
                     "p_spearman": (1 + (ns <= s0 + 1e-12).sum(0)) / (1 + N_PERM),
                     "null_rho_near_combined_95": np.percentile(nc, [2.5, 97.5])}
    return out


def groups_single_linkage(bl, d):
    n = len(bl)
    par = list(range(n))

    def f(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    for i, j in itertools.combinations(range(n), 2):
        if edge_gap(bl[i]["origin"], bl[j]["origin"], bl[i]["N"]) < d:
            par[f(i)] = f(j)
    roots = sorted({f(i) for i in range(n)})
    return [[bl[i] for i in range(n) if f(i) == r] for r in roots]


def thin_greedy(bl, d):
    """Remove the cutout with most neighbours below d (ties: smallest neighbour distance, then highest position)."""
    keep = list(bl)
    while True:
        nb = {b["pos"]: [c for c in keep if c is not b and edge_gap(b["origin"], c["origin"], b["N"]) < d] for b in keep}
        mx = max(len(v) for v in nb.values())
        if mx == 0:
            return keep
        cand = [b for b in keep if len(nb[b["pos"]]) == mx]
        cand.sort(key=lambda b: (min(edge_gap(b["origin"], c["origin"], b["N"]) for c in nb[b["pos"]]), -b["pos"]))
        keep.remove(cand[0])


def merge(bl):
    out = {}
    for key in ("base", "ref_gaussian"):
        out[key] = {"n": sum(b["r"][key]["n"] for b in bl), "s": sum(b["r"][key]["s"] for b in bl),
                    "a": sum(b["r"][key]["a"] for b in bl)}
    out["all"] = {"n": sum(b["r"]["all"]["n"] for b in bl), "ss": sum(b["r"]["all"]["ss"] for b in bl),
                  "w2": sum(b["r"]["all"]["w2"] for b in bl)}
    return out


def units_for(blocks, wins, rule):
    """Units (super-blocks) for a window under a rule: cutout | group13 | thin13 | thin26.
    Groups and thinning are formed within each snapshot set."""
    units, labels = [], []
    for s in sorted({SET_OF[w] for w in wins}):
        bl = [b for b in blocks if b["window"] in wins and SET_OF[b["window"]] == s]
        if rule == "cutout":
            parts = [[b] for b in bl]
        elif rule == "group13":
            parts = groups_single_linkage(bl, D_STAR)
        elif rule == "thin13":
            parts = [[b] for b in thin_greedy(bl, D_STAR)]
        elif rule == "thin26":
            parts = [[b] for b in thin_greedy(bl, D_STAR_2)]
        else:
            raise ValueError(rule)
        for p in parts:
            units.append(merge(p))
            labels.append([b["window"] + "_c%02d" % b["k"] for b in p])
    return units, labels


def robustness(blocks, rng):
    res = {}
    for rule in ("cutout", "group13", "thin13", "thin26"):
        res[rule] = {}
        pvals, keys = [], []
        for wkey, wins in WINDOWS.items():
            units, labels = units_for(blocks, wins, rule)
            n = len(units)
            full = DV._quantities(units, np.arange(n))
            boot = [DV._quantities(units, rng.integers(0, n, n)) for _ in range(N_BOOT)]
            d = {"n_units": n, "units": labels}
            for key, v in full.items():
                arr = np.array([x[key] for x in boot])
                lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
                d[key] = {"est": v, "lo": lo, "hi": hi}
            if wkey.startswith("holdout"):
                sf = SX.sign_flip([u["base"] for u in units])
                d["sign_flip"] = {"p": sf["p_signflip"], "k_units": sf["k_blocks"], "min_attainable_p": 2.0 ** -n}
                for i in range(3):
                    pvals.append(sf["p_signflip"][i]); keys.append((wkey, i))
            if wkey.startswith("pooled"):
                nm = d["norm"]
                crit = {"edges_sign_ci": bool(np.all((np.sign(nm["est"]) == SX.PREDICTED_SIGNS) &
                                                     ((nm["lo"] > 0) | (nm["hi"] < 0)))),
                        "N2_positive_ci": bool(d["N"]["lo"][1] > 0),
                        "ratios_exclude_1": bool(np.all((d["ratio_to_gaussian"]["lo"] > 1) |
                                                        (d["ratio_to_gaussian"]["hi"] < 1))),
                        "asymmetry_excludes_1": bool(d["asymmetry"]["lo"][0] > 1)}
                crit["all"] = all(crit.values())
                d["criteria"] = crit
            res[rule][wkey] = d
        adj = SX.holm(pvals)
        res[rule]["holm_holdout"] = {"%s_%s" % (k, I.EDGE_LABELS[i]): float(a) for (k, i), a in zip(keys, adj)}
    return res


# ---------------------------------------------------------------- B1: trimmed N2

def n2_estimates(arrs, srms):
    """arrs: list of (n_b, 3) arrays of J (edges 1<-2, 1<-3, 2<-3). Returns N2/S_rms for
    T0 (trimming on |n2|, n2 = -J12 + J23), V1 (edge-wise trimming as in Figure 6, combined) and
    V2 (removal of the largest |J23| from both edges), for the levels [0] + TRIMS + TRIMS_DESC."""
    X = np.concatenate(arrs)
    n2 = -X[:, 0] + X[:, 2]
    a, a12, a23 = np.abs(n2), np.abs(X[:, 0]), np.abs(X[:, 2])
    lev = TRIMS + TRIMS_DESC
    qs = [1 - q for q in lev]
    t0, t12, t23 = np.quantile(a, qs), np.quantile(a12, qs), np.quantile(a23, qs)
    full = n2.mean()
    T0 = [full] + [n2[a <= t].mean() for t in t0]
    V1 = [full] + [-X[a12 <= u, 0].mean() + X[a23 <= v, 2].mean() for u, v in zip(t12, t23)]
    V2 = [full] + [n2[a23 <= v].mean() for v in t23]
    return {"T0": np.array(T0) / srms, "V1": np.array(V1) / srms, "V2": np.array(V2) / srms}


def trimmed_n2(blocks, rng):
    res = {}
    for wkey in ("pooled_611", "pooled_1300", "holdout_611", "holdout_1300"):
        bl = [b for b in blocks if b["window"] in WINDOWS[wkey]]
        dns = [b["r"]["sub"].astype(np.float64) for b in bl]
        gau = []
        for b in bl:
            rg = np.random.default_rng(SEED_G + b["pos"])
            ss = b["r"]["all"]["ss"] / b["r"]["all"]["n"]
            w2 = b["r"]["all"]["w2"] / b["r"]["all"]["n"]
            SG, wG = DN.gaussian_reference(N_SUB_G, ss, w2, rg)
            JG, gG, _ = I.edge_currents(SG, wG)
            gau.append(JG[gG > DN.G_MIN])
        ssb = np.array([b["r"]["all"]["ss"] for b in bl])
        nb = np.array([b["r"]["all"]["n"] for b in bl], float)
        srms = math.sqrt(ssb.sum() / nb.sum())
        d = {"levels": [0.0] + TRIMS + TRIMS_DESC, "n_blocks": len(bl), "S_rms": srms,
             "DNS": n2_estimates(dns, srms), "Gaussian": n2_estimates(gau, srms)}
        if wkey.startswith("pooled"):
            bd, bg = [], []
            for _ in range(N_BOOT_TRIM):
                idx = rng.integers(0, len(bl), len(bl))
                s = math.sqrt(ssb[idx].sum() / nb[idx].sum())
                bd.append(n2_estimates([dns[i] for i in idx], s))
                bg.append(n2_estimates([gau[i] for i in idx], s))
            for tag, bb in (("DNS", bd), ("Gaussian", bg)):
                for var in ("T0", "V1", "V2"):
                    arr = np.array([x[var] for x in bb])
                    lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
                    d[tag + "_ci_" + var] = {"lo": lo, "hi": hi}
            nlev = 1 + len(TRIMS)
            d["criterion_T0"] = bool(np.all(d["DNS_ci_T0"]["lo"][:nlev] > 0))
            d["gaussian_T0_covers_zero"] = bool(np.all((d["Gaussian_ci_T0"]["lo"][:nlev] <= 0) &
                                                       (d["Gaussian_ci_T0"]["hi"][:nlev] >= 0)))
        res[wkey] = d
        print("B1", wkey, "done", flush=True)
    res["verdict"] = ("N2_ROBUST_UNDER_TRIMMING" if all(res[w]["criterion_T0"] for w in ("pooled_611", "pooled_1300"))
                      else "N2_NOT_ROBUST_UNDER_TRIMMING")
    return res


# ---------------------------------------------------------------- B2: characterisation of the 2<-3 tail

def _hist(x):
    return np.histogram(np.clip(x, LOG_BINS[0], LOG_BINS[-1] - 1e-9), bins=LOG_BINS)[0]


def _median_from_hist(h):
    c = np.cumsum(h)
    if c[-1] == 0:
        return None
    i = int(np.searchsorted(c, c[-1] / 2.0))
    return float(0.5 * (LOG_BINS[i] + LOG_BINS[i + 1]))


def _new_acc():
    return {"n": 0, "sum_J23": 0.0, "n_neg": 0, "n_chi_gt1": 0, "n_r_gt4": 0, "n_g23_lt005": 0, "n_c_clipped": 0,
            "sum_ln_r": 0.0, "sum_ln_c": 0.0, "sum_ln_d": 0.0, "sum_ln_b": 0.0, "sum_ln_rot": 0.0,
            "quad_n": {q: 0 for q in ("Q+R-", "Q+R+", "Q-R+", "Q-R-")},
            "quad_J23": {q: 0.0 for q in ("Q+R-", "Q+R+", "Q-R+", "Q-R-")},
            "h_ln_r": np.zeros(len(LOG_BINS) - 1, np.int64), "h_ln_g23": np.zeros(len(LOG_BINS) - 1, np.int64),
            "h_ln_d": np.zeros(len(LOG_BINS) - 1, np.int64), "h_ln_c": np.zeros(len(LOG_BINS) - 1, np.int64),
            "h_ln_chi": np.zeros(len(LOG_BINS) - 1, np.int64)}


def _add(acc, m, J23, lnr, lnc, lnd, lnb, lnrot, chi, r, g23, clipped, quads, lng23, lnchi):
    acc["n"] += int(m.sum())
    acc["sum_J23"] += float(J23[m].sum())
    acc["n_neg"] += int((J23[m] < 0).sum())
    acc["n_chi_gt1"] += int((chi[m] > 1).sum())
    acc["n_r_gt4"] += int((r[m] > 4).sum())
    acc["n_g23_lt005"] += int((g23[m] < 0.05).sum())
    acc["n_c_clipped"] += int(clipped[m].sum())
    acc["sum_ln_r"] += float(lnr[m].sum())
    acc["sum_ln_c"] += float(lnc[m].sum())
    acc["sum_ln_d"] += float(lnd[m].sum())
    acc["sum_ln_b"] += float(lnb * m.sum())
    acc["sum_ln_rot"] += float(lnrot[m].sum())
    for q, qm in quads.items():
        acc["quad_n"][q] += int((m & qm).sum())
        acc["quad_J23"][q] += float(J23[m & qm].sum())
    acc["h_ln_r"] += _hist(lnr[m]); acc["h_ln_g23"] += _hist(lng23[m]); acc["h_ln_d"] += _hist(lnd[m])
    acc["h_ln_c"] += _hist(lnc[m]); acc["h_ln_chi"] += _hist(lnchi[m])


def _summary(acc, total_J23):
    n = max(acc["n"], 1)
    return {"n": acc["n"], "share_of_total_J23": acc["sum_J23"] / total_J23, "frac_negative": acc["n_neg"] / n,
            "frac_chi23_gt1": acc["n_chi_gt1"] / n, "frac_enstrophy_gt4": acc["n_r_gt4"] / n,
            "frac_g23_lt005": acc["n_g23_lt005"] / n, "frac_c_clipped": acc["n_c_clipped"] / n,
            "mean_ln_r": acc["sum_ln_r"] / n, "mean_ln_c2c3": acc["sum_ln_c"] / n,
            "mean_ln_gap23_over_Srms": acc["sum_ln_d"] / n, "mean_ln_block": acc["sum_ln_b"] / n,
            "mean_ln_rot": acc["sum_ln_rot"] / n,
            "median_r": math.exp(_median_from_hist(acc["h_ln_r"])),
            "median_g23": math.exp(_median_from_hist(acc["h_ln_g23"])),
            "median_gap23_over_Srms": math.exp(_median_from_hist(acc["h_ln_d"])),
            "median_c2c3": math.exp(_median_from_hist(acc["h_ln_c"])),
            "median_chi23": math.exp(_median_from_hist(acc["h_ln_chi"])),
            "quadrant_frac_points": {q: v / n for q, v in acc["quad_n"].items()},
            "quadrant_share_J23": {q: v / acc["sum_J23"] if acc["sum_J23"] != 0 else None
                                   for q, v in acc["quad_J23"].items()}}


def tail_characterisation(blocks):
    res = {}
    for wkey in ("pooled_611", "pooled_1300"):
        bl = [b for b in blocks if b["window"] in WINDOWS[wkey]]
        sub = np.concatenate([b["r"]["sub"].astype(np.float64) for b in bl])
        a23 = np.abs(sub[:, 2])
        tau = {q: float(np.quantile(a23, 1 - q)) for q in TRIMS}
        acc_all = _new_acc()
        acc_t = {q: _new_acc() for q in TRIMS}
        for b in bl:
            u, dx, _ = D.load_cutout(D.cutout_path(DATA, b["window"], b["k"]))
            S, w = I.interior_points(u, dx, DN.MARGIN)
            J, g, st = I.edge_currents(S, w)
            om2 = st["om2"]
            om2_mean = om2.mean()
            srms_b = math.sqrt((S ** 2).sum((1, 2)).mean())
            sel = g > DN.G_MIN
            J23 = J[sel, 2]
            r = om2[sel] / om2_mean
            gap23 = st["gap"][sel, 2]
            g23 = st["g"][sel, 2]
            c = st["c"][sel]
            c23 = c[:, 1] * c[:, 2]
            clipped = c23 < 1e-12
            chi = om2[sel] / (4 * gap23 ** 2)
            S_, w_ = S[sel], w[sel]
            Om = np.zeros_like(S_)
            Om[:, 0, 1], Om[:, 0, 2], Om[:, 1, 2] = -0.5 * w_[:, 2], 0.5 * w_[:, 1], -0.5 * w_[:, 0]
            Om[:, 1, 0], Om[:, 2, 0], Om[:, 2, 1] = -Om[:, 0, 1], -Om[:, 0, 2], -Om[:, 1, 2]
            Q = 0.5 * (0.5 * om2[sel] - (S_ ** 2).sum((1, 2)))
            R = -np.linalg.det(S_ + Om)
            quads = {"Q+R-": (Q > 0) & (R <= 0), "Q+R+": (Q > 0) & (R > 0),
                     "Q-R+": (Q <= 0) & (R > 0), "Q-R-": (Q <= 0) & (R <= 0)}
            del Om, S_, w_, Q, R
            lnr = np.log(r)
            lnc = np.log(np.clip(c23, 1e-12, None))
            lnd = np.log(gap23 / srms_b)
            lnb = math.log(om2_mean / srms_b)
            # size of the rotation term |omega|^2 c2 c3 / (2 Delta23), computed directly; ln = lnr + lnb + lnc - lnd - ln 2
            lnrot = np.log(om2[sel] * np.clip(c23, 1e-12, None) / (2.0 * gap23))
            lng23 = np.log(g23)
            lnchi = np.log(chi)
            everything = np.ones(len(J23), bool)
            _add(acc_all, everything, J23, lnr, lnc, lnd, lnb, lnrot, chi, r, g23, clipped, quads, lng23, lnchi)
            for q in TRIMS:
                _add(acc_t[q], np.abs(J23) > tau[q], J23, lnr, lnc, lnd, lnb, lnrot, chi, r, g23, clipped, quads,
                     lng23, lnchi)
            del S, w, J, g, st, u
            print("B2", wkey, b["window"], b["k"], flush=True)
        tot = acc_all["sum_J23"]
        d = {"tau": tau, "all": _summary(acc_all, tot)}
        for q in TRIMS:
            T = _summary(acc_t[q], tot)
            A = d["all"]
            E_om = T["mean_ln_r"] - A["mean_ln_r"]
            E_c = T["mean_ln_c2c3"] - A["mean_ln_c2c3"]
            E_gap = -(T["mean_ln_gap23_over_Srms"] - A["mean_ln_gap23_over_Srms"])
            E_b = T["mean_ln_block"] - A["mean_ln_block"]
            dl = T["mean_ln_rot"] - A["mean_ln_rot"]
            if E_om >= 2 * E_gap and E_om > 0:
                cls = "ENSTROPHY-DOMINATED"
            elif E_gap >= 2 * E_om and E_gap > 0:
                cls = "GAP-DOMINATED"
            else:
                cls = "MIXED"
            T["log_excess"] = {"E_enstrophy": E_om, "E_alignment": E_c, "E_gap": E_gap, "E_block": E_b,
                               "total": dl, "identity_residual": dl - (E_om + E_c + E_gap + E_b)}
            T["class"] = cls
            d["tail_%s" % q] = T
        res[wkey] = d
    cl = [res[w]["tail_%s" % q]["class"] for w in ("pooled_611", "pooled_1300") for q in TRIMS]
    res["verdict"] = cl[0] if len(set(cl)) == 1 else "MIXED_ACROSS_WINDOWS"
    return res


# ---------------------------------------------------------------- main

def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    blocks = load_blocks()
    meta = {"script": os.path.basename(__file__), "analysis_plan": "prereg/PREREG_5_revision_analyses.md",
            "n_blocks": len(blocks), "D_STAR": D_STAR, "D_STAR_2": D_STAR_2, "L_INT": L_INT, "N_PERM": N_PERM,
            "N_BOOT": N_BOOT, "N_BOOT_TRIM": N_BOOT_TRIM, "TRIMS": TRIMS, "TRIMS_DESC": TRIMS_DESC,
            "N_SUB_G": N_SUB_G, "seeds": {"G": SEED_G, "boot": SEED_BOOT, "perm": SEED_PERM},
            "cache_check_max_abs_diff": check_cache(blocks)}
    print("cache check", meta["cache_check_max_abs_diff"], flush=True)
    if meta["cache_check_max_abs_diff"] > 1e-9:
        raise SystemExit("The cache does not reproduce the published values; stopping.")

    rng_perm = np.random.default_rng(SEED_PERM)
    rng_boot = np.random.default_rng(SEED_BOOT)
    A = {"geometry": {}, "dependence": {}}
    pv, keys = [], []
    for s, wins in SETS.items():
        bl = [b for b in blocks if b["window"] in wins]
        geo, gaps = geometry(bl, L_INT[s])
        A["geometry"][s] = geo
        dep = dependence(bl, gaps, rng_perm)
        A["dependence"][s] = dep
        pv.append(dep["raw"]["p_rho_near_combined"]); keys.append(s)
    adj = SX.holm(pv)
    A["holm_primary"] = {k: float(a) for k, a in zip(keys, adj)}
    A["verdict_dependence"] = "DEPENDENCE_DETECTED" if any(a < 0.05 for a in adj) else "NO_DEPENDENCE_DETECTED"
    print("A dependence done", flush=True)
    A["robustness"] = robustness(blocks, rng_boot)
    ok = all(A["robustness"][r][w]["criteria"]["all"] for r in ("group13", "thin13") for w in ("pooled_611", "pooled_1300"))
    A["verdict_robustness"] = "ROBUST_TO_SPATIAL_GROUPING" if ok else "NOT_ROBUST_TO_SPATIAL_GROUPING"
    print("A robustness done", flush=True)

    B1 = trimmed_n2(blocks, rng_boot)
    B2 = tail_characterisation(blocks)
    meta["runtime_s"] = time.time() - t0
    res = {"meta": meta, "A": A, "B1": B1, "B2": B2}
    with open(os.path.join(out_dir, "revision_analyses.json"), "w", encoding="utf-8") as f:
        json.dump(jsonable(res), f, indent=1, ensure_ascii=False)
    print("done", round(meta["runtime_s"]), "s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=RESULTS, help="directory with gates.json, dense_results.json and cache/")
    ap.add_argument("--data-dir", default=DATA, help="directory with the cutout files (used in B2)")
    ap.add_argument("--out", default=RESULTS, help="output directory")
    a = ap.parse_args()
    RESULTS, CACHE, DATA = a.results, os.path.join(a.results, "cache"), a.data_dir
    main(a.out)
