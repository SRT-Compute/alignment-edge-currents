"""Window-level results from the dense per-cutout statistics: significance,
reference comparisons, resolution, topology, tails, domain of validity."""
import json
import math

import numpy as np

import dense as DN
import instrument as I
import stats as SX

TRIMS = [0.001, 0.01, 0.05, 0.10]


def pooled_windows(blocks):
    """blocks: dict window -> list of per-cutout results (index order)."""
    W = {"holdout_611": blocks.get("holdout_611", []), "holdout_1300": blocks.get("holdout_1300", []),
         "extension_433": blocks.get("extension_433", []),
         "pooled_611": blocks.get("stage1_611", []) + blocks.get("holdout_611", []),
         "pooled_1300": blocks.get("stage1_1300", []) + blocks.get("holdout_1300", [])}
    return {k: v for k, v in W.items() if v}


def significance(W):
    primary = [k for k in ("holdout_611", "holdout_1300", "extension_433") if k in W]
    main = [k for k in ("pooled_611", "pooled_1300") if k in W]
    tests, p, keys = {}, [], []
    for k in primary:
        s = SX.sign_flip([b["base"] for b in W[k]])
        tests[k] = s
        for i in range(3):
            p.append(s["p_signflip"][i]); keys.append((k, i))
    adj = SX.holm(p)
    holm = {}
    for (k, i), pa in zip(keys, adj):
        holm.setdefault(k, [None] * 3)[i] = {"p_holm": float(pa), "significant": bool(pa < 0.01)}
    pooled = {k: SX.sign_flip([b["base"] for b in W[k]]) for k in main}
    ci99, norm = {}, {}
    for k in W:
        est, boot = SX.pooled_bootstrap([b["base"] for b in W[k]])
        ci99[k] = {"mean": est, "lo99": np.percentile(boot, 0.5, axis=0), "hi99": np.percentile(boot, 99.5, axis=0)}
        srms = math.sqrt(sum(b["all"]["ss"] for b in W[k]) / sum(b["all"]["n"] for b in W[k]))
        norm[k] = {"S_rms": srms, "mean_norm": est / srms, "se_norm": boot.std(0) / srms}
    het = {}
    hk = [k for k in ("extension_433", "pooled_611", "pooled_1300") if k in norm]
    if len(hk) >= 2:
        m = np.array([norm[k]["mean_norm"] for k in hk])
        se = np.array([norm[k]["se_norm"] for k in hk])
        wts = 1 / se ** 2
        mbar = (wts * m).sum(0) / wts.sum(0)
        Q = (wts * (m - mbar) ** 2).sum(0)
        het = {"windows": hk, "Q": Q, "df": len(hk) - 1, "I2": np.maximum(0, (Q - (len(hk) - 1)) / np.maximum(Q, 1e-12))}
    corr = {}
    for k in W:
        rs = [b["autocorr_12"] for b in W[k] if "autocorr_12" in b]
        if rs:
            rm = np.mean(rs, 0)
            z = np.where(rm < 0)[0]
            cut = z[0] if len(z) else len(rm)
            ell = float(rm[:cut].sum() * 2 - 1)
            corr[k] = {"rho": rm, "ell_grid": ell, "n_eff_per_block": (128 - 2 * DN.MARGIN) ** 3 / max(ell, 1) ** 3}
    allsig = all(c["significant"] for k in holm for c in holm[k])
    return {"primary": tests, "holm": holm, "pooled": pooled, "ci99": ci99, "normalised": norm,
            "heterogeneity": het, "autocorrelation": corr, "verdict": "SIGNIFICANT" if allsig else "NOT_ALL_SIGNIFICANT"}


def references(W):
    out = {}
    for k in W:
        real = [b["base"] for b in W[k]]
        d = {}
        for name in ("ref_gaussian", "ref_random_direction", "ref_permuted_enstrophy"):
            ref = [b[name] for b in W[k]]
            est, boot = SX.pooled_bootstrap(real, minus=ref)
            lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
            d[name] = {"reference": SX.cell(ref), "difference": est, "lo": lo, "hi": hi, "excludes_zero": (lo > 0) | (hi < 0)}
        d["verdict"] = ["TURBULENCE-SPECIFIC" if d["ref_gaussian"]["excludes_zero"][i] else "KINEMATIC" for i in range(3)]
        out[k] = d
    return out


def resolution(W):
    out, robust = {}, True
    for k in W:
        out[k] = {v: SX.cell([b[v] for b in W[k]]) for v in ("base", "FD4", "FD8", "DX2")}
        if k in ("pooled_611", "pooled_1300"):
            for v in ("FD4", "FD8", "DX2"):
                if out[k][v] is None or not out[k][v]["pass"].all():
                    robust = False
    out["verdict"] = "ROBUST" if robust else "RESOLUTION-SENSITIVE"
    return out


def topology(W):
    out = {}
    for k in W:
        tot = sum(b["base"]["n"] for b in W[k])
        tot_s = np.array([b["base"]["s"] for b in W[k]]).sum(0)
        out[k] = {}
        for q in ("Q+R-", "Q+R+", "Q-R+", "Q-R-"):
            bl = [b["qr"][q] for b in W[k]]
            n = sum(x["n"] for x in bl)
            out[k][q] = {"cell": SX.cell(bl) if n >= 500 else None, "fraction_of_points": n / tot,
                         "share_of_total": np.array([x["s"] for x in bl]).sum(0) / tot_s}
    main = [k for k in ("pooled_611", "pooled_1300") if k in W]
    out["verdict"] = ["QR-ROBUST" if all(out[k][q]["cell"] is not None and out[k][q]["cell"]["pass"][i]
                                         for k in main for q in out[k]) else "QR-PARTIAL" for i in range(3)]
    return out


def tails(W):
    out = {}
    for k in W:
        sub = np.concatenate([b["sub"] for b in W[k]]).astype(np.float64)
        d = {"mean_subsample": sub.mean(0), "median": np.median(sub, 0),
             "median_of_block_means": np.median([b["base"]["s"] / b["base"]["n"] for b in W[k]], 0), "trimmed": {}}
        for tr in TRIMS:
            tm = []
            for i in range(3):
                a = np.abs(sub[:, i])
                tm.append(sub[a <= np.quantile(a, 1 - tr), i].mean())
            d["trimmed"][str(tr)] = np.array(tm)
        verdict = []
        for i in range(3):
            flip = any(np.sign(d["trimmed"][str(tr)][i]) != np.sign(d["mean_subsample"][i]) for tr in TRIMS if tr <= 0.01)
            verdict.append("TAIL-CARRIED" if flip or np.sign(d["median"][i]) != np.sign(d["mean_subsample"][i]) else "BULK")
        d["verdict"] = verdict
        out[k] = d
    return out


def domain(W):
    out = {"gmin": {}, "filter": {}, "enstrophy": {}, "domain": {}}

    def contiguous(levels, cells, ref, i):
        ok = [c is not None and bool(c["pass"][i]) for c in cells]
        j = levels.index(ref)
        if not ok[j]:
            return []
        lo = hi = j
        while lo - 1 >= 0 and ok[lo - 1]:
            lo -= 1
        while hi + 1 < len(levels) and ok[hi + 1]:
            hi += 1
        return levels[lo:hi + 1]

    for k in W:
        out["gmin"][k] = {str(g): SX.cell([b["gmin"][str(g)] for b in W[k]]) for g in DN.G_MIN_LEVELS}
        out["filter"][k] = {str(D): (SX.cell([b["filter"][str(D)] for b in W[k]])
                                     if all(b["filter"][str(D)] is not None for b in W[k]) else None)
                            for D in DN.FILTER_LEVELS}
        out["enstrophy"][k] = [SX.cell([b["enstrophy_bins"][j] for b in W[k]]) for j in range(len(DN.ENSTROPHY_BINS) - 1)]
        out["domain"][k] = {I.EDGE_LABELS[i]: {
            "gmin": contiguous(DN.G_MIN_LEVELS, [out["gmin"][k][str(g)] for g in DN.G_MIN_LEVELS], 0.01, i),
            "filter_Delta_over_eta": contiguous(DN.FILTER_LEVELS, [out["filter"][k][str(D)] for D in DN.FILTER_LEVELS], 0, i),
            "enstrophy_bins_passing": [j for j, c in enumerate(out["enstrophy"][k]) if c is not None and c["pass"][i]]}
            for i in range(3)}
    return out


def to_json(obj):
    def default(o):
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.ndarray): return o.tolist()
        return str(o)
    return json.loads(json.dumps(obj, default=default))
