"""Reproduce all result files and figures.

    python run_all.py --data-dir DATA --out-dir results

DATA must contain the files written by fetch_data.py. Per-cutout statistics are
cached in OUT/cache, so an interrupted run can be resumed.
"""
import argparse
import json
import os
import pickle

import numpy as np

import channel as CH
import confirmatory as CF
import datasets as D
import dense as DN
import derived as DV
import figures as FG
import gates as G
import re_check as RE
import windows as WS


def dump(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(WS.to_json(obj), f, indent=1)


def main(data_dir, out_dir):
    os.makedirs(os.path.join(out_dir, "cache"), exist_ok=True)

    # 1. quality gates
    gates = {}
    for window, w in D.WINDOWS.items():
        res = []
        for k in range(len(w["origins"])):
            path = D.cutout_path(data_dir, window, k)
            if not os.path.exists(path):
                continue
            u, dx, _ = D.load_cutout(path)
            res.append(dict(G.isotropic(u, dx), cutout=k))
        nfail = sum(not r["pass"] for r in res)
        gates[window] = {"cutouts": res, "n_failed": nfail,
                         "verdict": "GATE_FAILURE" if nfail > G.MAX_FAILED else "GATE_PASS"}
    dump(gates, os.path.join(out_dir, "gates.json"))
    passing = {window: [r["cutout"] for r in g["cutouts"] if r["pass"]]
               for window, g in gates.items() if g["verdict"] == "GATE_PASS"}

    # 2. confirmatory sign and organisation test
    conf = {}
    for stage, wins in (("stage1", ("stage1_611", "stage1_1300")), ("holdout", ("holdout_611", "holdout_1300"))):
        conf[stage] = {}
        for window in wins:
            blocks, chis = [], []
            for k in passing.get(window, []):
                u, dx, _ = D.load_cutout(D.cutout_path(data_dir, window, k))
                J, c, x = CF.cutout_currents(u, dx)
                blocks.append(J); chis.append(x)
            conf[stage][window] = CF.window_result(blocks, chis) if blocks else {"window_pass": False, "note": "gate failure"}
        conf[stage]["verdict"] = "REPLICATED" if all(conf[stage][w]["window_pass"] for w in wins) else "NOT_REPLICATED"
    dump(conf, os.path.join(out_dir, "confirmatory.json"))

    # 3. dense per-cutout statistics
    blocks, position = {}, 0
    for window in D.PROCESSING_ORDER:
        for k in passing.get(window, []):
            cache = os.path.join(out_dir, "cache", "%s_c%02d.pkl" % (window, k))
            if os.path.exists(cache):
                r = pickle.load(open(cache, "rb"))
            else:
                u, dx, attrs = D.load_cutout(D.cutout_path(data_dir, window, k))
                r = DN.analyse_cutout(u, dx, D.WINDOWS[window]["eta"], DN.block_seed(window, position))
                pickle.dump(r, open(cache, "wb"))
                print("dense", window, k, flush=True)
            blocks.setdefault(window, []).append(r)
            position += 1

    # 4. window results
    W = WS.pooled_windows(blocks)
    inv = [b["invariance"] for v in blocks.values() for b in v if "invariance" in b]
    results = {"significance": WS.significance(W), "references": WS.references(W), "resolution": WS.resolution(W),
               "topology": WS.topology(W), "tails": WS.tails(W), "domain": WS.domain(W),
               "invariance": {"max_relative_change": max(inv), "pass": bool(max(inv) < 1e-10)}}

    # 5. channel flow
    results["channel_stage1"] = CH.stage1(D.channel_path(data_dir, "stage1"))
    dump(CH.stage2(D.channel_path(data_dir, "stage2")), os.path.join(out_dir, "channel_stage2.json"))
    dump(results, os.path.join(out_dir, "dense_results.json"))

    # 6. derived quantities
    derived = DV.edge_balance(W)
    dump(derived, os.path.join(out_dir, "edge_balance.json"))
    files = {"pooled_611": [D.cutout_path(data_dir, w, k) for w in ("stage1_611", "holdout_611") for k in passing[w]],
             "pooled_1300": [D.cutout_path(data_dir, w, k) for w in ("stage1_1300", "holdout_1300") for k in passing[w]]}
    dump(DV.stretching_efficiency(files, D.load_cutout), os.path.join(out_dir, "stretching_efficiency.json"))

    # 7. restricted-Euler check and figures
    dump(RE.run(), os.path.join(out_dir, "restricted_euler_check.json"))
    results_j = json.load(open(os.path.join(out_dir, "dense_results.json"), encoding="utf-8"))
    summary = FG.make_all(blocks, results_j, derived, data_dir, os.path.join(out_dir, "figures"))
    dump(summary, os.path.join(out_dir, "figure_values.json"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args()
    main(a.data_dir, a.out_dir)
