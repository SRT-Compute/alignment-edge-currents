"""Tables 6 and 7 of the revised manuscript (Appendix B), from the result files written by run_all.py.

Table 6: normalised mean currents <J>/S_rms of the pooled windows under the numerical variations (order of the central
differences, and differentiation on every second grid point), with the relative change of each mean with respect to
the sixth-order reference.
Table 7: relative divergence <|div u|>/<|grad u|> of the cutouts (gate G2) per window: median, range and number of
failed cutouts.

    cd code
    python appendix_b_tables.py --results ../results
"""
import argparse
import json
import math
import os
import statistics

import datasets as D

VARIANTS = (("base", "sixth order (reference)"), ("FD4", "fourth order"), ("FD8", "eighth order"),
            ("DX2", "every second grid point"))
POOLED = (("pooled_611", "stage1_611"), ("pooled_1300", "stage1_1300"))
WINDOWS7 = ("stage1_611", "holdout_611", "stage1_1300", "holdout_1300", "extension_433")


def dx_over_eta(window):
    w = D.WINDOWS[window]
    return 2 * math.pi / w["N"] / w["eta"]


def pct(x):
    r = round(x, 1)
    return "0.0" if r == 0 else "%+.1f" % r


def main(results):
    d = json.load(open(os.path.join(results, "dense_results.json"), encoding="utf-8"))
    g = json.load(open(os.path.join(results, "gates.json"), encoding="utf-8"))
    print("Table 6: <J>/S_rms for the edges 1<-2, 1<-3, 2<-3; relative change to the reference in parentheses")
    for w, src in POOLED:
        srms = d["significance"]["normalised"][w]["S_rms"]
        base = d["resolution"][w]["base"]["mean"]
        for v, name in VARIANTS:
            m = d["resolution"][w][v]["mean"]
            dxe = dx_over_eta(src) * (2 if v == "DX2" else 1)
            cells = ["%.4f" % (m[i] / srms) + ("" if v == "base" else " (%s %%)" % pct(100 * (m[i] / base[i] - 1)))
                     for i in range(3)]
            print("  %-11s %-26s dx/eta %.2f   %s" % (w, name, dxe, "   ".join(cells)))
    print("Table 7: relative divergence of the cutouts (gate limit 0.01)")
    for k in WINDOWS7:
        x = sorted(c["G2"] for c in g[k]["cutouts"])
        print("  %-14s dx/eta %.2f   median %.4f   range %.4f-%.4f   failed %d/%d"
              % (k, dx_over_eta(k), statistics.median(x), x[0], x[-1], g[k]["n_failed"], len(x)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join("..", "results"))
    main(ap.parse_args().results)
