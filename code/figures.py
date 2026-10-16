"""Figures 1-8 of the manuscript and supplementary figure S1, from the result
files and the per-cutout statistics."""
import glob
import json
import math
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import datasets as D
import instrument as I

LAB = [r"$1\leftarrow2$", r"$1\leftarrow3$", r"$2\leftarrow3$"]
WN = {"pooled_611": r"$Re_\lambda=611$", "pooled_1300": r"$Re_\lambda=1300$"}
COL = {"pooled_611": "#1f5f9e", "pooled_1300": "#c0392b"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def make_all(blocks, results, derived, data_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(5)
    VW = {"pooled_611": blocks["stage1_611"] + blocks["holdout_611"],
          "pooled_1300": blocks["stage1_1300"] + blocks["holdout_1300"]}
    SR = {k: results["significance"]["normalised"][k]["S_rms"] for k in VW}
    summary = {}

    def pooled(bl, key):
        return np.array([b[key]["s"] for b in bl]), np.array([b[key]["n"] for b in bl], float)

    def boot(bl, keys, fn, B=10000):
        arrs = [pooled(bl, k) for k in keys]
        idx = rng.integers(0, len(bl), (B, len(bl)))
        est = fn(*[a[0].sum(0) / a[1].sum() for a in arrs])
        bs = np.array([fn(*[a[0][i].sum(0) / a[1][i].sum() for a in arrs]) for i in idx])
        return est, np.percentile(bs, [2.5, 97.5], axis=0)

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, name + ".png"), dpi=300)
        fig.savefig(os.path.join(out_dir, name + ".pdf"))
        plt.close(fig)

    # Figure 2: DNS and the three reference ensembles
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.2))
    names = [("base", "DNS"), ("ref_gaussian", "Gaussian"), ("ref_random_direction", "random direction"),
             ("ref_permuted_enstrophy", r"shuffled $|\omega|^2$")]
    for i in range(3):
        for wi, k in enumerate(VW):
            for j, (key, lab) in enumerate(names):
                est, ci = boot(VW[k], [key], lambda m: m / SR[k], B=2000)
                x = wi * 5 + j
                ax[i].bar(x, est[i], color=["#333333", "#9aa5b1", "#d4a017", "#5b9e5b"][j], label=lab if (i == 0 and wi == 0) else None)
                ax[i].errorbar(x, est[i], yerr=[[est[i] - ci[0][i]], [ci[1][i] - est[i]]], color="k", lw=0.8, capsize=2)
        ax[i].axhline(0, color="k", lw=0.6)
        ax[i].set_xticks([1.5, 6.5]); ax[i].set_xticklabels([WN[k] for k in VW])
        ax[i].set_title(LAB[i]); ax[i].set_ylabel(r"$\langle J\rangle / S_{\rm rms}$")
    ax[0].legend(fontsize=7, frameon=False)
    save(fig, "figure2_references")

    # Figure 3: ratio DNS / Gaussian reference
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    summary["ratio_to_gaussian"] = {}
    for wi, k in enumerate(VW):
        est, ci = boot(VW[k], ["base", "ref_gaussian"], lambda a, g: a / g)
        summary["ratio_to_gaussian"][k] = {"est": est.tolist(), "ci": ci.tolist()}
        ax.errorbar(np.arange(3) + 0.12 * (wi - 0.5), est, yerr=[est - ci[0], ci[1] - est], fmt="o", color=COL[k], label=WN[k], capsize=3)
    ax.axhline(1, color="k", lw=0.6, ls="--"); ax.set_yscale("log")
    ax.set_xticks(range(3)); ax.set_xticklabels(LAB)
    ax.set_ylabel(r"$\langle J\rangle_{\rm DNS} / \langle J\rangle_{\rm Gaussian}$"); ax.legend(frameon=False, fontsize=8)
    save(fig, "figure3_ratio_to_gaussian")

    # Figure 1: distributions of chi (first four stage-1 cutouts per Reynolds number, every second interior point)
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    bins = np.logspace(-4, 4, 81)
    for k, window in (("pooled_611", "stage1_611"), ("pooled_1300", "stage1_1300")):
        chis, chisG = [], []
        for c in range(4):
            u, dx, _ = D.load_cutout(D.cutout_path(data_dir, window, c))
            A = I.gradient(u, dx)[:, :, 8:-8:2, 8:-8:2, 8:-8:2].reshape(3, 3, -1)
            S, w = I.strain_vorticity(A)
            Sp, wv = S.transpose(2, 0, 1), w.T
            st = I.alignment_state(Sp, wv)
            ok = (st["g"] > 0.01).all(1)
            chis.append(I.chi(st)[ok])
            X = rng.standard_normal((len(Sp), 3, 3))
            SG = 0.5 * (X + X.transpose(0, 2, 1))
            SG -= (np.trace(SG, axis1=1, axis2=2) / 3)[:, None, None] * np.eye(3)
            SG *= math.sqrt((Sp ** 2).sum((1, 2)).mean() / (SG ** 2).sum((1, 2)).mean())
            WG = rng.standard_normal((len(Sp), 3)) * math.sqrt(st["om2"].mean() / 3)
            stG = I.alignment_state(SG, WG)
            okG = (stG["g"] > 0.01).all(1)
            chisG.append(I.chi(stG)[okG])
        C, CG = np.concatenate(chis), np.concatenate(chisG)
        summary.setdefault("chi", {})[k] = {"median_DNS": np.median(C, 0).tolist(), "median_gaussian": np.median(CG, 0).tolist(),
                                            "fraction_above_1_DNS": (C > 1).mean(0).tolist(),
                                            "fraction_above_1_gaussian": (CG > 1).mean(0).tolist()}
        for i in range(3):
            h, _ = np.histogram(C[:, i], bins, density=True)
            hg, _ = np.histogram(CG[:, i], bins, density=True)
            xc = np.sqrt(bins[1:] * bins[:-1])
            ax[i].plot(xc, h * xc, color=COL[k], label="DNS " + WN[k])
            if k == "pooled_611":
                ax[i].plot(xc, hg * xc, color="#9aa5b1", ls="--", label="Gaussian")
    for i in range(3):
        ax[i].axvline(1, color="k", lw=0.6); ax[i].set_xscale("log"); ax[i].set_title(LAB[i])
        ax[i].set_xlabel(r"$\chi_{ij}=|\omega|^2/4(\lambda_i-\lambda_j)^2$"); ax[i].set_ylabel(r"$\chi\,p(\chi)$")
    ax[0].legend(fontsize=7, frameon=False)
    save(fig, "figure1_chi_distributions")

    # Figure 5: (Q, R) quadrants
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    qs = ["Q+R-", "Q+R+", "Q-R+", "Q-R-"]
    for i in range(3):
        for wi, k in enumerate(VW):
            sh = [results["topology"][k][q]["share_of_total"][i] for q in qs]
            ax[i].bar(np.arange(4) + 0.38 * (wi - 0.5), sh, 0.38, color=COL[k], label=WN[k])
        ax[i].axhline(0, color="k", lw=0.6); ax[i].set_xticks(range(4))
        ax[i].set_xticklabels(["Q>0,R<0", "Q>0,R>0", "Q<0,R>0", "Q<0,R<0"], fontsize=7)
        ax[i].set_title(LAB[i]); ax[i].set_ylabel("share of total mean")
    ax[0].legend(fontsize=7, frameon=False)
    save(fig, "figure5_QR_quadrants")

    # Figure 7: filter width
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    for k in VW:
        d = results["domain"]["filter"][k]
        Ds = [int(x) for x in d if d[x] is not None]
        for i in range(3):
            m = np.array([d[str(x)]["mean"][i] for x in Ds]) / SR[k]
            lo = np.array([d[str(x)]["lo"][i] for x in Ds]) / SR[k]
            hi = np.array([d[str(x)]["hi"][i] for x in Ds]) / SR[k]
            xs = [max(x, 1) for x in Ds]
            ax[i].plot(xs, m, "o-", color=COL[k], label=WN[k], ms=3)
            ax[i].fill_between(xs, lo, hi, color=COL[k], alpha=0.2)
    for i in range(3):
        ax[i].axhline(0, color="k", lw=0.6); ax[i].set_xscale("log"); ax[i].set_title(LAB[i])
        ax[i].set_xlabel(r"filter width $\Delta/\eta$ (1 = unfiltered)")
        ax[i].set_ylabel(r"$\langle J\rangle/S_{\rm rms}$ (unfiltered $S_{\rm rms}$)")
    ax[0].legend(fontsize=7, frameon=False)
    save(fig, "figure7_filter_width")

    # Figure 8: channel flow, normalised by the local S_rms (stage 1)
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    summary["channel"] = {}
    with h5py.File(D.channel_path(data_dir, "stage1"), "r") as ch:
        for b, yp in enumerate(D.YPLUS):
            A = np.concatenate([ch["b%d_k%02d" % (b, kk)][:].reshape(-1, 3, 3) for kk in range(12)])
            S = 0.5 * (A + A.transpose(0, 2, 1))
            srms = math.sqrt((S ** 2).sum((1, 2)).mean())
            c = results["channel_stage1"][str(yp)]["full"]
            summary["channel"][yp] = {"S_rms": srms, "J_norm": (np.array(c["mean"]) / srms).tolist()}
            for i in range(3):
                ax[i].errorbar(yp, c["mean"][i] / srms, yerr=[[(c["mean"][i] - c["lo"][i]) / srms], [(c["hi"][i] - c["mean"][i]) / srms]],
                               fmt="o", color="k", capsize=2)
    for i in range(3):
        for k in VW:
            n = results["significance"]["normalised"][k]
            ax[i].axhspan(n["mean_norm"][i] - 2 * n["se_norm"][i], n["mean_norm"][i] + 2 * n["se_norm"][i],
                          color=COL[k], alpha=0.15, label="isotropic " + WN[k])
        ax[i].axhline(0, color="k", lw=0.6); ax[i].set_xscale("log"); ax[i].set_xlabel(r"$y^+$"); ax[i].set_title(LAB[i])
        ax[i].set_ylabel(r"$\langle J\rangle/S_{\rm rms}$ (local)")
    ax[0].legend(fontsize=7, frameon=False)
    save(fig, "figure8_channel")

    # Figure 6: tails
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for k in VW:
        d = results["tails"][k]
        tr = [0] + [float(t) for t in d["trimmed"]]
        for i, mk in zip(range(3), ("o", "s", "^")):
            v = [d["mean_subsample"][i]] + [d["trimmed"][t][i] for t in d["trimmed"]]
            ax.plot([max(t, 1e-4) for t in tr], np.array(v) / abs(d["mean_subsample"][i]), mk + "-", color=COL[k], ms=3,
                    label=LAB[i] + " " + WN[k])
    ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log")
    ax.set_xlabel("fraction of largest |J| removed (1e-4 = none)"); ax.set_ylabel("trimmed mean / |full mean|")
    ax.legend(fontsize=6, frameon=False, ncol=2)
    save(fig, "figure6_tails")

    # Supplementary figure S1: cutout-level consistency
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.8))
    for i in range(3):
        for wi, w in enumerate(("holdout_611", "holdout_1300", "stage1_611", "stage1_1300")):
            k = "pooled_611" if w.endswith("611") else "pooled_1300"
            m = [b["base"]["s"][i] / b["base"]["n"] / SR[k] for b in blocks[w]]
            ax[i].plot(np.full(len(m), wi) + rng.uniform(-0.15, 0.15, len(m)), m, "o", ms=3, color=COL[k])
        ax[i].axhline(0, color="k", lw=0.6); ax[i].set_xticks(range(4))
        ax[i].set_xticklabels(["holdout 611", "holdout 1300", "stage 1 611", "stage 1 1300"], fontsize=7)
        ax[i].set_title(LAB[i]); ax[i].set_ylabel("cutout mean / $S_{rms}$")
    save(fig, "figureS1_cutout_consistency")

    # Figure 4: net drift into each alignment weight
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3), sharey=True)
    for a, (w, lab) in zip(ax, (("pooled_611", WN["pooled_611"]), ("pooled_1300", WN["pooled_1300"]))):
        for j, (key, nm, col) in enumerate((("N_gaussian", "Gaussian reference", "#9aa5b1"), ("N", "DNS", "#333333"))):
            est = np.array(derived[w][key]["est"]); lo = np.array(derived[w][key]["lo"]); hi = np.array(derived[w][key]["hi"])
            x = np.arange(3) + 0.38 * (j - 0.5)
            a.bar(x, est, 0.38, color=col, label=nm)
            a.errorbar(x, est, yerr=[est - lo, hi - est], fmt="none", color="k", lw=0.8, capsize=2)
        a.axhline(0, color="k", lw=0.6); a.set_xticks(range(3))
        a.set_xticklabels([r"$N_1$ (into $e_1$)", r"$N_2$ (into $e_2$)", r"$N_3$ (into $e_3$)"]); a.set_title(lab)
    ax[0].set_ylabel(r"mean restricted-Euler drift $N_i/S_{\rm rms}$"); ax[0].legend(frameon=False, fontsize=8)
    save(fig, "figure4_net_drift")
    return summary
