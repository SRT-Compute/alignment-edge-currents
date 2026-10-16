"""Block statistics: pooled ratios, block bootstrap, sign-flip test, Holm correction."""
import numpy as np
from scipy import stats as st

PREDICTED_SIGNS = np.array([-1, +1, -1])
O_MIN = 0.05
N_BOOT = 10000


def pooled_bootstrap(blocks, n_boot=N_BOOT, seed=7, minus=None):
    """Pooled mean sum(s)/sum(n) over blocks and its block-bootstrap replicates.
    With `minus`, the paired difference to a second list of blocks."""
    s = np.array([b["s"] for b in blocks])
    n = np.array([b["n"] for b in blocks], float)
    est = s.sum(0) / n.sum()
    idx = np.random.default_rng(seed).integers(0, len(blocks), (n_boot, len(blocks)))
    boot = s[idx].sum(1) / n[idx].sum(1)[:, None]
    if minus is not None:
        s2 = np.array([b["s"] for b in minus])
        n2 = np.array([b["n"] for b in minus], float)
        est = est - s2.sum(0) / n2.sum()
        boot = boot - s2[idx].sum(1) / n2[idx].sum(1)[:, None]
    return est, boot


def cell(blocks):
    """Pooled mean, 95 % interval, organisation and pass flag per edge
    (pass = predicted sign, interval excluding zero, O >= 0.05)."""
    blocks = [b for b in blocks if b is not None and b["n"] > 0]
    if len(blocks) < 3 or sum(b["n"] for b in blocks) < 100:
        return None
    est, boot = pooled_bootstrap(blocks)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    ntot = sum(b["n"] for b in blocks)
    O = np.abs(est) / (np.array([b["a"] for b in blocks]).sum(0) / ntot)
    ok = (np.sign(est) == PREDICTED_SIGNS) & ((lo > 0) | (hi < 0)) & (O >= O_MIN)
    return {"mean": est, "lo": lo, "hi": hi, "O": O, "pass": ok, "n": ntot, "n_blocks": len(blocks)}


def sign_flip(blocks):
    """One-sided sign-flip randomisation test of the pooled mean in the
    predicted direction (exact for up to 20 blocks), with binomial and t tests."""
    m = np.array([b["s"] / b["n"] for b in blocks])
    w = np.array([b["n"] for b in blocks], float)
    w /= w.sum()
    nb = len(m)
    T = (w[:, None] * m).sum(0)
    if 2 ** nb <= 1_000_000:
        signs = ((np.arange(2 ** nb)[:, None] >> np.arange(nb)) & 1) * 2 - 1
    else:
        signs = np.random.default_rng(3).choice([-1, 1], size=(200000, nb))
    Ts = signs.astype(float) @ (w[:, None] * m)
    p = (PREDICTED_SIGNS * Ts >= PREDICTED_SIGNS * T - 1e-15).sum(0) / len(Ts)
    k = (np.sign(m) == PREDICTED_SIGNS).sum(0)
    p_binom = np.array([st.binomtest(int(k[i]), nb, 0.5, alternative="greater").pvalue for i in range(3)])
    p_t = np.array([st.ttest_1samp(PREDICTED_SIGNS[i] * m[:, i], 0, alternative="greater").pvalue for i in range(3)])
    return {"T": T, "p_signflip": p, "k_blocks": k, "n_blocks": nb, "p_binom": p_binom, "p_t": p_t,
            "exact": bool(2 ** nb <= 1_000_000)}


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        adj[i] = run
    return adj
