"""Quality gates applied to every cutout or channel block before analysis."""
import numpy as np

import instrument as I

G1_RANGE = (0.45, 0.55)        # <S:S> / <|omega|^2>
G2_MAX_ISOTROPIC = 1.0e-2      # <|div u|> / <|grad u|>
G3_RANGE = (-0.70, -0.30)      # mean longitudinal derivative skewness
G2_MAX_CHANNEL = 0.05
MAX_FAILED = 2                 # a window (band) fails if more than this many of 12 blocks fail


def isotropic(u, dx):
    m = 4
    A = I.gradient(u, dx)[:, :, m:-m, m:-m, m:-m]
    S, w = I.strain_vorticity(A)
    g1 = float((S ** 2).sum((0, 1)).mean() / (w ** 2).sum(0).mean())
    g2 = float(np.abs(A[0, 0] + A[1, 1] + A[2, 2]).mean() / np.sqrt((A ** 2).sum((0, 1))).mean())
    g3 = float(np.mean([(A[i, i] ** 3).mean() / (A[i, i] ** 2).mean() ** 1.5 for i in range(3)]))
    ok = G1_RANGE[0] <= g1 <= G1_RANGE[1] and g2 <= G2_MAX_ISOTROPIC and G3_RANGE[0] <= g3 <= G3_RANGE[1]
    return {"G1": g1, "G2": g2, "G3": g3, "pass": bool(ok)}


def channel(A):
    """A: (n, 3, 3) velocity gradients of one block."""
    g2 = float(np.abs(np.trace(A, axis1=1, axis2=2)).mean() / np.sqrt((A ** 2).sum((1, 2))).mean())
    return {"G2": g2, "pass": bool(g2 <= G2_MAX_CHANNEL)}
