"""Velocity-gradient instrument: finite differences, strain eigenframe,
alignment weights and pairwise alignment currents.

Conventions
-----------
A = grad u with A[i, j] = d u_i / d x_j;  S = (A + A^T)/2;  omega = curl u.
Eigenvalues of S are ordered lambda_1 > lambda_2 > lambda_3 with eigenvectors e_i.
Alignment weights c_i = (omega_hat . e_i)^2, sum_i c_i = 1.
Edges are ordered (1<-2, 1<-3, 2<-3). J[i<-j] is the current into weight i from
weight j, so that dc_i/dt|_RE = sum_{j != i} J[i<-j] and J[j<-i] = -J[i<-j]:

    J[i<-j] = 2 c_i c_j (lambda_i - lambda_j) (1 - chi_ij),
    chi_ij  = |omega|^2 / (4 (lambda_i - lambda_j)^2).
"""
import numpy as np

STENCILS = {
    4: ((2 / 3.0, 1), (-1 / 12.0, 2)),
    6: ((3 / 4.0, 1), (-3 / 20.0, 2), (1 / 60.0, 3)),
    8: ((4 / 5.0, 1), (-1 / 5.0, 2), (4 / 105.0, 3), (-1 / 280.0, 4)),
}
EDGES = ((0, 1), (0, 2), (1, 2))
EDGE_LABELS = ("1<-2", "1<-3", "2<-3")


def derivative(f, axis, dx, order=6):
    """Central finite difference along one axis (periodic indexing; callers
    use only points far enough from the faces for the stencil to be exact)."""
    out = np.zeros_like(f)
    for coef, shift in STENCILS[order]:
        out += coef * (np.roll(f, -shift, axis=axis) - np.roll(f, shift, axis=axis))
    return out / dx


def gradient(u, dx, order=6):
    """A[i, j] = d u_i / d x_j for a velocity field u of shape (3, n0, n1, n2)
    whose component i is the velocity along array axis i."""
    return np.stack([np.stack([derivative(u[i], j, dx, order) for j in range(3)]) for i in range(3)])


def strain_vorticity(A):
    """Strain tensor and vorticity vector from A of shape (3, 3, ...)."""
    S = 0.5 * (A + np.swapaxes(A, 0, 1))
    w = np.stack([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]])
    return S, w


def interior_points(u, dx, margin, order=6):
    """Strain (n, 3, 3) and vorticity (n, 3) at all points at least `margin`
    grid points from the faces of the cube."""
    m = margin
    A = gradient(u, dx, order)[:, :, m:-m, m:-m, m:-m].reshape(3, 3, -1)
    S, w = strain_vorticity(A)
    return np.ascontiguousarray(S.transpose(2, 0, 1)), np.ascontiguousarray(w.T)


def alignment_state(S, w):
    """Eigenframe and alignment quantities for point samples S (n, 3, 3), w (n, 3)."""
    ev, E = np.linalg.eigh(S)
    lam, E = ev[:, ::-1], E[:, :, ::-1]
    norm = np.clip(np.linalg.norm(w, axis=1), 1e-300, None)
    xi = w / norm[:, None]
    a = np.einsum("ni,nij->nj", xi, E)
    S_frob = np.sqrt((S ** 2).sum((1, 2)))
    gap = np.stack([lam[:, 0] - lam[:, 1], lam[:, 0] - lam[:, 2], lam[:, 1] - lam[:, 2]], 1)
    return {"lam": lam, "E": E, "xi": xi, "c": a ** 2, "om2": norm ** 2,
            "gap": gap, "g": gap / np.clip(S_frob, 1e-300, None)[:, None]}


def edge_currents(S, w):
    """Edge currents J (n, 3), the smallest relative gap over the edges (n,),
    and the alignment state."""
    st = alignment_state(S, w)
    c, lam = st["c"], st["lam"]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = st["om2"][:, None] / (4 * st["gap"] ** 2)
        J = np.stack([2 * c[:, i] * c[:, j] * (lam[:, i] - lam[:, j]) for i, j in EDGES], 1) * (1 - chi)
    return J, st["g"].min(1), st


def chi(st):
    """chi_ij on the three edges for an alignment state."""
    return st["om2"][:, None] / (4 * st["gap"] ** 2)
