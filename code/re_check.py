"""Independent check of the pairwise current formula against direct
restricted-Euler evolution, dA/dt = -A^2 + (1/3) tr(A^2) I.

For random traceless A the alignment weights are computed before and after a
small time step (fourth-order Runge-Kutta, centred difference) and compared
with sum_j J[i<-j]."""
import numpy as np


def _rhs(A):
    A2 = A @ A
    return -A2 + np.trace(A2) / 3 * np.eye(3)


def _step(A, h):
    k1 = _rhs(A); k2 = _rhs(A + h / 2 * k1); k3 = _rhs(A + h / 2 * k2); k4 = _rhs(A + h * k3)
    return A + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def _weights(A):
    S = 0.5 * (A + A.T)
    w = np.array([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]])
    lam, E = np.linalg.eigh(S)
    lam, E = lam[::-1], E[:, ::-1]
    return lam, ((w / np.linalg.norm(w)) @ E) ** 2, w


def run(n=2000, h=1e-5, seed=1):
    rng = np.random.default_rng(seed)
    worst, used = 0.0, 0
    for _ in range(n):
        X = rng.standard_normal((3, 3))
        A = X - np.trace(X) / 3 * np.eye(3)
        lam, c, w = _weights(A)
        if min(lam[0] - lam[1], lam[1] - lam[2]) < 0.05:
            continue
        om2 = w @ w
        J = {}
        for i, j in ((0, 1), (0, 2), (1, 2)):
            d = lam[i] - lam[j]
            J[(i, j)] = 2 * c[i] * c[j] * d * (1 - om2 / (4 * d * d))
        pred = np.array([J[(0, 1)] + J[(0, 2)], -J[(0, 1)] + J[(1, 2)], -J[(0, 2)] - J[(1, 2)]])
        num = (_weights(_step(A, h))[1] - _weights(_step(A, -h))[1]) / (2 * h)
        worst = max(worst, float(np.abs(num - pred).max() / (np.abs(pred).max() + 1e-12)))
        used += 1
    return {"n_samples": used, "max_relative_deviation": worst}


if __name__ == "__main__":
    print(run())
