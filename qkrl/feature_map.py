"""Exact (noise-free, shot-free) classical reference for the quantum kernel.

    U(x)     = U_phi(x) H^{(x)9} U_phi(x) H^{(x)9}                       (reps = 2, H acts first)
    U_phi(x) = prod_k exp(i a x_k Z_k) * prod_{l<m} exp(i a (pi-x_l)(pi-x_m) Y_l Y_m)
               ORDERED: all Z terms first, then the 36 YY terms in lexicographic order
    k(x, z)  = |<0| U(z)^dag U(x) |0>|^2

Pure numpy; independent of Qiskit so it can cross-check the circuit construction.
Adapted from the historical ER-01A implementation (same authors); ORDERED mode only, which is the
mode that was executed on hardware.
"""
import numpy as np
from . import M, ALPHA, PAIRS


def _rz_layer(S, angles):
    """exp(i * angle_q * Z_q) on every qubit q. Qubit 0 is the most significant axis."""
    for q in range(M):
        T = S.reshape(2 ** q, 2, 2 ** (M - q - 1))
        p = np.exp(1j * angles[q])
        T[:, 0, :] *= p
        T[:, 1, :] *= np.conj(p)
    return S


def _yy(S, l, m, theta):
    """exp(i theta Y_l Y_m) = cos(theta) I + i sin(theta) Y_l Y_m."""
    T = S.reshape(2 ** l, 2, 2 ** (m - l - 1), 2, 2 ** (M - m - 1))
    c, s = np.cos(theta), 1j * np.sin(theta)
    s00, s01 = T[:, 0, :, 0, :].copy(), T[:, 0, :, 1, :].copy()
    s10, s11 = T[:, 1, :, 0, :].copy(), T[:, 1, :, 1, :].copy()
    T[:, 0, :, 0, :] = c * s00 - s * s11
    T[:, 0, :, 1, :] = c * s01 + s * s10
    T[:, 1, :, 0, :] = c * s10 + s * s01
    T[:, 1, :, 1, :] = c * s11 - s * s00
    return S


def _hadamard_all(S):
    for q in range(M):
        T = S.reshape(2 ** q, 2, 2 ** (M - q - 1))
        u, v = T[:, 0, :].copy(), T[:, 1, :].copy()
        T[:, 0, :], T[:, 1, :] = u + v, u - v
    return S * 2.0 ** (-M / 2.0)


def angles(x, alpha=ALPHA):
    """(single-qubit angles [9], pair angles [36]) already multiplied by alpha."""
    x = np.asarray(x, dtype=float)
    single = alpha * x
    pair = np.array([alpha * (np.pi - x[l]) * (np.pi - x[m]) for l, m in PAIRS])
    return single, pair


def statevector(x, alpha=ALPHA):
    """|psi(x)> = U(x)|0...0> as a length-512 complex vector (qubit 0 = most significant)."""
    single, pair = angles(x, alpha)
    S = np.zeros(2 ** M, dtype=complex)
    S[0] = 1.0
    for _ in range(2):
        S = _hadamard_all(S)
        S = _rz_layer(S, single)
        for t, (l, m) in enumerate(PAIRS):
            S = _yy(S, l, m, pair[t])
    return S


def kernel(x, z, alpha=ALPHA):
    """Exact fidelity kernel k(x, z) = |<psi(z)|psi(x)>|^2."""
    return float(abs(np.vdot(statevector(z, alpha), statevector(x, alpha))) ** 2)


def kernel_matrix(X, alpha=ALPHA):
    P = np.array([statevector(x, alpha) for x in X])
    return np.abs(P.conj() @ P.T) ** 2
