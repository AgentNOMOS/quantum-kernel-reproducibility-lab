"""Qiskit circuit constructions for the kernel-overlap workload.

Three representations of the SAME kernel entry k(x, z) = P(all-zero) of U(z)^dag U(x)|0>:

1. logical_overlap_circuit   -- the circuit executed on ibm_marrakesh (2026-08-23), built exactly
                                as the historical ER-01B code built it: 144 RYY, 4 H layers.
2. fused_overlap_circuit     -- exact algebraic simplification (ER-01B2): the two adjacent YY
                                layers (and Z layers) in the middle fuse, leaving 3 YY blocks
                                = 108 RYY. For x == z everything cancels to the empty circuit.
3. ionq_zz_overlap_circuit   -- the fused circuit with each YY block rotated into a ZZ block and
                                every pairwise phase emitted as ONE IonQ-native ZZ(theta) gate
                                (qiskit_ionq.ZZGate, angle in turns). Single-qubit gates are later
                                lowered to GPI/GPI2 by ionq_compile.py.
"""
import numpy as np
from qiskit import QuantumCircuit
from . import M, ALPHA, PAIRS


# --------------------------------------------------------------- 1. historical logical circuit
def _u_phi(qc, x):
    for k in range(M):
        qc.rz(-2.0 * ALPHA * float(x[k]), k)                       # exp(+i a x_k Z_k)
    for (l, m) in PAIRS:
        th = ALPHA * (np.pi - float(x[l])) * (np.pi - float(x[m]))
        qc.ryy(-2.0 * th, l, m)                                    # exp(+i th Y_l Y_m)


def _u_map(x):
    qc = QuantumCircuit(M, name="U")
    for _ in range(2):
        qc.h(range(M))
        _u_phi(qc, x)
    return qc


def logical_overlap_circuit(x, z, measure=True):
    """U(z)^dag U(x) followed by a 9-bit measurement (historical ER-01B construction)."""
    qc = QuantumCircuit(M, M)
    qc.compose(_u_map(x), inplace=True)
    qc.compose(_u_map(z).inverse(), inplace=True)
    if measure:
        qc.barrier()
        qc.measure(range(M), range(M))
    return qc


# --------------------------------------------------------------- 2. exact layer fusion
def _z_angles(w, s=1.0):
    return np.array([s * ALPHA * float(w[k]) for k in range(M)])


def _yy_angles(w, s=1.0):
    return np.array([s * ALPHA * (np.pi - float(w[l])) * (np.pi - float(w[m])) for (l, m) in PAIRS])


def build_layers(x, z):
    """Layer list of U(z)^dag U(x); ('RZ', a) means prod exp(i a_q Z_q), ('RYY', a) likewise."""
    L = [("H", None), ("RZ", _z_angles(x)), ("RYY", _yy_angles(x)),
         ("H", None), ("RZ", _z_angles(x)), ("RYY", _yy_angles(x))]
    L += [("RYY", _yy_angles(z, -1.0)), ("RZ", _z_angles(z, -1.0)), ("H", None),
          ("RYY", _yy_angles(z, -1.0)), ("RZ", _z_angles(z, -1.0)), ("H", None)]
    return L


def simplify(L):
    """Exact rewrites to fixpoint: fuse adjacent RZ/RYY layers (generators commute), cancel H.H,
    drop layers whose angles are all EXACTLY zero. No tolerance anywhere."""
    changed = True
    while changed:
        changed, out, i = False, [], 0
        while i < len(L):
            k, a = L[i]
            if i + 1 < len(L):
                k2, a2 = L[i + 1]
                if k in ("RZ", "RYY") and k2 == k:
                    out.append((k, a + a2)); i += 2; changed = True; continue
                if k == "H" and k2 == "H":
                    i += 2; changed = True; continue
            if k in ("RZ", "RYY") and np.all(a == 0.0):
                i += 1; changed = True; continue
            out.append((k, a)); i += 1
        L = out
    return L


def fused_overlap_circuit(x, z, measure=True):
    qc = QuantumCircuit(M, M)
    for k, a in simplify(build_layers(x, z)):
        if k == "H":
            qc.h(range(M))
        elif k == "RZ":
            for q in range(M):
                if a[q] != 0.0:
                    qc.rz(-2.0 * a[q], q)
        else:
            for t, (l, m) in enumerate(PAIRS):
                if a[t] != 0.0:
                    qc.ryy(-2.0 * a[t], l, m)
    if measure:
        qc.barrier()
        qc.measure(range(M), range(M))
    return qc


# --------------------------------------------------------------- 3. IonQ-native ZZ form
def _v(qc, q):          # operator V = S.H, with V Z V^dag = Y
    qc.h(q); qc.s(q)


def _vdg(qc, q):
    qc.sdg(q); qc.h(q)


def ionq_zz_overlap_circuit(x, z, measure=True, fuse=True, fold_quarter_turn=False):
    """fuse=True : exact fused form, 108 ZZ (0 for x == z).
    fuse=False: layer-for-layer image of the historical circuit, 144 ZZ for every pair including
                the self-overlap control, i.e. the same logical workload that ran on IBM.
    """
    layers = build_layers(x, z)
    return ionq_zz_from_layers(simplify(layers) if fuse else layers, measure, fold_quarter_turn)


def ionq_zz_from_layers(layers, measure=True, fold_quarter_turn=False):
    """Each YY block exp(i sum th Y_l Y_m) = V^{(x)9} exp(i sum th Z_l Z_m) V^dag{(x)9}; each
    exp(i th Z_l Z_m) is exactly IonQ ZZ(t) = exp(-i pi t Z(x)Z) with t = -th/pi (turns).
    fold_quarter_turn=True keeps every ZZ angle in [-1/4, 1/4] turns using the exact identity
    ZZ(t) = ZZ(t - k/2) . (Z(x)Z)^k up to global phase; the Z(x)Z factor is two 1q Z gates.
    IonQ documents no ZZ angle range; this option exists in case one applies."""
    from qiskit_ionq import ZZGate
    qc = QuantumCircuit(M, M)
    for k, a in layers:
        if k == "H":
            qc.h(range(M))
        elif k == "RZ":
            for q in range(M):
                if a[q] != 0.0:
                    qc.rz(-2.0 * a[q], q)
        else:
            idx = [t for t in range(len(PAIRS)) if a[t] != 0.0]
            if not idx:
                continue
            for q in range(M):
                _vdg(qc, q)
            for t in idx:
                l, m = PAIRS[t]
                turns = -float(a[t]) / np.pi
                if fold_quarter_turn:
                    folded = (turns + 0.25) % 0.5 - 0.25
                    k_ = int(round((turns - folded) / 0.5))
                    turns = folded
                    if k_ % 2:
                        qc.z(l); qc.z(m)
                qc.append(ZZGate(turns), [l, m])
            for q in range(M):
                _v(qc, q)
    if measure:
        qc.barrier()
        qc.measure(range(M), range(M))
    return qc


def half_depth_self_mirror_layers(x):
    """Control only (not part of the historical workload): reps=1 self-overlap
    U1(x)^dag U1(x), U1 = U_phi H, kept unfused -> 72 ZZ, exact value 1."""
    return [("H", None), ("RZ", _z_angles(x)), ("RYY", _yy_angles(x)),
            ("RYY", _yy_angles(x, -1.0)), ("RZ", _z_angles(x, -1.0)), ("H", None)]


def count_2q(qc):
    return sum(1 for inst in qc.data if inst.operation.num_qubits == 2 and inst.operation.name != "barrier")


def depth_2q(qc):
    return qc.depth(lambda inst: inst.operation.num_qubits == 2 and inst.operation.name != "barrier")
