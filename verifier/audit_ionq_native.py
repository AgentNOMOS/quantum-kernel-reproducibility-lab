"""Independent audit of the IonQ native-gate constructions (no network, no hardware).

1. Gate conventions: qiskit_ionq.ZZGate / GPIGate / GPI2Gate matrices are compared with the
   matrices printed in IonQ's native-gate guide (docs.ionq.com, read 2026-10-14). Wrong readings of
   "turns" are shown to differ, so the check can actually fail.
2. Constructions: for every kernel entry and every arm (144-ZZ full depth, 108-ZZ fused, each with
   and without quarter-turn folding), the FULL 512x512 unitary is compared with the historical
   logical circuit up to global phase, before and after lowering to {GPI, GPI2, ZZ}.
3. Lowering preserves the two-qubit layer: same number of ZZ, same (pair, angle) multiset.
4. Payload semantics: every JSON payload in experiments/ionq_native_payloads/ is re-read and
   simulated with ONLY the documented matrices (pure numpy, no qiskit-ionq), and compared with the
   reference unitary. Phase/angle ranges are reported against the documented limits.
Run from the repo root: python verifier/audit_ionq_native.py  ->  reports/ionq_native_audit.json
"""
import glob, json, os, socket, sys
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
socket.socket.connect = lambda *a, **k: (_ for _ in ()).throw(RuntimeError("network attempted"))

import numpy as np
from qiskit.quantum_info import Operator
from qiskit_ionq import ZZGate, GPIGate, GPI2Gate
from qkrl import M
from qkrl.circuits import logical_overlap_circuit, ionq_zz_overlap_circuit, count_2q
from qkrl.ionq_compile import to_forte_native

Z = np.diag([1.0, -1.0]); ZZ = np.kron(Z, Z)
def doc_zz(t):  return np.diag([np.exp(-1j*np.pi*t), np.exp(1j*np.pi*t), np.exp(1j*np.pi*t), np.exp(-1j*np.pi*t)])
def doc_gpi(p): return np.array([[0, np.exp(-2j*np.pi*p)], [np.exp(2j*np.pi*p), 0]])
def doc_gpi2(p): return np.array([[1, -1j*np.exp(-2j*np.pi*p)], [-1j*np.exp(2j*np.pi*p), 1]]) / np.sqrt(2)
def expm_diag(c): return np.diag(np.exp(c * np.diag(ZZ)))          # exp(c * Z(x)Z), ZZ diagonal


def phase_dist(A, B):
    """min over global phase of max|A - e^{i g} B|; 0 means equal up to global phase."""
    k = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    g = A[k] / B[k]; g /= abs(g)
    return float(np.max(np.abs(A - g * B)))


def conventions():
    rng = np.random.default_rng(1); out = {}
    ts = np.concatenate([rng.uniform(-0.5, 0.5, 50), [0.25, -0.25, 0.3113, -0.3113]])
    out["zz_qiskit_ionq_vs_doc_matrix"] = max(np.max(np.abs(np.asarray(ZZGate(t)) - doc_zz(t))) for t in ts)
    out["zz_doc_vs_expm(-i*pi*t*ZZ)"] = max(np.max(np.abs(doc_zz(t) - expm_diag(-1j*np.pi*t))) for t in ts)
    out["zz_doc_vs_qiskit_RZZ(2*pi*t)_radians"] = max(
        np.max(np.abs(Operator(ZZGate(t).definition).data - doc_zz(t))) for t in ts)
    # deliberately wrong readings -- these MUST differ
    out["WRONG_reading_exp(-i*2pi*t*ZZ)_max_diff"] = max(np.max(np.abs(doc_zz(t) - expm_diag(-2j*np.pi*t))) for t in ts)
    out["WRONG_reading_exp(-i*t*ZZ)_radians_max_diff"] = max(np.max(np.abs(doc_zz(t) - expm_diag(-1j*t))) for t in ts)
    ps = rng.uniform(-1, 1, 50)
    out["gpi_qiskit_ionq_vs_doc"] = max(np.max(np.abs(np.asarray(GPIGate(p)) - doc_gpi(p))) for p in ps)
    out["gpi2_qiskit_ionq_vs_doc"] = max(np.max(np.abs(np.asarray(GPI2Gate(p)) - doc_gpi2(p))) for p in ps)
    out["gpi_period_1_turn"] = max(np.max(np.abs(doc_gpi(p) - doc_gpi(p + 1))) for p in ps)
    # quarter-turn fold identity: ZZ(t) = ZZ(t - k/2) . (Z(x)Z)^k up to global phase
    out["fold_identity_max_phase_dist"] = max(phase_dist(doc_zz(t), doc_zz(t - k/2) @ np.linalg.matrix_power(ZZ, k))
                                              for t in ts for k in (-1, 1, 2))
    ok = (out["zz_qiskit_ionq_vs_doc_matrix"] < 1e-14 and out["zz_doc_vs_expm(-i*pi*t*ZZ)"] < 1e-14
          and out["zz_doc_vs_qiskit_RZZ(2*pi*t)_radians"] < 1e-12 and out["gpi_qiskit_ionq_vs_doc"] < 1e-14
          and out["gpi2_qiskit_ionq_vs_doc"] < 1e-14 and out["fold_identity_max_phase_dist"] < 1e-12
          and out["WRONG_reading_exp(-i*2pi*t*ZZ)_max_diff"] > 0.1 and out["WRONG_reading_exp(-i*t*ZZ)_radians_max_diff"] > 0.1)
    return ok, {k: float(v) for k, v in out.items()}


def zz_multiset(qc):
    return sorted((tuple(sorted(qc.find_bit(q).index for q in i.qubits)), round(float(i.operation.params[0]), 12))
                  for i in qc.data if i.operation.name == "zz")


def simulate_payload(payload):
    """Unitary of an IonQ native JSON circuit using ONLY the documented matrices.
    Qubit q <-> bit q of the basis index (qiskit little-endian) so it can be compared with Operator."""
    U = np.eye(2 ** M, dtype=complex).reshape([2] * M + [2 ** M])          # axes: q8..q0, column
    def ax(q): return M - 1 - q
    for g in payload["circuit"]:
        if g["gate"] in ("gpi", "gpi2"):
            G = (doc_gpi if g["gate"] == "gpi" else doc_gpi2)(g["phase"])
            U = np.moveaxis(np.tensordot(G, U, axes=([1], [ax(g["target"])])), 0, ax(g["target"]))
        elif g["gate"] == "zz":
            a, b = g["targets"]; t = g["angle"]
            sa = np.array([1, -1]).reshape([2 if i == ax(a) else 1 for i in range(M)] + [1])
            sb = np.array([1, -1]).reshape([2 if i == ax(b) else 1 for i in range(M)] + [1])
            U = U * np.exp(-1j * np.pi * t * sa * sb)
        else:
            raise ValueError(g["gate"])
    return U.reshape(2 ** M, 2 ** M)


def main():
    man = json.load(open(os.path.join(ROOT, "data/kernel_pair_manifest.json")))
    P = {p["pair_id"]: (np.array(p["sample_a"]["feature_vector_scaled"]), np.array(p["sample_b"]["feature_vector_scaled"]),
                        p["K_exact_er01a"]) for p in man["pairs"]}
    rep = {}
    ok_conv, rep["conventions"] = conventions()
    print("conventions:", "PASS" if ok_conv else "FAIL")

    rows, ok_c = {}, True
    refU = {pid: Operator(logical_overlap_circuit(a, b, measure=False)).data for pid, (a, b, _) in P.items()}
    for arm, kw, expect in [("A_full_depth", dict(fuse=False), 144), ("B_fused", dict(fuse=True), 108),
                            ("A_full_depth_folded", dict(fuse=False, fold_quarter_turn=True), 144),
                            ("B_fused_folded", dict(fuse=True, fold_quarter_turn=True), 108)]:
        for pid, (a, b, K) in P.items():
            pre = ionq_zz_overlap_circuit(a, b, measure=False, **kw)
            post = to_forte_native(pre)
            e = 0 if (pid == "P0" and arm.startswith("B")) else expect
            d_pre = phase_dist(Operator(pre).data, refU[pid]); d_post = phase_dist(Operator(post).data, refU[pid])
            ang = [x[1] for x in zz_multiset(post)]
            r = dict(zz_pre=count_2q(pre), zz_post=count_2q(post), zz_multiset_preserved=zz_multiset(pre) == zz_multiset(post),
                     unitary_phase_dist_pre=d_pre, unitary_phase_dist_post=d_post,
                     zz_angle_range_turns=[min(ang), max(ang)] if ang else None,
                     z_compensations=int(pre.count_ops().get("z", 0)))
            r["pass"] = bool(r["zz_pre"] == e and r["zz_post"] == e and r["zz_multiset_preserved"] and d_pre < 1e-10 and d_post < 1e-10
                             and (not arm.endswith("folded") or not ang or max(abs(x) for x in ang) <= 0.25 + 1e-12))
            ok_c &= r["pass"]; rows[f"{arm}:{pid}"] = r
    rep["constructions"] = rows
    print("constructions:", "PASS" if ok_c else "FAIL", f"({sum(r['pass'] for r in rows.values())}/{len(rows)})")

    pay, ok_p = {}, True
    files = sorted(glob.glob(os.path.join(ROOT, "experiments/ionq_native_payloads/**/*.json"), recursive=True))
    for f in files:
        d = json.load(open(f)); name = os.path.relpath(f, ROOT)
        pid = "P0" if ("_P0" in name or "spam" in name or "half_depth" in name) else name.split("_")[-1][:2]
        U = simulate_payload(d)
        ph = [g["phase"] for g in d["circuit"] if g["gate"] != "zz"]; an = [g["angle"] for g in d["circuit"] if g["gate"] == "zz"]
        if "half_depth" in name:
            dist = float(abs(abs(U[0, 0]) - 1.0)); target = "self-mirror: |<0|U|0>| = 1"
        else:
            dist = phase_dist(U, refU[pid]); target = "historical logical unitary"
        r = dict(compared_with=target, phase_dist=dist, n_zz=len(an), n_1q=len(ph),
                 phase_range=[min(ph), max(ph)] if ph else None, angle_range=[min(an), max(an)] if an else None,
                 P_allzero=float(abs(U[0, 0]) ** 2))
        r["phases_within_documented_-1_to_1"] = all(-1 <= x <= 1 for x in ph)
        r["zz_angles_within_quarter_turn"] = all(abs(x) <= 0.25 + 1e-12 for x in an)
        r["pass"] = bool(dist < 1e-10 and r["phases_within_documented_-1_to_1"])
        ok_p &= r["pass"]; pay[name] = r
    rep["payloads"] = pay
    print("payloads:", "PASS" if ok_p else "FAIL", f"({sum(r['pass'] for r in pay.values())}/{len(pay)})")
    rep["unresolved_hardware_restrictions"] = [
        "ZZ angle range: not documented (MS documents 0..0.25 turns; nothing for ZZ). Folded payloads stay within +/-0.25.",
        "Negative ZZ angles: not documented as accepted or rejected. Folded payloads still contain negative angles; "
        "an exact rewrite to non-negative angles is possible (ZZ(-t) = ZZ(1/2 - t).(Z(x)Z) up to global phase) but is not applied.",
        "ZZ angle precision/resolution: not documented (MS: ~3 decimal places stated).",
        "Gate parallelism, physical pulse implementation, max gates per circuit: not documented.",
        "Whether the access path offered by Qollab (IonQ cloud direct, AWS Braket, Azure, other) accepts Forte native ZZ: unverified."]
    rep["verdict"] = "PASS" if (ok_conv and ok_c and ok_p) else "FAIL"
    os.makedirs(os.path.join(ROOT, "reports"), exist_ok=True)
    json.dump(rep, open(os.path.join(ROOT, "reports/ionq_native_audit.json"), "w"), indent=1)
    print("verdict:", rep["verdict"]); return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
