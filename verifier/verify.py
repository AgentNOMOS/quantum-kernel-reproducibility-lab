"""Independent end-to-end verifier for the Quantum Kernel Reproducibility Lab.

Run from the repository root:   python verifier/verify.py [--evidence-dir DIR]
Writes reports/verification_report.json and exits non-zero if any check FAILS.

The public release withholds the raw IBM job files and provider metadata (see
WITHHELD_EVIDENCE.md). Checks that need them report WITHHELD, not PASS. A holder of the withheld
files can pass --evidence-dir to run them; the files are matched against the sealed SHA256 values
listed in WITHHELD_EVIDENCE.md.

Every number is recomputed from the files in this repository; the historical derived records in
evidence/historical_derived/ are only used as comparison targets. A socket guard makes any
network attempt an immediate failure: this verifier never contacts a quantum provider.
"""
import argparse, hashlib, json, os, socket, sys
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def _no_network(*a, **k):
    raise RuntimeError("network access attempted during verification")
socket.socket.connect = _no_network
socket.create_connection = _no_network

import numpy as np
import qiskit
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.qasm3 import dumps as q3dumps, loads as q3loads
from qiskit.transpiler import CouplingMap

from qkrl import M, PAIRS
from qkrl import feature_map as fm
from qkrl import metrics as mt
from qkrl.circuits import (logical_overlap_circuit, fused_overlap_circuit, ionq_zz_overlap_circuit,
                           count_2q, depth_2q)

TOL = 1e-12
IDS = ["P0", "P1", "P2", "P3", "P4"]
results = []


def rp(*p): return os.path.join(ROOT, *p)
def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()
def load(path): return json.load(open(rp(path)))


def check(cid, name, ok, **detail):
    if isinstance(ok, (bool, np.bool_)):
        ok = bool(ok)
    status = "PASS" if ok is True else ("FAIL" if ok is False else ok)
    results.append(dict(id=cid, name=name, status=status, **detail))
    print(f"[{status:^16}] {cid} {name}")


def pairs():
    man = load("data/kernel_pair_manifest.json")
    return [(p["pair_id"], np.array(p["sample_a"]["feature_vector_scaled"]),
             np.array(p["sample_b"]["feature_vector_scaled"]), p) for p in man["pairs"]]


def p_allzero_used(qc):
    """Exact P(all used qubits = 0) of a (possibly 156-qubit, physically laid out) circuit,
    after restricting it to the qubits it touches. Every used qubit is measured here."""
    used = sorted({qc.find_bit(q).index for inst in qc.data for q in inst.qubits})
    idx = {p: i for i, p in enumerate(used)}
    r = QuantumCircuit(len(used))
    measured = set()
    for inst in qc.data:
        n = inst.operation.name
        if n == "barrier":
            continue
        qs = [idx[qc.find_bit(q).index] for q in inst.qubits]
        if n == "measure":
            measured.add(qs[0]); continue
        r.append(inst.operation, qs)
    assert measured == set(range(len(used))), "not all used qubits are measured"
    return float(abs(Statevector.from_instruction(r).data[0]) ** 2), used


EVID = rp("evidence")


def ev(*p): return os.path.join(EVID, *p)


def main():
    global EVID
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence-dir", default=rp("evidence"))
    EVID = os.path.abspath(ap.parse_args().evidence_dir)
    P = pairs()
    K = {pid: p["K_exact_er01a"] for pid, _, _, p in P}

    # C01 -- input integrity
    dig = lambda v: hashlib.sha256(np.ascontiguousarray(v, dtype=np.float64).tobytes()).hexdigest()
    ok = all(dig(a) == p["sample_a"]["feature_vector_sha256"] and dig(b) == p["sample_b"]["feature_vector_sha256"]
             for _, a, b, p in P)
    distinct = {p["sample_a"]["sample_id"] for *_, p in P} | {p["sample_b"]["sample_id"] for *_, p in P}
    check("C01", "feature-vector digests match manifest; 9 distinct vectors", ok and len(distinct) == 9,
          distinct_vectors=len(distinct))

    # C02 -- exact kernel, pure numpy
    d = {pid: abs(fm.kernel(a, b) - K[pid]) for pid, a, b, _ in P}
    check("C02", "exact kernel (numpy statevector) reproduces sealed K_exact", max(d.values()) < TOL,
          max_abs_diff=max(d.values()), K_exact=K)

    # C03 -- exact kernel via the Qiskit logical circuit (independent construction)
    d = {pid: abs(abs(Statevector.from_instruction(logical_overlap_circuit(a, b, measure=False)).data[0]) ** 2 - K[pid])
         for pid, a, b, _ in P}
    check("C03", "Qiskit logical circuit P(all-zero) == K_exact", max(d.values()) < TOL, max_abs_diff=max(d.values()))

    # C04 -- logical circuit construction is byte-identical to the historical artefact
    same, ops = {}, {}
    for pid, a, b, _ in P:
        qc = logical_overlap_circuit(a, b)
        same[pid] = hashlib.sha256(q3dumps(qc).encode()).hexdigest() == sha(rp(f"circuits/historical/logical_{pid}.qasm3"))
        ops[pid] = {k: int(v) for k, v in qc.count_ops().items()}
    struct_ok = all(o.get("ryy") == 144 and o.get("rz") == 36 and o.get("h") == 36 for o in ops.values())
    if all(same.values()):
        status = struct_ok
    else:
        status = "VERSION_DEPENDENT" if qiskit.__version__ != "2.5.2" and struct_ok else False
    check("C04", "rebuilt logical circuits byte-identical to sealed QASM3 (qiskit 2.5.2)", status,
          byte_identical=same, logical_ops=ops["P1"], qiskit=qiskit.__version__)

    # C05 -- the executed IBM artefacts: 660 CZ / depth 1668 / 9 physical qubits, and re-transpilation
    tgt_path = ev("ibm_marrakesh_target_subset.json")
    cfg = dict(basis_gates=["cz", "rz", "sx", "x", "measure"], optimization_level=0, layout_method="sabre",
               routing_method="sabre", translation_method="translator", seed_transpiler=20260823)
    rec, retr = {}, {}
    for pid, a, b, _ in P:
        path = rp(f"circuits/historical/ibm_marrakesh_compiled_{pid}.qasm3")
        sealed = q3loads(open(path).read())
        p0, used = p_allzero_used(sealed)
        rec[pid] = dict(cz=int(sealed.count_ops().get("cz", 0)), depth=int(sealed.depth()), physical_qubits=used,
                        noiseless_allzero=p0, abs_diff=abs(p0 - K[pid]))
        if os.path.exists(tgt_path):
            cm = CouplingMap([tuple(e) for e in json.load(open(tgt_path))["coupling_edges"]])
            t = transpile(logical_overlap_circuit(a, b), coupling_map=cm, **cfg)
            retr[pid] = hashlib.sha256(q3dumps(t).encode()).hexdigest() == sha(path)
    ok = all(r["cz"] == 660 and r["depth"] == 1668 and len(r["physical_qubits"]) == 9 and r["abs_diff"] < TOL
             for r in rec.values())
    check("C05", "sealed IBM circuits: 660 CZ, depth 1668, 9 qubits, noiseless P(0)=K_exact", ok, per_pair=rec)
    if not os.path.exists(tgt_path):
        check("C06", "frozen IBM transpilation re-run offline reproduces sealed bytes", "WITHHELD",
              reason="needs the IBM coupling map (provider metadata), withheld from the public release")
    else:
        check("C06", "frozen IBM transpilation re-run offline reproduces sealed bytes",
              all(retr.values()) or ("VERSION_DEPENDENT" if qiskit.__version__ != "2.5.2" else False),
              byte_identical=retr)

    # C07 -- physical measurement evidence (IBM job da591ju1vhnc73flpn0g)
    wman = rp("WITHHELD_EVIDENCE.json")
    if os.path.exists(wman) and os.path.isdir(ev("ibm_marrakesh_job_da591ju1vhnc73flpn0g")):
        exp = json.load(open(wman))["files"]
        got = {f["path"]: (sha(ev(f["path"])) == f["sha256"]) if os.path.exists(ev(f["path"])) else None for f in exp}
        check("C00", "supplied withheld evidence matches the digests published in WITHHELD_EVIDENCE.json",
              all(v is True for v in got.values()), per_file=got)

    term_path = ev("ibm_marrakesh_job_da591ju1vhnc73flpn0g", "IBM_02_terminal.json")
    hist_path = ev("historical_derived", "er01b_stats_addendum.json")
    if os.path.exists(term_path) and os.path.exists(hist_path):
        hardware_checks(json.load(open(term_path)), json.load(open(hist_path)), K)
    else:
        for cid, name in [("C07", "raw IBM evidence: job id, backend, 20480 shots, PUB->artefact binding"),
                          ("C08", "kernel-contrast metrics recomputed from raw counts")]:
            check(cid, name, "WITHHELD", reason="raw IBM job output withheld from the public release; "
                  "run with --evidence-dir to verify (see WITHHELD_EVIDENCE.md)")
    offline_checks(P, K)


def hardware_checks(term, hist_doc, K):
    hw = []
    for r in term["pub_results"]:
        c = r["counts_from_sdk_get_counts"]
        hw.append(dict(pair_id=r["pair_id"], K_exact=K[r["pair_id"]], allzero=c.get("0" * 9, 0),
                       shots=sum(c.values()), distinct=len(c)))
    bind_ok = all(b["artifact_sha256"] == sha(rp(f"circuits/historical/ibm_marrakesh_compiled_{b['pair_id']}.qasm3"))
                  for b in term["pub_binding"])
    total = sum(h["shots"] for h in hw)
    check("C07", "raw IBM evidence: job id, backend, 5 x 4096 = 20480 shots, PUB->artefact digests bind",
          term["runtime_job_id"] == "da591ju1vhnc73flpn0g" and term["backend_from_job"] == "ibm_marrakesh"
          and term["final_status"] == "DONE" and total == 20480 and bind_ok and term["error_mitigation"] == "NONE",
          total_shots=total, allzero_counts={h["pair_id"]: h["allzero"] for h in hw},
          qpu_usage_seconds=term["usage_from_job"])

    s = mt.summarize(hw, f_max=0.05)
    hist = hist_doc["depolarizing_model"]
    ok = (abs(s["self_overlap_P0"] - 28 / 4096) < 1e-15 and abs(s["dynamic_range_retained"] - 0.003177162216662544) < 1e-12
          and abs(s["depolarising"]["f_mle"] - hist["f_mle"]) < 1e-12 and abs(s["mae_all"] - 0.515279794540962) < 1e-12
          and abs(s["spearman_p1_p4"] - 0.8) < 1e-12 and not any(r["exact_inside_ci"] for r in s["rows"]))
    check("C08", "kernel-contrast metrics recomputed from raw counts match sealed analysis", ok, metrics=s)


def offline_checks(P, K):

    # C09 -- offline-optimised IBM-style artefacts (NEVER executed on hardware)
    off = {}
    for fam, key in [("heavyhex_stock_transpiler_best", "ER-01B2 stock transpiler"),
                     ("heavyhex_exact_algebraic", "ER-01B2 exact algebraic"), ("heavyhex_paritynet", "ER-01B2.1"),
                     ("k9_paritynet", "ER-01B2.1 theoretical K9")]:
        for pid in IDS[1:]:
            qc = q3loads(open(rp(f"circuits/historical/offline_optimized/{fam}_{pid}.qasm3")).read())
            p0, used = p_allzero_used(qc)
            off[f"{fam}_{pid}"] = dict(two_qubit=count_2q(qc), ops={k: int(v) for k, v in qc.count_ops().items()},
                                       abs_diff=abs(p0 - K[pid]), source=key)
    ok = all(v["abs_diff"] < TOL for v in off.values())
    med = lambda fam: float(np.median([off[f"{fam}_{p}"]["two_qubit"] for p in IDS[1:]]))
    check("C09", "offline-optimised artefacts are exact and have the recorded 2Q counts",
          ok and med("heavyhex_stock_transpiler_best") == 496 and med("heavyhex_exact_algebraic") == 382
          and off["heavyhex_exact_algebraic_P1"]["two_qubit"] == 495
          and med("heavyhex_paritynet") == 300 and med("k9_paritynet") == 135,
          medians=dict(stock_transpiler_496=med("heavyhex_stock_transpiler_best"),
                       exact_algebraic_382=med("heavyhex_exact_algebraic"), heavyhex_paritynet_300=med("heavyhex_paritynet"),
                       k9_paritynet_135=med("k9_paritynet")), per_circuit=off)

    # C10 -- exact fusion: 144 -> 108 RYY (216 CNOT if each RYY costs 2), P0 -> empty
    fz = {}
    for pid, a, b, _ in P:
        qc = fused_overlap_circuit(a, b)
        p0 = abs(Statevector.from_instruction(qc.remove_final_measurements(inplace=False)).data[0]) ** 2
        fz[pid] = dict(ryy=int(qc.count_ops().get("ryy", 0)), abs_diff=abs(p0 - K[pid]))
    check("C10", "exact layer fusion: 108 RYY for P1-P4, 0 for P0, kernel unchanged",
          all(fz[p]["ryy"] == 108 for p in IDS[1:]) and fz["P0"]["ryy"] == 0 and max(v["abs_diff"] for v in fz.values()) < TOL,
          per_pair=fz)

    # C11 -- IonQ Forte-native lowering (local qiskit-ionq gate set; no API)
    from qkrl.ionq_compile import to_forte_native
    nat = {}
    for arm, kw in [("A_full_depth_144", dict(fuse=False)), ("B_fused_108", dict(fuse=True)),
                    ("A_full_depth_144_folded", dict(fuse=False, fold_quarter_turn=True)),
                    ("B_fused_108_folded", dict(fuse=True, fold_quarter_turn=True))]:
        for pid, a, b, _ in P:
            t = to_forte_native(ionq_zz_overlap_circuit(a, b, **kw))
            p0 = abs(Statevector.from_instruction(t.remove_final_measurements(inplace=False)).data[0]) ** 2
            ang = [float(i.operation.params[0]) for i in t.data if i.operation.name == "zz"]
            nat[f"{arm}:{pid}"] = dict(ops={k: int(v) for k, v in t.count_ops().items()}, zz=count_2q(t),
                                       zz_depth=depth_2q(t), abs_diff=abs(p0 - K[pid]),
                                       zz_turns_range=[min(ang), max(ang)] if ang else None)
    stock = to_forte_native(logical_overlap_circuit(*P[1][1:3]))
    ok = (all(nat[f"A_full_depth_144{s_}:{p}"]["zz"] == 144 for p in IDS for s_ in ("", "_folded"))
          and all(nat[f"B_fused_108{s_}:{p}"]["zz"] == 108 for p in IDS[1:] for s_ in ("", "_folded"))
          and nat["B_fused_108:P0"]["zz"] == 0 and max(v["abs_diff"] for v in nat.values()) < 1e-10
          and all(abs(x) <= 0.25 + 1e-12 for k_, v in nat.items() if "folded" in k_ and v["zz_turns_range"] for x in v["zz_turns_range"]))
    check("C11", "IonQ-native {GPI,GPI2,ZZ}: 144 ZZ (full depth) / 108 ZZ (fused), exact", ok, per_circuit=nat,
          stock_transpile_of_historical_circuit=dict(ops={k: int(v) for k, v in stock.count_ops().items()}, zz=count_2q(stock)))

    finish()


def finish():
    rep = dict(verifier="qkrl verifier v2", evidence_dir_supplied=os.path.exists(os.path.join(EVID, "ibm_marrakesh_job_da591ju1vhnc73flpn0g")), qiskit=qiskit.__version__, numpy=np.__version__, tolerance=TOL,
               network_guard="socket.connect disabled for the whole run", checks=results,
               summary={s_: sum(1 for r in results if r["status"] == s_) for s_ in ("PASS", "FAIL", "VERSION_DEPENDENT", "WITHHELD")})
    os.makedirs(rp("reports"), exist_ok=True)
    json.dump(rep, open(rp("reports/verification_report.json"), "w"), indent=1, default=float)
    print("\nsummary:", rep["summary"])
    sys.exit(1 if rep["summary"]["FAIL"] else 0)


if __name__ == "__main__":
    main()
