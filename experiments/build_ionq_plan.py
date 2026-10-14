"""Build the proposed IonQ Forte circuit set, native payloads and cost table. LOCAL ONLY.

Nothing here submits a job: payloads are written to experiments/ionq_native_payloads/ in IonQ's
documented native-gate JSON shape ({"gate": "gpi"|"gpi2"|"zz", ...}, angles in turns) so a
reviewer can see exactly what would run. Run from the repo root: python experiments/build_ionq_plan.py
"""
import json, os, sys
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
from qiskit.quantum_info import Statevector
from qkrl.circuits import (logical_overlap_circuit, ionq_zz_overlap_circuit, ionq_zz_from_layers, half_depth_self_mirror_layers,
                           count_2q, depth_2q)
from qkrl.ionq_compile import to_forte_native

PRICING = {   # retrieved 2026-10-14; see docs/IONQ_HARDWARE_PLAN.md for sources and caveats
    "aws_braket_ionq_forte": dict(per_task=0.30, per_shot=0.08,
                                  source="https://aws.amazon.com/braket/pricing/ (retrieved 2026-10-14)"),
    "azure_quantum_ionq_forte1": dict(per_1q_gate_shot=0.0001645, per_2q_gate_shot=0.001121,
                                      min_per_program_em_off=25.7899, min_per_program_em_on=168.195,
                                      source="https://learn.microsoft.com/en-us/azure/quantum/pricing (updated 2026-04-23, retrieved 2026-10-14)"),
}
BUDGETS = {
    "MINIMUM_PILOT": dict(sets=["A_full_depth"], spam=True, ladder=False, shots=1000, replicates=1),
    "RECOMMENDED_REPRODUCTION": dict(sets=["A_full_depth", "B_fused"], spam=True, ladder=True, shots=1000, replicates=2),
    "EXTENDED_ROBUSTNESS": dict(sets=["A_full_depth", "B_fused", "E_stock_fixed_angle"], spam=True, ladder=True,
                                shots=2000, replicates=3),
}


def native_payload(t, variant):
    """IonQ API v0.3 native JSON. Phases are wrapped into [-1, 1] turns (documented accepted range;
    GPI/GPI2 are 1-turn periodic, so this is exact)."""
    ops = []
    for inst in t.data:
        n = inst.operation.name; qs = [t.find_bit(q).index for q in inst.qubits]
        if n in ("gpi", "gpi2"):
            ph = float(inst.operation.params[0])
            ph = ph - np.round(ph) if abs(ph) > 1 else ph
            assert -1 <= ph <= 1
            ops.append({"gate": n, "target": qs[0], "phase": ph})
        elif n == "zz":
            ops.append({"gate": "zz", "targets": qs, "angle": float(inst.operation.params[0])})
    return {"format": "IonQ native-gate circuit (gateset=native), generated locally, NEVER SUBMITTED",
            "variant": variant, "gateset": "native", "qubits": 9, "circuit": ops}


def reorder_blocks(qc, seed):
    """Exact control: permute ZZ gates inside each commuting block (between 1q layers)."""
    rng = np.random.default_rng(seed); out = qc.copy_empty_like(); block = []
    def flush():
        for i in rng.permutation(len(block)):
            out.append(block[i])
        block.clear()
    for inst in qc.data:
        if inst.operation.name == "zz":
            block.append(inst)
        else:
            flush(); out.append(inst)
    flush(); return out


def main():
    man = json.load(open(os.path.join(ROOT, "data/kernel_pair_manifest.json")))
    P = {p["pair_id"]: (np.array(p["sample_a"]["feature_vector_scaled"]), np.array(p["sample_b"]["feature_vector_scaled"]),
                        p["K_exact_er01a"]) for p in man["pairs"]}
    circuits = {}
    def add(cid, role, builder, K):
        """builder(fold) -> circuit. Primary payload = folded (all |ZZ angle| <= 1/4 turn, safe if a
        range limit exists); unfolded alternative kept in payloads/unfolded/. Costs use the folded counts."""
        alt = to_forte_native(builder(False))
        os.makedirs(os.path.join(ROOT, "experiments/ionq_native_payloads/unfolded"), exist_ok=True)
        json.dump(native_payload(alt, "unfolded: ZZ angles as constructed"),
                  open(os.path.join(ROOT, f"experiments/ionq_native_payloads/unfolded/{cid}.json"), "w"), indent=0)
        t = to_forte_native(builder(True))
        p0 = float(abs(Statevector.from_instruction(t.remove_final_measurements(inplace=False)).data[0]) ** 2)
        assert abs(p0 - K) < 1e-10, cid
        ops = {k: int(v) for k, v in t.count_ops().items()}
        circuits[cid] = dict(role=role, K_exact=K, native_ops=ops, n_2q=count_2q(t),
                             n_1q=ops.get("gpi", 0) + ops.get("gpi2", 0), zz_depth=depth_2q(t), ideal_allzero=p0)
        os.makedirs(os.path.join(ROOT, "experiments/ionq_native_payloads"), exist_ok=True)
        json.dump(native_payload(t, "folded: all ZZ angles in [-1/4, 1/4] turns (primary)"), open(os.path.join(ROOT, f"experiments/ionq_native_payloads/{cid}.json"), "w"), indent=0)
    for pid, (a, b, K) in P.items():
        add(f"A_full_depth_{pid}", "primary: same logical circuit as the IBM run, arbitrary-angle ZZ, 144 ZZ",
            lambda f, a=a, b=b: ionq_zz_overlap_circuit(a, b, fuse=False, fold_quarter_turn=f), K)
        add(f"A_full_depth_reordered_{pid}", "optional control: A with ZZ order permuted inside commuting blocks (exact)",
            lambda f, a=a, b=b: reorder_blocks(ionq_zz_overlap_circuit(a, b, fuse=False, fold_quarter_turn=f), 20261014), K)
        add(f"E_stock_fixed_angle_{pid}", "extended control: stock qiskit-ionq lowering of the historical circuit "
            "(each RYY -> 2 fixed ZZ(1/4)); all-to-all but fixed-angle, isolates connectivity from entangler type",
            lambda f, a=a, b=b: logical_overlap_circuit(a, b), K)
        if pid != "P0":
            add(f"B_fused_{pid}", "secondary: exact fused form, 108 ZZ (best-known exact native cost)",
                lambda f, a=a, b=b: ionq_zz_overlap_circuit(a, b, fuse=True, fold_quarter_turn=f), K)
    a0, b0 = P["P0"][:2]
    add("C_spam_P0_fused", "control: P0 fused = empty circuit, measures state-prep + readout only",
        lambda f: ionq_zz_overlap_circuit(a0, b0, fuse=True, fold_quarter_turn=f), 1.0)
    add("D_half_depth_P0", "control: reps=1 self-mirror of P0 (72 ZZ), depth-decay point",
        lambda f: ionq_zz_from_layers(half_depth_self_mirror_layers(a0), fold_quarter_turn=f), 1.0)

    plans = {}
    for name, B in BUDGETS.items():
        ids = [c for c in circuits if any(c.startswith(s + "_P") for s in B["sets"])]
        ids += ["C_spam_P0_fused"] if B["spam"] else []
        ids += ["D_half_depth_P0"] if B["ladder"] else []
        tasks = len(ids) * B["replicates"]; shots = tasks * B["shots"]
        br = PRICING["aws_braket_ionq_forte"]
        az = PRICING["azure_quantum_ionq_forte1"]
        az_cost = 0.0
        for c in ids:
            g = B["shots"] * (circuits[c]["n_1q"] * az["per_1q_gate_shot"] + circuits[c]["n_2q"] * az["per_2q_gate_shot"])
            az_cost += B["replicates"] * max(g, az["min_per_program_em_off"])
        plans[name] = dict(circuits=ids, n_circuits=len(ids), shots_per_circuit=B["shots"], replicates=B["replicates"],
                           tasks=tasks, total_shots=shots,
                           kernel_entries=sorted({c.split("_")[-1] for c in ids if c[0] in "AB"}),
                           max_2q_per_circuit=max(circuits[c]["n_2q"] for c in ids),
                           usd_aws_braket=round(tasks * br["per_task"] + shots * br["per_shot"], 2),
                           usd_azure_quantum_em_off=round(az_cost, 2))
    out = dict(backend_family="IonQ Forte (native 2Q gate ZZ, arbitrary angle, all-to-all)", qubits=9,
               pricing=PRICING, circuits=circuits, budgets=plans,
               qollab_credit_conversion="UNKNOWN - Qollab has not published a credit accounting model")
    json.dump(out, open(os.path.join(ROOT, "experiments/ionq_plan.json"), "w"), indent=1)
    for n, p in plans.items():
        print(n, {k: v for k, v in p.items() if k != "circuits"})
    for c, v in circuits.items():
        print(f"{c:28s} 2Q={v['n_2q']:4d} 1Q={v['n_1q']:4d} zz_depth={v['zz_depth']:3d} K={v['K_exact']:.6f}")


if __name__ == "__main__":
    main()
