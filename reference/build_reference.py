"""Exact classical reference: the five historical kernel entries plus the full 9x9 Gram matrix of
the nine distinct input vectors (the Gram matrix is NEW in this package; it was not measured on
hardware). Run from the repo root: python reference/build_reference.py"""
import json, os, sys
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
from qkrl import feature_map as fm

man = json.load(open(os.path.join(ROOT, "data/kernel_pair_manifest.json")))
vec = {}
for p in man["pairs"]:
    for s in ("sample_a", "sample_b"):
        vec[p[s]["sample_id"]] = p[s]["feature_vector_scaled"]
ids = sorted(vec)
G = fm.kernel_matrix([vec[i] for i in ids])
entries = {p["pair_id"]: dict(a=p["sample_a"]["sample_id"], b=p["sample_b"]["sample_id"],
                              K_exact=fm.kernel(p["sample_a"]["feature_vector_scaled"], p["sample_b"]["feature_vector_scaled"]),
                              K_exact_sealed=p["K_exact_er01a"]) for p in man["pairs"]}
out = dict(definition="k(x,z)=|<0|U(z)^dag U(x)|0>|^2, U = U_phi H U_phi H, alpha=0.1, ORDERED, 9 qubits",
           method="exact numpy statevector (qkrl.feature_map), float64",
           historical_entries=entries, distinct_vectors=ids, gram_matrix=G.round(15).tolist(),
           gram_min_eigenvalue=float(np.linalg.eigvalsh(G).min()))
json.dump(out, open(os.path.join(ROOT, "reference/exact_kernel_reference.json"), "w"), indent=1)
for k, v in entries.items():
    print(k, v["a"], v["b"], f"{v['K_exact']:.15f}", f"|d|={abs(v['K_exact'] - v['K_exact_sealed']):.1e}")
