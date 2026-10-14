import hashlib, json, os, re
import numpy as np
import pytest
from qiskit.quantum_info import Statevector

from qkrl import feature_map as fm, metrics as mt
from qkrl.circuits import logical_overlap_circuit, fused_overlap_circuit, ionq_zz_overlap_circuit, count_2q

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = json.load(open(os.path.join(ROOT, "data/kernel_pair_manifest.json")))
PAIRS = [(p["pair_id"], np.array(p["sample_a"]["feature_vector_scaled"]),
          np.array(p["sample_b"]["feature_vector_scaled"]), p["K_exact_er01a"]) for p in MAN["pairs"]]
_TERM = os.path.join(ROOT, "evidence/ibm_marrakesh_job_da591ju1vhnc73flpn0g/IBM_02_terminal.json")
TERM = json.load(open(_TERM)) if os.path.exists(_TERM) else None
needs_raw = pytest.mark.skipif(TERM is None, reason="raw IBM job evidence withheld from the public release (WITHHELD_EVIDENCE.md)")


@pytest.mark.parametrize("pid,a,b,K", PAIRS)
def test_exact_kernel_numpy_and_qiskit_agree(pid, a, b, K):
    assert abs(fm.kernel(a, b) - K) < 1e-12
    sv = Statevector.from_instruction(logical_overlap_circuit(a, b, measure=False))
    assert abs(abs(sv.data[0]) ** 2 - K) < 1e-12


def test_kernel_matrix_is_symmetric_psd_unit_diagonal():
    X = [a for _, a, _, _ in PAIRS] + [b for _, _, b, _ in PAIRS]
    G = fm.kernel_matrix(X)
    assert np.allclose(G, G.T) and np.allclose(np.diag(G), 1.0)
    assert np.linalg.eigvalsh(G).min() > -1e-10


@needs_raw
def test_hardware_counts_total_20480_and_allzero_counts():
    counts = {r["pair_id"]: r["counts_from_sdk_get_counts"] for r in TERM["pub_results"]}
    assert sum(sum(c.values()) for c in counts.values()) == 20480
    assert [counts[p].get("000000000", 0) for p in ["P0", "P1", "P2", "P3", "P4"]] == [28, 27, 14, 17, 13]


@needs_raw
def test_contrast_metrics_reproduce_historical_values():
    K = {p: k for p, _, _, k in PAIRS}
    hw = [dict(pair_id=r["pair_id"], K_exact=K[r["pair_id"]], shots=4096,
               allzero=r["counts_from_sdk_get_counts"].get("000000000", 0)) for r in TERM["pub_results"]]
    s = mt.summarize(hw, f_max=0.05)
    assert s["self_overlap_P0"] == pytest.approx(0.0068359375, abs=1e-15)
    assert s["dynamic_range_retained"] == pytest.approx(0.003177162216662544, abs=1e-12)
    assert s["depolarising"]["f_mle"] == pytest.approx(0.004794928089640448, abs=1e-12)
    assert s["spearman_p1_p4"] == pytest.approx(0.8)
    assert not any(r["exact_inside_ci"] for r in s["rows"])


def test_fusion_gives_108_and_p0_empty():
    for pid, a, b, K in PAIRS:
        qc = fused_overlap_circuit(a, b)
        assert qc.count_ops().get("ryy", 0) == (0 if pid == "P0" else 108)
        assert logical_overlap_circuit(a, b).count_ops()["ryy"] == 144


@pytest.mark.parametrize("fuse,expected", [(False, 144), (True, 108)])
def test_ionq_native_counts_and_exactness(fuse, expected):
    from qkrl.ionq_compile import to_forte_native
    pid, a, b, K = PAIRS[3]
    t = to_forte_native(ionq_zz_overlap_circuit(a, b, fuse=fuse, fold_quarter_turn=True))
    assert set(t.count_ops()) <= {"gpi", "gpi2", "zz", "measure", "barrier"}
    assert count_2q(t) == expected
    assert max(abs(float(i.operation.params[0])) for i in t.data if i.operation.name == "zz") <= 0.25 + 1e-12
    sv = Statevector.from_instruction(t.remove_final_measurements(inplace=False))
    assert abs(abs(sv.data[0]) ** 2 - K) < 1e-10


def test_public_files_contain_no_account_identifiers_or_private_paths():
    # fragments are joined at runtime so this file does not match its own pattern
    pat = re.compile("|".join(["crn" + ":v1", "/" + "root/", "blue" + "mix",
                               r"api[_-]?key\s*[:=]\s*['\"][A-Za-z0-9]{16,}"]))
    for dp, dirs, fs in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in (".git", "__pycache__", ".pytest_cache")]
        for f in fs:
            if f.endswith((".json", ".md", ".py", ".qasm3", ".txt", ".toml")):
                assert not pat.search(open(os.path.join(dp, f), errors="ignore").read()), f
