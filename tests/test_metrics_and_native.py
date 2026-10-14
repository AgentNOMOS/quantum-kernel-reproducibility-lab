"""Tests that need no experimental data: metric definitions and IonQ native-gate conventions."""
import glob, json, os
import numpy as np
import pytest
from qkrl import metrics as mt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_wilson_matches_closed_form():
    # k = 0, n = 100: closed-form Wilson upper bound z^2 / (n + z^2)
    z = 1.959963984540054
    lo, hi = mt.wilson(0, 100)
    assert lo == pytest.approx(0.0, abs=1e-15) and hi == pytest.approx(z * z / (100 + z * z), rel=1e-12)


def test_spearman_limits():
    assert mt.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert mt.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)


def test_depolarising_fit_is_exact_on_its_own_expectation():
    # analytic consistency of the estimator (expected counts, not a measurement)
    K = np.array([1.0, 0.9996, 0.0006, 0.5, 0.0948]); n = 10 ** 7
    for f in (0.0048, 0.3):
        k = (f * K + (1 - f) / 512) * n
        assert mt.depolarising_fit(K, k, n, grid=200001)["f_mle"] == pytest.approx(f, abs=1e-5)


def test_ionq_zz_matches_documented_matrix_and_rejects_wrong_units():
    from qiskit_ionq import ZZGate
    for t in (-0.31, -0.25, 0.1, 0.25):
        doc = np.diag(np.exp(-1j * np.pi * t * np.array([1, -1, -1, 1])))     # docs.ionq.com native-gate guide
        assert np.allclose(np.asarray(ZZGate(t)), doc, atol=1e-15)
        wrong = np.diag(np.exp(-2j * np.pi * t * np.array([1, -1, -1, 1])))   # "exponent in turns" misreading
        assert not np.allclose(doc, wrong, atol=1e-3)


@pytest.mark.parametrize("path", sorted(glob.glob(os.path.join(ROOT, "experiments/ionq_native_payloads/*.json"))))
def test_primary_payloads_respect_documented_ranges(path):
    d = json.load(open(path))
    assert d["gateset"] == "native" and "NEVER SUBMITTED" in d["format"]
    for g in d["circuit"]:
        assert g["gate"] in ("gpi", "gpi2", "zz")
        if g["gate"] == "zz":
            assert abs(g["angle"]) <= 0.25 + 1e-12 and len(set(g["targets"])) == 2
        else:
            assert -1 <= g["phase"] <= 1
