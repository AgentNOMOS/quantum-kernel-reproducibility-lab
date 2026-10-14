"""Local Monte-Carlo power analysis for the preregistered decision rule (no hardware, no network).

Model: P(all-zero) = f * K_exact + (1 - f) / 512 for the five Arm-A circuits, binomial shots.
Decision rule (docs/IONQ_HARDWARE_PLAN.md section 6):
  IMPROVED                 if lower 95 % profile bound of f >= 0.05
  NO_MEANINGFUL_IMPROVEMENT if upper 95 % profile bound of f <  0.05
  INCONCLUSIVE             if the interval contains 0.05 (or, on hardware, a control fails)
f is not "fidelity"; it is the fraction of ideal kernel contrast surviving under this model.
"""
import json, os, sys
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
from qkrl.metrics import depolarising_fit

THRESHOLD = 0.05
K = np.array([p["K_exact_er01a"] for p in json.load(open(os.path.join(ROOT, "data/kernel_pair_manifest.json")))["pairs"]])


def verdict(fit):
    lo, hi = fit["f_ci95"]
    return "IMPROVED" if lo >= THRESHOLD else ("NO_MEANINGFUL_IMPROVEMENT" if hi < THRESHOLD else "INCONCLUSIVE")


def main(trials=400, seed=20261014):
    rng = np.random.default_rng(seed)
    out = {}
    for shots, reps in [(1000, 1), (1000, 2), (2000, 3)]:
        n = shots * reps
        for f in [0.0048, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.1, 0.2]:
            p = f * K + (1 - f) / 512
            tally = {"IMPROVED": 0, "NO_MEANINGFUL_IMPROVEMENT": 0, "INCONCLUSIVE": 0}
            for _ in range(trials):
                tally[verdict(depolarising_fit(K, rng.binomial(n, p), n, grid=20000))] += 1
            out[f"shots={shots}x{reps} f={f}"] = {k: v / trials for k, v in tally.items()}
            print(f"shots/circuit={n:5d}  f_true={f:<6}  " + "  ".join(f"{k}={v/trials:.3f}" for k, v in tally.items()))
    json.dump(dict(model="global depolarising, binomial", threshold=THRESHOLD, trials=trials, seed=seed, results=out),
              open(os.path.join(ROOT, "experiments/power_analysis.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
