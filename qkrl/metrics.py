"""Kernel-contrast metrics, exactly as defined for the 2026-08-23 IBM run (ER-01B).

K_hw          = count('000000000') / shots         (all-zero probability of U(z)^dag U(x)|0>)
dynamic range = K(P1) - K(P2)                       (max-kernel pair minus min-kernel pair)
retained      = DR_hw / DR_exact
self-overlap  = K_hw(P0), exact value 1
depolarising  P_obs = f * K_exact + (1 - f) / 512   (binomial MLE of the single parameter f,
                                                    95 % profile-likelihood interval)
"""
import numpy as np

UNIFORM = 1.0 / 512


def k_hw(counts, shots=None):
    shots = shots or sum(counts.values())
    return counts.get("0" * 9, 0) / shots


def wilson(k, n, z=1.959963984540054):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))


def depolarising_fit(K_exact, allzero_counts, shots, grid=200000, f_max=1.0):
    """Grid MLE of f on (0, f_max]; 95 % profile-likelihood CI (delta NLL = 1.920729)."""
    K = np.asarray(K_exact, float); k = np.asarray(allzero_counts, float)
    fs = np.linspace(1e-6, f_max, grid)
    p = np.clip(fs[:, None] * K[None, :] + (1 - fs[:, None]) * UNIFORM, 1e-12, 1 - 1e-12)
    nll = -(k * np.log(p) + (shots - k) * np.log(1 - p)).sum(axis=1)
    i = int(np.argmin(nll)); ok = nll <= nll[i] + 1.920729
    pred = fs[i] * K + (1 - fs[i]) * UNIFORM
    chi2 = float(np.sum((k - pred * shots) ** 2 / (pred * shots * (1 - pred))))
    return dict(f_mle=float(fs[i]), f_ci95=[float(fs[np.argmax(ok)]), float(fs[len(fs) - 1 - np.argmax(ok[::-1])])],
                chi2=chi2, dof=len(K) - 1, predicted=pred.tolist())


def summarize(pairs, f_max=1.0):
    """pairs: list of dict(pair_id, K_exact, allzero, shots). Returns the ER-01B metric set.
    f_max=0.05 reproduces the historical ER-01B fitting grid bit-for-bit; use 1.0 for new data."""
    by = {p["pair_id"]: p for p in pairs}
    rows = []
    for p in pairs:
        kh = p["allzero"] / p["shots"]; lo, hi = wilson(p["allzero"], p["shots"])
        rows.append(dict(pair_id=p["pair_id"], K_exact=p["K_exact"], allzero=p["allzero"], shots=p["shots"],
                         K_hw=kh, signed_error=kh - p["K_exact"], wilson95=[lo, hi],
                         exact_inside_ci=bool(lo <= p["K_exact"] <= hi)))
    r = {x["pair_id"]: x for x in rows}
    ids = ["P1", "P2", "P3", "P4"]
    dr_exact = by["P1"]["K_exact"] - by["P2"]["K_exact"]
    dr_hw = r["P1"]["K_hw"] - r["P2"]["K_hw"]
    return dict(
        rows=rows,
        mae_all=float(np.mean([abs(x["signed_error"]) for x in rows])),
        mae_p1_p4=float(np.mean([abs(r[i]["signed_error"]) for i in ids])),
        spearman_p1_p4=spearman([by[i]["K_exact"] for i in ids], [r[i]["K_hw"] for i in ids]),
        self_overlap_P0=r["P0"]["K_hw"],
        dynamic_range_exact=dr_exact, dynamic_range_hw=dr_hw, dynamic_range_retained=dr_hw / dr_exact,
        depolarising=depolarising_fit([p["K_exact"] for p in pairs], [p["allzero"] for p in pairs],
                                      pairs[0]["shots"], f_max=f_max),
    )
