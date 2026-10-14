"""Quantum Kernel Reproducibility Lab (QKRL).

Public reproduction package for a 9-qubit quantum-kernel workload that was executed once on
IBM Quantum hardware (ibm_marrakesh, job da591ju1vhnc73flpn0g, 2026-08-23) and analysed offline
afterwards. See README.md and docs/PROVENANCE.md for what is historical and what is new.
"""
M = 9                 # qubits = features
ALPHA = 0.1           # literal Eq. (7) coefficient of Duckett et al., PRD 109, 052002 (2024)
REPS = 2
PAIRS = [(l, m) for l in range(M) for m in range(l + 1, M)]   # 36 qubit pairs, lexicographic
__version__ = "0.1.0"
