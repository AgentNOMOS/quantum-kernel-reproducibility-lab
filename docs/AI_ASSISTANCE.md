# AI-assistance disclosure

**Historical research (ER-01A … ER-01 closure, August 2026).** This work was developed under the
project owner's direction with ChatGPT and AI-assisted coding workflows (owner's statement,
2026-10-14). The sealed research records do not say which tool produced which file. The owner
cannot say from memory which historical files were written independently and which may have
been adapted from public papers, examples or repositories. We therefore make no blanket claim of
independent authorship.

**This reproducibility package (October 2026).** The code and documentation in this repository
were produced with Claude Code (Anthropic Claude) under the owner's direction. They adapt the
historical code listed in `PUBLICATION_RIGHTS_MATRIX.md`.

How we handle the uncertainty, instead of asserting ownership:

1. **Provenance evidence per file.** `PUBLICATION_RIGHTS_MATRIX.md` classifies every file by
   origin and lists the evidence.
2. **Upstream similarity audit** (`verifier/upstream_similarity.py`, rerunnable). Every Python
   file here, and each historical file it adapts, was compared token by token with the Python sources of the libraries available locally: Qiskit,
   qiskit-ionq, qiskit-qasm3-import, qiskit-ibm-runtime, qiskit-aer, NumPy, SciPy,
   scikit-learn and pytest (3,240 files).
   - Settings: exact 25-token shingles, plus 45-token shingles with identifiers renamed.
   - Positive controls were detected: a verbatim upstream file 100 %, an identifier-renamed copy
     100 %, an embedded upstream class 43 %.
   - Result: **0 shared sequences in every release file.** The historical ancestors had two
     trivial hits, a literal Hadamard matrix and a run of tuple punctuation.
   - The audit **cannot** rule out adaptation from sources that are not installed locally (blogs,
     tutorials, other repositories). That residual risk is stated here and is not resolved.
3. **No historical code headers, licence notices or upstream URLs** were found in any of the 244
   historical Python files. The only URLs are API endpoints and the Zenodo dataset link.
4. Every numeric claim in `README.md` and `docs/` is produced or checked by code in this
   repository, except for the hardware values, whose raw evidence is withheld
   (`WITHHELD_EVIDENCE.md`).
5. Disagreements found while assembling this package were corrected and disclosed. For example,
   a sealed file named `best_G_*` holds the 496-CZ stock-transpiler circuits, not the 382-CZ ones.

Responsibility for the scientific claims rests with the human project owner.
