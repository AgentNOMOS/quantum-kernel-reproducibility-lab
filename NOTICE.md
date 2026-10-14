# NOTICE — licences and rights

**Licence scope.** The MIT licence in `LICENSE` (Copyright (c) 2026 Murat Keskin, AgentNOMOS
Research Lab) applies only to the project files marked "MIT" in `PUBLICATION_RIGHTS_MATRIX.md`,
and only to the extent of the rights the licensor holds in them. Those files were produced with
AI assistance (`docs/AI_ASSISTANCE.md`). MIT does **not** override or relicense third-party
material. In particular, the TrackML-derived input vectors remain under CC BY 4.0 (below).

The complete file-by-file classification, with provenance evidence, is in
`PUBLICATION_RIGHTS_MATRIX.md`. This notice summarises it.

| Material | Location | Licence / status |
|---|---|---|
| Source code, verifier, tests, experiment builders | `qkrl/`, `verifier/`, `tests/`, `experiments/*.py`, `reference/*.py` | MIT (see `LICENSE`) |
| Documentation | `README.md`, `docs/` | MIT |
| Circuits we generated (logical, compiled, native payloads) | `circuits/`, `experiments/ionq_native_payloads/` | MIT. These encode feature values derived from TrackML data, so the CC-BY-4.0 attribution below also applies |
| Input feature vectors (9 vectors, 9 features each) | `data/kernel_pair_manifest.json` | **CC-BY-4.0**. Derived from the TrackML Particle Tracking Challenge dataset, Zenodo DOI [10.5281/zenodo.4730167](https://doi.org/10.5281/zenodo.4730167), © the TrackML challenge organisers (Salzburger, Innocente, Vlimant, Rousseau, Gligorov, Estrade, Basara, Calafiura, Farell, Gray, Golling, Kiehn, Amrouche, Ustyuzhanin, Hushchyn, Moyse, Germain, Guyon). Changes: hit triplets were turned into 9 features, MinMax-scaled to [0,1], and 9 vectors were selected by a frozen rule. No raw TrackML file is redistributed |
| Raw output of IBM job `da591ju1vhnc73flpn0g` and provider metadata (coupling map, calibration values, job/usage records) | **not included** (`WITHHELD_EVIDENCE.md`) | **Withheld.** The redistribution terms could not be established from public IBM documentation. Never covered by MIT |
| Summary results of that job (5 all-zero counts and statistics derived from them) | `README.md`, `docs/` | Reported as the authors' research findings, with IBM's recommended acknowledgement. Not relicensed as data |
| Circuits we compiled for `ibm_marrakesh` (contain physical qubit indices) | `circuits/historical/ibm_marrakesh_compiled_*.qasm3`, `circuits/historical/offline_optimized/heavyhex_*` | Our compiler output, produced with open-source Qiskit. MIT. They are not provider files |
| Feature-map definition | `qkrl/feature_map.py` docstring | A scientific method from Duckett et al., *Phys. Rev. D* **109**, 052002 (2024), arXiv:2212.07279. Cited, not copied. The paper's PDF and text are **not** redistributed |

Runtime dependencies are installed from PyPI and are **not** vendored: qiskit (Apache-2.0),
qiskit-ionq (Apache-2.0), qiskit-qasm3-import (Apache-2.0), openqasm3 (Apache-2.0), rustworkx
(Apache-2.0), numpy (BSD-3-Clause and others), scipy (BSD-3-Clause), pytest (MIT).
`qkrl/ionq_compile.py` rebuilds the IonQ Forte native `Target` in the same way qiskit-ionq's
`IonQBackend._make_target()` does, using qiskit-ionq's public gate classes. No qiskit-ionq source
code is copied.
