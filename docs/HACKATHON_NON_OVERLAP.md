# Relationship to our Qollab × IonQ hackathon entries

## Chronology

1. **2026-08-22 → 2026-08-23** (verified: sealed local checkpoints). The IBM quantum-kernel
   research (ER-01) was done and closed. Its closure lists "apply for research credits" as a
   future option.
2. **Research-credit application** (owner-attested). After initially missing hackathon
   registration, the team contacted Qollab, learned about the separate research-credit programme,
   and applied with the existing ER-01 IBM kernel research as "Quantum Kernel Reproducibility
   Lab". The application was submitted **before** the hackathon entries. Its exact date is not
   recorded in this repository and is not stated here.
3. **2026-10-09 23:00 UTC → 2026-10-11 23:00 UTC** (verified from the hackathon repositories'
   own documentation and commit dates, first commit 2026-10-10 19:02 UTC). Five separate
   hackathon entries were built and submitted.
4. **2026-10** (owner-attested). Qollab replied to the earlier research-credit application. This
   repository is part of the answer.

The research-credit proposal comes from ER-01. It does **not** come from any hackathon entry.

## Audit method (2026-10-14)

All five public repositories were cloned read-only at the HEAD commits below. They were searched
(file contents and commit messages) for:
* ER-01 identifiers, the IBM job ID, backend names, TrackML, Duckett et al., QSVM, "quantum kernel",
  "feature map", RYY
* the distinctive historical numbers (660, 1668, 382, 135, 108, 0.0068, 0.318, 20480)
* the SHA256 prefixes of the 9 input vectors and 5 executed circuits
* ER-01 function names (`_apply_yy`, `overlap_circuit`, `u_phi`, `kernel_matrix`,
  `parity_network`, `exact_simplifier`, `native_mapping`)

No code from those repositories was executed.

| repository | HEAD | content | code reuse | result reuse | same experiment | funding / credits claimed |
|---|---|---|---|---|---|---|
| quantum-impossible-9of9 | `28bd72af` | Mermin–Peres magic-square game, 4-qubit circuits | none | none | no | none; no IonQ hardware used |
| quantum-allocation-lab | `2b68919a` | 8-qubit XY-mixer QAOA allocation vs classical baselines, synthetic data | none | none | no | none; its "credits" are a synthetic unit |
| quantum-research-governor | `289a8c21` | [[4,2,2]] error-detection run governor | none | none | no | none |
| quantum-claim-court | `32d1357d` | Claims about quantum cloning judged against literature, mathematics and a 3-qubit cloner | none | none | no | none; hardware use is refused in code |
| quantum-alliance | `0cfdb570` | Four-qubit CHSH/GHZ game | none | none | no | none; hardware use is refused in code |

Search hits were all false positives and are disclosed for completeness:
* "grant" in MIT licence text and API names
* the numbers 108/135/382/660 in QAOA benchmark tables
* `0.0068359375` (= 28/4096) in a quantum-allocation-lab local Aer simulation result for an
  11-qubit QAOA circuit. It is numerically equal to the IBM P0 self-overlap but unrelated.
* an unrelated *PhysRevD* reference in a literature record

## Shared generic methodology (disclosed, not scientific overlap)

The hackathon entries and this repository use the same engineering habits:
* SHA256 manifests
* standalone verifiers that recompute exported evidence
* backend policy gates (one job per run, hardware refused unless authorised)
* claim ceilings
* provenance and AI-use disclosures

These are general methods for working with evidence. No verifier, gate or manifest code is shared.

## No double funding

| item | earlier funding / credits | used in a hackathon entry? | requested from this grant? |
|---|---|---|---|
| ER-01A TrackML reproduction (simulation only) | none (self-funded, USD 0) | no | no |
| ER-01B IBM job `da591ju1vhnc73flpn0g` (20,480 shots) | IBM free Open plan, 9 QPU-seconds, USD 0; no Qollab/IonQ resource | no | no. Reported as the baseline only |
| ER-01B2–B2.2 offline compilation analysis | none (USD 0) | no | no |
| TrackML-derived input vectors (9) | public CC-BY-4.0 data | no | no |
| Circuits (logical, IBM-compiled, IonQ-native payloads) | none | no | **the IonQ circuits are the subject of the request**. They have never been executed |
| This reproducibility package | none (self-funded) | no | no |
| Proposed IonQ runs (arms A/B/C/D, optionally E) | none yet | no | **yes, this is the only item requested** |
| Five hackathon entries | hackathon participation only; no IonQ hardware used, no credits claimed | (they are the entries) | no |

No experiment, dataset, circuit or result appears on more than one funding line. A future
Qollab-funded IonQ run would be new work and would not be reported as a hackathon result.
