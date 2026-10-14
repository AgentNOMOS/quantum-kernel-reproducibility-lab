# Evidence (public release)

The raw IBM job output and provider metadata are withheld: see `../WITHHELD_EVIDENCE.md`.

`historical_derived/` holds the historical records that contain only our own exact or offline
computations. They are kept for comparison; the verifier recomputes the values independently:

* `local_reference.json`: exact kernel values frozen before the IBM run
* `finite_shot_reference.json`: shot-noise-only reference frozen before the IBM run
* `er01b21_reconstruct_216.json`: derivation of the 216 decomposition cost
* `er01b22_native_mapping.json`: the original 108-ZZ native mapping record
