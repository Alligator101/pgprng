# Changelog

All notable changes to pgprng are recorded here. This project doesn't
yet use git (that's the natural next step), so this file is the
version record in the meantime — each entry says what changed and,
where relevant, how it was verified not to change generator output.

## 1.1 — 2026-09-14

**Object-model refactor — internal representation only, no output change.**

- Added a `State` class (`pgprng_common.py`): each ensemble
  component's fixed `increment` and current running `state` are now a
  single object instead of two separate parallel-list entries.
- `ensemble_1`/`ensemble_2` replace the old `incs_11`/`states_11` and
  `incs_13`/`states_13` naming in `FilteredSelectionPRNG` and
  `SelectionCombinedPRNG`. Each is a plain list of pointers to `State`
  objects. `RotatingEnsemble` (used by `CombinedPRNG`) got the same
  treatment with a single `self.ensemble` list.
- Removed cached ensemble-size attributes (`self.e11`, `self.e13`,
  `self.e`) entirely. Ensemble size is now always read live via
  `len(ensemble)`, since ensemble sizes are meant to be mutable rather
  than baked into the code — which the old `_11`/`_13` naming did
  literally, since it hadn't matched the real ensemble sizes (101/103)
  since the project outgrew its original e=11/e=13 toy scale.
- Constructor parameter names updated to match: `incs_1`/`incs_2`/
  `seeds_1`/`seeds_2` (was `incs_11`/`incs_13`/`seeds_11`/`seeds_13`).
  `CombinedPRNG`'s two `RotatingEnsemble` attributes renamed
  `ens1`/`ens2` (was `ens11`/`ens13`).
- Verified behavior-preserving: all three generator classes
  (`FilteredSelectionPRNG`, `CombinedPRNG`, `SelectionCombinedPRNG`)
  produce byte-for-bit identical output against a fixed-seed baseline
  captured before this refactor.
- `pgprng_ensemble_construction.py`'s validation functions
  (`min_pairwise_gap`, `uniformity_chi2`, `smoke_test_ensemble`)
  deliberately left unchanged — they validate increment-set properties
  before any `State` objects exist, so the object-model change doesn't
  apply there.
- Corrected terminology throughout code, comments, and documentation:
  the per-component values were previously called "multiplier(s)" in
  places, which was wrong — they're used purely additively
  (`state = state + increment mod 2^64`), with no multiplicative term
  anywhere in that recurrence. `pgprng_multipliers_101_103.csv` was
  renamed to `pgprng_increments_101_103.csv` (header `ensemble,m` →
  `ensemble,increment`). Two genuinely multiplicative concepts were
  kept and clarified rather than removed: `mix64`'s actual multiply
  operations, and the "multiplying periods together" description of
  how the two ensembles' periods combine in `CombinedPRNG`.

## 1.0 — 2026-09-13

Initial public release.

- Two generator designs: `CombinedPRNG` (clock-stepping, exact period
  `2^64 × e1 × e2`) and `FilteredSelectionPRNG` (mutual selection with
  mix-then-combine — the recommended design, closing a trivial
  inversion flaw present in the earlier combine-then-mix ordering).
- Hand-built NIST SP 800-22-style statistical battery (12 tests, 500
  sequences × 1,000,000 bits each) — see `final_report.md`.
- Real PractRand 0.86 run against `FilteredSelectionPRNG`'s actual
  output stream: 8 GiB streamed clean, zero anomalies.
- Toy-scale cryptanalytic brute-force attempt (`pgprng_crack_attempt.py`)
  demonstrating the specific flaw the mix-then-combine redesign closes.
- MIT License.
