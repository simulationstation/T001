# Next Actions

This is the immediate owner-facing action list for the branch-aware export-failure upgrade.

## 1. Make the strict live rerun plus branch-aware scoring the new working baseline

Use the strict detector outputs, not the old gallery-first ledger, as the authoritative input for future ranking work.

Immediate artifacts to carry forward:
- strict live rerun in `outputs/live_difference_strict_20260311_0421/`
- strict benchmark in `outputs/benchmark_validation_rawresid_20260311_0416/`
- branch-aware outputs in `outputs/live_difference_strict_20260311_0421/branch_aware_v4/`

## 2. Harden empirical registration further

This remains the main detector weakness.

Priority tasks:
- tighten pair rejection or penalties for poor residual alignment
- characterize centroid migration more explicitly
- add targeted dipole regression examples to the test battery

## 3. Re-score and prune historical candidates under the strict branch-aware layer

The old historical gallery is still useful, but it is no longer the epistemic center of the project.

Required:
- compare historical `PASS` / `REVIEW` objects against strict branch-aware cluster rankings
- explicitly note which old objects survive, weaken, or flip to systematic/dust competitors

## 4. Automate the benchmark guardrail

The blind benchmark should become a standing precondition for trusting ranking changes.

Required:
- run benchmark regression whenever strict detector behavior or ranking logic changes materially
- record benchmark recovery fraction in branch-aware summaries
- block “serious lead” escalation when benchmark performance regresses

## 5. Improve provenance around same-location forced measurements

Current branch-aware logic should expose typed measurement outcomes.

Still useful next steps:
- separate chip-gap / footprint / masked-pixel causes more finely when possible
- preserve failed-pair context at cluster level
- carry explicit coverage fractions into every human-readable packet

## 6. Add stronger follow-up links for the top branch-ranked clusters

The branch-aware layer is triage, not the last word.

For the top clusters:
- run deeper all-epoch same-location follow-up
- add multi-band interpretation where coverage exists
- carry SED and dust-survivor competitors explicitly into the packet narrative

## 7. Define graduation rules clearly

A candidate should only become a serious lead when all of the following are true:
- benchmark guardrail passes
- branch class is `STRONG_EXPORT_FAILURE_LIKE` or a very strong `INTERMEDIATE_BRANCH_LIKE`
- registration confidence is acceptable
- artifact penalty is modest
- same-location forced measurements remain coherent
- dust-survivor competition is not dominating
- coverage completeness is not severely limited

## 8. Keep the right project posture

The upgraded mentality should remain:
- detector realism first
- candidate interpretation second
- publication claims last

That means the correct next step after any exciting candidate is usually not announcement. It is stricter measurement and stronger competitor elimination.
