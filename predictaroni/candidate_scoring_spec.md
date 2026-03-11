# Candidate Scoring Specification

This document defines the branch-aware ranking layer that sits on top of the strict `difference_upgrade` detector. It is intentionally factorized so a maintainer can audit how a candidate rose or fell.

## Scope

The scoring layer does not replace the detector. It interprets survivors from a strict difference-image search under an assumption-conditioned export-failure / branch-selection framework.

The detector remains responsible for:
- empirical registration
- residual-first detection
- pair-level artifact rejection
- provenance capture
- blind benchmark recovery

The scoring layer is responsible for:
- clustering same-location residuals into candidate sites
- forced same-location residual measurements across scanned pairs
- factorized interpretation scores
- class labels and rank ordering

## Input Artifacts

Expected strict-run inputs:
- `pair_queue.parquet`
- `pair_summary.parquet`
- `fade_candidates.parquet`
- `detection_measurements.parquet`
- benchmark summary from a validated strict benchmark run

Expected optional context:
- historical closure outputs
- extension outputs
- all-epoch follow-up products

## Candidate Unit

The preferred unit for branch-aware interpretation is not a single pairwise fade row. It is a same-location cluster of residual detections within one host galaxy.

Why:
- a single pair can be compelling but still artifact-prone
- export-failure coherence is a multi-epoch / multi-pair concept
- same-position forced measurements across pairs are more informative than one pairwise fade significance

## Score Components

All score components are expected to live on `[0, 1]` and to be derived from transparent heuristics, not hidden models.

### `detector_confidence`

Meaning:
- how convincing the underlying strict detections are before higher-level interpretation

Typical inputs:
- strongest `diff_sigma_min` or similar robust residual significance
- number of independent support pairs
- fraction of member detections that were already `PASS`
- pre-image S/N
- strong-fade member fraction

Expected range:
- `0.0-0.3`: weak or single-pair support
- `0.3-0.6`: credible but limited support
- `0.6-1.0`: strong strict-detector support

### `registration_confidence`

Meaning:
- confidence that the candidate survives empirical alignment quality checks

Typical inputs:
- median registration residual
- number of anchor matches
- cluster centroid scatter across member detections
- consistency of registration quality across forced-measurement pairs

Strong penalties:
- large residuals
- low anchor counts
- member centroid migration across pairs

### `subtraction_cleanliness`

Meaning:
- how clean the residual behavior looks once registration is done

Typical inputs:
- pair stability fraction
- pair extreme-fade background fraction
- warning flags such as blend/crowding/edge
- dispersion of same-location forced ratios across measured pairs

Strong penalties:
- field-wide instability
- dipole-like behavior
- edge/mask adjacency
- unstable residual sign

### `coverage_completeness`

Meaning:
- how much usable pair-level evidence exists for the candidate site

Typical inputs:
- measured-pair fraction among scanned pairs in the host
- filter diversity
- baseline span
- galaxy-level pair failure fraction

Interpretation:
- low coverage should lower confidence, not force a negative class

### `forced_photometry_coherence`

Meaning:
- how consistent the same-location residual measurements are across scanned pairs

Typical inputs:
- fraction of measured pairs showing positive fade-consistent residuals
- median post/pre ratio at the fixed location
- spread of post/pre ratios
- number of measured pairs
- absence of obvious return signatures

This is the core term that moves the search from “interesting fade” to “coherent export-failure candidate under control measurements.”

### `fade_persistence`

Meaning:
- whether the candidate stays suppressed over long baselines rather than appearing as a one-off transient

Typical inputs:
- severe fade fraction
- fade fraction in the latest available measured pairs
- baseline span
- low late-time return fraction

### `rebrightening_penalty`

Meaning:
- penalty for evidence that the source comes back or changes sign in a way more consistent with ordinary variability or transient return

Typical inputs:
- fraction of measured pairs with high post/pre ratio
- fraction of measured pairs with brightening-sign residuals
- late-time return fraction

### `artifact_penalty`

Meaning:
- direct summary of anti-self-deception evidence

Typical inputs:
- warning flags
- blend risk
- crowding
- edge ambiguity
- registration weakness
- galaxy-level pair failure fraction
- cluster centroid scatter

Interpretation:
- high artifact penalty should suppress rank strongly even if fade amplitude is large

### `dust_survivor_score`

Meaning:
- how strongly a candidate looks like an obscured survivor instead of an export-failure collapse

Typical inputs:
- strong optical suppression with weak or absent IR suppression when IR coverage exists
- rebrightening or partial return behavior
- lack of multi-band export support

Important limit:
- if IR coverage does not exist, this score must stay conservative rather than pretending to know the answer

### `export_failure_score`

Meaning:
- top-level score for strong branch-consistent export-failure interpretation

Should reward:
- strong strict-detector evidence
- clean registration and subtraction
- same-location forced coherence
- persistent suppression
- adequate coverage
- multi-band support when available

Should penalize:
- rebrightening
- artifact risk
- dust-survivor evidence

### `systematic_risk`

Meaning:
- probability-like artifact concern summary, not a detector significance

Typical inputs:
- artifact penalty
- registration weakness
- subtraction weakness
- low detector confidence

### `unresolved_variability_score`

Meaning:
- score for “not enough to claim export failure, not enough to force dust/systematic either”

Typical inputs:
- incomplete coverage
- weak forced coherence
- weak persistence
- return behavior that does not look decisively systematic

## Intermediate Branch Logic

The ranking layer must not collapse everything into “disappeared” versus “did not disappear.”

An intermediate branch-like candidate is expected to show some combination of:
- repeated same-location suppression
- partial rather than total loss
- heterogeneous suppression across filters or baselines
- low artifact burden
- limited but nontrivial return evidence

This is captured by `partial_branch_score`, which should reward:
- partial-suppression fraction
- coherence
- coverage
- registration quality

And should penalize:
- strong artifact indicators
- strong dust-survivor indicators
- strong rebrightening

## Aggregation Logic

Recommended aggregation pattern:
1. Build same-location clusters from strict detections.
2. Force-measure those cluster centers across all scanned pairs in the same galaxy.
3. Compute raw cluster support metrics.
4. Normalize raw metrics conservatively.
5. Compute factor scores.
6. Compute top-level interpretation scores.
7. Classify.
8. Rank.

Recommended class logic:
- `STRONG_EXPORT_FAILURE_LIKE`:
  requires high export-failure score plus adequate registration, coherence, and low artifact burden
- `INTERMEDIATE_BRANCH_LIKE`:
  requires moderate partial-branch score plus acceptable coherence and modest artifact burden
- `DUST_SURVIVOR_LIKE`:
  requires dust score to beat export interpretation by a real margin
- `SYSTEMATIC_LIKE`:
  requires systematic score to dominate
- `VARIABLE_OR_UNRESOLVED`:
  fallback when coverage or coherence is insufficient

## Benchmark Constraints

The blind known-supernova benchmark constrains tuning.

Required rule:
- do not accept a scoring change as operationally good if the strict benchmark regresses materially

Minimum expected checks:
- recovery fraction remains stable or improves
- localization error does not worsen materially
- candidate rank inflation does not come from weaker provenance or weaker artifact penalties

Recommended reporting:
- carry benchmark recovery fraction into each scoring run summary
- indicate whether the guardrail is currently passing
- do not mark branch-ranked leads as escalation-ready if benchmark performance is below threshold

## Interpretation of Final Ranking

`branch_rank_score` is not a claim of astrophysical truth.

It is a triage score meaning:
- how well the candidate fits the export-failure / branch-aware interpretation
- after strict detector controls and provenance penalties are applied

Recommended use:
- top branch-ranked clusters become follow-up targets
- they do not become publication claims automatically

## Required Output Columns

At minimum, the scored cluster table should carry:
- `cluster_id`
- `branch_rank`
- `branch_rank_score`
- `branch_class`
- `branch_confidence`
- `detector_confidence`
- `registration_confidence`
- `subtraction_cleanliness`
- `coverage_completeness`
- `forced_photometry_coherence`
- `fade_persistence`
- `rebrightening_penalty`
- `artifact_penalty`
- `dust_survivor_score`
- `partial_branch_score`
- `export_failure_score`
- `systematic_risk`
- `unresolved_variability_score`
- reason codes
- provenance counts

## Audit Questions

When reviewing a high-ranked candidate, ask:
1. Is the location coherent across more than one pair?
2. Is the registration actually good?
3. Are the forced same-location measurements consistent?
4. Is the site well covered or merely dramatic in one pair?
5. Could dust explain the wavelength pattern?
6. Is the artifact burden low enough to trust the ranking?
7. Is the benchmark currently healthy enough that this run should be trusted at all?
