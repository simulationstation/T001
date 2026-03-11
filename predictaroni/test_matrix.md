# Branch-Aware Test Matrix

This test battery is designed for the strict detector plus branch-aware interpretation layer. The detector must earn trust before the interpretation layer is allowed to matter.

## Status Categories

- `PASS`: acceptable for scientific use
- `WARNING`: usable with explicit caution
- `FAIL`: blocks escalation or invalidates the run

## 1. Detector Core Tests

### 1.1 Difference-first detection

Purpose:
- ensure surviving candidates come from residual-image detections, not only source-list priors

Pass:
- detector outputs residual detections from difference images directly

Warning:
- auxiliary source-list seeds are used only as a fallback with explicit labeling

Fail:
- pipeline silently depends on source catalogs for primary survivor creation

### 1.2 Empirical registration diagnostics

Purpose:
- ensure star-based registration is actually happening and recorded

Pass:
- pair outputs contain shift, residual, and anchor-count diagnostics
- survivors are penalized or rejected when registration residual is poor

Warning:
- diagnostics exist but are not yet used downstream

Fail:
- detector relies on WCS reprojection only

### 1.3 Subtraction overlap / coverage provenance

Purpose:
- make no-coverage and overlap failure explicit

Pass:
- pair-level outputs distinguish scanned, failed, skipped, and preflight-rejected cases

Fail:
- missing detections can be confused with missing coverage

## 2. Registration-Specific Tests

### 2.1 Low-anchor failure case

Inject or identify pairs with very few anchor stars.

Pass:
- registration confidence drops or the pair is explicitly downgraded

Fail:
- low-anchor pairs rank candidates normally

### 2.2 Large residual / dipole case

Use a pair known to produce dipole-like residuals.

Pass:
- artifact penalty and systematic risk increase
- candidate class trends toward `SYSTEMATIC_LIKE` or `VARIABLE_OR_UNRESOLVED`

Fail:
- candidate ranks highly on fade significance alone

## 3. Subtraction / Residual Tests

### 3.1 Edge and mask ambiguity

Purpose:
- verify edge and mask adjacency are treated as provenance, not silence

Pass:
- same-location forced rows produce explicit `OUTSIDE_FOOTPRINT` or `MASKED_OR_EDGE`
- coverage completeness declines

Fail:
- candidate is treated as a genuine non-detection without provenance penalty

### 3.2 Pair instability

Purpose:
- reject residual fields where many ordinary stars look unstable

Pass:
- pair stability fraction affects subtraction cleanliness or artifact penalty

Fail:
- field instability is recorded but ignored in ranking

## 4. Forced-Photometry Coherence Tests

### 4.1 Same-location persistence

Purpose:
- verify cluster centers are force-measured across scanned pairs

Pass:
- output table exists with one row per cluster per scanned pair
- measured, outside-footprint, and masked states are distinct

Fail:
- forced photometry is only done for pre-selected handpicked candidates

### 4.2 Non-detection support

Purpose:
- ensure absence of signal is represented numerically rather than dropped

Pass:
- low-S/N or near-zero residual rows remain in the forced-measurement table

Fail:
- only high-significance measurements are retained

### 4.3 Rebrightening detection

Purpose:
- ensure returns or sign flips are captured and penalized

Pass:
- rebrightening penalty rises when late-time ratios become large or residual sign flips

Fail:
- return behavior is invisible to ranking

## 5. Branch-Aware Interpretation Tests

### 5.1 Strong export-like archetype

Construct or identify a cluster with:
- multiple support pairs
- good registration
- persistent suppression
- low artifact burden

Pass:
- class is `STRONG_EXPORT_FAILURE_LIKE`

### 5.2 Intermediate branch-like archetype

Construct or identify a cluster with:
- repeated same-location suppression
- partial rather than total loss
- modest artifact burden

Pass:
- class is `INTERMEDIATE_BRANCH_LIKE`

Fail:
- candidate is forced into either “strong export” or “unresolved” with no intermediate route

### 5.3 Dust competitor

Use a cluster with strong short-wavelength suppression but weak long-wavelength suppression when IR coverage exists.

Pass:
- dust-survivor score rises materially
- class can flip to `DUST_SURVIVOR_LIKE`

### 5.4 Systematic competitor

Use a cluster with centroid migration, warning flags, poor registration, or edge ambiguity.

Pass:
- systematic risk dominates

## 6. Benchmark Regression Tests

### 6.1 Standing blind benchmark

Purpose:
- detector trust anchor

Pass:
- known-supernova blind recovery remains stable or improves versus the current validated run

Current validated reference:
- `12/12` benchmark pairs scanned
- `9/12` truth groups recovered
- recovery fraction `0.75`
- median localization error about `0.68 arcsec`

Fail:
- benchmark recovery regresses materially without an explicit accepted reason

### 6.2 Benchmark isolation

Purpose:
- detector must not read hidden truth during scoring or detection

Pass:
- benchmark truth remains in a separate tree and is consumed only by evaluation

Fail:
- hidden truth enters the detector configuration or ranking logic

## 7. Output Integrity Tests

### 7.1 Factorized explanation output

Pass:
- every scored cluster has explicit factor columns and reason codes

Fail:
- only a single rank number is emitted

### 7.2 Packet completeness

Pass:
- each cluster packet has machine-readable summary plus human-readable decision text

### 7.3 Provenance completeness

Pass:
- outputs distinguish:
  - measured
  - outside footprint
  - masked or edge-limited
  - pair prep failed

Fail:
- missingness is not typed

## 8. Promotion Rules

A candidate can graduate from “interesting” to “serious lead” only if:
- benchmark guardrail is currently passing
- branch class is export-like or intermediate branch-like
- registration confidence is acceptable
- artifact penalty is not high
- coverage completeness is not severely limited
- same-location forced measurements support persistence rather than return

Any failure in those conditions is at least a `WARNING`, and some are absolute `FAIL`s.
