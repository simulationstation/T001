# Observational Closure Test

## Goal
Adjudicate every disappearing-star candidate against four competing explanations using the already-built SUPERNOVA products.

## Competing Classes
- `EXPORT_FAILURE_LIKE`: strong permanent fade with low systematic and survivor evidence.
- `DUST_SURVIVOR_LIKE`: fade is accompanied by mid-IR or dusty survivor evidence.
- `SYSTEMATIC_LIKE`: crowding, blending, local coherence, or reduction disagreement dominate.
- `VARIABLE_OR_UNRESOLVED`: still interesting, but current evidence does not break the tie.

## Input Products
- candidate ledger
- scene-model metrics
- neighborhood/systematics metrics
- mid-IR crossmatch
- all-epoch follow-up summaries and photometry

## Scoring Logic
Each candidate gets:
- `fade_evidence_score`: how extreme and repeatable the fade looks in pairwise and scene-level products.
- `uniqueness_evidence_score`: whether the source is locally unique instead of following field-wide behavior.
- `followup_support_score`: whether the all-epoch packet supports a persistent multi-filter disappearance.
- `systematic_penalty_score`: how strongly crowding, blending, disagreement, residuals, or local trends argue for an artifact.
- `survivor_evidence_score`: whether mid-IR or dusty colors support an obscured surviving star.

The exported hypothesis score is a weighted combination of the first three terms minus the last two. Dust and systematic scores are computed separately so the test can force a competing explanation when warranted.

## Decision Rules
1. `DUST_SURVIVOR_LIKE` if survivor evidence is dominant.
2. `SYSTEMATIC_LIKE` if systematic evidence is dominant.
3. `EXPORT_FAILURE_LIKE` if export evidence is high and both competing penalties are low.
4. Otherwise `VARIABLE_OR_UNRESOLVED`.

## Closure Success
The test succeeds if it either:
- isolates candidates that strongly prefer the export-failure interpretation, or
- kills the sample cleanly enough to define an upper limit on that channel.
