# Extended Plan

Date: 2026-03-10

## Goal

Extend the current archival failed-supernova pipeline with additional search modes that can recover candidates missed by the pairwise fading scan.

These extensions are meant to answer a different question than the original discovery pass:

`What else could find real disappearing-star candidates that pairwise subtraction and threshold cuts might miss?`

## Current baseline

The existing pipeline already provides:

- nearby-galaxy selection
- HST/JWST archive coverage matrix
- matched epoch-pair discovery
- fading-source candidate ledger
- per-candidate follow-up photometry and SED products for the current `PASS` objects

The extensions below sit on top of that baseline rather than replacing it.

## Extension methods

### 1. All-epoch scene modeling

Idea:

- stop treating each candidate as only a pre/post pair
- model the source flux at a fixed sky position across every usable epoch
- use the full time series to test whether the fade is coherent and unique

Deliverables:

- per-candidate all-epoch scene photometry table
- scene-model variability metrics
- scene-model plots

Execution target:

- all current `PASS` and `REVIEW` candidates, prioritized by current ledger score

### 2. Deep-reference forced photometry

Idea:

- build a deep reference stack for a host galaxy
- detect the source list once on the deep image
- force photometry for that fixed catalog across all selected epochs

Deliverables:

- per-galaxy deep reference catalog
- all-epoch forced-photometry matrix for bright stars
- bright-star fade ranking table

Execution target:

- top host galaxies from the current candidate ledger, starting with `NGC 5861`, `MESSIER 101`, and `MESSIER 077`

### 3. Mid-IR survivor / dust branch

Idea:

- crossmatch optical/NIR fading candidates against mid-IR surveys and products
- flag candidates consistent with an obscured surviving star instead of a true disappearance

Deliverables:

- candidate-level WISE / NEOWISE / IRSA crossmatch table
- mid-IR survivor flags
- dust-survivor triage summary

Execution target:

- all current ledger candidates

### 4. Self-supervised anomaly ranking

Idea:

- build compact feature vectors from cutout sequences and candidate metrics
- rank objects that look unlike the normal candidate population
- use reconstruction error and distance-in-feature-space as additional novelty signals

Deliverables:

- anomaly feature table
- anomaly scores
- combined anomaly-ranked candidate list

Execution target:

- all current ledger candidates

### 5. Population-prior search

Idea:

- prioritize candidates that look astrophysically plausible as failed-supernova progenitors
- combine host-galaxy properties, local environment, fade strength, and discovery context into a population prior

Deliverables:

- candidate population-prior table
- updated combined ranking

Execution target:

- all current ledger candidates

### 6. Injection-recovery optimization

Idea:

- inject fake disappearing stars into representative real image pairs
- rerun the discovery logic
- measure what kinds of fades are recovered or missed

Deliverables:

- injection trial table
- recovery summary by flux ratio, brightness, and crowding
- recommended threshold adjustments

Execution target:

- representative high-value hosts and filters, starting with `NGC 5861`, `MESSIER 101`, and `MESSIER 077`

### 7. Local-neighborhood empirical systematics

Idea:

- compare each candidate against nearby stars with similar brightness and background
- down-rank candidates whose signal looks like a local systematic trend
- up-rank candidates that are locally unique

Deliverables:

- local-neighborhood comparison table
- local uniqueness / systematic-risk scores
- candidate-level neighborhood verdicts

Execution target:

- all current ledger candidates

## Execution rule

All seven methods should be implemented as reproducible pipeline stages under `src/supernova_pipeline/` and written to `extensions/`.

The first execution pass should:

- run all candidate-level stages on the current ledger
- run galaxy-level stages on the top host galaxies
- leave machine-readable summaries for every stage

## Success condition

The extension pass is successful if it produces:

- runnable code for all seven methods
- written outputs for each method under `extensions/`
- refreshed rankings or vetoes that materially change how candidates are triaged
