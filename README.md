# SUPERNOVA

Archival failed-supernova / disappearing-star search workspace built around Ryan's JWST proposal and extended into a live HST/JWST discovery pipeline.

## What this repo contains

- The science framing in the local, gitignored JWST proposal PDF plus `super_plan.md` and `super_report.md`
- A runnable Python package in `src/supernova_pipeline/`
- Derived catalog, archive, candidate, follow-up, and strict benchmark products from the current live run
- Branch-aware export-failure ranking outputs plus operator docs in `predictaroni/`
- Visual review artifacts for every currently surviving candidate

## Current status

The project is past the original archive-mining pilot and into the stricter detector-validation phase Ryan actually cares about.

- Nearby-galaxy master catalog: `12,419` search-eligible systems
- Expanded archive run: `130` galaxies queried, `3,930` observations, `120,995` matched epoch pairs
- Historical candidate ledger from the earlier live search: `18` surviving objects, `3` `PASS`, `15` `REVIEW`
- Per-candidate packets live in `candidate_packets/`; gallery assets live in `figures/candidates/`
- Completed follow-up products now exist for all `3` current `PASS` objects in `sed/`
- Strict blinded benchmark after detector fixes: `12/12` truth-centered benchmark pairs scanned, `9/12` known supernova truth groups recovered, recovery fraction `0.75`, median localization error `0.68 arcsec`
- Strict live rerun complete in `outputs/live_difference_strict_20260311_0421/`: `60` pairs completed, `49` scanned, `11` failed, `1201` fade detections, `65` raw `PASS`, `1136` raw `REVIEW`
- Branch-aware rescoring complete in `outputs/live_difference_strict_20260311_0421/branch_aware_v4/`: `901` same-location clusters, `30` `STRONG_EXPORT_FAILURE_LIKE`, `58` `INTERMEDIATE_BRANCH_LIKE`, `21` `DUST_SURVIVOR_LIKE`, `792` `VARIABLE_OR_UNRESOLVED`

The main state change is that the detector is no longer just producing interesting candidates. It now passes a nontrivial blinded supernova recovery test, completes a strict live rerun, and supports a second-stage branch-aware interpretation layer built on same-location forced measurements. The gallery below is still historical/provisional, but the repo now has a stricter authoritative ranking path than it did before.

## Ryan acceptance tests

Ryan's feedback narrowed the real acceptance criteria. The project is only scientifically mature if it can do these things:

1. detect on difference images rather than relying on source catalogs alone
2. register images empirically from stars instead of trusting nominal WCS alignment
3. do forced photometry at the same residual location across all epochs, including non-detections
4. carry explicit provenance for chip gaps, masks, low-S/N, and no-coverage cases
5. recover known targeted supernovae in a blinded benchmark before trusting disappearing-star claims

## Progress against those tests

- `Difference-first detection`: implemented in the stricter `difference_upgrade` path rather than only the older pre-seeded candidate workflow.
- `Empirical registration`: partially improved, but this is still the main weak point. The detector is much better than the original WCS-only version, but image registration remains the first thing to harden further.
- `Forced photometry on residual locations`: now part of the stricter benchmark / residual-measurement path and the branch-aware rescoring stage.
- `Coverage and provenance tracking`: pair preflight, overlap checks, restricted-rights skips, and structured checkpoint outputs now exist.
- `Blind recovery test`: passed in the useful sense of being operational and nontrivial. The corrected benchmark in `outputs/benchmark_validation_rawresid_20260311_0416/` recovers `9` of `12` truth groups.

## What the results currently mean

- The old candidate set is still useful as a discovery ledger and follow-up target list.
- The closure layer currently classifies `3` objects as `EXPORT_FAILURE_LIKE`, `1` as `DUST_SURVIVOR_LIKE`, `1` as `SYSTEMATIC_LIKE`, and `13` as `VARIABLE_OR_UNRESOLVED`.
- The strict branch-aware layer is now the preferred interpretation path for the current live rerun, not the old gallery-first ledger.
- The current branch-aware result is intentionally selective: `30` strong export-like clusters, `58` intermediate branch-like clusters, `21` dust-like competitors, and the rest unresolved rather than overclaimed.
- The project is still not at a point where the gallery below should be treated as publication-grade failed-supernova claims.

## What must be done next

1. Treat `outputs/live_difference_strict_20260311_0421/branch_aware_v4/` as the current strict ranking baseline, not the old gallery-first ledger.
2. Improve empirical star-based registration to reduce subtraction dipoles and residual astrometric artifacts.
3. Compare the historical `PASS` / `REVIEW` objects against the branch-aware strict clusters instead of carrying them forward unchallenged.
4. Use known targeted supernova recovery as the standing regression test for every future detector or ranking change.
5. Do deeper all-epoch follow-up only on the top branch-ranked strict clusters, especially the strongest `MESSIER 061` and `MESSIER 066` sites.

## Current interpretation of the best objects

The strongest historical `PASS` objects are two long-baseline `HST/WFC3-IR F160W` candidates in `NGC 5861` and one very long-baseline `HST/ACS F555W` candidate in `MESSIER 101`. Based on direct expert feedback, `NGC 5861` likely contains at least one real supernova-like transient, while `MESSIER 101` remains highly suspect for artifact / subtraction issues until it survives the strict rerun.

## Branch-Aware Strict Outputs

The new strict interpretation outputs live in `outputs/live_difference_strict_20260311_0421/branch_aware_v4/`.

Most useful artifacts:
- `branch_cluster_scores.parquet`
- `branch_detection_scores.parquet`
- `forced_measurements.parquet`
- `branch_summary.json`
- `branch_report.md`
- `packets/`

These are the outputs to use if you want the repository's current best attempt at ranking survivors by export-failure coherence under detector controls rather than by raw fade amplitude alone.

## Hidden web board

The repo now also contains a deployable hidden review page for the current strict branch-aware candidates:

- Source payload: `site/uber_nova_search/`
- Build command: `supernova-mission build-uber-nova-page --branch-output-dir outputs/live_difference_strict_20260311_0421/branch_aware_v4 --output-dir site/uber_nova_search`
- Intended deploy target: `/uber_nova_search` under `/var/www/quasardipolephenomenon/`

This board is not meant to be linked from the homepage. It is a Ryan-facing review surface that shows the current non-unresolved branch-aware candidates with slow registered blink GIFs, triptych panels, filter context, and score decomposition.
It now includes client-side filters for:

- `Matched to blind`: same-galaxy association to a benchmark truth position within the page's operator radius
- `Highest Value`: the top-quality branch-aware survivors under stricter support-depth and systematic-risk cuts

## Prime `PASS` objects

These are the top surviving objects from the earlier cleaned ledger. They are kept here because they remain the main follow-up targets, but they should now be read as provisional historical candidates pending the strict live rerun described above. Each one is shown with the same visual treatment as the earlier single-object example: animated blink, overview panel, and packet cutout. Where follow-up has finished, the light curve and SED fit are also included.

### `NGC 5861` `PASS` object A

Candidate: `NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023`. Filter pair: `F160W -> F160W`. Baseline: `1494.5` days. Fade: `97.49%`. Significance: `111.45 sigma`. Priority: `91.48`.

![NGC 5861 PASS A blink](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_blink.gif)

![NGC 5861 PASS A overview](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_overview.png)

![NGC 5861 PASS A packet](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_packet.png)

![NGC 5861 PASS A light curve](sed/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023/lightcurve.png)

![NGC 5861 PASS A SED fit](sed/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023/sed_fit.png)

### `NGC 5861` `PASS` object B

Candidate: `NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025`. Filter pair: `F160W -> F160W`. Baseline: `1456.6` days. Fade: `97.06%`. Significance: `82.81 sigma`. Priority: `71.27`.

![NGC 5861 PASS B blink](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_blink.gif)

![NGC 5861 PASS B overview](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_overview.png)

![NGC 5861 PASS B packet](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_packet.png)

![NGC 5861 PASS B light curve](sed/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025/lightcurve.png)

![NGC 5861 PASS B SED fit](sed/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025/sed_fit.png)

### `MESSIER 101` `PASS` object

Candidate: `MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `99.99%`. Significance: `12.62 sigma`. Priority: `21.59`.

![MESSIER 101 PASS blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_blink.gif)

![MESSIER 101 PASS overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_overview.png)

![MESSIER 101 PASS packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_packet.png)

![MESSIER 101 PASS light curve](sed/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003/lightcurve.png)

![MESSIER 101 PASS SED fit](sed/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003/sed_fit.png)

## Candidate gallery

Every surviving ledger entry is shown below with the same three artifacts: animated blink, overview panel, and packet cutout.

### PASS candidates

<details>
<summary><code>NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023</code> | PASS | NGC 5861 | F160W -> F160W | 1494.5 days | fade 97.49% | 111.45 sigma | priority 91.48</summary>

Galaxy: `NGC 5861`. Filter pair: `F160W -> F160W`. Baseline: `1494.5` days. Fade: `97.49%`. Significance: `111.45 sigma`. Priority: `91.48`.

![NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023 blink](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_blink.gif)

![NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023 overview](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_overview.png)

![NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023 packet](figures/candidates/NGC_5861_hst_15145_74_wfc3_ir_f160w_idgg74_hst_16250_14_wfc3_ir_f160w_iebs14_0023_packet.png)
</details>

<details>
<summary><code>NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025</code> | PASS | NGC 5861 | F160W -> F160W | 1456.6 days | fade 97.06% | 82.81 sigma | priority 71.27</summary>

Galaxy: `NGC 5861`. Filter pair: `F160W -> F160W`. Baseline: `1456.6` days. Fade: `97.06%`. Significance: `82.81 sigma`. Priority: `71.27`.

![NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025 blink](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_blink.gif)

![NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025 overview](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_overview.png)

![NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025 packet](figures/candidates/NGC_5861_hst_15145_78_wfc3_ir_f160w_idgg78_hst_16250_14_wfc3_ir_f160w_iebs14_0025_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003</code> | PASS | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 99.99% | 12.62 sigma | priority 21.59</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `99.99%`. Significance: `12.62 sigma`. Priority: `21.59`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0003_packet.png)
</details>

### REVIEW candidates

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 67.15% | 65.65 sigma | priority 61.10</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `67.15%`. Significance: `65.65 sigma`. Priority: `61.10`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1397_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 90.88% | 52.14 sigma | priority 57.37</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `90.88%`. Significance: `52.14 sigma`. Priority: `57.37`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1314_packet.png)
</details>

<details>
<summary><code>MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257</code> | REVIEW | MESSIER 077 | F814W -> F814W | 1904.6 days | fade 64.53% | 69.69 sigma | priority 45.23</summary>

Galaxy: `MESSIER 077`. Filter pair: `F814W -> F814W`. Baseline: `1904.6` days. Fade: `64.53%`. Significance: `69.69 sigma`. Priority: `45.23`.

![MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257 blink](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257_blink.gif)

![MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257 overview](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257_overview.png)

![MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257 packet](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f814w_idm104_hst_17502_01_wfc3_uvis_f814w_ifab01_0257_packet.png)
</details>

<details>
<summary><code>MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217</code> | REVIEW | MESSIER 066 | F438W -> F438W | 2203.7 days | fade 54.32% | 58.02 sigma | priority 21.90</summary>

Galaxy: `MESSIER 066`. Filter pair: `F438W -> F438W`. Baseline: `2203.7` days. Fade: `54.32%`. Significance: `58.02 sigma`. Priority: `21.90`.

![MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217 blink](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217_blink.gif)

![MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217 overview](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217_overview.png)

![MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217 packet](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f438w_ic0h01_hst_15654_22_wfc3_uvis_f438w_idxr22_0217_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 68.11% | 45.82 sigma | priority 19.27</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `68.11%`. Significance: `45.82 sigma`. Priority: `19.27`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0449_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 87.09% | 35.36 sigma | priority 15.28</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `87.09%`. Significance: `35.36 sigma`. Priority: `15.28`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0225_packet.png)
</details>

<details>
<summary><code>MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091</code> | REVIEW | MESSIER 066 | F814W -> F814W | 2203.6 days | fade 68.01% | 33.14 sigma | priority 13.37</summary>

Galaxy: `MESSIER 066`. Filter pair: `F814W -> F814W`. Baseline: `2203.6` days. Fade: `68.01%`. Significance: `33.14 sigma`. Priority: `13.37`.

![MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091 blink](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091_blink.gif)

![MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091 overview](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091_overview.png)

![MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091 packet](figures/candidates/MESSIER_066_hst_12968_01_wfc3_uvis_f814w_ic0h01_hst_15654_22_wfc3_uvis_f814w_idxr22_0091_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 100.12% | 983.42 sigma | priority 7.50</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `100.12%`. Significance: `983.42 sigma`. Priority: `7.50`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_0008_packet.png)
</details>

<details>
<summary><code>6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061</code> | REVIEW | 6dF J0245599-073443 | F814W -> F814W | 2139.9 days | fade 59.91% | 10.53 sigma | priority 7.44</summary>

Galaxy: `6dF J0245599-073443`. Filter pair: `F814W -> F814W`. Baseline: `2139.9` days. Fade: `59.91%`. Significance: `10.53 sigma`. Priority: `7.44`.

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061 blink](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061_blink.gif)

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061 overview](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061_overview.png)

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061 packet](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0061_packet.png)
</details>

<details>
<summary><code>6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002</code> | REVIEW | 6dF J0245599-073443 | F814W -> F814W | 2139.9 days | fade 58.29% | 11.36 sigma | priority 5.72</summary>

Galaxy: `6dF J0245599-073443`. Filter pair: `F814W -> F814W`. Baseline: `2139.9` days. Fade: `58.29%`. Significance: `11.36 sigma`. Priority: `5.72`.

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002 blink](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002_blink.gif)

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002 overview](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002_overview.png)

![6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002 packet](figures/candidates/6dF_J0245599-073443_hst_11575_02_acs_wfc_f814w_jb4u02_hst_14226_01_acs_wfc_f814w_jcuc01_0002_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 117.68% | 11.97 sigma | priority 5.42</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `117.68%`. Significance: `11.97 sigma`. Priority: `5.42`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1266_packet.png)
</details>

<details>
<summary><code>MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482</code> | REVIEW | MESSIER 101 | F555W -> F555W | 5328.0 days | fade 66.55% | 7.55 sigma | priority 2.41</summary>

Galaxy: `MESSIER 101`. Filter pair: `F555W -> F555W`. Baseline: `5328.0` days. Fade: `66.55%`. Significance: `7.55 sigma`. Priority: `2.41`.

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482 blink](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482_blink.gif)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482 overview](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482_overview.png)

![MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482 packet](figures/candidates/MESSIER_101_hst_9490_02_acs_wfc_f555w_j8d602_hst_14678_01_acs_wfc_f555w_jd6v01_1482_packet.png)
</details>

<details>
<summary><code>MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138</code> | REVIEW | MESSIER 077 | F555W -> F555W | 1904.6 days | fade 56.57% | 15.17 sigma | priority 2.37</summary>

Galaxy: `MESSIER 077`. Filter pair: `F555W -> F555W`. Baseline: `1904.6` days. Fade: `56.57%`. Significance: `15.17 sigma`. Priority: `2.37`.

![MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138 blink](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138_blink.gif)

![MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138 overview](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138_overview.png)

![MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138 packet](figures/candidates/MESSIER_077_hst_15151_04_wfc3_uvis_f555w_idm104_hst_17502_01_wfc3_uvis_f555w_ifab01_0138_packet.png)
</details>

<details>
<summary><code>MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033</code> | REVIEW | MESSIER 061 | F814W -> F814W | 3368.4 days | fade 64.34% | 5.88 sigma | priority 2.08</summary>

Galaxy: `MESSIER 061`. Filter pair: `F814W -> F814W`. Baseline: `3368.4` days. Fade: `64.34%`. Significance: `5.88 sigma`. Priority: `2.08`.

![MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033 blink](figures/candidates/MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033_blink.gif)

![MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033 overview](figures/candidates/MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033_overview.png)

![MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033 packet](figures/candidates/MESSIER_061_hst_skycell-p1458x13y04_acs_wfc_f814w_all_hst_12574_01_acs_wfc_f814w_jbs801_0033_packet.png)
</details>

<details>
<summary><code>NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137</code> | REVIEW | NGC 5728 | F160W -> F160W | 104.5 days | fade 43.78% | 5.38 sigma | priority 0.03</summary>

Galaxy: `NGC 5728`. Filter pair: `F160W -> F160W`. Baseline: `104.5` days. Fade: `43.78%`. Significance: `5.38 sigma`. Priority: `0.03`.

![NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137 blink](figures/candidates/NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137_blink.gif)

![NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137 overview](figures/candidates/NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137_overview.png)

![NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137 packet](figures/candidates/NGC_5728_hst_15640_21_wfc3_ir_f160w_idxb21_hst_15640_2b_wfc3_ir_f160w_idxb2b_0137_packet.png)
</details>

## Reproduce the pipeline

```bash
python3 -m venv .venv
./.venv/bin/pip install -e .
./.venv/bin/supernova-mission run-pilot --top-n 130 --archive-top-n 25
```

For the current expanded discovery stage:

```bash
./.venv/bin/supernova-mission build-archive \
  --archive-top-n 130 \
  --radius-deg 0.05 \
  --min-baseline-days 30

./.venv/bin/supernova-mission run-pixel-search \
  --epoch-pairs-path archive/epoch_pairs.parquet \
  --observation-matrix-path archive/observation_matrix.parquet \
  --max-pairs 113 \
  --per-galaxy 3 \
  --compatibilities exact \
  --same-collection-only \
  --max-workers 4

./.venv/bin/supernova-mission run-followup \
  --candidate-ledger-path candidates/candidate_ledger.parquet \
  --observation-matrix-path archive/observation_matrix.parquet \
  --galaxy-master-path catalogs/galaxy_master.parquet \
  --statuses PASS \
  --max-workers 4

PYTHONPATH=src .venv/bin/python -m supernova_pipeline.cli run-branch-scoring \
  --strict-output-dir outputs/live_difference_strict_20260311_0421 \
  --benchmark-dir outputs/benchmark_validation_rawresid_20260311_0416 \
  --precomputed-dir outputs/live_difference_strict_20260311_0421/branch_aware_v2 \
  --output-dir outputs/live_difference_strict_20260311_0421/branch_aware_v4
```

## Key outputs

- `catalogs/galaxy_master.parquet`
- `catalogs/galaxy_scores.parquet`
- `catalogs/galaxy_exclusions.csv.gz`
- `archive/observation_matrix.parquet`
- `archive/epoch_pairs.parquet`
- `archive/archive_failures.csv`
- `candidates/candidate_ledger.parquet`
- `candidates/review_queue.csv`
- `candidate_packets/`
- `figures/candidates/`
- `sed/`
- `outputs/live_difference_strict_20260311_0421/branch_aware_v4/`
- `predictaroni/`

## Notes

- The raw `data/cache/` download cache is intentionally not versioned because it is large and mostly reproducible from MAST.
- `catalogs/galaxy_exclusions.csv` is preserved in compressed form as `catalogs/galaxy_exclusions.csv.gz` to stay within normal GitHub file-size limits.
- The ledger is the authoritative current candidate list for the gallery above; follow-up products are still being written by the active `run-followup` stage.
