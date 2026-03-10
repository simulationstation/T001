# Supernova Mission Plan

Date: 2026-03-10
Primary inputs:
- Ryan's email in the task prompt
- `JWST_Failed_SNe_Cycle_5.pdf`
- Local `/home/primary` project inventory

## Mission

Build a reproducible search for failed supernovae and other extreme fading massive stars in nearby galaxies using existing HST/JWST data first, then package the best candidates for follow-up and paper-ready analysis.

## What Ryan Is Actually Asking For

Ryan's message is a nine-step pipeline:

1. Build the nearby-galaxy target list.
2. Find real HST/JWST science images covering those galaxies.
3. Match multi-epoch observations in compatible filters with useful time baselines.
4. Reduce the data and build matched photometric catalogs.
5. Do subtraction / difference imaging.
6. Find stars that faded strongly.
7. Pull every other HST/JWST epoch for those positions.
8. Fit pre/post SEDs with extinction.
9. Produce cutouts, dates, positions, and candidate packets.

The proposal PDF adds the critical scientific frame:

- The goal is not just "find a transient." It is to determine which massive stars explode and which directly collapse.
- The intended sample is large: about 130 nearby star-forming galaxies.
- The proposal expected archival leverage from existing JWST and HST imaging, but also depended on new epochs to become definitive.

## Why The Straightforward Version Can Stall

The original plan is scientifically sound, but operationally fragile.

Main failure modes:

- Too much manual target vetting.
- Too much full image reduction before ranking candidates.
- HST/JWST filter, pointing, and camera heterogeneity.
- Crowding and host-background systematics dominating formal photometric errors.
- Finicky reduction wrappers (`hst123`) with little protection against bad runs.
- No explicit "blend / hideability / subtraction-failure" model.
- No strong candidate ledger with reason codes and re-review packets.
- Archival-only searches can be interesting but not always constraining if the baselines or overlap are poor.

## Core Idea: A Negative-Space Failed-SN Engine

Do not treat this as a pure image-reduction project.

Treat it as a three-layer search engine:

1. `Galaxy layer`
   Build the best possible nearby-galaxy universe and score each galaxy by detectability, not just distance.
2. `Coverage layer`
   Build a machine-readable matrix of all useful HST/JWST/Spitzer epochs and filter compatibility before downloading or reducing everything.
3. `Candidate layer`
   Use a DESI-style anomaly ledger for disappearing stars: robust score, leave-one-out stability, crowding penalties, and PASS/REVIEW/FAIL reason codes.

This is the main novel change from the proposal workflow.

Instead of reducing everything and then seeing what is interesting, we first rank:

- which galaxies are worth the effort,
- which footprints actually support a disappearance search,
- which stars are plausible candidates,
- which candidates survive negative-evidence checks.

## Data Universe

### A. Galaxy catalogs

Use multiple catalogs, not one.

Primary:

- `NED Local Volume Sample (NED-LVS)` for a vetted nearby-galaxy starting set.
- `GLADE+` for broad completeness and cross-identifiers.
- `NASA 50 Mpc Galaxy Catalog (50MGC)` for a modern local-universe galaxy census to 50 Mpc.

Secondary:

- HyperLEDA and EDD as adjudication layers when distances, morphologies, or photometry disagree.
- The Updated Nearby Galaxy Catalog for very local systems.

Hard exclusions:

- M31
- Magellanic Clouds
- Any galaxy where crowding makes a clean disappearing-star search unrealistic for the available data

### B. Space-based imaging and catalogs

Primary:

- MAST HST archive
- MAST JWST archive
- Hubble Source Catalog / HSCv3
- Hubble Catalog of Variables / HCV

High-value nearby-galaxy HLSPs:

- PHANGS-HST
- PHANGS-JWST
- LEGUS

Important complementary imaging:

- Spitzer Heritage Archive / IRSA

### C. Time-domain and veto data

Use these to reject impostors and classify survivors:

- ZTF public data at IRSA
- ATLAS
- ASAS-SN
- TESS / TESSCut where useful
- NEOWISE / unWISE time-domain
- Gaia DR3 + SIMBAD + NED object context
- HEASARC / Chandra / XMM / Swift for AGN or X-ray source context
- Radio context where available

## Target-Selection Model

Build a `galaxy_master.parquet` table and score galaxies by a `searchability_score`.

The score should combine:

- distance
- galaxy inclination
- Galactic extinction
- angular size
- local star-formation proxy
- host surface brightness
- existing HST/JWST epoch count
- overlap area between epochs
- availability of similar filters
- crowding proxy

This should replace simple hard cuts as the first ordering.

Recommended outputs:

- `catalogs/galaxy_master.parquet`
- `catalogs/galaxy_aliases.parquet`
- `catalogs/galaxy_scores.parquet`
- `catalogs/galaxy_exclusions.csv`

## Archive-Matrix Stage

Build a full observation ledger before heavy downloads.

For each galaxy:

- query MAST metadata for HST and JWST science products,
- reject calibrations and engineering products,
- compute footprint overlap with the galaxy,
- group by instrument, detector, filter, date, and depth,
- identify same-filter or transformable-filter epoch pairs,
- record baseline in days.

Also query:

- HSC/HCV coverage
- PHANGS/LEGUS/Spitzer products where present

Recommended outputs:

- `archive/observation_matrix.parquet`
- `archive/epoch_pairs.parquet`
- `archive/coverage_footprints/`
- `archive/filter_compatibility.csv`

## Candidate-Discovery Stage

Run discovery in two passes.

### Pass 1: Cheap prefilter

Use existing catalogs and high-level products to find likely interesting locations before custom reduction.

Examples:

- HSC/HCV stars with strong fading or disappearance signatures
- sources present in deep pre-imaging and absent in post stacks
- PHANGS / LEGUS stars with strong epoch mismatch
- locations where post-image depth is comfortably below pre-image flux

### Pass 2: Targeted image-level work

Only for shortlisted epoch pairs and candidate regions:

- align images
- perform forced photometry
- build difference images
- run fake-star disappearance injections to calibrate significance

## The DESI Transfer: Candidate Ledger Instead Of Loose Notes

Borrow the best pattern from `DESI-BH-CANDIDATE-SEARCH`.

Each candidate should get:

- a `status`: `PASS`, `REVIEW`, `FAIL`
- a `priority_score`
- a flat list of reason codes
- a flat list of warning flags
- a packet with cutouts and provenance

Recommended metrics:

- `fade_sigma`
- `fade_fraction`
- `S_raw`
- `S_min_LOO`
- `S_robust`
- `d_max`
- `depth_margin_post`
- `crowding_index`
- `blend_risk`
- `host_bg_penalty`
- `astrometric_residual`
- `cross_reduction_agreement`

Recommended reason codes:

- `single_epoch_only`
- `fragile_after_loo`
- `high_crowding`
- `possible_blend_sink`
- `filter_mismatch`
- `bad_overlap`
- `dust_survivor_possible`
- `known_variable_type`
- `possible_cluster_confusion`
- `poor_subtraction_quality`
- `insufficient_post_depth`

Recommended outputs:

- `candidates/candidate_ledger.parquet`
- `candidates/candidate_summary.csv`
- `candidates/review_queue.csv`
- `candidate_packets/<candidate_id>/`

## Reduction Strategy

Do not make `hst123` the only path.

Preferred order:

1. Catalog-first / HLSP-first screening
2. Targeted custom photometry for candidates
3. Full reduction only when justified

Reduction stack:

- `astroquery` + MAST APIs for discovery and metadata
- `photutils`, `regions`, `reproject`, `astropy`
- `drizzlepac` for HST where needed
- `jwst` / `stcal` for JWST where needed
- `space_phot` as a serious alternative path
- `hst123` only as one reduction option, not the pipeline truth

Every key candidate should, when possible, be checked with two independent reductions.

## SED And Physical Inference

For each surviving candidate:

- gather every HST/JWST/Spitzer epoch at that position,
- convert photometry to flux units,
- fit pre/post SEDs,
- include Milky Way extinction and broad internal-extinction priors,
- explicitly allow non-detections,
- test whether a dusty survivor fits better than true disappearance.

Outputs:

- `sed/<candidate_id>/photometry.csv`
- `sed/<candidate_id>/fit_summary.json`
- `sed/<candidate_id>/pre_post_sed.png`

## Local Assets To Reuse Immediately

### Already on disk

- Processed `GLADE+` HEALPix indices:
  `../Dark_Siren_Ladder_Audit/data/processed/galaxies/gladeplus/`
- Large astroquery cache:
  `~/.astropy/cache/astroquery`
- Existing MAST download cache:
  `~/.lightkurve/cache/mastDownload`
- Active Globus CLI login:
  `~/.venv_globuscli/bin/globus`

### Reusable local code patterns

- `../Dark_Siren_Ladder_Audit/scripts/fetch_gladeplus.py`
- `../Dark_Siren_Ladder_Audit/scripts/build_gladeplus_index.py`
- `../DESI-BH-CANDIDATE-SEARCH/scripts/master_verifier.py`
- `../DESI-BH-CANDIDATE-SEARCH/build_priority_packet.py`
- `../DESI-BH-CANDIDATE-SEARCH/scripts/exhaustive_archival_search.py`
- `../DESI-BH-CANDIDATE-SEARCH/scripts/ztf_long_baseline.py`
- `../DESI-BH-CANDIDATE-SEARCH/verify_candidates.py`

### Current environment gaps

Across the inspected environments, `astroquery` exists in some venvs, but I did not find local installs of:

- `photutils`
- `jwst`
- `drizzlepac`
- `space_phot`
- `hst123`

That is a real blocker for Ryan steps 4-8.

## Suggested Project Structure

```text
SUPERNOVA/
  JWST_Failed_SNe_Cycle_5.pdf
  super_plan.md
  super_report.md
  catalogs/
  archive/
  candidates/
  candidate_packets/
  sed/
  figures/
  notes/
```

## First Two Weeks

### Week 1

- Create dedicated supernova venv.
- Build merged nearby-galaxy master table from NED-LVS + GLADE+ + 50MGC.
- Exclude M31, Magellanic Clouds, and obviously hopeless systems.
- Build archive metadata matrix for the top-ranked galaxies.
- Produce first table of matched epoch pairs.

### Week 2

- Run catalog-first screening with HSC/HCV and HLSP products.
- Choose 5-10 pilot galaxies.
- Build first difference-imaging / forced-photometry prototypes.
- Define final candidate metrics and reason codes.
- Produce first candidate packet template.

## Risks And Mitigations

| Risk | Why it matters | Mitigation |
|---|---|---|
| JWST/HST heterogeneity | Different filters, footprints, and resolutions | Precompute compatibility matrix before reduction |
| Crowding | False disappearances or hidden survivors | Mandatory crowding/blend model plus injection tests |
| Correlated subtraction noise | Formal errors will be too optimistic | Empirical local-noise and fake-source calibration |
| Dusty survivor / merger confusion | Could mimic a failed SN | Use Spitzer/NEOWISE/JWST MIR plus SED fits |
| Data volume | Full archive pulls are expensive | Metadata-first ranking, targeted downloads, cache reuse |
| Pipeline fragility | `hst123` can fail badly | Maintain two reduction paths |
| Incomplete public access | Some data may still be in EAP | Work public-first, maintain watchlist for later releases |

## Definition Of Success

This project is working if it produces:

- a defensible ranked nearby-galaxy search universe,
- a full machine-readable HST/JWST epoch-pair matrix,
- a candidate ledger with reproducible scores and reason codes,
- packetized disappearing-star candidates,
- pre/post SED fits with extinction,
- enough audit trail that a paper cannot drift away from the pipeline outputs.

## External Sources

- NED Local Volume Sample: https://ned.ipac.caltech.edu/NED::LVS/
- GLADE+: https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.2503D/abstract
- NASA 50 Mpc Galaxy Catalog: https://heasarc.gsfc.nasa.gov/W3Browse/all/nasa50mwc.html
- MAST missions via astroquery: https://astroquery.readthedocs.io/en/latest/mast/mast_missions.html
- MAST cloud-hosted public datasets: https://github.com/spacetelescope/mast_notebooks/blob/main/notebooks/Cloud_access_to_public_astronomy_data/Cloud_access_to_public_astronomy_data.ipynb
- Hubble Source Catalog v3: https://archive.stsci.edu/hst/hsc/
- Hubble Catalog of Variables: https://archive.stsci.edu/hst/hsc/help/HCV/
- PHANGS-HST: https://archive.stsci.edu/hlsp/phangs-hst
- PHANGS-JWST: https://archive.stsci.edu/hlsp/phangs-jwst
- LEGUS: https://archive.stsci.edu/hlsp/legus
- JWST Cycle 5 APT and exclusive access notes: https://www.stsci.edu/jwst/science-execution/jwst-apt-form-information/jwst-apt-cycle-5-form-information
- Spitzer Heritage Archive note: https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzerhelpspitzer/dataanalysistools/spitzerheritagearchive/
- ZTF public data at IRSA: https://irsa.ipac.caltech.edu/Missions/ztf.html
- NEOWISE single-exposure database: https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/expsup/
