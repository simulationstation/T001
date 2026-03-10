# Supernova Report

Date: 2026-03-10

## Executive Summary

Ryan is trying to turn a hard but clean question into a systematic survey:

`Which nearby massive stars disappear without producing a normal luminous supernova?`

The attached proposal is not a failed idea. It is a good science case that depends on three things being true at the same time:

- the galaxy sample must be chosen well,
- the multi-epoch archive coverage must be real and usable,
- the candidate triage must be rigorous enough to survive crowding, filter mismatch, and dust-survivor impostors.

My main conclusion is:

`The archival-only version is still worth doing, but it needs a different operating model than the proposal text implies.`

The best route is a `negative-space disappearing-star engine`:

- build the nearby-galaxy universe first,
- build a full archive-coverage matrix second,
- use catalog-level and HLSP-level screening before expensive reductions,
- then run a DESI-style candidate ledger with reason codes, packets, and robust ranking.

## What The Proposal Says

The proposal is a Cycle 5 + Cycle 9 JWST/HST program aimed at a definitive nearby-galaxy failed-SN search.

Headline numbers from the PDF:

- about `130` nearby star-forming galaxies
- about `92` with previous JWST imaging
- about `58` with previous HST WFC3/IR imaging
- expected by the full program: about `6` failed SNe and about `36` SN progenitors

The important nuance is that the proposal is strongest when new epochs are added. Ryan explicitly says the existing observations are still interesting, but maybe not fully constraining. That is consistent with the proposal itself.

## What Exists Locally Right Now

`/home/primary/SUPERNOVA` contains only the proposal PDF. There is no supernova-specific codebase yet.

The useful local infrastructure is elsewhere.

### Strong local assets

- Processed `GLADE+` index already exists in:
  `../Dark_Siren_Ladder_Audit/data/processed/galaxies/gladeplus/`
- Astroquery cache exists and is large:
  `~/.astropy/cache/astroquery` is about `4.1G`
- Existing MAST download cache exists:
  `~/.lightkurve/cache/mastDownload` is about `9.9G`
- Globus CLI is live:
  `~/.venv_globuscli/bin/globus whoami` returned the logged-in account

### Strong local methodology to reuse

The best transferable patterns are in `../DESI-BH-CANDIDATE-SEARCH/`.

What is valuable there:

- robust anomaly scores
- leave-one-out diagnostics
- high-leverage / fragile-candidate flags
- PASS / REVIEW / FAIL reason-code ledgers
- candidate packets and summary markdown
- crossmatch-driven veto logic

The key conceptual transfer is this:

`Do not only ask "is the signal large?" Ask "can a boring failure mode explain it?"`

That is exactly the right mindset for disappearing stars.

### Environment gap

The current inspected Python environments do not contain the real HST/JWST reduction stack you need for Ryan steps 4-8.

I found no local installs of:

- `photutils`
- `jwst`
- `drizzlepac`
- `space_phot`
- `hst123`

So the machine is ready for catalog work, archive discovery, and data engineering, but not yet for serious HST/JWST reduction.

## What I Think Probably Went Wrong Before

Not a single fatal flaw. More likely a stack of practical ones.

Most likely issues:

- Archive heterogeneity was underestimated.
- The pipeline leaned too heavily on full reduction before triage.
- `hst123` is useful but too fragile to be the only reduction path.
- There was no strong candidate ledger with negative-evidence accounting.
- The proposal-scale science case depended on future epochs, while archival-only work has patchy leverage.

That means the fix is not "work harder on the same pipeline."

The fix is:

- better galaxy ranking,
- better archive ranking,
- earlier use of catalog-level information,
- much stricter candidate packaging and veto logic.

## The New Strategy

### 1. Build the galaxy universe like a survey problem

Use multiple catalogs together:

- NED-LVS
- GLADE+
- 50MGC
- secondary arbitration from HyperLEDA / EDD when needed

Then score galaxies by actual searchability:

- distance
- inclination
- extinction
- surface brightness
- crowding
- archival epoch availability
- filter compatibility

### 2. Build the archive universe like a database problem

Before downloading pixels, build a full HST/JWST observation matrix:

- science-only products
- footprint overlap
- epoch baseline
- filter comparability
- depth
- HLSP availability

This is where the current project should become machine-like instead of notebook-like.

### 3. Build the candidate universe like a forensic anomaly problem

This is the DESI lesson.

Each candidate needs:

- a robust fade score
- leave-one-out stability
- crowding and blend penalties
- subtraction-quality flags
- a dusty-survivor test
- a packet with cutouts and provenance

I would explicitly maintain both:

- a `high-confidence queue`
- a `discovery salvage queue`

That second queue matters because the best failed-SN candidates may be visually dramatic but formally fragile.

## Datasets That Matter Most

### Highest priority

- MAST HST archive
- MAST JWST archive
- HSCv3 and HCV
- PHANGS-HST
- PHANGS-JWST
- LEGUS
- Spitzer Heritage Archive

### Best veto / context layers

- ZTF
- Gaia DR3
- SIMBAD
- NED
- NEOWISE / unWISE
- HEASARC high-energy archives

### Why this matters

The science is not just "disappearance."

The real question is whether the star:

- truly vanished,
- is hidden by dust,
- merged,
- is a variable in a confused region,
- or is an artifact of bad subtraction or incomplete depth.

That means mid-IR, time-domain optical, and environment/context data are part of the core search, not optional extras.

## Globus Assessment

Globus is configured and logged in locally, which is useful for large external pulls.

A quick endpoint check found obvious value for:

- NOIRLab
- SDSS
- unWISE / related large catalog movement

I did not find an obvious dedicated STScI / MAST endpoint from a quick search, so I would not assume Globus solves the HST/JWST side. For MAST, I would plan around the official APIs, archive downloads, and cloud-hosted public data options.

## Bottom Line

Ryan's mission is still viable.

But the version that has a real chance to work is not:

- "query a catalog, reduce everything, hope a few stars disappear."

It is:

- "assemble the nearby-galaxy universe,"
- "score it,"
- "assemble the archive universe,"
- "score it,"
- "run a negative-space candidate engine,"
- "reduce only where the evidence justifies it,"
- "package every survivor like a forensic case."

That is the novel path I would take.

## Execution Status

The plan has now been executed as a live pilot in this workspace.

Artifacts created by the runnable pipeline include:

- `catalogs/galaxy_master.parquet`
- `archive/observation_matrix.parquet`
- `archive/epoch_pairs.parquet`
- `archive/pixel_pair_queue.parquet`
- `candidates/pixel_pair_summary.parquet`
- `candidates/source_measurements.parquet`
- `candidates/candidate_ledger.parquet`
- `candidate_packets/`

Current live results from the tightened pixel search:

- `11` exact-filter pilot pairs selected
- `9` scanned successfully
- `2` failed for known reasons
- `1` surviving `REVIEW` candidate
- `0` surviving `PASS` candidates

Important pair-level findings:

- `MESSIER 061` `JWST/JWST F2100W` is photometrically stable enough to trust at the pair level, but yielded no strong fading source.
- `MESSIER 077` `HST/HST WFC3 F555W` is the cleanest HST archival pair so far and produced one surviving `REVIEW` object.
- Several other archival pairs that initially looked promising failed the bulk-stability test, exactly the failure mode Ryan warned about when "most of the stars" do not remain constant.

Best surviving candidate at the moment:

- Galaxy: `MESSIER 077`
- Pair: `hst_15151_04_wfc3_uvis_f555w_idm104` to `hst_17502_01_wfc3_uvis_f555w_ifab01`
- Baseline: about `1904.61` days
- Position: `RA 40.6696134 deg`, `Dec -0.0132900 deg`
- Fade fraction: about `0.565`
- Fade significance: about `15.16 sigma`
- Status: `REVIEW`
- Warnings: `CROWDED`, `BLEND_RISK`

That means the current answer is:

- archival-only discovery is possible,
- most naive candidates are not trustworthy,
- pair-stability gating is essential,
- and the only object that currently survives those cuts is interesting enough to inspect manually but not yet strong enough to call a failed-supernova candidate.

## Most Important Next Actions

1. Create the dedicated supernova environment and install the missing HST/JWST reduction stack.
2. Build the merged nearby-galaxy master table from NED-LVS + GLADE+ + 50MGC.
3. Build the HST/JWST archive matrix and matched epoch-pair table.
4. Pilot the workflow on `5-10` galaxies before scaling.
5. Implement the candidate ledger and packet format before heavy reduction.

## Sources

- NED Local Volume Sample: https://ned.ipac.caltech.edu/NED::LVS/
- GLADE+: https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.2503D/abstract
- NASA 50 Mpc Galaxy Catalog: https://heasarc.gsfc.nasa.gov/W3Browse/all/nasa50mwc.html
- MAST missions via astroquery: https://astroquery.readthedocs.io/en/latest/mast/mast_missions.html
- Hubble Source Catalog v3: https://archive.stsci.edu/hst/hsc/
- Hubble Catalog of Variables: https://archive.stsci.edu/hst/hsc/help/HCV/
- PHANGS-HST: https://archive.stsci.edu/hlsp/phangs-hst
- PHANGS-JWST: https://archive.stsci.edu/hlsp/phangs-jwst
- LEGUS: https://archive.stsci.edu/hlsp/legus
- JWST Cycle 5 APT form information: https://www.stsci.edu/jwst/science-execution/jwst-apt-form-information/jwst-apt-cycle-5-form-information
- Spitzer Heritage Archive note: https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzerhelpspitzer/dataanalysistools/spitzerheritagearchive/
- ZTF public data at IRSA: https://irsa.ipac.caltech.edu/Missions/ztf.html
- NEOWISE single-exposure database: https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/expsup/
