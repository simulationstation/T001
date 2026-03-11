# Agent Instructions (SUPERNOVA)

This folder inherits the reusable runner and reporting rules from the neighboring
Braket-heavy workflow in `/home/primary/Hubble-Systematics-Review-Chain`.

`/home/primary/BAO_Ruler_Experiments` does not currently contain an `AGENT.md`
or `AGENTS.md`, so there were no folder-local agent rules to merge from that
project.

## Running jobs

- Keep runs reproducible, config-driven, and conservative about decision-grade claims.
- For long runs, create a dedicated output directory under `outputs/`.
- Write run metadata into each output directory:
  - `config.json`
  - `env.json`
  - `git_like.json`
  - `progress.json`
- For detached launches, use the established pattern:
  - write a `job.sh`
  - launch via `setsid` + `taskset`
  - persist `pid.txt`
  - stream logs to `run.log`
- Default to single-thread BLAS/OpenMP to avoid oversubscription:
  - `OMP_NUM_THREADS=1`
  - `MKL_NUM_THREADS=1`
  - `OPENBLAS_NUM_THREADS=1`
  - `NUMEXPR_NUM_THREADS=1`

## Monitoring policy

- Use micro-polling rather than coarse blocking waits.
- Surface results as soon as artifacts appear rather than waiting on long fixed sleeps.
- Prefer checking `progress.json`, stage files, logs, and output artifacts frequently.
- Maintain granular status tracking for any long or staged evaluation.

## Quantum experiment policy

- Always start with the smallest pilot likely to reveal the effect.
- Preflight locally before any hardware spend.
- Use the minimum qubit count, minimum circuit count, minimum cycle depth, and minimum task count needed to test whether the effect exists at all.
- Do not fan out into multi-task or batch hardware submissions until:
  - local preflight has passed
  - packaging has been validated
  - hard task caps have been checked
- Treat packaging and scale control as hard gates, not optimizations.
- Keep submission opt-in by default:
  - precheck and packaging can run automatically
  - actual queueing must remain explicit

## Execution priority

- Prioritize forward execution that produces new empirical outputs.
- Prefer the next measurable pilot or report over administrative workflow steps.
- Do not introduce freeze/tag/state-lock steps unless explicitly requested or required for a concrete run.

## Code style

- Prefer small, testable pure functions for core math and metric computation.
- Do not silently change defaults that affect inference or cost.
- Surface run-critical choices in configs or machine-readable outputs.
- Make staged pipelines resumable:
  - each stage writes checkpoint artifacts
  - reruns can skip completed stages

## Data acquisition and external systems

- For decision-grade work, new web downloads are allowed when they improve evidence quality or reproducibility.
- For decision-grade bulk data acquisition, use Globus connectors first when a public source exists.
- Always parallelize independent work when possible unless ordering constraints require serialization.

## Writing outputs

- All manuscript, README, report, and summary text must be reader-facing.
- Never include operator-directed meta language such as:
  - `as requested`
  - `you asked`
  - `per your message`
- Narrative outputs must stand alone for external readers without conversation context.
