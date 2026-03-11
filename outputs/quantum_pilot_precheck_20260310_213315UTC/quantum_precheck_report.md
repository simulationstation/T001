# Quantum Pilot Precheck Report

Created UTC: `2026-03-10T21:33:21.961590Z`

## Scope

This packet defines pilot-level dry-run designs for the six quantum tests in the current proposal.
Nothing in this packet was submitted to AWS Braket.

## Agent/Workflow Notes

- `/home/primary/Hubble-Systematics-Review-Chain` supplied the Braket runner pattern and pilot preflight model.
- `/home/primary/BAO_Ruler_Experiments` currently contains no `AGENT.md` or `AGENTS.md`, and no Braket submission code was found there.
- Credential source used for live device checks: `/home/primary/Hubble-Systematics-Review-Chain/.env`

## Live Device Inventory

- `IQM Garnet` in `eu-north-1`: status `ONLINE`; actions `braket.ir.openqasm.program, braket.ir.openqasm.program_set`
- `Rigetti Ankaa-3` in `us-west-1`: status `ONLINE`; actions `braket.ir.openqasm.program, braket.ir.openqasm.program_set`
- `QuEra Aquila` in `us-east-1`: status `ONLINE`; actions `braket.ir.ahs.program`

## Pilot Tests

### `core_shell_collapse_ladder`

- Question: Can the core stay organized while shell-mediated export fails sharply?
- Platform: Superconducting qubits
- Planned hardware footprint: `14` tasks at `1500` shots
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$34.65`

### `quantum_neutrinosphere_analog`

- Question: Can gross leakage survive while structured transport fails?
- Platform: Gate-model pilot on Rigetti as a low-cost proxy for the broader trapped-ion/cold-atom concept
- Planned hardware footprint: `10` tasks at `1200` shots
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$13.80`

### `shock_revival_threshold_simulator`

- Question: Does a small neutral-atom analog show a threshold-like revival split?
- Platform: QuEra Aquila analog Hamiltonian simulation
- Planned hardware footprint: `8` tasks at `100` shots
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$10.40`

### `shell_only_visibility_test`

- Question: Can shell readout geometry change visibility while the interior stays fixed?
- Platform: Superconducting qubits
- Planned hardware footprint: `8` tasks at `1200` shots
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$16.32`

### `cross_platform_universality_test`

- Question: Does the same reduced branch ordering package cleanly for two hardware families?
- Platform: IQM Garnet + Rigetti Ankaa-3
- Planned hardware footprint: `12` tasks total
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$17.70`

### `collapse_rate_atlas_microgrid`

- Question: Can a tiny pilot grid already separate more than one branch class?
- Platform: Gate-model microgrid on Rigetti
- Planned hardware footprint: `9` tasks at `800` shots
- Precheck verdict: `READY` (`preflight_pass`)
- Estimated pilot spend: `$9.18`

## Total Estimated Spend

Current total for all six pilot quantum tests: `$102.05`

Pricing assumptions were taken from the official AWS Braket pricing page on 2026-03-10:
https://aws.amazon.com/braket/pricing/

## Output Files

- [stage01_device_inventory.json](/home/primary/SUPERNOVA/outputs/quantum_pilot_precheck_20260310_213315UTC/stage01_device_inventory.json)
- [stage02_test_designs.json](/home/primary/SUPERNOVA/outputs/quantum_pilot_precheck_20260310_213315UTC/stage02_test_designs.json)
- [stage03_precheck_results.json](/home/primary/SUPERNOVA/outputs/quantum_pilot_precheck_20260310_213315UTC/stage03_precheck_results.json)
- [stage04_cost_summary.json](/home/primary/SUPERNOVA/outputs/quantum_pilot_precheck_20260310_213315UTC/stage04_cost_summary.json)
