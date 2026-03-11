# Quantum Pilot Result Review Update

Created UTC: `2026-03-11T02:54:43.151112+00:00`

## Live Status
- `collapse_rate_atlas_microgrid`: `COMPLETED` x `9`
- `core_shell_collapse_ladder`: `COMPLETED` x `7`
- `cross_platform_universality_test_iqm`: `COMPLETED` x `3`
- `cross_platform_universality_test_rigetti`: `COMPLETED` x `3`
- `quantum_neutrinosphere_analog`: `COMPLETED` x `10`
- `shell_only_visibility_test`: `COMPLETED` x `8`
- `shock_revival_threshold_simulator`: `COMPLETED` x `8`

## Shock Revival Update
- All `8/8` QuEra Aquila tasks completed with `100/100` successful shots per task.
- Edge occupancy spans `0.37` to `0.58`, while center occupancy spans `0.66` to `0.90`.
- Strongest suppression occurred at `scale=1.00`, `detuning=2.0e+06`, where edge occupancy fell to `0.37` and edge-loss proxy rose to `0.63`.
- Weakest suppression occurred at `scale=0.85`, `detuning=1.0e+06`, where edge occupancy was `0.58`.
- Interpretation: completed cleanly with full shot success; boundary-sensitive edge suppression is present, but the pilot does not yet show a crisp two-branch bifurcation.

## Prior Pilot Headlines
- `core_shell_collapse_ladder`: shield and open stayed strongly positive in `ZZ` (`0.927`, `0.927`), while wrong flipped negative (`-0.905`).
- `quantum_neutrinosphere_analog`: matched runs kept higher inner `ZZ` than mismatched (`0.882` vs `0.806`) and lower outer excitation (`0.132` vs `0.329`); blocked runs suppressed outer excitation to `0.014`.
- `shell_only_visibility_test`: shell readout moved strongly (`0.204` vs `0.035`) while inner `ZZ` changed only modestly (`0.931` vs `0.911`).
- `cross_platform_universality_test`: both IQM and Rigetti preserved `shield > open > wrong` ordering.
- `collapse_rate_atlas_microgrid`: the current grid stayed in one intermediate regime (`ZZ 0.740-0.792`, outer excitation `0.256-0.510`).
