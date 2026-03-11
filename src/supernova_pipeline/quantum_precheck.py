from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import boto3
import numpy as np
from braket.ahs import AnalogHamiltonianSimulation, AtomArrangement, DrivingField, Hamiltonian
from braket.aws import AwsDevice, AwsSession
from braket.circuits import Circuit
from braket.devices import LocalSimulator
from dotenv import dotenv_values


ROOT = Path("/home/primary/SUPERNOVA")
HUBBLE_ROOT = Path("/home/primary/Hubble-Systematics-Review-Chain")
PRICING_SOURCE_URL = "https://aws.amazon.com/braket/pricing/"


@dataclass(frozen=True)
class DevicePlan:
    key: str
    provider: str
    name: str
    region: str
    arn: str
    pricing_key: str
    program_type: str


PRICING = {
    "iqm_garnet": {
        "label": "IQM Garnet",
        "task_fee_usd": 0.30,
        "shot_fee_usd": 0.00145,
        "source_url": PRICING_SOURCE_URL,
        "pricing_checked_utc": "2026-03-10T00:00:00Z",
    },
    "rigetti_ankaa": {
        "label": "Rigetti Ankaa-3",
        "task_fee_usd": 0.30,
        "shot_fee_usd": 0.00090,
        "source_url": PRICING_SOURCE_URL,
        "pricing_checked_utc": "2026-03-10T00:00:00Z",
    },
    "quera_aquila": {
        "label": "QuEra Aquila",
        "task_fee_usd": 0.30,
        "shot_fee_usd": 0.01000,
        "source_url": PRICING_SOURCE_URL,
        "pricing_checked_utc": "2026-03-10T00:00:00Z",
    },
}

DEVICES = {
    "iqm_garnet": DevicePlan(
        key="iqm_garnet",
        provider="IQM",
        name="Garnet",
        region="eu-north-1",
        arn="arn:aws:braket:eu-north-1::device/qpu/iqm/Garnet",
        pricing_key="iqm_garnet",
        program_type="openqasm",
    ),
    "rigetti_ankaa": DevicePlan(
        key="rigetti_ankaa",
        provider="Rigetti",
        name="Ankaa-3",
        region="us-west-1",
        arn="arn:aws:braket:us-west-1::device/qpu/rigetti/Ankaa-3",
        pricing_key="rigetti_ankaa",
        program_type="openqasm",
    ),
    "quera_aquila": DevicePlan(
        key="quera_aquila",
        provider="QuEra",
        name="Aquila",
        region="us-east-1",
        arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila",
        pricing_key="quera_aquila",
        program_type="ahs",
    ),
}

_GATE_SIM = LocalSimulator("braket_sv")
_AHS_SIM = LocalSimulator("braket_ahs")


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _stamp() -> str:
    return _utc_now().strftime("%Y%m%d_%H%M%SUTC")


def _iso() -> str:
    return _utc_now().isoformat().replace("+00:00", "Z")


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _git_like(root: Path) -> dict[str, Any]:
    def run(args: list[str]) -> str | None:
        try:
            out = subprocess.check_output(args, cwd=root, stderr=subprocess.DEVNULL)
            return out.decode("utf-8", errors="replace").strip()
        except Exception:
            return None

    return {
        "head": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "status_porcelain": run(["git", "status", "--porcelain"]),
    }


def _env_like() -> dict[str, Any]:
    return {
        "timestamp_utc": _iso(),
        "python": sys.version,
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "mkl_num_threads": os.environ.get("MKL_NUM_THREADS"),
        "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        "numexpr_num_threads": os.environ.get("NUMEXPR_NUM_THREADS"),
    }


def _write_progress(path: Path, stage: str, status: str, extra: dict[str, Any] | None = None) -> None:
    payload = {"timestamp_utc": _iso(), "stage": stage, "status": status}
    if extra:
        payload.update(extra)
    _write_json(path, payload)


def _load_aws_credentials(root_dir: Path) -> tuple[dict[str, str] | None, dict[str, Any]]:
    search_paths = [
        root_dir / ".env",
        HUBBLE_ROOT / ".env",
    ]
    env_key = os.environ.get("AWS_ACCESS_KEY_ID") or os.environ.get("AWS_ACCESS_KEY")
    env_secret = os.environ.get("AWS_SECRET_ACCESS_KEY") or os.environ.get("AWS_ACCESS_SECRET")
    if env_key and env_secret:
        return (
            {"aws_access_key_id": env_key, "aws_secret_access_key": env_secret},
            {"source": "environment", "search_paths": [str(p) for p in search_paths]},
        )
    for candidate in search_paths:
        if not candidate.exists():
            continue
        vals = dotenv_values(candidate)
        key = vals.get("AWS_ACCESS_KEY") or vals.get("AWS_ACCESS_KEY_ID")
        secret = vals.get("AWS_ACCESS_SECRET") or vals.get("AWS_SECRET_ACCESS_KEY")
        if key and secret:
            return (
                {"aws_access_key_id": key, "aws_secret_access_key": secret},
                {"source": str(candidate), "search_paths": [str(p) for p in search_paths]},
            )
    return None, {"source": None, "search_paths": [str(p) for p in search_paths]}


def _aws_session(creds: dict[str, str], region: str) -> AwsSession:
    return AwsSession(
        boto_session=boto3.Session(
            aws_access_key_id=creds["aws_access_key_id"],
            aws_secret_access_key=creds["aws_secret_access_key"],
            region_name=region,
        )
    )


def _find_chain(connectivity_graph: dict[str, list[str]], length: int) -> list[int] | None:
    graph = {int(k): sorted(int(v) for v in values) for k, values in connectivity_graph.items()}

    def dfs(path: list[int]) -> list[int] | None:
        if len(path) == length:
            return path
        current = path[-1]
        for nxt in graph.get(current, []):
            if nxt in path:
                continue
            found = dfs(path + [nxt])
            if found:
                return found
        return None

    for node in sorted(graph):
        found = dfs([node])
        if found:
            return found
    return None


def _device_inventory(creds: dict[str, str] | None) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, spec in DEVICES.items():
        info: dict[str, Any] = {
            "provider": spec.provider,
            "name": spec.name,
            "region": spec.region,
            "arn": spec.arn,
            "pricing": PRICING[spec.pricing_key],
            "program_type": spec.program_type,
        }
        if creds is None:
            info.update({"status": "UNKNOWN", "device_check_ok": False, "reason": "missing_credentials"})
            payload[key] = info
            continue
        try:
            session = _aws_session(creds, spec.region)
            device = AwsDevice(spec.arn, aws_session=session)
            info["status"] = str(device.status)
            info["device_check_ok"] = True
            paradigm = getattr(device.properties, "paradigm", None)
            info["qubit_count"] = int(getattr(paradigm, "qubitCount", 0) or 0)
            action_keys = [str(k.value if hasattr(k, "value") else k) for k in getattr(device.properties, "action", {}).keys()]
            info["action_keys"] = sorted(action_keys)
            connectivity = getattr(paradigm, "connectivity", None)
            graph = getattr(connectivity, "connectivityGraph", {}) or {}
            if graph:
                info["chain3"] = _find_chain(graph, 3)
                info["chain4"] = _find_chain(graph, 4)
                info["connectivity_nodes"] = len(graph)
            else:
                info["chain3"] = None
                info["chain4"] = None
                info["connectivity_nodes"] = 0
            payload[key] = info
        except Exception as exc:
            info.update({"status": "ERROR", "device_check_ok": False, "reason": f"{type(exc).__name__}: {exc}"})
            payload[key] = info
    return payload


def _cost(pricing_key: str, *, shots: int, task_count: int) -> dict[str, Any]:
    model = PRICING[pricing_key]
    per_task = float(model["task_fee_usd"] + model["shot_fee_usd"] * shots)
    total = float(per_task * task_count)
    return {
        "pricing_label": model["label"],
        "source_url": model["source_url"],
        "pricing_checked_utc": model["pricing_checked_utc"],
        "task_fee_usd": float(model["task_fee_usd"]),
        "shot_fee_usd": float(model["shot_fee_usd"]),
        "shots_per_task": int(shots),
        "task_count": int(task_count),
        "per_task_estimate_usd": per_task,
        "total_estimate_usd": total,
    }


def _instruction_names(circuit: Circuit) -> list[str]:
    names: list[str] = []
    for instruction in circuit.instructions:
        operator = getattr(instruction, "operator", None)
        names.append(type(operator).__name__)
    return names


def _probabilities(circuit: Circuit, targets: list[int]) -> np.ndarray:
    c = circuit.copy()
    c.probability(target=targets)
    return np.asarray(_GATE_SIM.run(c, shots=0).result().values[0], dtype=float)


def _bell_proxy(circuit: Circuit, qa: int = 0, qb: int = 1) -> dict[str, float]:
    zz_probs = _probabilities(circuit, [qa, qb])
    zz_corr = float((zz_probs[0] + zz_probs[3]) - (zz_probs[1] + zz_probs[2]))
    xx_circuit = circuit.copy()
    xx_circuit.h(qa)
    xx_circuit.h(qb)
    xx_probs = _probabilities(xx_circuit, [qa, qb])
    xx_corr = float((xx_probs[0] + xx_probs[3]) - (xx_probs[1] + xx_probs[2]))
    return {
        "zz_corr": zz_corr,
        "xx_corr": xx_corr,
        "bell_proxy": 0.5 * (zz_corr + xx_corr),
    }


def _shell_excitation(circuit: Circuit, target: int) -> float:
    probs = _probabilities(circuit, [target])
    return float(probs[1])


def _outer_excitation(circuit: Circuit, targets: list[int]) -> float:
    probs = _probabilities(circuit, targets)
    gross = 0.0
    for idx, prob in enumerate(probs):
        gross += prob * sum((idx >> shift) & 1 for shift in range(len(targets) - 1, -1, -1))
    return float(gross)


def _core_shell_circuit(shell_label: str, family: str) -> Circuit:
    circuit = Circuit().h(0).cnot(0, 1)
    if shell_label == "plus":
        circuit.h(2)
    if family != "uncoupled":
        repeats = {"shield": 2, "open": 1, "wrong": 1}[family]
        for _ in range(repeats):
            circuit.cz(2, 1)
        if family == "wrong":
            circuit.x(1)
    return circuit


def _neutrinosphere_circuit(theta: float, family: str) -> Circuit:
    circuit = Circuit().h(0).cnot(0, 1).ry(2, theta).ry(3, theta / 2.0)
    if family == "matched":
        circuit.cz(1, 2).cz(1, 2).cz(1, 3).cz(1, 3)
    elif family == "mismatched":
        circuit.cz(1, 2).cz(1, 3)
    elif family == "blocked":
        circuit.ry(2, -theta).ry(3, -theta / 2.0)
    return circuit


def _shell_visibility_circuit(angle: float, mode: str) -> Circuit:
    shell_angle = float(angle if mode == "bright" else angle / 3.0)
    return Circuit().h(0).cnot(0, 1).ry(2, shell_angle)


def _collapse_atlas_circuit(theta: float, phase: float) -> Circuit:
    return (
        Circuit()
        .h(0)
        .cnot(0, 1)
        .ry(2, theta)
        .ry(3, theta / 2.0)
        .rz(1, phase)
        .cz(1, 2)
        .cz(1, 3)
    )


def _shock_revival_program(scale: float, detuning_peak: float) -> AnalogHamiltonianSimulation:
    register = AtomArrangement()
    spacing = 6.1e-6 * scale
    for idx in range(4):
        register.add((idx * spacing, 0.0))
    times = [0.0, 1.0e-6, 2.0e-6, 3.0e-6]
    amplitudes = [0.0, 1.5e6, 1.5e6, 0.0]
    detunings = [-3.0e6, -1.0e6, detuning_peak, -2.0e6]
    phases = [0.0, 0.0, 0.0, 0.0]
    hamiltonian = Hamiltonian()
    hamiltonian += DrivingField.from_lists(
        times=times,
        amplitudes=amplitudes,
        detunings=detunings,
        phases=phases,
    )
    return AnalogHamiltonianSimulation(register=register, hamiltonian=hamiltonian)


def _ahs_edge_loss_fraction(program: AnalogHamiltonianSimulation, shots: int = 64) -> dict[str, float]:
    result = _AHS_SIM.run(program, shots=shots).result()
    successful = [shot for shot in result.measurements if getattr(shot, "status", None).value == "Success"]
    if not successful:
        return {"edge_loss_fraction": 0.0, "mean_post_sequence": 0.0}
    post = np.asarray([shot.post_sequence for shot in successful], dtype=float)
    edge_mean = float(np.mean(post[:, -1]))
    return {
        "edge_loss_fraction": float(1.0 - edge_mean),
        "mean_post_sequence": edge_mean,
    }


def _core_shell_precheck(device_inventory: dict[str, Any]) -> dict[str, Any]:
    device = device_inventory["iqm_garnet"]
    shots = 1500
    manifest: list[dict[str, Any]] = []
    for family in ["shield", "open", "wrong"]:
        for shell_label in ["zero", "plus"]:
            for basis in ["ZZ", "XX"]:
                manifest.append({"family": family, "shell_label": shell_label, "basis": basis})
    for basis in ["ZZ", "XX"]:
        manifest.append({"family": "uncoupled", "shell_label": "plus", "basis": basis})

    summary: dict[str, Any] = {"families": {}, "task_count": len(manifest), "shots": shots}
    for family in ["shield", "open", "wrong", "uncoupled"]:
        zero = _bell_proxy(_core_shell_circuit("zero", family))
        plus = _bell_proxy(_core_shell_circuit("plus", family))
        summary["families"][family] = {
            "bell_zero": zero["bell_proxy"],
            "bell_plus": plus["bell_proxy"],
            "shell_sensitivity": abs(plus["bell_proxy"] - zero["bell_proxy"]),
            "zz_plus": plus["zz_corr"],
            "xx_plus": plus["xx_corr"],
        }
    ordering_ok = (
        summary["families"]["shield"]["bell_plus"]
        > summary["families"]["open"]["bell_plus"]
        > summary["families"]["wrong"]["bell_plus"]
    )
    submit_ready = bool(
        ordering_ok
        and summary["families"]["shield"]["shell_sensitivity"] < 0.05
        and summary["families"]["open"]["shell_sensitivity"] > 0.30
        and summary["families"]["wrong"]["bell_plus"] < 0.0
        and device.get("status") == "ONLINE"
        and device.get("chain3")
    )
    return {
        "label": "core_shell_collapse_ladder",
        "question": "Can the core stay organized while shell-mediated export fails sharply?",
        "platform": "Superconducting qubits",
        "device": device,
        "pilot_plan": {
            "device_key": "iqm_garnet",
            "task_count": len(manifest),
            "shots": shots,
            "hard_task_cap": 16,
            "chosen_chain": device.get("chain3"),
            "submit_mode": "off_by_default",
        },
        "local_precheck": {
            "ordering_ok": ordering_ok,
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
            "summary": summary,
        },
        "cost_estimate_usd": _cost("iqm_garnet", shots=shots, task_count=len(manifest)),
        "manifest_preview": manifest[:6],
        "sample_circuit_instructions": _instruction_names(_core_shell_circuit("plus", "open")),
    }


def _neutrinosphere_precheck(device_inventory: dict[str, Any]) -> dict[str, Any]:
    device = device_inventory["rigetti_ankaa"]
    shots = 1200
    theta_values = [0.4, 0.6, 0.8, 1.0]
    manifest = [{"family": family, "theta": theta} for family in ["matched", "mismatched"] for theta in theta_values]
    manifest.extend([{"family": "blocked", "theta": theta} for theta in [0.6, 1.0]])
    summary: dict[str, Any] = {"families": {}}
    for family in ["matched", "mismatched", "blocked"]:
        rows = []
        for theta in theta_values if family != "blocked" else [0.6, 1.0]:
            circuit = _neutrinosphere_circuit(theta, family)
            bell = _bell_proxy(circuit)["bell_proxy"]
            gross = _outer_excitation(circuit, [2, 3])
            rows.append({"theta": theta, "inner_bell_proxy": bell, "gross_outer_excitation": gross})
        summary["families"][family] = rows
    matched = {float(row["theta"]): row for row in summary["families"]["matched"]}
    mismatched = {float(row["theta"]): row for row in summary["families"]["mismatched"]}
    bell_gaps = [matched[theta]["inner_bell_proxy"] - mismatched[theta]["inner_bell_proxy"] for theta in matched]
    gross_gaps = [abs(matched[theta]["gross_outer_excitation"] - mismatched[theta]["gross_outer_excitation"]) for theta in matched]
    submit_ready = bool(
        max(bell_gaps) >= 0.15
        and max(gross_gaps) <= 1e-6 + 0.02
        and device.get("status") == "ONLINE"
        and device.get("chain4")
    )
    return {
        "label": "quantum_neutrinosphere_analog",
        "question": "Can gross leakage survive while structured transport fails?",
        "platform": "Gate-model pilot on Rigetti as a low-cost proxy for the broader trapped-ion/cold-atom concept",
        "device": device,
        "pilot_plan": {
            "device_key": "rigetti_ankaa",
            "task_count": len(manifest),
            "shots": shots,
            "hard_task_cap": 12,
            "chosen_chain": device.get("chain4"),
            "submit_mode": "off_by_default",
        },
        "local_precheck": {
            "max_bell_gap": float(max(bell_gaps)),
            "max_gross_gap": float(max(gross_gaps)),
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
            "summary": summary,
        },
        "cost_estimate_usd": _cost("rigetti_ankaa", shots=shots, task_count=len(manifest)),
        "manifest_preview": manifest[:6],
        "sample_circuit_instructions": _instruction_names(_neutrinosphere_circuit(0.8, "mismatched")),
    }


def _shock_revival_precheck(device_inventory: dict[str, Any], creds: dict[str, str] | None) -> dict[str, Any]:
    device = device_inventory["quera_aquila"]
    shots = 100
    detuning_values = [1.0e6, 2.0e6, 3.0e6, 4.0e6]
    scale_values = [0.85, 1.00]
    manifest = [{"scale": scale, "detuning_peak": det} for scale in scale_values for det in detuning_values]
    rows = []
    discretized_ok = False
    if creds is not None and device.get("status") == "ONLINE":
        session = _aws_session(creds, DEVICES["quera_aquila"].region)
        aquila = AwsDevice(DEVICES["quera_aquila"].arn, aws_session=session)
    else:
        aquila = None
    for item in manifest:
        program = _shock_revival_program(item["scale"], item["detuning_peak"])
        if aquila is not None:
            program = program.discretize(aquila)
            discretized_ok = True
        metrics = _ahs_edge_loss_fraction(program, shots=64)
        rows.append(
            {
                "scale": item["scale"],
                "detuning_peak": item["detuning_peak"],
                **metrics,
                "atom_count": len(program.register),
            }
        )
    edge_losses = [row["edge_loss_fraction"] for row in rows]
    submit_ready = bool(
        (max(edge_losses) - min(edge_losses)) >= 0.12
        and device.get("status") == "ONLINE"
        and "braket.ir.ahs.program" in device.get("action_keys", [])
    )
    return {
        "label": "shock_revival_threshold_simulator",
        "question": "Does a small neutral-atom analog show a threshold-like revival split?",
        "platform": "QuEra Aquila analog Hamiltonian simulation",
        "device": device,
        "pilot_plan": {
            "device_key": "quera_aquila",
            "task_count": len(manifest),
            "shots": shots,
            "hard_task_cap": 10,
            "atom_count": 4,
            "submit_mode": "off_by_default",
            "discretized_against_device": discretized_ok,
        },
        "local_precheck": {
            "edge_loss_span": float(max(edge_losses) - min(edge_losses)),
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
            "summary": rows,
        },
        "cost_estimate_usd": _cost("quera_aquila", shots=shots, task_count=len(manifest)),
        "manifest_preview": manifest[:6],
        "sample_program_preview": {
            "times_s": [0.0, 1.0e-6, 2.0e-6, 3.0e-6],
            "amplitudes_rad_per_s": [0.0, 1.5e6, 1.5e6, 0.0],
            "detunings_rad_per_s": [-3.0e6, -1.0e6, detuning_values[0], -2.0e6],
        },
    }


def _shell_visibility_precheck(device_inventory: dict[str, Any]) -> dict[str, Any]:
    device = device_inventory["iqm_garnet"]
    shots = 1200
    angles = [0.3, 0.7, 1.1, 1.4]
    manifest = [{"mode": mode, "angle": angle} for mode in ["bright", "quiet"] for angle in angles]
    summary: dict[str, list[dict[str, float]]] = {"bright": [], "quiet": []}
    for mode in ["bright", "quiet"]:
        for angle in angles:
            circuit = _shell_visibility_circuit(angle, mode)
            summary[mode].append(
                {
                    "angle": angle,
                    "inner_bell_proxy": _bell_proxy(circuit)["bell_proxy"],
                    "shell_excitation": _shell_excitation(circuit, 2),
                }
            )
    mean_bright = float(np.mean([row["shell_excitation"] for row in summary["bright"]]))
    mean_quiet = float(np.mean([row["shell_excitation"] for row in summary["quiet"]]))
    bell_drift = abs(
        float(np.mean([row["inner_bell_proxy"] for row in summary["bright"]]))
        - float(np.mean([row["inner_bell_proxy"] for row in summary["quiet"]]))
    )
    submit_ready = bool(
        (mean_bright - mean_quiet) >= 0.15
        and bell_drift <= 0.02
        and device.get("status") == "ONLINE"
        and device.get("chain3")
    )
    return {
        "label": "shell_only_visibility_test",
        "question": "Can shell readout geometry change visibility while the interior stays fixed?",
        "platform": "Superconducting qubits",
        "device": device,
        "pilot_plan": {
            "device_key": "iqm_garnet",
            "task_count": len(manifest),
            "shots": shots,
            "hard_task_cap": 10,
            "chosen_chain": device.get("chain3"),
            "submit_mode": "off_by_default",
        },
        "local_precheck": {
            "mean_bright_visibility": mean_bright,
            "mean_quiet_visibility": mean_quiet,
            "inner_bell_drift": bell_drift,
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
            "summary": summary,
        },
        "cost_estimate_usd": _cost("iqm_garnet", shots=shots, task_count=len(manifest)),
        "manifest_preview": manifest[:6],
        "sample_circuit_instructions": _instruction_names(_shell_visibility_circuit(1.1, "bright")),
    }


def _cross_platform_precheck(device_inventory: dict[str, Any]) -> dict[str, Any]:
    manifest = [{"family": family, "shell_label": "plus"} for family in ["shield", "open", "wrong"]]
    manifest = [dict(item, basis=basis) for item in manifest for basis in ["ZZ", "XX"]]
    local_metrics = {
        family: _bell_proxy(_core_shell_circuit("plus", family))["bell_proxy"]
        for family in ["shield", "open", "wrong"]
    }
    ordering_ok = local_metrics["shield"] > local_metrics["open"] > local_metrics["wrong"]
    device_runs = []
    total_cost = 0.0
    for key in ["iqm_garnet", "rigetti_ankaa"]:
        device = device_inventory[key]
        shots = 1000
        cost = _cost(key, shots=shots, task_count=len(manifest))
        total_cost += cost["total_estimate_usd"]
        device_runs.append(
            {
                "device_key": key,
                "device": device,
                "shots": shots,
                "task_count": len(manifest),
                "chosen_chain": device.get("chain3"),
                "cost_estimate_usd": cost,
            }
        )
    submit_ready = bool(
        ordering_ok
        and all(run["device"].get("status") == "ONLINE" and run["chosen_chain"] for run in device_runs)
    )
    return {
        "label": "cross_platform_universality_test",
        "question": "Does the same reduced branch ordering package cleanly for two hardware families?",
        "platform": "IQM Garnet + Rigetti Ankaa-3",
        "devices": device_runs,
        "pilot_plan": {
            "total_task_count": len(manifest) * len(device_runs),
            "submit_mode": "off_by_default",
            "hard_task_cap_per_device": 8,
        },
        "local_precheck": {
            "ordering_ok": ordering_ok,
            "local_metrics": local_metrics,
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
        },
        "cost_estimate_usd": {
            "source_url": PRICING_SOURCE_URL,
            "total_estimate_usd": float(total_cost),
        },
        "manifest_preview": manifest[:6],
    }


def _collapse_rate_classification(bell: float, gross: float) -> str:
    if bell >= 0.82 and gross <= 0.15:
        return "luminous_export_like"
    if bell < 0.55 and gross >= 0.05:
        return "silent_collapse_like"
    return "partial_fallback_like"


def _collapse_rate_atlas_precheck(device_inventory: dict[str, Any]) -> dict[str, Any]:
    device = device_inventory["rigetti_ankaa"]
    shots = 800
    theta_values = [0.4, 0.7, 1.0]
    phase_values = [0.0, 0.8, 1.6]
    manifest = [{"theta": theta, "phase": phase} for theta in theta_values for phase in phase_values]
    rows = []
    for item in manifest:
        circuit = _collapse_atlas_circuit(item["theta"], item["phase"])
        bell = _bell_proxy(circuit)["bell_proxy"]
        gross = _outer_excitation(circuit, [2, 3])
        rows.append(
            {
                "theta": item["theta"],
                "phase": item["phase"],
                "inner_bell_proxy": bell,
                "gross_outer_excitation": gross,
                "classification": _collapse_rate_classification(bell, gross),
            }
        )
    class_counts: dict[str, int] = {}
    for row in rows:
        class_counts[row["classification"]] = class_counts.get(row["classification"], 0) + 1
    submit_ready = bool(len(class_counts) >= 2 and device.get("status") == "ONLINE" and device.get("chain4"))
    return {
        "label": "collapse_rate_atlas_microgrid",
        "question": "Can a tiny pilot grid already separate more than one branch class?",
        "platform": "Gate-model microgrid on Rigetti",
        "device": device,
        "pilot_plan": {
            "device_key": "rigetti_ankaa",
            "task_count": len(manifest),
            "shots": shots,
            "hard_task_cap": 10,
            "chosen_chain": device.get("chain4"),
            "submit_mode": "off_by_default",
        },
        "local_precheck": {
            "classification_counts": class_counts,
            "submit_ready": submit_ready,
            "submit_reason": "preflight_pass" if submit_ready else "preflight_fail",
            "summary": rows,
        },
        "cost_estimate_usd": _cost("rigetti_ankaa", shots=shots, task_count=len(manifest)),
        "manifest_preview": manifest[:6],
        "sample_circuit_instructions": _instruction_names(_collapse_atlas_circuit(0.7, 0.8)),
    }


def _build_report(
    *,
    output_dir: Path,
    credential_meta: dict[str, Any],
    device_inventory: dict[str, Any],
    tests: list[dict[str, Any]],
    total_cost: float,
) -> str:
    lines = [
        "# Quantum Pilot Precheck Report",
        "",
        f"Created UTC: `{_iso()}`",
        "",
        "## Scope",
        "",
        "This packet defines pilot-level dry-run designs for the six quantum tests in the current proposal.",
        "Nothing in this packet was submitted to AWS Braket.",
        "",
        "## Agent/Workflow Notes",
        "",
        "- `/home/primary/Hubble-Systematics-Review-Chain` supplied the Braket runner pattern and pilot preflight model.",
        "- `/home/primary/BAO_Ruler_Experiments` currently contains no `AGENT.md` or `AGENTS.md`, and no Braket submission code was found there.",
        f"- Credential source used for live device checks: `{credential_meta.get('source')}`",
        "",
        "## Live Device Inventory",
        "",
    ]
    for key, info in device_inventory.items():
        lines.append(f"- `{info['provider']} {info['name']}` in `{info['region']}`: status `{info.get('status')}`; actions `{', '.join(info.get('action_keys', [])) or 'n/a'}`")
    lines.extend(["", "## Pilot Tests", ""])
    for test in tests:
        lines.append(f"### `{test['label']}`")
        lines.append("")
        lines.append(f"- Question: {test['question']}")
        lines.append(f"- Platform: {test['platform']}")
        pilot = test["pilot_plan"]
        if "task_count" in pilot:
            lines.append(f"- Planned hardware footprint: `{pilot['task_count']}` tasks at `{pilot['shots']}` shots")
        else:
            lines.append(f"- Planned hardware footprint: `{pilot['total_task_count']}` tasks total")
        submit = test["local_precheck"]["submit_ready"]
        lines.append(f"- Precheck verdict: `{'READY' if submit else 'NOT_READY'}` (`{test['local_precheck']['submit_reason']}`)")
        cost = test["cost_estimate_usd"]["total_estimate_usd"]
        lines.append(f"- Estimated pilot spend: `${cost:.2f}`")
        lines.append("")
    lines.extend(
        [
            "## Total Estimated Spend",
            "",
            f"Current total for all six pilot quantum tests: `${total_cost:.2f}`",
            "",
            "Pricing assumptions were taken from the official AWS Braket pricing page on 2026-03-10:",
            PRICING_SOURCE_URL,
            "",
            "## Output Files",
            "",
            f"- [stage01_device_inventory.json]({output_dir / 'stage01_device_inventory.json'})",
            f"- [stage02_test_designs.json]({output_dir / 'stage02_test_designs.json'})",
            f"- [stage03_precheck_results.json]({output_dir / 'stage03_precheck_results.json'})",
            f"- [stage04_cost_summary.json]({output_dir / 'stage04_cost_summary.json'})",
        ]
    )
    return "\n".join(lines) + "\n"


def run_quantum_precheck(root_dir: Path, output_dir: Path | None = None) -> dict[str, Path]:
    out_dir = output_dir or (root_dir / "outputs" / f"quantum_pilot_precheck_{_stamp()}")
    _ensure_dir(out_dir)
    progress = out_dir / "progress.json"

    _write_json(out_dir / "env.json", _env_like())
    _write_json(out_dir / "git_like.json", _git_like(root_dir))

    creds, credential_meta = _load_aws_credentials(root_dir)
    config = {
        "run_label": out_dir.name,
        "created_utc": _iso(),
        "pricing_source_url": PRICING_SOURCE_URL,
        "credential_discovery": credential_meta,
        "tests": [
            "core_shell_collapse_ladder",
            "quantum_neutrinosphere_analog",
            "shock_revival_threshold_simulator",
            "shell_only_visibility_test",
            "cross_platform_universality_test",
            "collapse_rate_atlas_microgrid",
        ],
    }
    _write_json(out_dir / "config.json", config)

    _write_progress(progress, "stage01_device_inventory", "running")
    inventory = _device_inventory(creds)
    _write_json(out_dir / "stage01_device_inventory.json", inventory)

    _write_progress(progress, "stage02_test_designs", "running")
    designs = {
        "core_shell_collapse_ladder": {
            "device_key": "iqm_garnet",
            "planned_tasks": 14,
            "planned_shots": 1500,
            "goal": "Recover shield > open > wrong ordering with low shell sensitivity for the shield branch.",
        },
        "quantum_neutrinosphere_analog": {
            "device_key": "rigetti_ankaa",
            "planned_tasks": 10,
            "planned_shots": 1200,
            "goal": "Hold gross outer leakage roughly fixed while degrading inner Bell structure in the mismatched branch.",
        },
        "shock_revival_threshold_simulator": {
            "device_key": "quera_aquila",
            "planned_tasks": 8,
            "planned_shots": 100,
            "goal": "Find a small AHS pulse family with a nontrivial edge-loss span across a narrow parameter sweep.",
        },
        "shell_only_visibility_test": {
            "device_key": "iqm_garnet",
            "planned_tasks": 8,
            "planned_shots": 1200,
            "goal": "Change shell visibility strongly while keeping the interior Bell proxy stable.",
        },
        "cross_platform_universality_test": {
            "device_keys": ["iqm_garnet", "rigetti_ankaa"],
            "planned_tasks_per_device": 6,
            "planned_shots_per_device": 1000,
            "goal": "Package the same reduced branch-ordering circuit family on two hardware classes.",
        },
        "collapse_rate_atlas_microgrid": {
            "device_key": "rigetti_ankaa",
            "planned_tasks": 9,
            "planned_shots": 800,
            "goal": "Ensure a tiny hardware grid already separates more than one branch class.",
        },
    }
    _write_json(out_dir / "stage02_test_designs.json", designs)

    _write_progress(progress, "stage03_precheck_results", "running")
    tests = [
        _core_shell_precheck(inventory),
        _neutrinosphere_precheck(inventory),
        _shock_revival_precheck(inventory, creds),
        _shell_visibility_precheck(inventory),
        _cross_platform_precheck(inventory),
        _collapse_rate_atlas_precheck(inventory),
    ]
    precheck_results = {test["label"]: test for test in tests}
    _write_json(out_dir / "stage03_precheck_results.json", precheck_results)

    _write_progress(progress, "stage04_cost_summary", "running")
    total_cost = float(sum(test["cost_estimate_usd"]["total_estimate_usd"] for test in tests))
    cost_summary = {
        "pricing_source_url": PRICING_SOURCE_URL,
        "pricing_checked_utc": "2026-03-10T00:00:00Z",
        "tests": {test["label"]: test["cost_estimate_usd"] for test in tests},
        "total_estimate_usd": total_cost,
    }
    _write_json(out_dir / "stage04_cost_summary.json", cost_summary)

    report_text = _build_report(
        output_dir=out_dir,
        credential_meta=credential_meta,
        device_inventory=inventory,
        tests=tests,
        total_cost=total_cost,
    )
    _write_text(out_dir / "quantum_precheck_report.md", report_text)

    _write_progress(progress, "complete", "done", {"total_estimate_usd": total_cost, "test_count": len(tests)})
    return {
        "config": out_dir / "config.json",
        "device_inventory": out_dir / "stage01_device_inventory.json",
        "test_designs": out_dir / "stage02_test_designs.json",
        "precheck_results": out_dir / "stage03_precheck_results.json",
        "cost_summary": out_dir / "stage04_cost_summary.json",
        "report": out_dir / "quantum_precheck_report.md",
        "progress": progress,
    }
