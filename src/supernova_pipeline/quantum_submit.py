from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import boto3
from braket.aws import AwsDevice, AwsSession
from braket.circuits import Circuit

from .quantum_precheck import (
    DEVICES,
    ROOT,
    _collapse_atlas_circuit,
    _core_shell_circuit,
    _ensure_dir,
    _git_like,
    _iso,
    _load_aws_credentials,
    _neutrinosphere_circuit,
    _shell_visibility_circuit,
    _shock_revival_program,
    _stamp,
    _utc_now,
    _write_json,
    _write_progress,
)


def _aws_session(creds: dict[str, str], region: str) -> AwsSession:
    return AwsSession(
        boto_session=boto3.Session(
            aws_access_key_id=creds["aws_access_key_id"],
            aws_secret_access_key=creds["aws_secret_access_key"],
            region_name=region,
        )
    )


def _aws_boto_session(creds: dict[str, str], region: str) -> boto3.Session:
    return boto3.Session(
        aws_access_key_id=creds["aws_access_key_id"],
        aws_secret_access_key=creds["aws_secret_access_key"],
        region_name=region,
    )


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _env_like() -> dict[str, Any]:
    return {
        "timestamp_utc": _iso(),
        "submission_mode": "live_hardware",
    }


def _ensure_bucket(creds: dict[str, str], region: str) -> dict[str, Any]:
    session = _aws_boto_session(creds, region)
    sts = session.client("sts")
    account_id = sts.get_caller_identity()["Account"]
    bucket = f"amazon-braket-{account_id}-{region}"
    s3 = session.client("s3")
    existing = {bucket_info["Name"] for bucket_info in s3.list_buckets().get("Buckets", [])}
    created = False
    if bucket not in existing:
        kwargs: dict[str, Any] = {"Bucket": bucket}
        if region != "us-east-1":
            kwargs["CreateBucketConfiguration"] = {"LocationConstraint": region}
        s3.create_bucket(**kwargs)
        created = True
    probe_key = f"supernova/bootstrap/{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.txt"
    s3.put_object(Bucket=bucket, Key=probe_key, Body=b"braket-ready")
    return {"account_id": account_id, "bucket": bucket, "region": region, "created": created, "probe_key": probe_key}


def _measured_circuit(base: Circuit, measured_qubits: list[int]) -> Circuit:
    circuit = base.copy()
    for q in measured_qubits:
        circuit.measure(q)
    return circuit


def _core_shell_manifest() -> tuple[list[dict[str, Any]], list[Circuit]]:
    manifest: list[dict[str, Any]] = []
    circuits: list[Circuit] = []
    for family in ["shield", "open", "wrong"]:
        for shell_label in ["zero", "plus"]:
            base = _core_shell_circuit(shell_label, family)
            circuits.append(_measured_circuit(base, [0, 1, 2]))
            manifest.append({"family": family, "shell_label": shell_label, "measured_qubits": [0, 1, 2]})
    base = _core_shell_circuit("plus", "uncoupled")
    circuits.append(_measured_circuit(base, [0, 1, 2]))
    manifest.append({"family": "uncoupled", "shell_label": "plus", "measured_qubits": [0, 1, 2]})
    return manifest, circuits


def _neutrinosphere_manifest() -> tuple[list[dict[str, Any]], list[Circuit]]:
    manifest: list[dict[str, Any]] = []
    circuits: list[Circuit] = []
    for family in ["matched", "mismatched"]:
        for theta in [0.4, 0.6, 0.8, 1.0]:
            base = _neutrinosphere_circuit(theta, family)
            circuits.append(_measured_circuit(base, [0, 1, 2, 3]))
            manifest.append({"family": family, "theta": theta, "measured_qubits": [0, 1, 2, 3]})
    for theta in [0.6, 1.0]:
        base = _neutrinosphere_circuit(theta, "blocked")
        circuits.append(_measured_circuit(base, [0, 1, 2, 3]))
        manifest.append({"family": "blocked", "theta": theta, "measured_qubits": [0, 1, 2, 3]})
    return manifest, circuits


def _shell_visibility_manifest() -> tuple[list[dict[str, Any]], list[Circuit]]:
    manifest: list[dict[str, Any]] = []
    circuits: list[Circuit] = []
    for mode in ["bright", "quiet"]:
        for angle in [0.3, 0.7, 1.1, 1.4]:
            base = _shell_visibility_circuit(angle, mode)
            circuits.append(_measured_circuit(base, [0, 1, 2]))
            manifest.append({"mode": mode, "angle": angle, "measured_qubits": [0, 1, 2]})
    return manifest, circuits


def _cross_platform_manifest() -> tuple[list[dict[str, Any]], list[Circuit]]:
    manifest: list[dict[str, Any]] = []
    circuits: list[Circuit] = []
    for family in ["shield", "open", "wrong"]:
        base = _core_shell_circuit("plus", family)
        circuits.append(_measured_circuit(base, [0, 1, 2]))
        manifest.append({"family": family, "shell_label": "plus", "measured_qubits": [0, 1, 2]})
    return manifest, circuits


def _collapse_atlas_manifest() -> tuple[list[dict[str, Any]], list[Circuit]]:
    manifest: list[dict[str, Any]] = []
    circuits: list[Circuit] = []
    for theta in [0.4, 0.7, 1.0]:
        for phase in [0.0, 0.8, 1.6]:
            base = _collapse_atlas_circuit(theta, phase)
            circuits.append(_measured_circuit(base, [0, 1, 2, 3]))
            manifest.append({"theta": theta, "phase": phase, "measured_qubits": [0, 1, 2, 3]})
    return manifest, circuits


def _shock_revival_manifest() -> tuple[list[dict[str, Any]], list[Any]]:
    manifest: list[dict[str, Any]] = []
    programs: list[Any] = []
    for scale in [0.85, 1.00]:
        for detuning_peak in [1.0e6, 2.0e6, 3.0e6, 4.0e6]:
            programs.append(_shock_revival_program(scale, detuning_peak))
            manifest.append({"scale": scale, "detuning_peak": detuning_peak, "atom_count": 4})
    return manifest, programs


def _task_snapshot(session: boto3.Session, task_id: str) -> dict[str, Any]:
    client = session.client("braket")
    meta = client.get_quantum_task(quantumTaskArn=task_id)
    return {
        "task_id": task_id,
        "status": meta.get("status"),
        "device_arn": meta.get("deviceArn"),
        "createdAt": meta.get("createdAt").isoformat() if meta.get("createdAt") else None,
        "queueInfo": meta.get("queueInfo"),
        "failureReason": meta.get("failureReason"),
    }


def _submit_gate_family(
    *,
    name: str,
    device_key: str,
    shots: int,
    manifest: list[dict[str, Any]],
    circuits: list[Circuit],
    creds: dict[str, str],
    output_dir: Path,
) -> dict[str, Any]:
    spec = DEVICES[device_key]
    region = spec.region
    bucket_info = _ensure_bucket(creds, region)
    aws_session = _aws_session(creds, region)
    boto_session = _aws_boto_session(creds, region)
    device = AwsDevice(spec.arn, aws_session=aws_session)
    family_dir = _ensure_dir(output_dir / name)
    _write_json(family_dir / "manifest.json", {"entries": manifest, "shots": shots, "device": spec.arn})
    task_ids: list[str] = []
    snapshots: list[dict[str, Any]] = []
    if len(circuits) != len(manifest):
        raise RuntimeError(f"{name}: manifest/circuit length mismatch")
    for idx, (entry, circuit) in enumerate(zip(manifest, circuits, strict=True), start=1):
        task = device.run(
            circuit,
            shots=shots,
            s3_destination_folder=(bucket_info["bucket"], f"supernova/{output_dir.name}/{name}/task_{idx:03d}"),
        )
        task_id = task.id
        task_ids.append(task_id)
        snapshots.append(
            {
                **entry,
                "task_id": task_id,
                "task_index": idx,
            }
        )
    for row in snapshots:
        row["initial_status"] = _task_snapshot(boto_session, row["task_id"])
    result = {
        "name": name,
        "device_key": device_key,
        "device_arn": spec.arn,
        "bucket": bucket_info,
        "shots": shots,
        "task_count": len(task_ids),
        "task_ids": task_ids,
        "tasks": snapshots,
        "submitted_utc": _iso(),
    }
    _write_json(family_dir / "submit_results.json", result)
    return result


def _submit_ahs_family(
    *,
    name: str,
    device_key: str,
    shots: int,
    manifest: list[dict[str, Any]],
    programs: list[Any],
    creds: dict[str, str],
    output_dir: Path,
) -> dict[str, Any]:
    spec = DEVICES[device_key]
    region = spec.region
    bucket_info = _ensure_bucket(creds, region)
    aws_session = _aws_session(creds, region)
    boto_session = _aws_boto_session(creds, region)
    device = AwsDevice(spec.arn, aws_session=aws_session)
    family_dir = _ensure_dir(output_dir / name)
    _write_json(family_dir / "manifest.json", {"entries": manifest, "shots": shots, "device": spec.arn})
    task_ids: list[str] = []
    snapshots: list[dict[str, Any]] = []
    for idx, (entry, program) in enumerate(zip(manifest, programs, strict=True), start=1):
        disc = program.discretize(device)
        task = device.run(
            disc,
            shots=shots,
            s3_destination_folder=(bucket_info["bucket"], f"supernova/{output_dir.name}/{name}/task_{idx:03d}"),
        )
        task_id = task.id
        task_ids.append(task_id)
        snapshots.append({**entry, "task_id": task_id, "task_index": idx})
    for row in snapshots:
        row["initial_status"] = _task_snapshot(boto_session, row["task_id"])
    result = {
        "name": name,
        "device_key": device_key,
        "device_arn": spec.arn,
        "bucket": bucket_info,
        "shots": shots,
        "task_count": len(task_ids),
        "task_ids": task_ids,
        "tasks": snapshots,
        "submitted_utc": _iso(),
    }
    _write_json(family_dir / "submit_results.json", result)
    return result


def _build_summary(output_dir: Path, submissions: list[dict[str, Any]]) -> str:
    total_tasks = sum(item["task_count"] for item in submissions)
    lines = [
        "# Quantum Pilot Submission Summary",
        "",
        f"Submitted UTC: `{_iso()}`",
        "",
        f"Submission packet: `{output_dir.name}`",
        "",
        f"Total tasks submitted: `{total_tasks}`",
        "",
        "## Families",
        "",
    ]
    for item in submissions:
        lines.append(f"### `{item['name']}`")
        lines.append("")
        lines.append(f"- Device: `{item['device_key']}`")
        lines.append(f"- Task count: `{item['task_count']}`")
        lines.append(f"- Bucket: `{item['bucket']['bucket']}`")
        lines.append(f"- First task: `{item['task_ids'][0]}`")
        statuses = {}
        for row in item["tasks"]:
            status = row["initial_status"]["status"]
            statuses[status] = statuses.get(status, 0) + 1
        lines.append(f"- Initial status counts: `{statuses}`")
        lines.append("")
    return "\n".join(lines) + "\n"


def submit_quantum_pilots(root_dir: Path, output_dir: Path | None = None) -> dict[str, Path]:
    out_dir = output_dir or (root_dir / "outputs" / f"quantum_pilot_submit_{_stamp()}")
    _ensure_dir(out_dir)
    progress = out_dir / "progress.json"
    _write_json(out_dir / "env.json", _env_like())
    _write_json(out_dir / "git_like.json", _git_like(root_dir))
    creds, credential_meta = _load_aws_credentials(root_dir)
    if creds is None:
        raise RuntimeError("AWS credentials not found; refusing live submission")
    _write_json(out_dir / "config.json", {"created_utc": _iso(), "credential_source": credential_meta, "mode": "submit"})

    submissions: list[dict[str, Any]] = []

    _write_progress(progress, "stage01_submit_core_shell", "running")
    manifest, circuits = _core_shell_manifest()
    submissions.append(
        _submit_gate_family(
            name="core_shell_collapse_ladder",
            device_key="iqm_garnet",
            shots=1500,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )

    _write_progress(progress, "stage02_submit_neutrinosphere", "running")
    manifest, circuits = _neutrinosphere_manifest()
    submissions.append(
        _submit_gate_family(
            name="quantum_neutrinosphere_analog",
            device_key="rigetti_ankaa",
            shots=1200,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )

    _write_progress(progress, "stage03_submit_shock_revival", "running")
    manifest, programs = _shock_revival_manifest()
    submissions.append(
        _submit_ahs_family(
            name="shock_revival_threshold_simulator",
            device_key="quera_aquila",
            shots=100,
            manifest=manifest,
            programs=programs,
            creds=creds,
            output_dir=out_dir,
        )
    )

    _write_progress(progress, "stage04_submit_shell_visibility", "running")
    manifest, circuits = _shell_visibility_manifest()
    submissions.append(
        _submit_gate_family(
            name="shell_only_visibility_test",
            device_key="iqm_garnet",
            shots=1200,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )

    _write_progress(progress, "stage05_submit_cross_platform", "running")
    manifest, circuits = _cross_platform_manifest()
    submissions.append(
        _submit_gate_family(
            name="cross_platform_universality_test_iqm",
            device_key="iqm_garnet",
            shots=1000,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )
    submissions.append(
        _submit_gate_family(
            name="cross_platform_universality_test_rigetti",
            device_key="rigetti_ankaa",
            shots=1000,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )

    _write_progress(progress, "stage06_submit_collapse_rate", "running")
    manifest, circuits = _collapse_atlas_manifest()
    submissions.append(
        _submit_gate_family(
            name="collapse_rate_atlas_microgrid",
            device_key="rigetti_ankaa",
            shots=800,
            manifest=manifest,
            circuits=circuits,
            creds=creds,
            output_dir=out_dir,
        )
    )

    summary_payload = {
        "created_utc": _iso(),
        "submission_count": len(submissions),
        "total_tasks": sum(item["task_count"] for item in submissions),
        "families": submissions,
    }
    _write_json(out_dir / "submit_summary.json", summary_payload)
    _write_text(out_dir / "submit_summary.md", _build_summary(out_dir, submissions))
    _write_progress(progress, "complete", "done", {"total_tasks": summary_payload["total_tasks"]})

    return {
        "config": out_dir / "config.json",
        "submit_summary": out_dir / "submit_summary.json",
        "submit_report": out_dir / "submit_summary.md",
        "progress": progress,
    }
