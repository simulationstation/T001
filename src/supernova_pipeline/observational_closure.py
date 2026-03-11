from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .utils import ensure_dir, normalize_01, utc_stamp, write_json


CLASS_EXPORT = "EXPORT_FAILURE_LIKE"
CLASS_DUST = "DUST_SURVIVOR_LIKE"
CLASS_SYSTEMATIC = "SYSTEMATIC_LIKE"
CLASS_UNRESOLVED = "VARIABLE_OR_UNRESOLVED"


def _safe_json_list(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, str) and value.strip():
        try:
            payload = json.loads(value)
        except json.JSONDecodeError:
            return []
        if isinstance(payload, list):
            return [str(item) for item in payload]
    return []


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def _weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not np.any(good):
        return float("nan"), float("nan")
    weights = 1.0 / np.square(errors[good])
    mean = float(np.sum(values[good] * weights) / np.sum(weights))
    err = float(np.sqrt(1.0 / np.sum(weights)))
    return mean, err


def _norm_series(
    series: pd.Series,
    *,
    log: bool = False,
    clip_lower: float | None = 0.0,
    upper_quantile: float = 0.95,
    default: float = 0.0,
) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    if clip_lower is not None:
        values = values.clip(lower=clip_lower)
    if log:
        values = np.log10(np.clip(values, 1e-9, np.inf))
    finite = np.isfinite(values)
    if not np.any(finite):
        return pd.Series(default, index=series.index, dtype=float)
    hi = float(np.nanquantile(values[finite], upper_quantile))
    lo = float(np.nanmin(values[finite]))
    clipped = values.clip(lower=lo, upper=hi)
    out = normalize_01(pd.Series(clipped, index=series.index, dtype=float))
    return out.fillna(default)


def _load_optional_csv(path: Path, *, rename: dict[str, str] | None = None, usecols: list[str] | None = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=usecols or [])
    data = pd.read_csv(path)
    if usecols is not None:
        missing = [column for column in usecols if column not in data.columns]
        for column in missing:
            data[column] = np.nan
        data = data[usecols]
    if rename:
        data = data.rename(columns=rename)
    return data


def _summarize_followup_photometry(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "followup_exists": False,
            "followup_measured_epochs": 0,
            "followup_measured_filters": 0,
            "followup_common_filters": 0,
            "followup_median_filter_ratio": float("nan"),
            "followup_min_filter_ratio": float("nan"),
            "followup_max_filter_ratio": float("nan"),
            "followup_persistent_fade_fraction": float("nan"),
            "followup_rebrightening_flag": False,
        }

    phot = pd.read_parquet(path)
    measured = phot[
        (phot["status"].astype(str) == "MEASURED")
        & np.isfinite(pd.to_numeric(phot["flux_jy"], errors="coerce"))
        & np.isfinite(pd.to_numeric(phot["err_jy"], errors="coerce"))
        & (pd.to_numeric(phot["err_jy"], errors="coerce") > 0)
    ].copy()
    if measured.empty:
        return {
            "followup_exists": True,
            "followup_measured_epochs": 0,
            "followup_measured_filters": 0,
            "followup_common_filters": 0,
            "followup_median_filter_ratio": float("nan"),
            "followup_min_filter_ratio": float("nan"),
            "followup_max_filter_ratio": float("nan"),
            "followup_persistent_fade_fraction": float("nan"),
            "followup_rebrightening_flag": False,
        }

    filter_rows: list[dict[str, Any]] = []
    for filter_name, filter_df in measured.groupby("filter_name", dropna=False):
        pre_mean, _ = _weighted_mean(
            filter_df.loc[filter_df["stage"].astype(str) == "pre", "flux_jy"].to_numpy(dtype=float),
            filter_df.loc[filter_df["stage"].astype(str) == "pre", "err_jy"].to_numpy(dtype=float),
        )
        post_mean, _ = _weighted_mean(
            filter_df.loc[filter_df["stage"].astype(str) == "post", "flux_jy"].to_numpy(dtype=float),
            filter_df.loc[filter_df["stage"].astype(str) == "post", "err_jy"].to_numpy(dtype=float),
        )
        ratio = float(post_mean / pre_mean) if math.isfinite(pre_mean) and abs(pre_mean) > 0 and math.isfinite(post_mean) else float("nan")
        filter_rows.append(
            {
                "filter_name": str(filter_name),
                "pre_mean": pre_mean,
                "post_mean": post_mean,
                "ratio": ratio,
            }
        )

    filter_table = pd.DataFrame(filter_rows)
    common = filter_table[
        np.isfinite(filter_table["pre_mean"])
        & np.isfinite(filter_table["post_mean"])
        & (np.abs(filter_table["pre_mean"]) > 0)
    ].copy()
    ratios = common["ratio"].to_numpy(dtype=float)
    persistent = float(np.mean(ratios <= 0.2)) if len(ratios) else float("nan")
    rebrightening = bool(np.any(ratios > 0.5)) if len(ratios) else False
    return {
        "followup_exists": True,
        "followup_measured_epochs": int(len(measured)),
        "followup_measured_filters": int(measured["filter_name"].nunique()),
        "followup_common_filters": int(len(common)),
        "followup_median_filter_ratio": float(np.nanmedian(ratios)) if len(ratios) else float("nan"),
        "followup_min_filter_ratio": float(np.nanmin(ratios)) if len(ratios) else float("nan"),
        "followup_max_filter_ratio": float(np.nanmax(ratios)) if len(ratios) else float("nan"),
        "followup_persistent_fade_fraction": persistent,
        "followup_rebrightening_flag": rebrightening,
    }


def _load_followup_table(root_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted((root_dir / "sed").glob("*/followup_summary.json")):
        payload = json.loads(path.read_text())
        phot_path = Path(str(payload.get("photometry_path", "")))
        metrics = _summarize_followup_photometry(phot_path) if phot_path else _summarize_followup_photometry(Path(""))
        rows.append(
            {
                "candidate_id": str(payload.get("candidate_id", path.parent.name)),
                "followup_created_utc": str(payload.get("created_utc", "")),
                "followup_lightcurve_path": str(payload.get("lightcurve_path", "")),
                "followup_sed_plot_path": str(payload.get("sed_plot_path", "")),
                "followup_summary_path": str(path),
                "followup_photometry_path": str(payload.get("photometry_path", "")),
                "followup_n_observations_considered": int(payload.get("n_observations_considered", 0) or 0),
                "followup_n_measurements": int(payload.get("n_measurements", 0) or 0),
                "followup_n_stage_points": int(payload.get("n_stage_points", 0) or 0),
                "followup_pre_fit_ok": bool(payload.get("fits", {}).get("pre", {}).get("fit_ok", False)),
                "followup_post_fit_ok": bool(payload.get("fits", {}).get("post", {}).get("fit_ok", False)),
                "followup_pre_fit_points": int(payload.get("fits", {}).get("pre", {}).get("n_points", 0) or 0),
                "followup_post_fit_points": int(payload.get("fits", {}).get("post", {}).get("n_points", 0) or 0),
                **metrics,
            }
        )
    return pd.DataFrame(rows)


def _build_design_markdown() -> str:
    return """# Observational Closure Test

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
"""


def _reason_codes(row: pd.Series) -> list[str]:
    reasons: list[str] = []
    if float(row.get("fade_fraction", 0.0)) >= 0.95:
        reasons.append("NEAR_TOTAL_DISAPPEARANCE")
    elif float(row.get("fade_fraction", 0.0)) >= 0.70:
        reasons.append("STRONG_FADE")
    if bool(row.get("has_post_lt_5pct", False)):
        reasons.append("SUB_5PCT_POST")
    if float(row.get("fade_sigma", 0.0)) >= 10.0:
        reasons.append("HIGH_SIGNIFICANCE_FADE")
    if float(row.get("baseline_days", 0.0)) >= 365.0:
        reasons.append("LONG_BASELINE")
    if bool(row.get("followup_exists", False)):
        reasons.append("ALL_EPOCH_FOLLOWUP")
    if float(row.get("followup_persistent_fade_fraction", 0.0)) >= 0.67:
        reasons.append("MULTIFILTER_PERSISTENCE")
    if float(row.get("local_unique_effective", 0.0)) >= 3.0:
        reasons.append("LOCALLY_UNIQUE")
    if bool(row.get("mid_ir_survivor_flag", False)):
        reasons.append("MID_IR_SURVIVOR")
    if bool(row.get("dusty_color_flag", False)):
        reasons.append("DUSTY_COLOR")
    if float(row.get("blend_risk", 0.0)) >= 0.5 or bool(row.get("has_blend_warning", False)):
        reasons.append("BLEND_RISK")
    if float(row.get("crowding_index", 0.0)) >= 5.0 or bool(row.get("has_crowded_warning", False)):
        reasons.append("CROWDED_FIELD")
    if float(row.get("cross_reduction_agreement", 1.0)) < 0.8:
        reasons.append("CROSS_REDUCTION_TENSION")
    if bool(row.get("followup_rebrightening_flag", False)):
        reasons.append("POST_EPOCH_REBRIGHTENING")
    if not bool(row.get("followup_exists", False)):
        reasons.append("NO_ALL_EPOCH_FOLLOWUP")
    return reasons


def _classify_candidate(row: pd.Series) -> tuple[str, str]:
    export_score = float(row["export_failure_score"])
    dust_score = float(row["dust_survivor_score"])
    systematic_score = float(row["systematic_like_score"])
    unresolved_score = float(row["unresolved_score"])
    top_class, top_score = max(
        [
            (CLASS_EXPORT, export_score),
            (CLASS_DUST, dust_score),
            (CLASS_SYSTEMATIC, systematic_score),
            (CLASS_UNRESOLVED, unresolved_score),
        ],
        key=lambda item: item[1],
    )
    ordered = sorted([export_score, dust_score, systematic_score, unresolved_score], reverse=True)
    margin = float(ordered[0] - ordered[1]) if len(ordered) > 1 else 0.0
    completeness = float(row["closure_data_completeness"])
    if margin >= 0.20 and completeness >= 0.55:
        confidence = "HIGH"
    elif margin >= 0.10 and completeness >= 0.35:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    if dust_score >= 0.60 and dust_score >= export_score + 0.05:
        return CLASS_DUST, confidence
    if systematic_score >= 0.60 and systematic_score >= export_score + 0.05:
        return CLASS_SYSTEMATIC, confidence
    if export_score >= 0.60 and float(row["survivor_evidence_score"]) <= 0.35 and float(row["systematic_penalty_score"]) <= 0.45:
        return CLASS_EXPORT, confidence
    if top_class == CLASS_EXPORT and export_score >= 0.50 and float(row["systematic_penalty_score"]) <= 0.55:
        return CLASS_EXPORT, confidence
    return top_class if top_class != CLASS_EXPORT else CLASS_UNRESOLVED, confidence


def _candidate_packet_markdown(row: pd.Series) -> str:
    lines = [
        f"# {row['candidate_id']}",
        "",
        f"- Galaxy: `{row['galaxy_name']}`",
        f"- Ledger status: `{row['status']}`",
        f"- Closure class: `{row['closure_class']}`",
        f"- Confidence: `{row['closure_confidence']}`",
        f"- Export score: `{row['export_failure_score']:.3f}`",
        f"- Dust-survivor score: `{row['dust_survivor_score']:.3f}`",
        f"- Systematic score: `{row['systematic_like_score']:.3f}`",
        f"- Unresolved score: `{row['unresolved_score']:.3f}`",
        "",
        "## Key Signals",
        f"- Fade fraction: `{row['fade_fraction']:.4f}`",
        f"- Fade sigma: `{row['fade_sigma']:.2f}`",
        f"- Baseline days: `{row['baseline_days']:.2f}`",
        f"- Cross-reduction agreement: `{row['cross_reduction_agreement']:.3f}`",
        f"- Blend risk: `{row['blend_risk']:.3f}`",
        f"- Crowding index: `{row['crowding_index']:.3f}`",
        f"- Mid-IR survivor flag: `{bool(row['mid_ir_survivor_flag'])}`",
        f"- Dusty color flag: `{bool(row['dusty_color_flag'])}`",
        f"- Follow-up exists: `{bool(row['followup_exists'])}`",
        f"- Follow-up persistent fade fraction: `{float(row['followup_persistent_fade_fraction']) if pd.notna(row['followup_persistent_fade_fraction']) else float('nan'):.3f}`",
        "",
        "## Reason Codes",
    ]
    for code in _safe_json_list(row["closure_reason_codes_json"]):
        lines.append(f"- `{code}`")
    return "\n".join(lines) + "\n"


def _prepare_frame(root_dir: Path, candidate_ledger_path: Path) -> pd.DataFrame:
    ledger = pd.read_parquet(candidate_ledger_path).copy()
    ext = _load_optional_csv(
        root_dir / "extensions" / "combined" / "candidate_extension_ranking.csv",
        usecols=["candidate_id", "extended_rank_score", "extension_verdict", "population_prior_score", "anomaly_score"],
    )
    mid_ir = _load_optional_csv(
        root_dir / "extensions" / "mid_ir" / "mid_ir_crossmatch.csv",
        rename={"status": "mid_ir_match_status"},
        usecols=["candidate_id", "status", "match_sep_arcsec", "mid_ir_survivor_flag", "dusty_color_flag", "w1_w2"],
    )
    neighborhood = _load_optional_csv(
        root_dir / "extensions" / "neighborhood" / "neighborhood_scores.csv",
        usecols=["candidate_id", "local_unique_sigma", "local_systematic_risk", "local_coherent_trend", "local_neighbor_count", "local_verdict"],
    )
    scene = _load_optional_csv(
        root_dir / "extensions" / "scene_model" / "scene_model_metrics.csv",
        usecols=[
            "candidate_id",
            "scene_mode",
            "n_scene_epochs",
            "scene_fade_fraction",
            "scene_fade_sigma",
            "scene_unique_fade_sigma",
            "scene_local_comparators",
            "scene_monotonic_score",
            "scene_residual_center",
            "scene_residual_std",
        ],
    )
    followup = _load_followup_table(root_dir)

    merged = ledger.merge(ext, on="candidate_id", how="left")
    merged = merged.merge(mid_ir, on="candidate_id", how="left")
    merged = merged.merge(neighborhood, on="candidate_id", how="left")
    merged = merged.merge(scene, on="candidate_id", how="left")
    merged = merged.merge(followup, on="candidate_id", how="left")

    merged["reason_codes"] = merged["reason_codes_json"].map(_safe_json_list)
    merged["warning_flags"] = merged["warning_flags_json"].map(_safe_json_list)
    merged["has_post_lt_5pct"] = merged["reason_codes"].map(lambda values: "POST_LT_5PCT" in values)
    merged["has_robust_fade_sig"] = merged["reason_codes"].map(lambda values: "ROBUST_FADE_SIG" in values)
    merged["has_blend_warning"] = merged["warning_flags"].map(lambda values: "BLEND_RISK" in values)
    merged["has_crowded_warning"] = merged["warning_flags"].map(lambda values: "CROWDED" in values)
    merged["has_edge_warning"] = merged["warning_flags"].map(lambda values: "EDGE_NEAR_MASK" in values)
    merged["baseline_days"] = (pd.to_numeric(merged["post_mjd"], errors="coerce") - pd.to_numeric(merged["pre_mjd"], errors="coerce")).abs()
    merged["mid_ir_survivor_flag"] = merged["mid_ir_survivor_flag"].map(_as_bool)
    merged["dusty_color_flag"] = merged["dusty_color_flag"].map(_as_bool)
    merged["followup_exists"] = merged["followup_exists"].map(_as_bool)
    merged["followup_rebrightening_flag"] = merged["followup_rebrightening_flag"].map(_as_bool)
    merged["scene_residual_snr"] = np.divide(
        np.abs(pd.to_numeric(merged["scene_residual_center"], errors="coerce")),
        np.clip(pd.to_numeric(merged["scene_residual_std"], errors="coerce"), 1e-6, np.inf),
        out=np.full(len(merged), np.nan, dtype=float),
        where=np.isfinite(pd.to_numeric(merged["scene_residual_std"], errors="coerce")),
    )
    merged["local_unique_effective"] = pd.to_numeric(merged["local_unique_sigma"], errors="coerce").combine_first(
        pd.to_numeric(merged["scene_unique_fade_sigma"], errors="coerce")
    )
    merged["local_systematic_effective"] = pd.to_numeric(merged["local_systematic_risk"], errors="coerce")
    fallback_systematic = 1.0 / (1.0 + pd.to_numeric(merged["local_unique_effective"], errors="coerce").clip(lower=0.0))
    merged["local_systematic_effective"] = merged["local_systematic_effective"].combine_first(fallback_systematic)

    merged["fade_frac_score"] = normalize_01(pd.to_numeric(merged["fade_fraction"], errors="coerce").clip(lower=0.0, upper=1.5)).fillna(0.0)
    merged["fade_sigma_score"] = _norm_series(merged["fade_sigma"], log=True, default=0.0)
    merged["scene_fade_frac_score"] = normalize_01(pd.to_numeric(merged["scene_fade_fraction"], errors="coerce").clip(lower=0.0, upper=1.5)).fillna(0.0)
    merged["scene_fade_sigma_score"] = _norm_series(merged["scene_fade_sigma"], log=True, default=0.0)
    merged["baseline_score"] = _norm_series(merged["baseline_days"], log=True, default=0.0)
    merged["agreement_score"] = pd.to_numeric(merged["cross_reduction_agreement"], errors="coerce").clip(lower=0.0, upper=1.0).fillna(0.0)
    merged["local_unique_score"] = _norm_series(merged["local_unique_effective"], log=True, default=0.5)
    merged["local_systematic_penalty"] = pd.to_numeric(merged["local_systematic_effective"], errors="coerce").clip(lower=0.0, upper=1.0).fillna(0.30)
    merged["scene_comparator_support"] = _norm_series(merged["scene_local_comparators"], default=0.40)
    monotonic = pd.to_numeric(merged["scene_monotonic_score"], errors="coerce")
    monotonic_default = pd.Series(np.where(merged["followup_exists"], 0.50, 0.30), index=merged.index, dtype=float)
    merged["monotonic_support"] = ((monotonic + 1.0) / 2.0).clip(lower=0.0, upper=1.0).fillna(monotonic_default)
    merged["followup_measurement_score"] = _norm_series(merged["followup_n_measurements"], log=True, default=0.0)
    merged["followup_stage_score"] = _norm_series(merged["followup_n_stage_points"], default=0.0)
    merged["followup_common_filter_score"] = _norm_series(merged["followup_common_filters"], default=0.0)
    merged["followup_persistent_fade_score"] = pd.to_numeric(merged["followup_persistent_fade_fraction"], errors="coerce").clip(lower=0.0, upper=1.0).fillna(0.0)
    merged["crowding_penalty"] = (pd.to_numeric(merged["crowding_index"], errors="coerce") / 10.0).clip(lower=0.0, upper=1.0).fillna(0.0)
    merged["host_bg_penalty_norm"] = _norm_series(merged["host_bg_penalty"], default=0.0)
    merged["astrometric_penalty"] = normalize_01(pd.to_numeric(merged["astrometric_residual"], errors="coerce").clip(lower=0.0)).fillna(0.0)
    merged["agreement_penalty"] = 1.0 - merged["agreement_score"]
    merged["blend_penalty"] = pd.to_numeric(merged["blend_risk"], errors="coerce").clip(lower=0.0, upper=1.0).fillna(0.0)
    merged["scene_residual_penalty"] = _norm_series(merged["scene_residual_snr"], log=True, default=0.30)
    merged["warning_penalty"] = (
        0.35 * merged["has_blend_warning"].astype(float)
        + 0.25 * merged["has_crowded_warning"].astype(float)
        + 0.15 * merged["has_edge_warning"].astype(float)
    ).clip(upper=1.0)
    merged["survivor_evidence_score"] = (
        0.80 * merged["mid_ir_survivor_flag"].astype(float)
        + 0.20 * merged["dusty_color_flag"].astype(float)
    ).clip(upper=1.0)

    merged["fade_evidence_score"] = (
        0.23 * merged["fade_frac_score"]
        + 0.20 * merged["fade_sigma_score"]
        + 0.12 * merged["scene_fade_frac_score"]
        + 0.12 * merged["scene_fade_sigma_score"]
        + 0.11 * merged["agreement_score"]
        + 0.08 * merged["has_post_lt_5pct"].astype(float)
        + 0.05 * merged["has_robust_fade_sig"].astype(float)
        + 0.05 * merged["baseline_score"]
        + 0.04 * merged["monotonic_support"]
    ).clip(0.0, 1.0)
    merged["uniqueness_evidence_score"] = (
        0.45 * merged["local_unique_score"]
        + 0.30 * (1.0 - merged["local_systematic_penalty"])
        + 0.15 * merged["scene_comparator_support"]
        + 0.10 * merged["monotonic_support"]
    ).clip(0.0, 1.0)
    merged["followup_support_score"] = np.where(
        merged["followup_exists"],
        (
            0.30 * merged["followup_measurement_score"]
            + 0.20 * merged["followup_stage_score"]
            + 0.20 * merged["followup_common_filter_score"]
            + 0.20 * merged["followup_persistent_fade_score"]
            + 0.10 * (1.0 - merged["followup_rebrightening_flag"].astype(float))
        ).clip(0.0, 1.0),
        0.0,
    )
    merged["systematic_penalty_score"] = (
        0.28 * merged["blend_penalty"]
        + 0.18 * merged["crowding_penalty"]
        + 0.10 * merged["host_bg_penalty_norm"]
        + 0.14 * merged["astrometric_penalty"]
        + 0.12 * merged["local_systematic_penalty"]
        + 0.10 * merged["scene_residual_penalty"]
        + 0.08 * merged["agreement_penalty"]
        + 0.10 * merged["warning_penalty"]
    ).clip(0.0, 1.0)
    merged["closure_data_completeness"] = (
        0.35 * merged["followup_exists"].astype(float)
        + 0.20 * merged["baseline_score"]
        + 0.20 * merged["scene_comparator_support"]
        + 0.15 * np.isfinite(pd.to_numeric(merged["scene_fade_fraction"], errors="coerce")).astype(float)
        + 0.10 * np.isfinite(pd.to_numeric(merged["local_unique_effective"], errors="coerce")).astype(float)
    ).clip(0.0, 1.0)

    merged["export_failure_score"] = (
        0.42 * merged["fade_evidence_score"]
        + 0.22 * merged["uniqueness_evidence_score"]
        + 0.16 * merged["followup_support_score"]
        + 0.10 * merged["baseline_score"]
        + 0.10 * merged["agreement_score"]
        - 0.35 * merged["systematic_penalty_score"]
        - 0.55 * merged["survivor_evidence_score"]
    ).clip(0.0, 1.0)
    merged["dust_survivor_score"] = (
        0.70 * merged["survivor_evidence_score"]
        + 0.15 * (1.0 - merged["has_post_lt_5pct"].astype(float))
        + 0.15 * merged["followup_rebrightening_flag"].astype(float)
    ).clip(0.0, 1.0)
    merged["systematic_like_score"] = (
        0.60 * merged["systematic_penalty_score"]
        + 0.15 * merged["agreement_penalty"]
        + 0.15 * merged["warning_penalty"]
        + 0.10 * (1.0 - merged["followup_support_score"])
    ).clip(0.0, 1.0)
    merged["unresolved_score"] = (
        0.40 * (1.0 - merged["followup_support_score"])
        + 0.25 * (1.0 - merged["uniqueness_evidence_score"])
        + 0.20 * (1.0 - merged["closure_data_completeness"])
        + 0.15 * (1.0 - merged["fade_evidence_score"])
    ).clip(0.0, 1.0)
    return merged


def run_observational_closure(
    *,
    root_dir: Path,
    candidate_ledger_path: Path,
) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "closure")
    packets_dir = ensure_dir(out_dir / "packets")

    closure = _prepare_frame(root_dir=root_dir, candidate_ledger_path=candidate_ledger_path)
    closure["closure_reason_codes_json"] = ""
    closure["closure_class"] = ""
    closure["closure_confidence"] = ""

    for index, row in closure.iterrows():
        candidate_class, confidence = _classify_candidate(row)
        closure.at[index, "closure_class"] = candidate_class
        closure.at[index, "closure_confidence"] = confidence
        closure.at[index, "closure_reason_codes_json"] = json.dumps(_reason_codes(closure.loc[index]))

    closure = closure.sort_values(
        by=["closure_class", "export_failure_score", "fade_sigma"],
        ascending=[True, False, False],
    ).reset_index(drop=True)
    closure["closure_rank"] = np.arange(1, len(closure) + 1, dtype=int)

    full_path = out_dir / "closure_all_features.parquet"
    full_csv = out_dir / "closure_all_features.csv"
    class_path = out_dir / "closure_classification.parquet"
    class_csv = out_dir / "closure_classification.csv"
    summary_path = out_dir / "closure_summary.json"
    method_path = out_dir / "closure_method.md"
    report_path = out_dir / "closure_report.md"
    manifest_path = out_dir / "closure_manifest.json"

    closure.to_parquet(full_path, index=False)
    closure.to_csv(full_csv, index=False)
    keep_cols = [
        "closure_rank",
        "candidate_id",
        "galaxy_name",
        "status",
        "closure_class",
        "closure_confidence",
        "export_failure_score",
        "dust_survivor_score",
        "systematic_like_score",
        "unresolved_score",
        "fade_fraction",
        "fade_sigma",
        "baseline_days",
        "mid_ir_survivor_flag",
        "dusty_color_flag",
        "followup_exists",
        "followup_persistent_fade_fraction",
        "closure_reason_codes_json",
    ]
    closure[keep_cols].to_parquet(class_path, index=False)
    closure[keep_cols].to_csv(class_csv, index=False)

    class_counts = closure["closure_class"].value_counts(dropna=False).to_dict()
    summary = {
        "created_utc": utc_stamp(),
        "n_candidates": int(len(closure)),
        "class_counts": {str(key): int(value) for key, value in class_counts.items()},
        "n_export_failure_like": int((closure["closure_class"] == CLASS_EXPORT).sum()),
        "n_dust_survivor_like": int((closure["closure_class"] == CLASS_DUST).sum()),
        "n_systematic_like": int((closure["closure_class"] == CLASS_SYSTEMATIC).sum()),
        "n_variable_or_unresolved": int((closure["closure_class"] == CLASS_UNRESOLVED).sum()),
        "top_export_candidates": closure.loc[closure["closure_class"] == CLASS_EXPORT, ["candidate_id", "galaxy_name", "export_failure_score"]]
        .sort_values("export_failure_score", ascending=False)
        .head(10)
        .to_dict(orient="records"),
        "top_dust_candidates": closure.loc[closure["closure_class"] == CLASS_DUST, ["candidate_id", "galaxy_name", "dust_survivor_score"]]
        .sort_values("dust_survivor_score", ascending=False)
        .head(10)
        .to_dict(orient="records"),
        "top_systematic_candidates": closure.loc[closure["closure_class"] == CLASS_SYSTEMATIC, ["candidate_id", "galaxy_name", "systematic_like_score"]]
        .sort_values("systematic_like_score", ascending=False)
        .head(10)
        .to_dict(orient="records"),
    }
    write_json(summary_path, summary)
    method_path.write_text(_build_design_markdown())

    report_lines = [
        "# Observational Closure Report",
        "",
        "## Summary",
        f"- Candidates adjudicated: `{len(closure)}`",
        f"- `EXPORT_FAILURE_LIKE`: `{summary['n_export_failure_like']}`",
        f"- `DUST_SURVIVOR_LIKE`: `{summary['n_dust_survivor_like']}`",
        f"- `SYSTEMATIC_LIKE`: `{summary['n_systematic_like']}`",
        f"- `VARIABLE_OR_UNRESOLVED`: `{summary['n_variable_or_unresolved']}`",
        "",
        "## Top Export-Failure-Like Candidates",
    ]
    export_rows = closure.loc[closure["closure_class"] == CLASS_EXPORT, ["candidate_id", "galaxy_name", "export_failure_score", "closure_confidence"]].sort_values("export_failure_score", ascending=False)
    if export_rows.empty:
        report_lines.append("- None under the current rubric.")
    else:
        for _, row in export_rows.head(10).iterrows():
            report_lines.append(
                f"- `{row['candidate_id']}` in `{row['galaxy_name']}` with score `{row['export_failure_score']:.3f}` and confidence `{row['closure_confidence']}`"
            )
    report_lines.extend(
        [
            "",
            "## Top Competing Explanations",
            "### Dust-Survivor-Like",
        ]
    )
    dust_rows = closure.loc[closure["closure_class"] == CLASS_DUST, ["candidate_id", "galaxy_name", "dust_survivor_score", "closure_confidence"]].sort_values("dust_survivor_score", ascending=False)
    if dust_rows.empty:
        report_lines.append("- None under the current rubric.")
    else:
        for _, row in dust_rows.head(10).iterrows():
            report_lines.append(
                f"- `{row['candidate_id']}` in `{row['galaxy_name']}` with score `{row['dust_survivor_score']:.3f}` and confidence `{row['closure_confidence']}`"
            )
    report_lines.extend(["", "### Systematic-Like"])
    sys_rows = closure.loc[closure["closure_class"] == CLASS_SYSTEMATIC, ["candidate_id", "galaxy_name", "systematic_like_score", "closure_confidence"]].sort_values("systematic_like_score", ascending=False)
    if sys_rows.empty:
        report_lines.append("- None under the current rubric.")
    else:
        for _, row in sys_rows.head(10).iterrows():
            report_lines.append(
                f"- `{row['candidate_id']}` in `{row['galaxy_name']}` with score `{row['systematic_like_score']:.3f}` and confidence `{row['closure_confidence']}`"
            )
    report_path.write_text("\n".join(report_lines) + "\n")

    packet_paths: dict[str, str] = {}
    for _, row in closure.iterrows():
        candidate_dir = ensure_dir(packets_dir / str(row["candidate_id"]))
        decision_path = candidate_dir / "decision.md"
        json_path = candidate_dir / "decision.json"
        decision_path.write_text(_candidate_packet_markdown(row))
        packet_payload = {
            "candidate_id": str(row["candidate_id"]),
            "galaxy_name": str(row["galaxy_name"]),
            "closure_class": str(row["closure_class"]),
            "closure_confidence": str(row["closure_confidence"]),
            "export_failure_score": float(row["export_failure_score"]),
            "dust_survivor_score": float(row["dust_survivor_score"]),
            "systematic_like_score": float(row["systematic_like_score"]),
            "unresolved_score": float(row["unresolved_score"]),
            "closure_reason_codes": _safe_json_list(row["closure_reason_codes_json"]),
            "followup_summary_path": str(row.get("followup_summary_path", "")),
            "followup_lightcurve_path": str(row.get("followup_lightcurve_path", "")),
            "followup_sed_plot_path": str(row.get("followup_sed_plot_path", "")),
        }
        write_json(json_path, packet_payload)
        packet_paths[str(row["candidate_id"])] = str(json_path)

    write_json(
        manifest_path,
        {
            "created_utc": utc_stamp(),
            "method_path": str(method_path),
            "report_path": str(report_path),
            "summary_path": str(summary_path),
            "classification_csv": str(class_csv),
            "packets": packet_paths,
        },
    )
    return {
        "closure_all_features": full_path,
        "closure_all_features_csv": full_csv,
        "closure_classification": class_path,
        "closure_classification_csv": class_csv,
        "closure_summary": summary_path,
        "closure_method": method_path,
        "closure_report": report_path,
        "closure_manifest": manifest_path,
    }
