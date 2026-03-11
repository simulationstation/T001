from __future__ import annotations

import json
import math
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import center_of_mass, distance_transform_edt, gaussian_filter, label, maximum_filter, shift as ndi_shift
from scipy.spatial import cKDTree

from .candidate_followup import _candidate_cutout
from .pixel_search import (
    ScienceImage,
    _aperture_fluxes,
    _background_stats,
    _detect_sources,
    _estimate_scale_factor,
    _finite_mask,
    _load_science_image,
    _pair_stability_metrics,
    _product_query,
    _save_candidate_figure,
    _select_best_product,
    _select_pixel_pairs,
    _download_product,
    _cached_products,
)
from .utils import append_log_line, ensure_dir, env_like, git_like, json_list, read_json, utc_stamp, write_json, write_progress, write_text


DETECTION_COLUMNS = [
    "candidate_id",
    "pair_id",
    "galaxy_name",
    "ra_deg",
    "dec_deg",
    "event_sign",
    "status",
    "diff_sigma",
    "diff_sigma_min",
    "post_to_pre_ratio",
    "pre_flux",
    "post_flux",
    "pre_snr",
    "post_snr",
    "blend_risk",
    "crowding_index",
    "edge_distance_px",
    "registration_dx_px",
    "registration_dy_px",
    "registration_residual_px",
    "registration_n_matches",
    "pair_stable_fraction",
    "pair_extreme_fade_fraction",
    "pre_obsid",
    "post_obsid",
    "pre_obs_id",
    "post_obs_id",
    "pre_obs_collection",
    "post_obs_collection",
    "pre_filter",
    "post_filter",
    "pre_instrument",
    "post_instrument",
    "pre_mjd",
    "post_mjd",
    "cutout_path",
    "reason_codes_json",
    "warning_flags_json",
    "provenance_json",
]

MEASUREMENT_COLUMNS = [
    "x",
    "y",
    "ra_deg",
    "dec_deg",
    "flux",
    "peak",
    "event_sign",
    "diff_flux",
    "diff_flux_std",
    "diff_sigma",
    "diff_sigma_min",
    "diff_sigma_max",
    "flux_pre",
    "flux_post",
    "err_pre",
    "err_post",
    "post_to_pre_ratio",
    "pre_snr",
    "post_snr",
    "crowding_index",
    "blend_risk",
    "edge_distance_px",
    "pair_is_stable",
    "registration_residual_px",
    "registration_dx_px",
    "registration_dy_px",
    "registration_n_matches",
    "pair_stable_fraction",
    "pair_extreme_fade_fraction",
    "pair_id",
    "galaxy_name",
    "obs_id_1",
    "obs_id_2",
    "filter_1",
    "filter_2",
    "pre_mjd",
    "post_mjd",
]


@dataclass(slots=True)
class ScienceFootprint:
    path: Path
    wcs: WCS
    shape: tuple[int, int]
    pixel_scale_arcsec: float
    extname: str


def _iso() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _log(out_dir: Path, message: str) -> None:
    append_log_line(out_dir / "run.log", f"[{_iso()}] {message}")


def _checkpoint_dirs(out_dir: Path) -> dict[str, Path]:
    root = ensure_dir(out_dir / "checkpoints")
    return {
        "root": root,
        "pair_summaries": ensure_dir(root / "pair_summaries"),
        "candidates": ensure_dir(root / "candidates"),
        "measurements": ensure_dir(root / "measurements"),
    }


def _pair_checkpoint_paths(out_dir: Path, pair_id: str) -> dict[str, Path]:
    slug = _slugify(pair_id)
    dirs = _checkpoint_dirs(out_dir)
    return {
        "pair_summary": dirs["pair_summaries"] / f"{slug}.json",
        "candidates": dirs["candidates"] / f"{slug}.parquet",
        "measurements": dirs["measurements"] / f"{slug}.parquet",
        "packet_root": out_dir / "packets" / slug,
    }


def _empty_detections_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=DETECTION_COLUMNS)


def _empty_measurements_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=MEASUREMENT_COLUMNS)


def _write_run_metadata(out_dir: Path, root_dir: Path, config: dict[str, Any], *, resume: bool) -> None:
    if resume and (out_dir / "config.json").exists():
        return
    write_json(out_dir / "config.json", config)
    write_json(out_dir / "env.json", env_like())
    write_json(out_dir / "git_like.json", git_like(root_dir))


def _write_pair_checkpoint(
    *,
    out_dir: Path,
    pair_summary: dict[str, Any],
    candidates: list[dict[str, Any]],
    detection_table: pd.DataFrame,
) -> None:
    paths = _pair_checkpoint_paths(out_dir, str(pair_summary["pair_id"]))
    candidate_df = pd.DataFrame(candidates)
    if candidate_df.empty:
        candidate_df = _empty_detections_frame()
    if detection_table.empty and len(detection_table.columns) == 0:
        detection_table = _empty_measurements_frame()
    candidate_df.to_parquet(paths["candidates"], index=False)
    detection_table.to_parquet(paths["measurements"], index=False)
    write_json(paths["pair_summary"], pair_summary)


def _load_existing_checkpoint_rows(out_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int, set[str]]:
    dirs = _checkpoint_dirs(out_dir)
    pair_rows: list[dict[str, Any]] = []
    detections: list[dict[str, Any]] = []
    measurement_count = 0
    completed_pair_ids: set[str] = set()
    for path in sorted(dirs["pair_summaries"].glob("*.json")):
        payload = read_json(path)
        pair_rows.append(payload)
        if payload.get("status") in {"SCANNED", "FAILED", "SKIPPED"}:
            completed_pair_ids.add(str(payload["pair_id"]))
        cand_path = dirs["candidates"] / f"{path.stem}.parquet"
        if cand_path.exists():
            cand_df = pd.read_parquet(cand_path)
            if not cand_df.empty:
                detections.extend(cand_df.to_dict(orient="records"))
        meas_path = dirs["measurements"] / f"{path.stem}.parquet"
        if meas_path.exists():
            measurement_count += int(len(pd.read_parquet(meas_path)))
    return pair_rows, detections, measurement_count, completed_pair_ids


def _build_detection_outputs(
    *,
    out_dir: Path,
    pair_rows: list[dict[str, Any]],
    detections: list[dict[str, Any]],
    n_detection_measurements: int,
    sign_mode: str,
    total_pairs: int,
) -> dict[str, Any]:
    pair_summary_df = pd.DataFrame(pair_rows)
    detections_df = pd.DataFrame(detections)
    if detections_df.empty:
        detections_df = _empty_detections_frame()
    else:
        detections_df = detections_df.sort_values(
            ["event_sign", "status", "diff_sigma_min"],
            ascending=[True, True, False],
        ).reset_index(drop=True)
    fade_df = detections_df[detections_df["event_sign"].astype(str).eq("fade")].copy()
    pair_summary_df.to_parquet(out_dir / "pair_summary.parquet", index=False)
    pair_summary_df.to_csv(out_dir / "pair_summary.csv", index=False)
    detections_df.to_parquet(out_dir / "detections.parquet", index=False)
    detections_df.to_csv(out_dir / "detections.csv", index=False)
    fade_df.to_parquet(out_dir / "fade_candidates.parquet", index=False)
    fade_df.to_csv(out_dir / "fade_candidates.csv", index=False)
    summary = {
        "created_utc": utc_stamp(),
        "sign_mode": sign_mode,
        "n_pairs_selected": int(total_pairs),
        "n_pairs_completed": int(len(pair_rows)),
        "n_pairs_scanned": int((pair_summary_df["status"] == "SCANNED").sum()) if not pair_summary_df.empty else 0,
        "n_pairs_failed": int((pair_summary_df["status"] == "FAILED").sum()) if not pair_summary_df.empty else 0,
        "n_pairs_skipped": int((pair_summary_df["status"] == "SKIPPED").sum()) if not pair_summary_df.empty else 0,
        "n_detection_measurements": int(n_detection_measurements),
        "n_detections": int(len(detections_df)),
        "n_fade_candidates": int(len(fade_df)),
        "n_pass": int((fade_df["status"] == "PASS").sum()) if not fade_df.empty else 0,
        "n_review": int((fade_df["status"] == "REVIEW").sum()) if not fade_df.empty else 0,
        "n_benchmark_signals": int((detections_df["status"] == "BENCHMARK_SIGNAL").sum()) if not detections_df.empty else 0,
    }
    write_json(out_dir / "summary.json", summary)
    return summary


def _write_detection_measurements(out_dir: Path) -> int:
    frames: list[pd.DataFrame] = []
    for path in sorted((_checkpoint_dirs(out_dir)["measurements"]).glob("*.parquet")):
        frames.append(pd.read_parquet(path))
    measurement_df = pd.concat(frames, ignore_index=True) if frames else _empty_measurements_frame()
    measurement_df.to_parquet(out_dir / "detection_measurements.parquet", index=False)
    measurement_df.to_csv(out_dir / "detection_measurements.csv", index=False)
    return int(len(measurement_df))


def _curate_usable_pairs(
    candidate_pairs: pd.DataFrame,
    *,
    obs_meta: dict[str, dict[str, Any]],
    root_dir: Path,
    target_pairs: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]] = {}
    usable_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    for candidate_rank, (_, pair) in enumerate(candidate_pairs.iterrows(), start=1):
        issue = _pair_preflight_issue(pair, obs_meta, root_dir=root_dir, footprint_cache=footprint_cache)
        if issue is None:
            usable_rows.append(pair.to_dict())
            if target_pairs > 0 and len(usable_rows) >= target_pairs:
                break
            continue
        rejected_rows.append(
            {
                "candidate_rank": int(candidate_rank),
                "pair_id": str(pair["pair_id"]),
                "galaxy_name": str(pair["galaxy_name"]),
                "obs_id_1": str(pair["obs_id_1"]),
                "obs_id_2": str(pair["obs_id_2"]),
                "compatibility": str(pair.get("compatibility", "")),
                "baseline_days": float(pair.get("baseline_days", 0.0) or 0.0),
                "pair_score": float(pair.get("pair_score", 0.0) or 0.0),
                "preflight_reason": str(issue),
            }
        )
    usable = pd.DataFrame(usable_rows).reset_index(drop=True) if usable_rows else candidate_pairs.head(0).copy()
    rejected = pd.DataFrame(rejected_rows).reset_index(drop=True)
    return usable, rejected


def _safe_divide(num: np.ndarray | float, den: np.ndarray | float, *, default: float = np.nan) -> np.ndarray:
    num_arr = np.asarray(num, dtype=float)
    den_arr = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num_arr, den_arr).shape, default, dtype=float)
    good = np.isfinite(num_arr) & np.isfinite(den_arr) & (np.abs(den_arr) > 0)
    out[good] = num_arr[good] / den_arr[good]
    return out


def _slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return text.strip("_") or "item"


def _imaging_like_mask(frame: pd.DataFrame, *, filter_col: str = "filters", instrument_col: str = "instrument_name") -> pd.Series:
    filters = frame[filter_col].astype(str)
    instruments = frame[instrument_col].astype(str)
    return (
        ~filters.str.contains(";", regex=False, na=False)
        & ~filters.str.contains("^G\\d+", case=False, na=False, regex=True)
        & ~instruments.str.contains("NIRSPEC|STIS|COS|FOC|PHOTOMETER", case=False, na=False)
        & ~instruments.str.contains("MIRI/IFU", case=False, na=False)
    )


def _shift_image(data: np.ndarray, *, dx: float, dy: float, order: int = 1) -> np.ndarray:
    mask = np.isfinite(data).astype(float)
    filled = np.nan_to_num(np.asarray(data, dtype=float), nan=0.0)
    shifted_data = ndi_shift(filled, shift=(dy, dx), order=order, mode="constant", cval=0.0, prefilter=False)
    shifted_mask = ndi_shift(mask, shift=(dy, dx), order=1 if order > 0 else 0, mode="constant", cval=0.0, prefilter=False)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = shifted_data / np.clip(shifted_mask, 1e-6, np.inf)
    out[shifted_mask < 0.5] = np.nan
    return out


def _estimate_registration_shift(
    pre_data: np.ndarray,
    post_data: np.ndarray,
    valid_mask: np.ndarray,
    *,
    pixel_scale_arcsec: float,
    max_sources: int = 300,
) -> dict[str, float]:
    pre_src = _detect_sources(pre_data, valid_mask, pixel_scale_arcsec=pixel_scale_arcsec, max_sources=max_sources)
    post_src = _detect_sources(post_data, valid_mask, pixel_scale_arcsec=pixel_scale_arcsec, max_sources=max_sources)
    if pre_src.empty or post_src.empty:
        return {"dx_px": 0.0, "dy_px": 0.0, "residual_px": float("nan"), "n_matches": 0}

    pre_src = pre_src.sort_values("flux", ascending=False).head(max_sources).reset_index(drop=True)
    post_src = post_src.sort_values("flux", ascending=False).head(max_sources).reset_index(drop=True)
    pre_xy = pre_src[["x", "y"]].to_numpy(dtype=float)
    post_xy = post_src[["x", "y"]].to_numpy(dtype=float)
    dx = 0.0
    dy = 0.0
    residual = float("nan")
    n_matches = 0

    for radius in (5.0, 3.0, 2.0, 1.5):
        shifted_post = post_xy + np.array([dx, dy], dtype=float)
        tree = cKDTree(shifted_post)
        dist, idx = tree.query(pre_xy, distance_upper_bound=radius)
        good = np.isfinite(dist) & (dist < radius) & (idx < len(post_xy))
        if good.sum() < 6:
            continue
        deltas = pre_xy[good] - post_xy[idx[good]]
        median = np.nanmedian(deltas, axis=0)
        mad = np.nanmedian(np.abs(deltas - median), axis=0)
        scatter = np.where(np.isfinite(mad) & (mad > 0), 1.4826 * mad, 0.5)
        keep = np.all(np.abs(deltas - median) <= np.maximum(3.0 * scatter, 0.5), axis=1)
        if keep.sum() >= 6:
            deltas = deltas[keep]
        dx = float(np.nanmedian(deltas[:, 0]))
        dy = float(np.nanmedian(deltas[:, 1]))
        residual = float(np.nanmedian(np.hypot(deltas[:, 0] - dx, deltas[:, 1] - dy)))
        n_matches = int(len(deltas))
    return {"dx_px": dx, "dy_px": dy, "residual_px": residual, "n_matches": n_matches}


def _measure_residual_positions(
    pre_data: np.ndarray,
    post_data: np.ndarray,
    diff_data: np.ndarray,
    valid_mask: np.ndarray,
    positions: np.ndarray,
    *,
    radius: float,
) -> pd.DataFrame:
    flux_pre, err_pre, bg_pre = _aperture_fluxes(pre_data, valid_mask, positions, radius=radius)
    flux_post, err_post, bg_post = _aperture_fluxes(post_data, valid_mask, positions, radius=radius)
    flux_diff, err_diff, bg_diff = _aperture_fluxes(diff_data, valid_mask, positions, radius=radius)
    return pd.DataFrame(
        {
            "radius": radius,
            "flux_pre": flux_pre,
            "err_pre": err_pre,
            "flux_post": flux_post,
            "err_post": err_post,
            "flux_diff": flux_diff,
            "err_diff": err_diff,
            "bg_pre": bg_pre,
            "bg_post": bg_post,
            "bg_diff": bg_diff,
            "diff_sigma": _safe_divide(flux_diff, err_diff, default=np.nan),
            "post_to_pre_ratio": _safe_divide(flux_post, flux_pre, default=np.nan),
        }
    )


def _deduplicate_residuals(detections: pd.DataFrame, *, min_sep_px: float = 2.0) -> pd.DataFrame:
    if detections.empty:
        return detections
    keep_rows: list[int] = []
    for sign, sign_df in detections.groupby("event_sign", dropna=False):
        ordered = sign_df.sort_values("flux", ascending=False)
        accepted: list[np.ndarray] = []
        for idx, row in ordered.iterrows():
            pos = np.array([float(row["x"]), float(row["y"])], dtype=float)
            if accepted and np.min([np.hypot(*(pos - ref)) for ref in accepted]) < min_sep_px:
                continue
            accepted.append(pos)
            keep_rows.append(idx)
    return detections.loc[sorted(keep_rows)].reset_index(drop=True)


def _fallback_residual_peaks(
    signed_image: np.ndarray,
    valid_mask: np.ndarray,
    *,
    max_sources: int = 250,
    min_separation_px: int = 5,
) -> pd.DataFrame:
    good = valid_mask & np.isfinite(signed_image)
    if good.sum() < 100:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    median, std = _background_stats(signed_image, good)
    if not math.isfinite(std) or std <= 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])

    work = np.where(good, signed_image - median, 0.0)
    smooth = gaussian_filter(work, sigma=1.0)
    threshold = max(4.5 * std, float(np.nanpercentile(smooth[good], 99.5)) if good.any() else 0.0)
    if not math.isfinite(threshold) or threshold <= 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])

    local_max = smooth == maximum_filter(smooth, size=min_separation_px, mode="nearest")
    peak_mask = good & local_max & (smooth >= threshold)
    ys, xs = np.where(peak_mask)
    if len(xs) == 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])

    rows: list[dict[str, float]] = []
    for y_idx, x_idx in zip(ys, xs, strict=True):
        rows.append(
            {
                "x": float(x_idx),
                "y": float(y_idx),
                "flux": float(work[y_idx, x_idx]),
                "peak": float(smooth[y_idx, x_idx]),
            }
        )
    frame = pd.DataFrame(rows)
    if len(frame) > max_sources:
        frame = frame.sort_values("peak", ascending=False).head(max_sources)
    return frame.reset_index(drop=True)


def _fallback_component_centroids(
    signed_image: np.ndarray,
    valid_mask: np.ndarray,
    *,
    max_sources: int = 150,
    focus_center: tuple[float, float] | None = None,
    focus_half_size: int | None = None,
) -> pd.DataFrame:
    x_offset = 0
    y_offset = 0
    work_image = signed_image
    work_mask = valid_mask
    if focus_center is not None and focus_half_size is not None:
        cx, cy = focus_center
        x0 = max(int(round(cx)) - int(focus_half_size), 0)
        x1 = min(int(round(cx)) + int(focus_half_size) + 1, signed_image.shape[1])
        y0 = max(int(round(cy)) - int(focus_half_size), 0)
        y1 = min(int(round(cy)) + int(focus_half_size) + 1, signed_image.shape[0])
        x_offset = x0
        y_offset = y0
        work_image = signed_image[y0:y1, x0:x1]
        work_mask = valid_mask[y0:y1, x0:x1]

    good = work_mask & np.isfinite(work_image)
    if good.sum() < 100:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    median, std = _background_stats(work_image, good)
    if not math.isfinite(std) or std <= 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    work = np.where(good, work_image - median, 0.0)
    smooth = gaussian_filter(work, sigma=1.2)
    threshold = max(3.5 * std, float(np.nanpercentile(smooth[good], 99.0)) if good.any() else 0.0)
    if not math.isfinite(threshold) or threshold <= 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    mask = good & (smooth >= threshold)
    labels, n_labels = label(mask)
    if n_labels == 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    rows: list[dict[str, float]] = []
    for component_id in range(1, n_labels + 1):
        comp = labels == component_id
        area = int(comp.sum())
        if area < 4 or area > 400:
            continue
        weights = np.where(comp, np.clip(work, 0.0, None), 0.0)
        if float(np.nansum(weights)) <= 0:
            weights = np.where(comp, np.clip(smooth, 0.0, None), 0.0)
        if float(np.nansum(weights)) <= 0:
            continue
        cy, cx = center_of_mass(weights)
        if not (math.isfinite(cx) and math.isfinite(cy)):
            continue
        peak = float(np.nanmax(smooth[comp]))
        yi = int(np.clip(round(cy), 0, work_image.shape[0] - 1))
        xi = int(np.clip(round(cx), 0, work_image.shape[1] - 1))
        rows.append(
            {
                "x": float(x_offset + cx),
                "y": float(y_offset + cy),
                "flux": float(work[yi, xi]),
                "peak": peak,
            }
        )
    if not rows:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])
    frame = pd.DataFrame(rows).sort_values("peak", ascending=False)
    if len(frame) > max_sources:
        frame = frame.head(max_sources)
    return frame.reset_index(drop=True)


def _refine_residual_centroids(
    detections: pd.DataFrame,
    diff_image: np.ndarray,
    valid_mask: np.ndarray,
    *,
    stamp_half_size: int = 5,
    max_shift_px: float = 3.0,
) -> pd.DataFrame:
    if detections.empty:
        return detections
    refined = detections.copy()
    for idx, row in refined.iterrows():
        sign = str(row.get("event_sign", "fade"))
        signed = diff_image if sign == "fade" else -diff_image
        x0 = int(round(float(row["x"])))
        y0 = int(round(float(row["y"])))
        xs = slice(max(x0 - stamp_half_size, 0), min(x0 + stamp_half_size + 1, signed.shape[1]))
        ys = slice(max(y0 - stamp_half_size, 0), min(y0 + stamp_half_size + 1, signed.shape[0]))
        stamp = np.asarray(signed[ys, xs], dtype=float)
        stamp_valid = np.asarray(valid_mask[ys, xs], dtype=bool) & np.isfinite(stamp)
        if stamp_valid.sum() < 9:
            continue
        bg = float(np.nanmedian(stamp[stamp_valid]))
        weights = np.where(stamp_valid, np.clip(stamp - bg, 0.0, None), 0.0)
        if not np.isfinite(weights).any() or float(np.nansum(weights)) <= 0:
            continue
        cy, cx = center_of_mass(weights)
        if not (math.isfinite(cx) and math.isfinite(cy)):
            continue
        x_new = float(xs.start + cx)
        y_new = float(ys.start + cy)
        shift = float(np.hypot(x_new - float(row["x"]), y_new - float(row["y"])))
        if shift > max_shift_px:
            continue
        refined.at[idx, "x"] = x_new
        refined.at[idx, "y"] = y_new
        xi = int(np.clip(round(x_new), 0, signed.shape[1] - 1))
        yi = int(np.clip(round(y_new), 0, signed.shape[0] - 1))
        refined.at[idx, "peak"] = float(signed[yi, xi]) if math.isfinite(float(signed[yi, xi])) else float(row.get("peak", np.nan))
    return refined


def _residual_status(
    *,
    event_sign: str,
    pair_is_stable: bool,
    diff_sigma: float,
    post_to_pre_ratio: float,
    pre_snr: float,
    post_snr: float,
    edge_distance_px: float,
    blend_risk: float,
    crowding_index: float,
    registration_residual_px: float,
) -> str | None:
    if not pair_is_stable:
        return None
    if event_sign == "fade":
        if (
            math.isfinite(diff_sigma)
            and diff_sigma >= 7.0
            and math.isfinite(post_to_pre_ratio)
            and post_to_pre_ratio <= 0.05
            and post_to_pre_ratio > -0.1
            and post_snr <= 3.0
            and edge_distance_px >= 10.0
            and blend_risk < 0.5
            and crowding_index < 4.0
            and (not math.isfinite(registration_residual_px) or registration_residual_px <= 0.75)
        ):
            return "PASS"
        if (
            math.isfinite(diff_sigma)
            and diff_sigma >= 5.0
            and math.isfinite(post_to_pre_ratio)
            and post_to_pre_ratio <= 0.50
            and edge_distance_px >= 6.0
            and (not math.isfinite(registration_residual_px) or registration_residual_px <= 1.0)
        ):
            return "REVIEW"
        return None
    if event_sign == "brighten":
        if math.isfinite(diff_sigma) and diff_sigma >= 5.0 and post_snr >= max(5.0, pre_snr + 2.0) and edge_distance_px >= 6.0:
            return "BENCHMARK_SIGNAL"
    return None


def _candidate_packet(
    *,
    packet_dir: Path,
    detection_id: str,
    title: str,
    pre_crop: np.ndarray,
    post_crop: np.ndarray,
    diff_crop: np.ndarray,
    x_pix: float,
    y_pix: float,
    payload: dict[str, Any],
) -> Path:
    ensure_dir(packet_dir)
    half = 18
    x0 = max(int(round(x_pix)) - half, 0)
    x1 = min(int(round(x_pix)) + half + 1, pre_crop.shape[1])
    y0 = max(int(round(y_pix)) - half, 0)
    y1 = min(int(round(y_pix)) + half + 1, pre_crop.shape[0])
    figure_path = _save_candidate_figure(
        packet_dir,
        detection_id,
        pre_crop[y0:y1, x0:x1],
        post_crop[y0:y1, x0:x1],
        diff_crop[y0:y1, x0:x1],
        title=title,
    )
    write_json(packet_dir / "summary.json", {**payload, "figure_path": str(figure_path)})
    return figure_path


def _resolve_observation_product(row: pd.Series, *, root_dir: Path) -> tuple[dict[str, Any], ScienceImage]:
    products = _cached_products(root_dir, obs_collection=str(row["obs_collection"]), obs_id=str(row["obs_id"]))
    if products.empty:
        products = _product_query(int(row["obsid"]))
    product = _select_best_product(products, obs_collection=str(row["obs_collection"]))
    if product is None:
        rights = sorted({str(value).upper() for value in products.get("dataRights", pd.Series(dtype=object)).dropna().tolist()})
        if rights and not any(value in {"", "PUBLIC"} for value in rights):
            raise RuntimeError(f"Restricted calibrated products for {row['obs_id']}: {','.join(rights)}")
        raise RuntimeError(f"No suitable calibrated FITS product for {row['obs_id']}")
    path = _download_product(product, root_dir=root_dir, obs_collection=str(row["obs_collection"]), obs_id=str(row["obs_id"]))
    image = _load_science_image(
        path,
        obs_collection=str(row["obs_collection"]),
        obs_id=str(row["obs_id"]),
        filter_name=str(row["filters"]),
        instrument_name=str(row["instrument_name"]),
    )
    return {
        "product_filename": str(product["productFilename"]),
        "product_path": str(path),
        "product_uri": str(product["dataURI"]),
        "product_priority": float(product["priority"]),
        "image_ext": image.extname,
    }, image


def _load_science_footprint(path: Path) -> ScienceFootprint:
    best_score = -1.0
    best_wcs: WCS | None = None
    best_shape: tuple[int, int] | None = None
    best_extname = "PRIMARY"
    with fits.open(path, memmap=False, do_not_scale_image_data=True) as hdul:
        for hdu in hdul:
            header = hdu.header
            naxis = int(header.get("NAXIS", 0) or 0)
            if naxis < 2:
                continue
            try:
                wcs = WCS(header).celestial
            except Exception:
                continue
            if not wcs.has_celestial:
                continue
            extname = (getattr(hdu, "name", "") or "PRIMARY").upper()
            if naxis == 2:
                shape = (int(header.get("NAXIS2", 0) or 0), int(header.get("NAXIS1", 0) or 0))
                score = float(shape[0] * shape[1])
            elif naxis == 3 and 1 <= int(header.get("NAXIS3", 0) or 0) <= 8:
                shape = (int(header.get("NAXIS2", 0) or 0), int(header.get("NAXIS1", 0) or 0))
                score = float(shape[0] * shape[1]) + 1e7
            else:
                continue
            if min(shape) <= 0:
                continue
            if extname == "SCI":
                score += 1e9
            elif extname == "PRIMARY":
                score += 5e8
            if score > best_score:
                best_score = score
                best_wcs = wcs
                best_shape = shape
                best_extname = extname
    if best_wcs is None or best_shape is None:
        raise RuntimeError(f"No celestial footprint found in {path}")
    pix_scale = float(np.nanmedian(proj_plane_pixel_scales(best_wcs)) * 3600.0)
    return ScienceFootprint(
        path=path,
        wcs=best_wcs,
        shape=best_shape,
        pixel_scale_arcsec=pix_scale if math.isfinite(pix_scale) and pix_scale > 0 else 0.05,
        extname=best_extname,
    )


def _resolve_observation_footprint(row: pd.Series, *, root_dir: Path) -> tuple[dict[str, Any], ScienceFootprint]:
    products = _cached_products(root_dir, obs_collection=str(row["obs_collection"]), obs_id=str(row["obs_id"]))
    if products.empty:
        products = _product_query(int(row["obsid"]))
    product = _select_best_product(products, obs_collection=str(row["obs_collection"]))
    if product is None:
        rights = sorted({str(value).upper() for value in products.get("dataRights", pd.Series(dtype=object)).dropna().tolist()})
        if rights and not any(value in {"", "PUBLIC"} for value in rights):
            raise RuntimeError(f"Restricted calibrated products for {row['obs_id']}: {','.join(rights)}")
        raise RuntimeError(f"No suitable calibrated FITS product for {row['obs_id']}")
    path = _download_product(product, root_dir=root_dir, obs_collection=str(row["obs_collection"]), obs_id=str(row["obs_id"]))
    footprint = _load_science_footprint(path)
    return {
        "product_filename": str(product["productFilename"]),
        "product_path": str(path),
        "product_uri": str(product["dataURI"]),
        "product_priority": float(product["priority"]),
        "image_ext": footprint.extname,
    }, footprint


def _resolve_cached_observation_image(
    obs_id: str,
    *,
    obs_meta: dict[str, dict[str, Any]],
    root_dir: Path,
    image_cache: dict[str, tuple[dict[str, Any], ScienceImage]],
) -> tuple[dict[str, Any], ScienceImage]:
    key = str(obs_id)
    cached = image_cache.get(key)
    if cached is not None:
        return cached
    row_payload = dict(obs_meta[key])
    row_payload.setdefault("obs_id", key)
    row = pd.Series(row_payload)
    cached = _resolve_observation_product(row, root_dir=root_dir)
    image_cache[key] = cached
    return cached


def _resolve_cached_observation_footprint(
    obs_id: str,
    *,
    obs_meta: dict[str, dict[str, Any]],
    root_dir: Path,
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]],
) -> tuple[dict[str, Any], ScienceFootprint]:
    key = str(obs_id)
    cached = footprint_cache.get(key)
    if cached is not None:
        return cached
    row_payload = dict(obs_meta[key])
    row_payload.setdefault("obs_id", key)
    row = pd.Series(row_payload)
    cached = _resolve_observation_footprint(row, root_dir=root_dir)
    footprint_cache[key] = cached
    return cached


def _prepare_registered_pair(
    ref_image: ScienceImage,
    cmp_image: ScienceImage,
    *,
    center_world: SkyCoord | None,
    cutout_arcsec: float,
) -> dict[str, Any]:
    if center_world is None:
        ref_cut = ref_image
        cmp_cut = cmp_image
    else:
        ref_cut = _candidate_cutout(ref_image, center_world, size_arcsec=cutout_arcsec, min_pixels=128)
        cmp_cut = _candidate_cutout(cmp_image, center_world, size_arcsec=cutout_arcsec, min_pixels=128)
    from reproject import reproject_interp

    cmp_reproj, footprint = reproject_interp((cmp_cut.data, cmp_cut.wcs), ref_cut.wcs, shape_out=ref_cut.data.shape)
    overlap = _finite_mask(ref_cut.data) & _finite_mask(cmp_reproj) & (footprint > 0)
    if overlap.sum() < 2500:
        raise RuntimeError("Insufficient overlap after WCS reprojection.")

    ref_bg, _ = _background_stats(ref_cut.data, overlap)
    cmp_bg, _ = _background_stats(cmp_reproj, overlap)
    ref_sub = ref_cut.data - ref_bg
    cmp_sub = cmp_reproj - cmp_bg

    reg = _estimate_registration_shift(ref_sub, cmp_sub, overlap, pixel_scale_arcsec=ref_cut.pixel_scale_arcsec)
    cmp_aligned = _shift_image(cmp_sub, dx=reg["dx_px"], dy=reg["dy_px"], order=1)
    footprint_aligned = _shift_image(np.asarray(footprint, dtype=float), dx=reg["dx_px"], dy=reg["dy_px"], order=0)
    valid = _finite_mask(ref_sub) & _finite_mask(cmp_aligned) & (footprint_aligned > 0.5)
    if valid.sum() < 2500:
        raise RuntimeError("Insufficient overlap after empirical registration.")

    anchors = _detect_sources(ref_sub, valid, pixel_scale_arcsec=ref_cut.pixel_scale_arcsec, max_sources=300)
    positions = anchors[["x", "y"]].to_numpy(dtype=float) if not anchors.empty else np.empty((0, 2), dtype=float)
    scale_factor, n_scale_sources = _estimate_scale_factor(ref_sub, cmp_aligned, valid, positions)
    cmp_scaled = cmp_aligned / scale_factor
    diff_image = ref_sub - cmp_scaled
    return {
        "ref_image": ref_cut,
        "cmp_image": cmp_cut,
        "ref_sub": ref_sub,
        "cmp_scaled": cmp_scaled,
        "diff_image": diff_image,
        "valid_mask": valid,
        "anchors": anchors,
        "registration": reg,
        "scale_factor": float(scale_factor),
        "n_scale_sources": int(n_scale_sources),
        "edge_distance_map": distance_transform_edt(valid),
    }


def _covers_center_with_margin(image: ScienceImage | ScienceFootprint, center_world: SkyCoord, *, cutout_arcsec: float) -> bool:
    x_pix, y_pix = image.wcs.world_to_pixel(center_world)
    if not (np.isfinite(x_pix) and np.isfinite(y_pix)):
        return False
    half_size_pix = max(64, int(0.5 * cutout_arcsec / max(image.pixel_scale_arcsec, 1e-3)))
    shape = image.data.shape if hasattr(image, "data") else image.shape
    return (
        half_size_pix <= float(x_pix) < shape[1] - half_size_pix
        and half_size_pix <= float(y_pix) < shape[0] - half_size_pix
    )


def _pair_has_explicit_center(pair: pd.Series) -> bool:
    ra = pd.to_numeric(pd.Series([pair.get("center_ra_deg")]), errors="coerce").iloc[0]
    dec = pd.to_numeric(pd.Series([pair.get("center_dec_deg")]), errors="coerce").iloc[0]
    return bool(math.isfinite(ra) and math.isfinite(dec))


def _pair_center_world(pair: pd.Series, obs_meta: dict[str, dict[str, Any]]) -> SkyCoord:
    meta_1 = obs_meta[str(pair["obs_id_1"])]
    center_ra = pd.to_numeric(pd.Series([pair.get("center_ra_deg", meta_1["galaxy_ra"])]), errors="coerce").iloc[0]
    center_dec = pd.to_numeric(pd.Series([pair.get("center_dec_deg", meta_1["galaxy_dec"])]), errors="coerce").iloc[0]
    if not math.isfinite(center_ra) or not math.isfinite(center_dec):
        raise RuntimeError("Non-finite pair center coordinates.")
    return SkyCoord(float(center_ra) * u.deg, float(center_dec) * u.deg, frame="icrs")


def _pair_cutout_arcsec(pair: pd.Series, obs_meta: dict[str, dict[str, Any]]) -> float:
    meta_1 = obs_meta[str(pair["obs_id_1"])]
    meta_2 = obs_meta[str(pair["obs_id_2"])]
    default_cutout = max(240.0, 2.2 * max(float(meta_1.get("search_radius_deg", 0.0) or 0.0), float(meta_2.get("search_radius_deg", 0.0) or 0.0)) * 3600.0)
    center_source = str(pair.get("center_source", "") or "")
    if center_source and center_source != "galaxy_center":
        return float(np.clip(default_cutout, 96.0, 160.0))
    return float(default_cutout)


def _compact_exception_message(exc: Exception) -> str:
    message = re.sub(r"\s+", " ", str(exc)).strip()
    return message[:160] if len(message) > 160 else message


def _benchmark_pair_covers_truth(
    pair: pd.Series,
    *,
    obs_meta: dict[str, dict[str, Any]],
    root_dir: Path,
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]],
) -> bool:
    return _pair_preflight_issue(pair, obs_meta, root_dir=root_dir, footprint_cache=footprint_cache) is None


def _pair_preflight_issue(
    pair: pd.Series,
    obs_meta: dict[str, dict[str, Any]],
    *,
    root_dir: Path | None = None,
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]] | None = None,
) -> str | None:
    for key in [str(pair["obs_id_1"]), str(pair["obs_id_2"])]:
        meta = obs_meta[key]
        rights = str(meta.get("dataRights", "") or "").upper()
        if rights not in {"", "PUBLIC"}:
            return f"RESTRICTED_DATA_RIGHTS:{key}:{rights}"
    if root_dir is None or footprint_cache is None:
        return None
    try:
        _, image_1 = _resolve_cached_observation_footprint(str(pair["obs_id_1"]), obs_meta=obs_meta, root_dir=root_dir, footprint_cache=footprint_cache)
        _, image_2 = _resolve_cached_observation_footprint(str(pair["obs_id_2"]), obs_meta=obs_meta, root_dir=root_dir, footprint_cache=footprint_cache)
    except Exception as exc:
        return f"PRODUCT_PRECHECK_FAILED:{_compact_exception_message(exc)}"
    if _pair_has_explicit_center(pair):
        center_world = _pair_center_world(pair, obs_meta)
        cutout_arcsec = _pair_cutout_arcsec(pair, obs_meta)
        if not _covers_center_with_margin(image_1, center_world, cutout_arcsec=cutout_arcsec) or not _covers_center_with_margin(image_2, center_world, cutout_arcsec=cutout_arcsec):
            return "TARGET_CENTER_OUTSIDE_FOOTPRINT"
    return None


def _scan_pair_difference(
    pair: pd.Series,
    obs_meta: dict[str, dict[str, Any]],
    *,
    root_dir: Path,
    out_dir: Path,
    sign_mode: str,
    image_cache: dict[str, tuple[dict[str, Any], ScienceImage]],
) -> tuple[dict[str, Any], list[dict[str, Any]], pd.DataFrame]:
    meta_1 = obs_meta[str(pair["obs_id_1"])]
    meta_2 = obs_meta[str(pair["obs_id_2"])]
    has_explicit_center = _pair_has_explicit_center(pair)
    product_1, image_1 = _resolve_cached_observation_image(str(pair["obs_id_1"]), obs_meta=obs_meta, root_dir=root_dir, image_cache=image_cache)
    product_2, image_2 = _resolve_cached_observation_image(str(pair["obs_id_2"]), obs_meta=obs_meta, root_dir=root_dir, image_cache=image_cache)
    product_info = {
        "obsid_1": int(pair["obsid_1"]),
        "obsid_2": int(pair["obsid_2"]),
        "obs_id_1": str(pair["obs_id_1"]),
        "obs_id_2": str(pair["obs_id_2"]),
        "product_filename_1": str(product_1["product_filename"]),
        "product_filename_2": str(product_2["product_filename"]),
        "product_uri_1": str(product_1["product_uri"]),
        "product_uri_2": str(product_2["product_uri"]),
        "product_path_1": str(product_1["product_path"]),
        "product_path_2": str(product_2["product_path"]),
        "product_priority_1": float(product_1["product_priority"]),
        "product_priority_2": float(product_2["product_priority"]),
        "image_ext_1": str(product_1["image_ext"]),
        "image_ext_2": str(product_2["image_ext"]),
    }
    if has_explicit_center:
        center_world = _pair_center_world(pair, obs_meta)
        cutout_arcsec = _pair_cutout_arcsec(pair, obs_meta)
        if not _covers_center_with_margin(image_1, center_world, cutout_arcsec=cutout_arcsec) or not _covers_center_with_margin(image_2, center_world, cutout_arcsec=cutout_arcsec):
            raise RuntimeError("Target center outside usable footprint before reprojection.")
        prepared = _prepare_registered_pair(image_1, image_2, center_world=center_world, cutout_arcsec=cutout_arcsec)
    else:
        prepared = _prepare_registered_pair(image_1, image_2, center_world=None, cutout_arcsec=0.0)

    ref_sub = prepared["ref_sub"]
    cmp_scaled = prepared["cmp_scaled"]
    diff_image = prepared["diff_image"]
    valid = prepared["valid_mask"]
    anchors = prepared["anchors"]
    edge_distance_map = prepared["edge_distance_map"]
    reg = prepared["registration"]

    if anchors.empty:
        raise RuntimeError("No anchor stars detected after registration.")

    anchor_positions = anchors[["x", "y"]].to_numpy(dtype=float)
    anchor_metrics = _measure_residual_positions(ref_sub, cmp_scaled, diff_image, valid, anchor_positions, radius=3.5)
    anchor_catalog = anchors.copy()
    anchor_catalog["flux_pre"] = anchor_metrics["flux_pre"].to_numpy()
    anchor_catalog["flux_post"] = anchor_metrics["flux_post"].to_numpy()
    anchor_catalog["flux_ratio"] = anchor_metrics["post_to_pre_ratio"].to_numpy()
    stable_fraction, extreme_fraction, pair_is_stable = _pair_stability_metrics(anchor_catalog)

    residual_frames: list[pd.DataFrame] = []
    if sign_mode in {"fade", "both"}:
        fade_resid = _detect_sources(diff_image, valid, pixel_scale_arcsec=prepared["ref_image"].pixel_scale_arcsec, max_sources=500)
        if not fade_resid.empty:
            fade_resid["event_sign"] = "fade"
            residual_frames.append(fade_resid)
        fade_peaks = _fallback_residual_peaks(diff_image, valid, max_sources=250)
        if not fade_peaks.empty:
            fade_peaks["event_sign"] = "fade"
            residual_frames.append(fade_peaks)
        if has_explicit_center:
            fade_components = _fallback_component_centroids(
                diff_image,
                valid,
                max_sources=80,
                focus_center=(0.5 * diff_image.shape[1], 0.5 * diff_image.shape[0]),
                focus_half_size=48,
            )
            if not fade_components.empty:
                fade_components["event_sign"] = "fade"
                residual_frames.append(fade_components)
    if sign_mode in {"brighten", "both"}:
        bright_resid = _detect_sources(-diff_image, valid, pixel_scale_arcsec=prepared["ref_image"].pixel_scale_arcsec, max_sources=500)
        if not bright_resid.empty:
            bright_resid["event_sign"] = "brighten"
            residual_frames.append(bright_resid)
        bright_peaks = _fallback_residual_peaks(-diff_image, valid, max_sources=250)
        if not bright_peaks.empty:
            bright_peaks["event_sign"] = "brighten"
            residual_frames.append(bright_peaks)
        if has_explicit_center:
            bright_components = _fallback_component_centroids(
                -diff_image,
                valid,
                max_sources=80,
                focus_center=(0.5 * diff_image.shape[1], 0.5 * diff_image.shape[0]),
                focus_half_size=48,
            )
            if not bright_components.empty:
                bright_components["event_sign"] = "brighten"
                residual_frames.append(bright_components)

    residual_sources = _deduplicate_residuals(pd.concat(residual_frames, ignore_index=True) if residual_frames else pd.DataFrame(columns=["x", "y", "flux", "peak", "event_sign"]))
    residual_sources = _refine_residual_centroids(residual_sources, diff_image, valid)
    residual_sources = _deduplicate_residuals(residual_sources)
    if residual_sources.empty:
        pair_summary = {
            "pair_id": str(pair["pair_id"]),
            "galaxy_name": str(pair["galaxy_name"]),
            "obs_id_1": str(pair["obs_id_1"]),
            "obs_id_2": str(pair["obs_id_2"]),
            "filter_1": str(pair["filter_1"]),
            "filter_2": str(pair["filter_2"]),
            "compatibility": str(pair["compatibility"]),
            "baseline_days": float(pair["baseline_days"]),
            "pair_score": float(pair["pair_score"]),
            "registration_dx_px": float(reg["dx_px"]),
            "registration_dy_px": float(reg["dy_px"]),
            "registration_residual_px": float(reg["residual_px"]) if math.isfinite(reg["residual_px"]) else np.nan,
            "registration_n_matches": int(reg["n_matches"]),
            "scale_factor_post_to_pre": float(prepared["scale_factor"]),
            "n_scale_sources": int(prepared["n_scale_sources"]),
            "pair_stable_fraction": stable_fraction,
            "pair_extreme_fade_fraction": extreme_fraction,
            "pair_is_stable": pair_is_stable,
            "n_residual_detections": 0,
            **product_info,
        }
        return pair_summary, [], residual_sources

    det_positions = residual_sources[["x", "y"]].to_numpy(dtype=float)
    measures = [_measure_residual_positions(ref_sub, cmp_scaled, diff_image, valid, det_positions, radius=radius) for radius in (2.5, 3.5, 4.5)]
    metrics = measures[1].copy()
    metrics["diff_sigma_min"] = np.nanmin(np.vstack([m["diff_sigma"].to_numpy() for m in measures]), axis=0)
    metrics["diff_sigma_max"] = np.nanmax(np.vstack([m["diff_sigma"].to_numpy() for m in measures]), axis=0)
    metrics["post_to_pre_ratio_med"] = np.nanmedian(np.vstack([m["post_to_pre_ratio"].to_numpy() for m in measures]), axis=0)
    metrics["diff_flux_med"] = np.nanmedian(np.vstack([m["flux_diff"].to_numpy() for m in measures]), axis=0)
    metrics["diff_flux_std"] = np.nanstd(np.vstack([m["flux_diff"].to_numpy() for m in measures]), axis=0)

    anchor_tree = cKDTree(anchor_positions)
    det_dist, _ = anchor_tree.query(det_positions, k=1)
    crowding_index = np.array([len(anchor_tree.query_ball_point(pos, r=10.0)) for pos in det_positions], dtype=float)
    blend_risk = np.clip(1.0 - det_dist / 10.0, 0.0, 1.0)
    xi = np.clip(np.round(det_positions[:, 0]).astype(int), 0, edge_distance_map.shape[1] - 1)
    yi = np.clip(np.round(det_positions[:, 1]).astype(int), 0, edge_distance_map.shape[0] - 1)
    edge_distance_px = edge_distance_map[yi, xi]

    detection_table = residual_sources.copy()
    detection_table["diff_flux"] = metrics["diff_flux_med"].to_numpy()
    detection_table["diff_flux_std"] = metrics["diff_flux_std"].to_numpy()
    detection_table["diff_sigma"] = metrics["diff_sigma"].to_numpy()
    detection_table["diff_sigma_min"] = metrics["diff_sigma_min"].to_numpy()
    detection_table["diff_sigma_max"] = metrics["diff_sigma_max"].to_numpy()
    detection_table["flux_pre"] = metrics["flux_pre"].to_numpy()
    detection_table["flux_post"] = metrics["flux_post"].to_numpy()
    detection_table["err_pre"] = metrics["err_pre"].to_numpy()
    detection_table["err_post"] = metrics["err_post"].to_numpy()
    detection_table["post_to_pre_ratio"] = metrics["post_to_pre_ratio_med"].to_numpy()
    detection_table["pre_snr"] = _safe_divide(np.abs(metrics["flux_pre"].to_numpy()), metrics["err_pre"].to_numpy(), default=np.nan)
    detection_table["post_snr"] = _safe_divide(np.abs(metrics["flux_post"].to_numpy()), metrics["err_post"].to_numpy(), default=np.nan)
    detection_table["crowding_index"] = crowding_index
    detection_table["blend_risk"] = blend_risk
    detection_table["edge_distance_px"] = edge_distance_px
    detection_table["pair_is_stable"] = bool(pair_is_stable)
    detection_table["registration_residual_px"] = float(reg["residual_px"]) if math.isfinite(reg["residual_px"]) else np.nan
    detection_table["registration_dx_px"] = float(reg["dx_px"])
    detection_table["registration_dy_px"] = float(reg["dy_px"])
    detection_table["registration_n_matches"] = int(reg["n_matches"])
    detection_table["pair_stable_fraction"] = float(stable_fraction) if math.isfinite(stable_fraction) else np.nan
    detection_table["pair_extreme_fade_fraction"] = float(extreme_fraction) if math.isfinite(extreme_fraction) else np.nan
    detection_table["pair_id"] = str(pair["pair_id"])
    detection_table["galaxy_name"] = str(pair["galaxy_name"])
    detection_table["obs_id_1"] = str(pair["obs_id_1"])
    detection_table["obs_id_2"] = str(pair["obs_id_2"])
    detection_table["filter_1"] = str(pair["filter_1"])
    detection_table["filter_2"] = str(pair["filter_2"])
    detection_table["pre_mjd"] = float(meta_1["t_min"])
    detection_table["post_mjd"] = float(meta_2["t_min"])
    det_world = prepared["ref_image"].wcs.pixel_to_world(
        detection_table["x"].to_numpy(dtype=float),
        detection_table["y"].to_numpy(dtype=float),
    )
    detection_table["ra_deg"] = np.asarray(det_world.ra.deg, dtype=float)
    detection_table["dec_deg"] = np.asarray(det_world.dec.deg, dtype=float)

    candidates: list[dict[str, Any]] = []
    pair_packet_root = _pair_checkpoint_paths(out_dir, str(pair["pair_id"]))["packet_root"]
    shutil.rmtree(pair_packet_root, ignore_errors=True)
    packet_root = ensure_dir(pair_packet_root)
    for idx, row in detection_table.iterrows():
        status = _residual_status(
            event_sign=str(row["event_sign"]),
            pair_is_stable=bool(row["pair_is_stable"]),
            diff_sigma=float(row["diff_sigma"]) if math.isfinite(row["diff_sigma"]) else float("nan"),
            post_to_pre_ratio=float(row["post_to_pre_ratio"]) if math.isfinite(row["post_to_pre_ratio"]) else float("nan"),
            pre_snr=float(row["pre_snr"]) if math.isfinite(row["pre_snr"]) else float("nan"),
            post_snr=float(row["post_snr"]) if math.isfinite(row["post_snr"]) else float("nan"),
            edge_distance_px=float(row["edge_distance_px"]),
            blend_risk=float(row["blend_risk"]),
            crowding_index=float(row["crowding_index"]),
            registration_residual_px=float(row["registration_residual_px"]) if math.isfinite(row["registration_residual_px"]) else float("nan"),
        )
        if status is None:
            continue

        world = prepared["ref_image"].wcs.pixel_to_world(float(row["x"]), float(row["y"]))
        detection_id = _slugify(f"{pair['galaxy_name']}_{pair['obs_id_1']}_{pair['obs_id_2']}_{row['event_sign']}_{idx:04d}")
        reason_codes: list[str] = []
        warning_flags: list[str] = []
        ratio = float(row["post_to_pre_ratio"]) if math.isfinite(row["post_to_pre_ratio"]) else float("nan")
        if str(row["event_sign"]) == "fade":
            if math.isfinite(ratio) and ratio <= 0.05:
                reason_codes.append("POST_LT_5PCT")
            if math.isfinite(float(row["diff_sigma_min"])) and float(row["diff_sigma_min"]) >= 5.0:
                reason_codes.append("ROBUST_DIFF_SIG")
            if math.isfinite(float(row["registration_residual_px"])) and float(row["registration_residual_px"]) <= 0.5:
                reason_codes.append("REGISTRATION_OK")
        else:
            reason_codes.append("BRIGHTENING_SIGNAL")
        if float(row["blend_risk"]) >= 0.5:
            warning_flags.append("BLEND_RISK")
        if float(row["crowding_index"]) >= 4.0:
            warning_flags.append("CROWDED")
        if float(row["edge_distance_px"]) < 10.0:
            warning_flags.append("EDGE_NEAR_MASK")
        if math.isfinite(float(row["registration_residual_px"])) and float(row["registration_residual_px"]) > 0.75:
            warning_flags.append("ASTROMETRY_RESIDUAL")
        if int(row["registration_n_matches"]) < 10:
            warning_flags.append("LOW_REGISTRATION_ANCHORS")

        packet_dir = ensure_dir(packet_root / detection_id)
        figure_path = _candidate_packet(
            packet_dir=packet_dir,
            detection_id=detection_id,
            title=f"{pair['galaxy_name']}  {pair['filter_1']}->{pair['filter_2']}  {row['event_sign']}",
            pre_crop=ref_sub,
            post_crop=cmp_scaled,
            diff_crop=diff_image,
            x_pix=float(row["x"]),
            y_pix=float(row["y"]),
            payload={
                "candidate_id": detection_id,
                "status": status,
                "event_sign": str(row["event_sign"]),
                "ra_deg": float(world.ra.deg),
                "dec_deg": float(world.dec.deg),
                "galaxy_name": str(pair["galaxy_name"]),
                "pair_id": str(pair["pair_id"]),
                "diff_sigma": float(row["diff_sigma"]),
                "diff_sigma_min": float(row["diff_sigma_min"]),
                "post_to_pre_ratio": ratio,
                "pre_flux": float(row["flux_pre"]),
                "post_flux": float(row["flux_post"]),
                "registration_dx_px": float(row["registration_dx_px"]),
                "registration_dy_px": float(row["registration_dy_px"]),
                "registration_residual_px": float(row["registration_residual_px"]) if math.isfinite(row["registration_residual_px"]) else np.nan,
                "registration_n_matches": int(row["registration_n_matches"]),
                "warning_flags": warning_flags,
                "reason_codes": reason_codes,
                "product_info": product_info,
            },
        )
        candidates.append(
            {
                "candidate_id": detection_id,
                "pair_id": str(pair["pair_id"]),
                "galaxy_name": str(pair["galaxy_name"]),
                "ra_deg": float(world.ra.deg),
                "dec_deg": float(world.dec.deg),
                "event_sign": str(row["event_sign"]),
                "status": status,
                "diff_sigma": float(row["diff_sigma"]),
                "diff_sigma_min": float(row["diff_sigma_min"]),
                "post_to_pre_ratio": ratio,
                "pre_flux": float(row["flux_pre"]),
                "post_flux": float(row["flux_post"]),
                "pre_snr": float(row["pre_snr"]) if math.isfinite(row["pre_snr"]) else np.nan,
                "post_snr": float(row["post_snr"]) if math.isfinite(row["post_snr"]) else np.nan,
                "blend_risk": float(row["blend_risk"]),
                "crowding_index": float(row["crowding_index"]),
                "edge_distance_px": float(row["edge_distance_px"]),
                "registration_dx_px": float(row["registration_dx_px"]),
                "registration_dy_px": float(row["registration_dy_px"]),
                "registration_residual_px": float(row["registration_residual_px"]) if math.isfinite(row["registration_residual_px"]) else np.nan,
                "registration_n_matches": int(row["registration_n_matches"]),
                "pair_stable_fraction": float(row["pair_stable_fraction"]) if math.isfinite(row["pair_stable_fraction"]) else np.nan,
                "pair_extreme_fade_fraction": float(row["pair_extreme_fade_fraction"]) if math.isfinite(row["pair_extreme_fade_fraction"]) else np.nan,
                "pre_obsid": int(pair["obsid_1"]),
                "post_obsid": int(pair["obsid_2"]),
                "pre_obs_id": str(pair["obs_id_1"]),
                "post_obs_id": str(pair["obs_id_2"]),
                "pre_obs_collection": str(pair["obs_collection_1"]),
                "post_obs_collection": str(pair["obs_collection_2"]),
                "pre_filter": str(pair["filter_1"]),
                "post_filter": str(pair["filter_2"]),
                "pre_instrument": str(pair["instrument_1"]),
                "post_instrument": str(pair["instrument_2"]),
                "pre_mjd": float(meta_1["t_min"]),
                "post_mjd": float(meta_2["t_min"]),
                "cutout_path": str(figure_path),
                "reason_codes_json": json_list(reason_codes),
                "warning_flags_json": json_list(warning_flags),
                "provenance_json": json.dumps(
                    {
                        "baseline_days": float(pair["baseline_days"]),
                        "pair_score": float(pair["pair_score"]),
                        "product_info": product_info,
                        "scale_factor_post_to_pre": float(prepared["scale_factor"]),
                        "n_scale_sources": int(prepared["n_scale_sources"]),
                    },
                    sort_keys=True,
                ),
            }
        )

    pair_summary = {
        "pair_id": str(pair["pair_id"]),
        "galaxy_name": str(pair["galaxy_name"]),
        "obs_id_1": str(pair["obs_id_1"]),
        "obs_id_2": str(pair["obs_id_2"]),
        "filter_1": str(pair["filter_1"]),
        "filter_2": str(pair["filter_2"]),
        "compatibility": str(pair["compatibility"]),
        "baseline_days": float(pair["baseline_days"]),
        "pair_score": float(pair["pair_score"]),
        "registration_dx_px": float(reg["dx_px"]),
        "registration_dy_px": float(reg["dy_px"]),
        "registration_residual_px": float(reg["residual_px"]) if math.isfinite(reg["residual_px"]) else np.nan,
        "registration_n_matches": int(reg["n_matches"]),
        "scale_factor_post_to_pre": float(prepared["scale_factor"]),
        "n_scale_sources": int(prepared["n_scale_sources"]),
        "pair_stable_fraction": stable_fraction,
        "pair_extreme_fade_fraction": extreme_fraction,
        "pair_is_stable": pair_is_stable,
        "n_residual_detections": int(len(residual_sources)),
        "n_survivors": int(len(candidates)),
        "n_fade_survivors": int(sum(c["event_sign"] == "fade" for c in candidates)),
        "n_brighten_survivors": int(sum(c["event_sign"] == "brighten" for c in candidates)),
        **product_info,
    }
    return pair_summary, candidates, detection_table


def _failed_pair_row(pair_row: pd.Series, exc: Exception) -> dict[str, Any]:
    return {
        "pair_id": str(pair_row["pair_id"]),
        "galaxy_name": str(pair_row["galaxy_name"]),
        "obs_id_1": str(pair_row["obs_id_1"]),
        "obs_id_2": str(pair_row["obs_id_2"]),
        "filter_1": str(pair_row["filter_1"]),
        "filter_2": str(pair_row["filter_2"]),
        "compatibility": str(pair_row["compatibility"]),
        "baseline_days": float(pair_row["baseline_days"]),
        "pair_score": float(pair_row["pair_score"]),
        "status": "FAILED",
        "error": str(exc),
        "n_residual_detections": 0,
        "n_survivors": 0,
    }


def _skipped_pair_row(pair_row: pd.Series, reason: str) -> dict[str, Any]:
    return {
        "pair_id": str(pair_row["pair_id"]),
        "galaxy_name": str(pair_row["galaxy_name"]),
        "obs_id_1": str(pair_row["obs_id_1"]),
        "obs_id_2": str(pair_row["obs_id_2"]),
        "filter_1": str(pair_row["filter_1"]),
        "filter_2": str(pair_row["filter_2"]),
        "compatibility": str(pair_row["compatibility"]),
        "baseline_days": float(pair_row["baseline_days"]),
        "pair_score": float(pair_row["pair_score"]),
        "status": "SKIPPED",
        "error": str(reason),
        "n_residual_detections": 0,
        "n_survivors": 0,
    }


def _load_pairs_for_search(
    *,
    epoch_pairs_path: Path,
    pair_subset_path: Path | None,
    max_pairs: int,
    per_galaxy: int,
    compatibilities: tuple[str, ...],
    include_cross_collection: bool,
) -> pd.DataFrame:
    if pair_subset_path is not None:
        subset = pd.read_parquet(pair_subset_path) if str(pair_subset_path).endswith(".parquet") else pd.read_csv(pair_subset_path)
        return subset.reset_index(drop=True)
    pairs = pd.read_parquet(epoch_pairs_path)
    return _select_pixel_pairs(
        pairs,
        max_pairs=max_pairs,
        per_galaxy=per_galaxy,
        compatibilities=set(compatibilities),
        include_cross_collection=include_cross_collection,
    )


def run_difference_upgrade(
    *,
    root_dir: Path,
    epoch_pairs_path: Path,
    observation_matrix_path: Path,
    output_dir: Path | None = None,
    pair_subset_path: Path | None = None,
    max_pairs: int = 60,
    per_galaxy: int = 2,
    compatibilities: tuple[str, ...] = ("exact", "very_similar"),
    include_cross_collection: bool = True,
    sign_mode: str = "fade",
    resume: bool = False,
) -> dict[str, Path]:
    out_dir = ensure_dir(output_dir or (root_dir / "difference_upgrade"))
    progress_path = out_dir / "progress.json"
    _write_run_metadata(
        out_dir,
        root_dir,
        {
            "created_utc": _iso(),
            "mode": "difference_upgrade",
            "resume": bool(resume),
            "epoch_pairs_path": str(epoch_pairs_path),
            "observation_matrix_path": str(observation_matrix_path),
            "pair_subset_path": str(pair_subset_path) if pair_subset_path is not None else None,
            "max_pairs": int(max_pairs),
            "per_galaxy": int(per_galaxy),
            "compatibilities": list(compatibilities),
            "include_cross_collection": bool(include_cross_collection),
            "sign_mode": str(sign_mode),
        },
        resume=resume,
    )
    _log(out_dir, f"run start mode=difference_upgrade resume={resume} sign_mode={sign_mode}")
    write_progress(progress_path, "stage01_prepare_queue", "running")
    candidate_max_pairs = max_pairs
    candidate_per_galaxy = per_galaxy
    if pair_subset_path is None:
        if candidate_max_pairs > 0:
            candidate_max_pairs = max(candidate_max_pairs * 8, candidate_max_pairs)
        if candidate_per_galaxy > 0:
            candidate_per_galaxy = max(candidate_per_galaxy * 4, candidate_per_galaxy)
    candidate_pairs = _load_pairs_for_search(
        epoch_pairs_path=epoch_pairs_path,
        pair_subset_path=pair_subset_path,
        max_pairs=candidate_max_pairs,
        per_galaxy=candidate_per_galaxy,
        compatibilities=compatibilities,
        include_cross_collection=include_cross_collection,
    )

    observations = pd.read_parquet(observation_matrix_path)
    obs_meta = observations.drop_duplicates("obs_id").set_index("obs_id").to_dict(orient="index")
    selected, preflight_rejected = _curate_usable_pairs(
        candidate_pairs,
        obs_meta=obs_meta,
        root_dir=root_dir,
        target_pairs=max_pairs,
    )
    candidate_pairs.to_parquet(out_dir / "candidate_pair_pool.parquet", index=False)
    candidate_pairs.to_csv(out_dir / "candidate_pair_pool.csv", index=False)
    preflight_rejected.to_parquet(out_dir / "pair_queue_preflight_skips.parquet", index=False)
    preflight_rejected.to_csv(out_dir / "pair_queue_preflight_skips.csv", index=False)
    selected.to_parquet(out_dir / "pair_queue.parquet", index=False)
    selected.to_csv(out_dir / "pair_queue.csv", index=False)
    write_progress(
        progress_path,
        "stage01_prepare_queue",
        "completed",
        {
            "n_candidate_pairs": int(len(candidate_pairs)),
            "n_pairs_selected": int(len(selected)),
            "n_preflight_rejected": int(len(preflight_rejected)),
        },
    )
    _log(
        out_dir,
        f"queue prepared candidate_pairs={len(candidate_pairs)} selected_pairs={len(selected)} preflight_rejected={len(preflight_rejected)}",
    )
    image_cache: dict[str, tuple[dict[str, Any], ScienceImage]] = {}
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]] = {}

    pair_rows, detections, n_detection_measurements, completed_pair_ids = _load_existing_checkpoint_rows(out_dir) if resume else ([], [], 0, set())
    if resume and completed_pair_ids:
        _log(out_dir, f"resume loaded completed_pairs={len(completed_pair_ids)} detections={len(detections)} measurements={n_detection_measurements}")
    completed_lookup = set(completed_pair_ids)
    total_pairs = int(len(selected))
    write_progress(
        progress_path,
        "stage02_scan_pairs",
        "running",
        {
            "n_pairs_selected": total_pairs,
            "n_pairs_completed": int(len(pair_rows)),
            "n_pairs_remaining": int(max(total_pairs - len(pair_rows), 0)),
            "n_detections": int(len(detections)),
            "n_detection_measurements": int(n_detection_measurements),
        },
    )

    for pair_index, (_, pair) in enumerate(selected.iterrows(), start=1):
        pair_id = str(pair["pair_id"])
        if pair_id in completed_lookup:
            _log(out_dir, f"pair skip resume pair_id={pair_id} index={pair_index}/{total_pairs}")
            continue
        preflight_issue = _pair_preflight_issue(pair, obs_meta, root_dir=root_dir, footprint_cache=footprint_cache)
        if preflight_issue is not None:
            skipped_row = _skipped_pair_row(pair, preflight_issue)
            _write_pair_checkpoint(
                out_dir=out_dir,
                pair_summary=skipped_row,
                candidates=[],
                detection_table=_empty_measurements_frame(),
            )
            pair_rows.append(skipped_row)
            _log(out_dir, f"pair skipped pair_id={pair_id} reason={preflight_issue}")
            summary = _build_detection_outputs(
                out_dir=out_dir,
                pair_rows=pair_rows,
                detections=detections,
                n_detection_measurements=n_detection_measurements,
                sign_mode=sign_mode,
                total_pairs=total_pairs,
            )
            write_progress(
                progress_path,
                "stage02_scan_pairs",
                "running",
                {
                    **summary,
                    "n_pairs_remaining": int(max(total_pairs - len(pair_rows), 0)),
                    "active_pair_id": pair_id,
                    "active_pair_index": int(pair_index),
                },
            )
            continue
        _log(
            out_dir,
            "pair start "
            f"pair_id={pair_id} index={pair_index}/{total_pairs} "
            f"galaxy={pair['galaxy_name']} filter={pair['filter_1']}->{pair['filter_2']} baseline_days={float(pair['baseline_days']):.2f}",
        )
        write_progress(
            progress_path,
            "stage02_scan_pairs",
            "running",
            {
                "n_pairs_selected": total_pairs,
                "n_pairs_completed": int(len(pair_rows)),
                "n_pairs_remaining": int(max(total_pairs - len(pair_rows), 0)),
                "active_pair_id": pair_id,
                "active_pair_index": int(pair_index),
                "active_galaxy_name": str(pair["galaxy_name"]),
                "active_filter_1": str(pair["filter_1"]),
                "active_filter_2": str(pair["filter_2"]),
                "active_baseline_days": float(pair["baseline_days"]),
                "n_detections": int(len(detections)),
                "n_detection_measurements": int(n_detection_measurements),
            },
        )
        try:
            pair_summary, pair_candidates, detection_table = _scan_pair_difference(
                pair,
                obs_meta,
                root_dir=root_dir,
                out_dir=out_dir,
                sign_mode=sign_mode,
                image_cache=image_cache,
            )
            pair_summary["status"] = "SCANNED"
            pair_summary["error"] = None
            _write_pair_checkpoint(
                out_dir=out_dir,
                pair_summary=pair_summary,
                candidates=pair_candidates,
                detection_table=detection_table,
            )
            pair_rows.append(pair_summary)
            detections.extend(pair_candidates)
            n_detection_measurements += int(len(detection_table))
            _log(
                out_dir,
                "pair complete "
                f"pair_id={pair_id} status=SCANNED residuals={int(pair_summary['n_residual_detections'])} survivors={int(pair_summary['n_survivors'])} "
                f"reg_resid_px={pair_summary.get('registration_residual_px')}",
            )
        except Exception as exc:
            failure_row = _failed_pair_row(pair, exc)
            _write_pair_checkpoint(
                out_dir=out_dir,
                pair_summary=failure_row,
                candidates=[],
                detection_table=pd.DataFrame(),
            )
            pair_rows.append(failure_row)
            _log(out_dir, f"pair failed pair_id={pair_id} error={exc}")

        summary = _build_detection_outputs(
            out_dir=out_dir,
            pair_rows=pair_rows,
            detections=detections,
            n_detection_measurements=n_detection_measurements,
            sign_mode=sign_mode,
            total_pairs=total_pairs,
        )
        write_progress(
            progress_path,
            "stage02_scan_pairs",
            "running",
            {
                **summary,
                "n_pairs_remaining": int(max(total_pairs - len(pair_rows), 0)),
                "active_pair_id": pair_id,
                "active_pair_index": int(pair_index),
            },
        )

    write_progress(progress_path, "stage03_aggregate_outputs", "running", {"n_pairs_completed": int(len(pair_rows))})
    n_detection_measurements = _write_detection_measurements(out_dir)
    summary = _build_detection_outputs(
        out_dir=out_dir,
        pair_rows=pair_rows,
        detections=detections,
        n_detection_measurements=n_detection_measurements,
        sign_mode=sign_mode,
        total_pairs=total_pairs,
    )
    write_progress(progress_path, "stage03_aggregate_outputs", "completed", summary)
    _log(out_dir, f"run completed scanned={summary['n_pairs_scanned']} failed={summary['n_pairs_failed']} detections={summary['n_detections']}")
    return {
        "pair_queue": out_dir / "pair_queue.parquet",
        "pair_summary": out_dir / "pair_summary.parquet",
        "pair_summary_csv": out_dir / "pair_summary.csv",
        "detection_measurements": out_dir / "detection_measurements.parquet",
        "detections": out_dir / "detections.parquet",
        "detections_csv": out_dir / "detections.csv",
        "fade_candidates": out_dir / "fade_candidates.parquet",
        "fade_candidates_csv": out_dir / "fade_candidates.csv",
        "summary": out_dir / "summary.json",
    }


def _candidate_stage(epoch_mjd: float, *, pre_mjd: float, post_mjd: float) -> str:
    if not math.isfinite(epoch_mjd):
        return "unknown"
    if epoch_mjd <= pre_mjd + 0.5:
        return "pre"
    if epoch_mjd >= post_mjd - 0.5:
        return "post"
    return "window"


def run_difference_followup(
    *,
    root_dir: Path,
    detections_path: Path,
    observation_matrix_path: Path,
    output_dir: Path | None = None,
    statuses: tuple[str, ...] = ("PASS", "REVIEW"),
    max_candidates: int = 0,
) -> dict[str, Path]:
    out_dir = ensure_dir(output_dir or (root_dir / "difference_upgrade" / "followup"))
    detections = pd.read_parquet(detections_path).copy()
    detections = detections[detections["event_sign"].astype(str).eq("fade")].copy()
    detections = detections[detections["status"].astype(str).isin(set(statuses))].copy()
    detections = detections.sort_values(["status", "diff_sigma_min"], ascending=[True, False]).reset_index(drop=True)
    if max_candidates > 0:
        detections = detections.head(max_candidates).copy()

    observations = pd.read_parquet(observation_matrix_path).copy()
    observations = observations[_imaging_like_mask(observations)].copy()
    observations = observations[observations["is_public"].fillna(False)].copy()

    rows: list[dict[str, Any]] = []
    index_rows: list[dict[str, Any]] = []
    image_cache: dict[str, ScienceImage] = {}

    for _, cand in detections.iterrows():
        world = SkyCoord(float(cand["ra_deg"]) * u.deg, float(cand["dec_deg"]) * u.deg, frame="icrs")
        galaxy_obs = observations[
            observations["galaxy_name"].astype(str).eq(str(cand["galaxy_name"]))
            & observations["filters"].astype(str).eq(str(cand["pre_filter"]))
        ].drop_duplicates("obs_id").copy()
        if galaxy_obs.empty:
            continue
        galaxy_obs = galaxy_obs.sort_values("t_min").reset_index(drop=True)
        ref_row = galaxy_obs[galaxy_obs["obs_id"].astype(str).eq(str(cand["pre_obs_id"]))]
        ref_row = ref_row.iloc[0] if not ref_row.empty else galaxy_obs.iloc[0]

        if str(ref_row["obs_id"]) not in image_cache:
            _, image_cache[str(ref_row["obs_id"])] = _resolve_observation_product(ref_row, root_dir=root_dir)
        ref_image = image_cache[str(ref_row["obs_id"])]

        candidate_dir = ensure_dir(out_dir / str(cand["candidate_id"]))
        candidate_rows: list[dict[str, Any]] = []
        failures = 0
        for _, obs in galaxy_obs.iterrows():
            row = {
                "candidate_id": str(cand["candidate_id"]),
                "galaxy_name": str(cand["galaxy_name"]),
                "reference_obs_id": str(ref_row["obs_id"]),
                "epoch_obs_id": str(obs["obs_id"]),
                "filter_name": str(obs["filters"]),
                "instrument_name": str(obs["instrument_name"]),
                "t_min": float(obs["t_min"]) if pd.notna(obs["t_min"]) else np.nan,
                "stage": _candidate_stage(float(obs["t_min"]) if pd.notna(obs["t_min"]) else float("nan"), pre_mjd=float(cand["pre_mjd"]), post_mjd=float(cand["post_mjd"])),
            }
            try:
                if str(obs["obs_id"]) not in image_cache:
                    _, image_cache[str(obs["obs_id"])] = _resolve_observation_product(obs, root_dir=root_dir)
                cmp_image = image_cache[str(obs["obs_id"])]
                prepared = _prepare_registered_pair(ref_image, cmp_image, center_world=world, cutout_arcsec=90.0)
                x_pix, y_pix = prepared["ref_image"].wcs.world_to_pixel(world)
                if not (np.isfinite(x_pix) and np.isfinite(y_pix)):
                    raise RuntimeError("Candidate position cannot be projected into reference cutout.")
                metrics = _measure_residual_positions(
                    prepared["ref_sub"],
                    prepared["cmp_scaled"],
                    prepared["diff_image"],
                    prepared["valid_mask"],
                    np.array([[float(x_pix), float(y_pix)]], dtype=float),
                    radius=3.5,
                ).iloc[0]
                row.update(
                    {
                        "status": "MEASURED",
                        "diff_flux": float(metrics["flux_diff"]),
                        "diff_err": float(metrics["err_diff"]),
                        "diff_sigma": float(metrics["diff_sigma"]),
                        "flux_ref": float(metrics["flux_pre"]),
                        "flux_epoch": float(metrics["flux_post"]),
                        "post_to_pre_ratio": float(metrics["post_to_pre_ratio"]) if math.isfinite(metrics["post_to_pre_ratio"]) else np.nan,
                        "registration_dx_px": float(prepared["registration"]["dx_px"]),
                        "registration_dy_px": float(prepared["registration"]["dy_px"]),
                        "registration_residual_px": float(prepared["registration"]["residual_px"]) if math.isfinite(prepared["registration"]["residual_px"]) else np.nan,
                        "registration_n_matches": int(prepared["registration"]["n_matches"]),
                        "scale_factor_post_to_pre": float(prepared["scale_factor"]),
                    }
                )
            except Exception as exc:
                failures += 1
                row.update({"status": "FAILED", "error": str(exc)})
            candidate_rows.append(row)
            rows.append(row)

        candidate_df = pd.DataFrame(candidate_rows)
        candidate_df.to_parquet(candidate_dir / "difference_lightcurve.parquet", index=False)
        candidate_df.to_csv(candidate_dir / "difference_lightcurve.csv", index=False)
        summary = {
            "candidate_id": str(cand["candidate_id"]),
            "galaxy_name": str(cand["galaxy_name"]),
            "reference_obs_id": str(ref_row["obs_id"]),
            "n_epochs_considered": int(len(candidate_df)),
            "n_epochs_measured": int((candidate_df["status"] == "MEASURED").sum()) if not candidate_df.empty else 0,
            "n_failures": int(failures),
            "measured_diff_sigma_max": float(candidate_df.loc[candidate_df["status"] == "MEASURED", "diff_sigma"].max()) if not candidate_df.empty and (candidate_df["status"] == "MEASURED").any() else np.nan,
            "measured_ratio_min": float(candidate_df.loc[candidate_df["status"] == "MEASURED", "post_to_pre_ratio"].min()) if not candidate_df.empty and (candidate_df["status"] == "MEASURED").any() else np.nan,
        }
        write_json(candidate_dir / "followup_summary.json", summary)
        index_rows.append(summary)

    all_rows = pd.DataFrame(rows)
    all_rows.to_parquet(out_dir / "followup_measurements.parquet", index=False)
    all_rows.to_csv(out_dir / "followup_measurements.csv", index=False)
    index_df = pd.DataFrame(index_rows)
    index_df.to_csv(out_dir / "followup_index.csv", index=False)
    write_json(
        out_dir / "followup_summary.json",
        {
            "created_utc": utc_stamp(),
            "n_candidates": int(len(index_df)),
            "n_epochs_measured": int((all_rows["status"] == "MEASURED").sum()) if not all_rows.empty else 0,
            "n_epoch_failures": int((all_rows["status"] == "FAILED").sum()) if not all_rows.empty else 0,
        },
    )
    return {
        "followup_measurements": out_dir / "followup_measurements.parquet",
        "followup_index": out_dir / "followup_index.csv",
        "summary": out_dir / "followup_summary.json",
    }


def _canonical_sn_name(text: str) -> str | None:
    text = str(text).upper().replace("_", "-")
    match = re.search(r"(SN-?\d{4}[A-Z0-9]+)", text)
    if not match:
        return None
    return match.group(1).replace("-", "")


def _build_hidden_supernova_truth(observations: pd.DataFrame) -> pd.DataFrame:
    work = observations.copy()
    work["sn_name"] = work["target_name"].map(_canonical_sn_name)
    work = work[work["sn_name"].notna()].copy()
    work = work[_imaging_like_mask(work)].copy()
    truth = (
        work.groupby(["galaxy_name", "sn_name"], as_index=False)
        .agg(
            truth_ra_deg=("s_ra", "median"),
            truth_dec_deg=("s_dec", "median"),
            n_observations=("obs_id", "nunique"),
            first_mjd=("t_min", "min"),
            last_mjd=("t_min", "max"),
        )
        .sort_values(["galaxy_name", "sn_name"])
        .reset_index(drop=True)
    )
    return truth


def _build_benchmark_pairs(observations: pd.DataFrame, epoch_pairs: pd.DataFrame, *, root_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = observations.copy()
    work["sn_name"] = work["target_name"].map(_canonical_sn_name)
    work = work[work["sn_name"].notna()].copy()
    work = work[_imaging_like_mask(work)].copy()
    truth = _build_hidden_supernova_truth(observations)
    obs_meta = observations.drop_duplicates("obs_id").set_index("obs_id").to_dict(orient="index")
    footprint_cache: dict[str, tuple[dict[str, Any], ScienceFootprint]] = {}
    sn_obs = work[["obs_id", "galaxy_name", "sn_name"]].drop_duplicates()

    benchmark_pairs = epoch_pairs.merge(sn_obs.rename(columns={"obs_id": "obs_id_1"}), on=["obs_id_1", "galaxy_name"], how="left")
    benchmark_pairs = benchmark_pairs.merge(sn_obs.rename(columns={"obs_id": "obs_id_2", "sn_name": "sn_name_2"}), on=["obs_id_2", "galaxy_name"], how="left")
    benchmark_pairs["sn_name"] = benchmark_pairs["sn_name"].combine_first(benchmark_pairs["sn_name_2"])
    benchmark_pairs = benchmark_pairs[benchmark_pairs["sn_name"].notna()].copy()
    benchmark_pairs = benchmark_pairs[benchmark_pairs["compatibility"].astype(str).eq("exact")].copy()
    benchmark_pairs = benchmark_pairs[benchmark_pairs["baseline_days"].astype(float) >= 30.0].copy()
    for col in ["instrument_1", "instrument_2"]:
        benchmark_pairs = benchmark_pairs[~benchmark_pairs[col].astype(str).str.contains("NIRSPEC|STIS|COS|FOC|PHOTOMETER|MIRI/IFU", case=False, na=False)].copy()
    for col in ["filter_1", "filter_2"]:
        benchmark_pairs = benchmark_pairs[~benchmark_pairs[col].astype(str).str.contains(";", regex=False, na=False)].copy()
    benchmark_pairs["rank_score"] = benchmark_pairs["pair_score"].astype(float) + 0.02 * benchmark_pairs["baseline_days"].astype(float)
    benchmark_pairs = benchmark_pairs.merge(truth, on=["galaxy_name", "sn_name"], how="left")
    benchmark_pairs["center_ra_deg"] = pd.to_numeric(benchmark_pairs["truth_ra_deg"], errors="coerce")
    benchmark_pairs["center_dec_deg"] = pd.to_numeric(benchmark_pairs["truth_dec_deg"], errors="coerce")
    benchmark_pairs["center_source"] = np.where(
        benchmark_pairs["center_ra_deg"].notna() & benchmark_pairs["center_dec_deg"].notna(),
        "truth_supernova",
        "galaxy_center",
    )
    benchmark_pairs = benchmark_pairs.sort_values(
        ["galaxy_name", "sn_name", "rank_score"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    selected_rows: list[pd.Series] = []
    for _, group in benchmark_pairs.groupby(["galaxy_name", "sn_name"], sort=False):
        chosen: pd.Series | None = None
        for _, row in group.iterrows():
            if _benchmark_pair_covers_truth(row, obs_meta=obs_meta, root_dir=root_dir, footprint_cache=footprint_cache):
                chosen = row
                break
        if chosen is not None:
            selected_rows.append(chosen)

    selected_pairs = pd.DataFrame(selected_rows).reset_index(drop=True) if selected_rows else benchmark_pairs.head(0).copy()
    if not selected_pairs.empty:
        truth = truth.merge(selected_pairs[["galaxy_name", "sn_name"]].drop_duplicates(), on=["galaxy_name", "sn_name"], how="inner")
    else:
        truth = truth.head(0).copy()
    return truth.reset_index(drop=True), selected_pairs.reset_index(drop=True)


def run_supernova_benchmark(
    *,
    root_dir: Path,
    epoch_pairs_path: Path,
    observation_matrix_path: Path,
    output_dir: Path | None = None,
    resume: bool = False,
) -> dict[str, Path]:
    out_dir = ensure_dir(output_dir or (root_dir / "difference_upgrade" / "benchmark"))
    progress_path = out_dir / "progress.json"
    _write_run_metadata(
        out_dir,
        root_dir,
        {
            "created_utc": _iso(),
            "mode": "supernova_benchmark",
            "resume": bool(resume),
            "epoch_pairs_path": str(epoch_pairs_path),
            "observation_matrix_path": str(observation_matrix_path),
        },
        resume=resume,
    )
    _log(out_dir, f"run start mode=supernova_benchmark resume={resume}")
    hidden_dir = ensure_dir(out_dir / "hidden_truth")
    truth_path = hidden_dir / "truth.parquet"
    benchmark_pairs_path = out_dir / "benchmark_pairs.parquet"

    if resume and truth_path.exists() and benchmark_pairs_path.exists():
        truth = pd.read_parquet(truth_path)
        benchmark_pairs = pd.read_parquet(benchmark_pairs_path)
        _log(out_dir, f"benchmark resume reused truth_groups={len(truth)} benchmark_pairs={len(benchmark_pairs)}")
        write_progress(progress_path, "stage01_prepare_truth", "completed", {"n_truth_groups": int(len(truth)), "n_benchmark_pairs": int(len(benchmark_pairs))})
    else:
        write_progress(progress_path, "stage01_prepare_truth", "running")
        observations = pd.read_parquet(observation_matrix_path)
        epoch_pairs = pd.read_parquet(epoch_pairs_path)
        truth, benchmark_pairs = _build_benchmark_pairs(observations, epoch_pairs, root_dir=root_dir)
        truth.to_parquet(truth_path, index=False)
        truth.to_csv(hidden_dir / "truth.csv", index=False)
        benchmark_pairs.to_parquet(benchmark_pairs_path, index=False)
        benchmark_pairs.to_csv(out_dir / "benchmark_pairs.csv", index=False)
        write_progress(progress_path, "stage01_prepare_truth", "completed", {"n_truth_groups": int(len(truth)), "n_benchmark_pairs": int(len(benchmark_pairs))})
        _log(out_dir, f"benchmark prepared truth_groups={len(truth)} benchmark_pairs={len(benchmark_pairs)}")

    write_progress(progress_path, "stage02_scan_pairs", "running", {"n_benchmark_pairs": int(len(benchmark_pairs))})
    search_paths = run_difference_upgrade(
        root_dir=root_dir,
        epoch_pairs_path=epoch_pairs_path,
        observation_matrix_path=observation_matrix_path,
        output_dir=ensure_dir(out_dir / "scan"),
        pair_subset_path=benchmark_pairs_path,
        max_pairs=0,
        per_galaxy=0,
        compatibilities=("exact",),
        include_cross_collection=True,
        sign_mode="both",
        resume=resume,
    )
    write_progress(progress_path, "stage02_scan_pairs", "completed", {"scan_summary_path": str(search_paths["summary"])})

    evaluation_path = out_dir / "benchmark_evaluation.parquet"
    summary_out = out_dir / "benchmark_summary.json"
    report_out = out_dir / "benchmark_report.md"
    if resume and evaluation_path.exists() and summary_out.exists() and report_out.exists():
        summary = read_json(summary_out)
        _log(out_dir, f"benchmark resume reused evaluation recovered={summary.get('n_recovered_truth_groups')}")
        write_progress(progress_path, "stage03_evaluate_recovery", "completed", summary)
    else:
        write_progress(progress_path, "stage03_evaluate_recovery", "running")
        residuals = pd.read_parquet(search_paths["detection_measurements"])
        if residuals.empty:
            detections = pd.DataFrame(columns=["galaxy_name", "ra_deg", "dec_deg", "event_sign", "pair_id", "diff_sigma"])
        else:
            residuals["ra_deg"] = pd.to_numeric(residuals["ra_deg"], errors="coerce")
            residuals["dec_deg"] = pd.to_numeric(residuals["dec_deg"], errors="coerce")
            residuals["diff_sigma"] = pd.to_numeric(residuals["diff_sigma"], errors="coerce")
            residuals["edge_distance_px"] = pd.to_numeric(residuals["edge_distance_px"], errors="coerce")
            detections = residuals[
                residuals["ra_deg"].notna()
                & residuals["dec_deg"].notna()
                & residuals["diff_sigma"].abs().ge(5.0)
                & residuals["edge_distance_px"].fillna(0.0).ge(6.0)
            ].copy()

        eval_rows: list[dict[str, Any]] = []
        for _, row in truth.iterrows():
            group_detections = detections[detections["galaxy_name"].astype(str).eq(str(row["galaxy_name"]))].copy()
            if group_detections.empty:
                eval_rows.append(
                    {
                        "galaxy_name": str(row["galaxy_name"]),
                        "sn_name": str(row["sn_name"]),
                        "recovered": False,
                        "best_sep_arcsec": np.nan,
                        "matched_candidate_id": None,
                        "matched_pair_id": None,
                        "matched_event_sign": None,
                        "matched_diff_sigma": np.nan,
                    }
                )
                continue
            dra = (group_detections["ra_deg"].to_numpy(dtype=float) - float(row["truth_ra_deg"])) * np.cos(np.deg2rad(float(row["truth_dec_deg"])))
            ddec = group_detections["dec_deg"].to_numpy(dtype=float) - float(row["truth_dec_deg"])
            sep_arcsec = np.hypot(dra, ddec) * 3600.0
            best_idx = int(np.nanargmin(sep_arcsec))
            eval_rows.append(
                {
                    "galaxy_name": str(row["galaxy_name"]),
                    "sn_name": str(row["sn_name"]),
                    "recovered": bool(sep_arcsec[best_idx] <= 1.5),
                    "best_sep_arcsec": float(sep_arcsec[best_idx]),
                    "matched_candidate_id": str(group_detections.iloc[best_idx].get("candidate_id", "")) or None,
                    "matched_pair_id": str(group_detections.iloc[best_idx].get("pair_id", "")) or None,
                    "matched_event_sign": str(group_detections.iloc[best_idx]["event_sign"]),
                    "matched_diff_sigma": float(group_detections.iloc[best_idx]["diff_sigma"]),
                }
            )

        evaluation = pd.DataFrame(eval_rows)
        evaluation.to_parquet(evaluation_path, index=False)
        evaluation.to_csv(out_dir / "benchmark_evaluation.csv", index=False)

        summary = {
            "created_utc": utc_stamp(),
            "n_truth_groups": int(len(truth)),
            "n_benchmark_pairs": int(len(benchmark_pairs)),
            "n_truth_groups_with_pairs": int(benchmark_pairs[["galaxy_name", "sn_name"]].drop_duplicates().shape[0]),
            "n_detected_signals": int(len(detections)),
            "n_recovered_truth_groups": int(evaluation["recovered"].sum()) if not evaluation.empty else 0,
            "recovery_fraction": float(evaluation["recovered"].mean()) if not evaluation.empty else 0.0,
            "median_recovered_sep_arcsec": float(evaluation.loc[evaluation["recovered"], "best_sep_arcsec"].median()) if not evaluation.empty and evaluation["recovered"].any() else np.nan,
        }
        write_json(summary_out, summary)
        report_lines = [
            "# Blind Supernova Benchmark",
            "",
            f"- Truth groups: `{summary['n_truth_groups']}`",
            f"- Truth groups with pairable benchmark data: `{summary['n_truth_groups_with_pairs']}`",
            f"- Benchmark pairs scanned: `{summary['n_benchmark_pairs']}`",
            f"- Detected residual signals: `{summary['n_detected_signals']}`",
            f"- Recovered truth groups within 1.5 arcsec: `{summary['n_recovered_truth_groups']}`",
            f"- Recovery fraction: `{summary['recovery_fraction']:.3f}`",
            f"- Median recovered localization error (arcsec): `{summary['median_recovered_sep_arcsec'] if math.isfinite(summary['median_recovered_sep_arcsec']) else float('nan'):.3f}`",
        ]
        write_text(report_out, "\n".join(report_lines) + "\n")
        write_progress(progress_path, "stage03_evaluate_recovery", "completed", summary)
        _log(out_dir, f"benchmark completed recovered={summary['n_recovered_truth_groups']} fraction={summary['recovery_fraction']:.3f}")
    return {
        "truth": truth_path,
        "benchmark_pairs": benchmark_pairs_path,
        "scan_summary": search_paths["summary"],
        "evaluation": evaluation_path,
        "summary": summary_out,
        "report": report_out,
    }
