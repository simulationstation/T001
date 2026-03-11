from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

from .difference_upgrade import _measure_residual_positions, _prepare_registered_pair
from .pixel_search import ScienceImage, _load_science_image
from .utils import (
    angular_sep_deg,
    append_log_line,
    ensure_dir,
    env_like,
    git_like,
    normalize_01,
    utc_stamp,
    write_json,
    write_progress,
    write_text,
)

CLASS_STRONG_EXPORT = "STRONG_EXPORT_FAILURE_LIKE"
CLASS_INTERMEDIATE = "INTERMEDIATE_BRANCH_LIKE"
CLASS_DUST = "DUST_SURVIVOR_LIKE"
CLASS_SYSTEMATIC = "SYSTEMATIC_LIKE"
CLASS_UNRESOLVED = "VARIABLE_OR_UNRESOLVED"


@dataclass(slots=True)
class ClusterGroup:
    cluster_id: str
    galaxy_name: str
    ra_deg: float
    dec_deg: float
    member_indices: list[int]


def _iso() -> str:
    return pd.Timestamp.utcnow().isoformat().replace("+00:00", "Z")


def _log(out_dir: Path, message: str) -> None:
    append_log_line(out_dir / "run.log", f"[{_iso()}] {message}")


def _write_run_metadata(out_dir: Path, root_dir: Path, config: dict[str, Any]) -> None:
    write_json(out_dir / "config.json", config)
    write_json(out_dir / "env.json", env_like())
    write_json(out_dir / "git_like.json", git_like(root_dir))


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


def _json_list(values: list[str]) -> str:
    return json.dumps(sorted({str(value) for value in values}))


def _clip01(values: Any) -> Any:
    return np.clip(values, 0.0, 1.0)


def _frac(numerator: float, denominator: float, *, default: float = 0.0) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator <= 0:
        return default
    return float(numerator / denominator)


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


def _weighted_ra_dec(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    weights = np.asarray(weights, dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)
    dec_rad = np.deg2rad(dec_deg)
    x = np.cos(dec_rad) * np.cos(np.deg2rad(ra_deg))
    y = np.cos(dec_rad) * np.sin(np.deg2rad(ra_deg))
    z = np.sin(dec_rad)
    vec = np.array(
        [
            np.sum(weights * x),
            np.sum(weights * y),
            np.sum(weights * z),
        ],
        dtype=float,
    )
    norm = float(np.linalg.norm(vec))
    if not math.isfinite(norm) or norm <= 0:
        return float(np.nanmedian(ra_deg)), float(np.nanmedian(dec_deg))
    vec /= norm
    ra = math.degrees(math.atan2(vec[1], vec[0])) % 360.0
    dec = math.degrees(math.asin(np.clip(vec[2], -1.0, 1.0)))
    return float(ra), float(dec)


def _slugify(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value)).strip("_") or "item"


def _band_group(filter_name: str, instrument_name: str) -> str:
    inst = str(instrument_name).upper()
    filt = str(filter_name).upper()
    if "MIRI" in inst:
        return "MIR"
    if "NIRCAM" in inst or "NIRISS" in inst or "WFC3/IR" in inst or "NICMOS" in inst:
        return "NIR"
    if "ACS" in inst or "WFC3/UVIS" in inst or "WFPC" in inst:
        return "OPTICAL"
    if filt.startswith("F") and any(code in filt for code in ["555", "606", "625", "775", "814", "438"]):
        return "OPTICAL"
    return "UNKNOWN"


def _cluster_single_galaxy(sub: pd.DataFrame, *, radius_arcsec: float) -> list[ClusterGroup]:
    if sub.empty:
        return []
    ra = sub["ra_deg"].to_numpy(dtype=float)
    dec = sub["dec_deg"].to_numpy(dtype=float)
    ra0 = float(np.nanmedian(ra))
    dec0 = float(np.nanmedian(dec))
    x = (ra - ra0) * 3600.0 * math.cos(math.radians(dec0))
    y = (dec - dec0) * 3600.0
    pts = np.column_stack([x, y])
    tree = cKDTree(pts)
    parent = np.arange(len(sub), dtype=int)

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a: int, b: int) -> None:
        ra_root = find(a)
        rb_root = find(b)
        if ra_root != rb_root:
            parent[rb_root] = ra_root

    for i, j in tree.query_pairs(radius_arcsec):
        union(i, j)

    components: dict[int, list[int]] = {}
    for idx in range(len(sub)):
        components.setdefault(find(idx), []).append(idx)

    clusters: list[ClusterGroup] = []
    galaxy_name = str(sub["galaxy_name"].iloc[0])
    for serial, members in enumerate(
        sorted(components.values(), key=lambda rows: (-len(rows), rows[0])),
        start=1,
    ):
        member_frame = sub.iloc[members]
        weights = pd.to_numeric(member_frame["diff_sigma_min"], errors="coerce").fillna(
            pd.to_numeric(member_frame["diff_sigma"], errors="coerce")
        )
        center_ra, center_dec = _weighted_ra_dec(
            member_frame["ra_deg"].to_numpy(dtype=float),
            member_frame["dec_deg"].to_numpy(dtype=float),
            np.clip(weights.to_numpy(dtype=float), 1.0, np.inf),
        )
        clusters.append(
            ClusterGroup(
                cluster_id=_slugify(f"{galaxy_name}_branch_{serial:04d}"),
                galaxy_name=galaxy_name,
                ra_deg=center_ra,
                dec_deg=center_dec,
                member_indices=list(member_frame.index),
            )
        )
    return clusters


def _build_cluster_tables(detections: pd.DataFrame, *, radius_arcsec: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    member_rows: list[dict[str, Any]] = []
    cluster_rows: list[dict[str, Any]] = []
    for _, sub in detections.groupby("galaxy_name", sort=True):
        clusters = _cluster_single_galaxy(sub.reset_index(), radius_arcsec=radius_arcsec)
        sub_reset = sub.reset_index()
        for cluster in clusters:
            members = sub_reset.iloc[cluster.member_indices].copy()
            reason_values = members["reason_codes_json"].map(_safe_json_list)
            warning_values = members["warning_flags_json"].map(_safe_json_list)
            reason_union = sorted({code for codes in reason_values for code in codes})
            warning_union = sorted({code for codes in warning_values for code in codes})
            weights = pd.to_numeric(members["diff_sigma_min"], errors="coerce").fillna(
                pd.to_numeric(members["diff_sigma"], errors="coerce")
            ).clip(lower=1.0)
            sep_arcsec = angular_sep_deg(
                cluster.ra_deg,
                cluster.dec_deg,
                members["ra_deg"].to_numpy(dtype=float),
                members["dec_deg"].to_numpy(dtype=float),
            ) * 3600.0
            for _, member in members.iterrows():
                member_rows.append(
                    {
                        "cluster_id": cluster.cluster_id,
                        "candidate_id": str(member["candidate_id"]),
                        "pair_id": str(member["pair_id"]),
                        "galaxy_name": str(member["galaxy_name"]),
                        "status": str(member["status"]),
                        "ra_deg": float(member["ra_deg"]),
                        "dec_deg": float(member["dec_deg"]),
                        "cluster_center_ra_deg": float(cluster.ra_deg),
                        "cluster_center_dec_deg": float(cluster.dec_deg),
                        "cluster_sep_arcsec": float(
                            angular_sep_deg(
                                cluster.ra_deg,
                                cluster.dec_deg,
                                np.asarray([float(member["ra_deg"])]),
                                np.asarray([float(member["dec_deg"])]),
                            )[0]
                            * 3600.0
                        ),
                        "reason_codes_json": str(member["reason_codes_json"]),
                        "warning_flags_json": str(member["warning_flags_json"]),
                        "cutout_path": str(member["cutout_path"]),
                    }
                )
            cluster_rows.append(
                {
                    "cluster_id": cluster.cluster_id,
                    "galaxy_name": cluster.galaxy_name,
                    "ra_deg": float(cluster.ra_deg),
                    "dec_deg": float(cluster.dec_deg),
                    "n_members": int(len(members)),
                    "n_pass_members": int((members["status"].astype(str) == "PASS").sum()),
                    "n_review_members": int((members["status"].astype(str) == "REVIEW").sum()),
                    "n_support_pairs": int(members["pair_id"].nunique()),
                    "n_support_filters": int(members["pre_filter"].nunique()),
                    "n_support_instruments": int(members["pre_instrument"].nunique()),
                    "detection_baseline_span_days": float(
                        pd.to_numeric(members["post_mjd"], errors="coerce").max()
                        - pd.to_numeric(members["pre_mjd"], errors="coerce").min()
                    ),
                    "member_diff_sigma_max": float(pd.to_numeric(members["diff_sigma"], errors="coerce").max()),
                    "member_diff_sigma_minmax": float(pd.to_numeric(members["diff_sigma_min"], errors="coerce").max()),
                    "member_diff_sigma_median": float(pd.to_numeric(members["diff_sigma_min"], errors="coerce").median()),
                    "member_ratio_median": float(pd.to_numeric(members["post_to_pre_ratio"], errors="coerce").median()),
                    "member_ratio_min": float(pd.to_numeric(members["post_to_pre_ratio"], errors="coerce").min()),
                    "member_ratio_max": float(pd.to_numeric(members["post_to_pre_ratio"], errors="coerce").max()),
                    "member_severe_fraction": float((pd.to_numeric(members["post_to_pre_ratio"], errors="coerce") <= 0.05).mean()),
                    "member_partial_fraction": float(
                        (
                            (pd.to_numeric(members["post_to_pre_ratio"], errors="coerce") > 0.05)
                            & (pd.to_numeric(members["post_to_pre_ratio"], errors="coerce") <= 0.50)
                        ).mean()
                    ),
                    "member_pre_snr_median": float(pd.to_numeric(members["pre_snr"], errors="coerce").median()),
                    "member_post_snr_median": float(pd.to_numeric(members["post_snr"], errors="coerce").median()),
                    "member_blend_risk_median": float(pd.to_numeric(members["blend_risk"], errors="coerce").median()),
                    "member_crowding_median": float(pd.to_numeric(members["crowding_index"], errors="coerce").median()),
                    "member_edge_distance_median": float(pd.to_numeric(members["edge_distance_px"], errors="coerce").median()),
                    "member_registration_residual_median": float(
                        pd.to_numeric(members["registration_residual_px"], errors="coerce").median()
                    ),
                    "member_registration_matches_median": float(
                        pd.to_numeric(members["registration_n_matches"], errors="coerce").median()
                    ),
                    "member_pair_stable_fraction_median": float(
                        pd.to_numeric(members["pair_stable_fraction"], errors="coerce").median()
                    ),
                    "member_pair_extreme_fade_fraction_median": float(
                        pd.to_numeric(members["pair_extreme_fade_fraction"], errors="coerce").median()
                    ),
                    "member_centroid_scatter_arcsec": float(np.nanmedian(sep_arcsec)) if len(sep_arcsec) else float("nan"),
                    "member_reason_codes_json": _json_list(reason_union),
                    "member_warning_flags_json": _json_list(warning_union),
                    "member_candidate_ids_json": _json_list(members["candidate_id"].astype(str).tolist()),
                    "member_pair_ids_json": _json_list(members["pair_id"].astype(str).tolist()),
                }
            )
    cluster_df = pd.DataFrame(cluster_rows).sort_values(
        ["n_pass_members", "n_support_pairs", "member_diff_sigma_minmax"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    members_df = pd.DataFrame(member_rows).sort_values(["cluster_id", "status", "cluster_sep_arcsec"]).reset_index(drop=True)
    return cluster_df, members_df


def _load_pair_catalog(strict_output_dir: Path) -> pd.DataFrame:
    pair_queue = pd.read_parquet(strict_output_dir / "pair_queue.parquet").copy()
    pair_summary = pd.read_parquet(strict_output_dir / "pair_summary.parquet").copy()
    pair_catalog = pair_queue.merge(
        pair_summary[
            [
                "pair_id",
                "status",
                "error",
                "registration_dx_px",
                "registration_dy_px",
                "registration_residual_px",
                "registration_n_matches",
                "scale_factor_post_to_pre",
                "n_scale_sources",
                "pair_stable_fraction",
                "pair_extreme_fade_fraction",
                "pair_is_stable",
                "n_residual_detections",
                "n_survivors",
                "product_filename_1",
                "product_filename_2",
                "product_path_1",
                "product_path_2",
                "image_ext_1",
                "image_ext_2",
            ]
        ],
        on="pair_id",
        how="left",
        suffixes=("", "_summary"),
    )
    return pair_catalog


def _load_pair_images(pair_row: pd.Series, image_cache: dict[str, ScienceImage]) -> tuple[ScienceImage, ScienceImage]:
    path_1 = Path(str(pair_row["product_path_1"]))
    path_2 = Path(str(pair_row["product_path_2"]))
    key_1 = f"{pair_row['obs_id_1']}::{path_1}"
    key_2 = f"{pair_row['obs_id_2']}::{path_2}"
    image_1 = image_cache.get(key_1)
    if image_1 is None:
        image_1 = _load_science_image(
            path_1,
            obs_collection=str(pair_row["obs_collection_1"]),
            obs_id=str(pair_row["obs_id_1"]),
            filter_name=str(pair_row["filter_1"]),
            instrument_name=str(pair_row["instrument_1"]),
        )
        image_cache[key_1] = image_1
    image_2 = image_cache.get(key_2)
    if image_2 is None:
        image_2 = _load_science_image(
            path_2,
            obs_collection=str(pair_row["obs_collection_2"]),
            obs_id=str(pair_row["obs_id_2"]),
            filter_name=str(pair_row["filter_2"]),
            instrument_name=str(pair_row["instrument_2"]),
        )
        image_cache[key_2] = image_2
    return image_1, image_2


def _force_measure_cluster_centers(
    *,
    root_dir: Path,
    strict_output_dir: Path,
    output_dir: Path,
    clusters: pd.DataFrame,
    pair_catalog: pd.DataFrame,
) -> pd.DataFrame:
    progress_path = output_dir / "progress.json"
    scanned_pairs = pair_catalog[pair_catalog["status"].astype(str) == "SCANNED"].copy().reset_index(drop=True)
    cluster_by_galaxy = {key: frame.reset_index(drop=True) for key, frame in clusters.groupby("galaxy_name", sort=False)}
    image_cache: dict[str, ScienceImage] = {}
    rows: list[dict[str, Any]] = []

    write_progress(
        progress_path,
        "stage01_forced_measurements",
        "running",
        {
            "n_clusters": int(len(clusters)),
            "n_scanned_pairs": int(len(scanned_pairs)),
            "n_forced_rows_written": 0,
        },
    )
    for pair_index, (_, pair) in enumerate(scanned_pairs.iterrows(), start=1):
        galaxy_clusters = cluster_by_galaxy.get(str(pair["galaxy_name"]))
        if galaxy_clusters is None or galaxy_clusters.empty:
            continue
        _log(
            output_dir,
            f"forced pair start pair_id={pair['pair_id']} index={pair_index}/{len(scanned_pairs)} galaxy={pair['galaxy_name']} clusters={len(galaxy_clusters)}",
        )
        try:
            image_1, image_2 = _load_pair_images(pair, image_cache)
            prepared = _prepare_registered_pair(image_1, image_2, center_world=None, cutout_arcsec=0.0)
            ref_image = prepared["ref_image"]
            x_pix, y_pix = ref_image.wcs.world_to_pixel(
                SkyCoord(
                    galaxy_clusters["ra_deg"].to_numpy(dtype=float) * u.deg,
                    galaxy_clusters["dec_deg"].to_numpy(dtype=float) * u.deg,
                    frame="icrs",
                )
            )
            width = ref_image.data.shape[1]
            height = ref_image.data.shape[0]
            inside = (
                np.isfinite(x_pix)
                & np.isfinite(y_pix)
                & (x_pix >= 0.0)
                & (y_pix >= 0.0)
                & (x_pix < width)
                & (y_pix < height)
            )
            xi = np.clip(np.round(x_pix).astype(int, copy=False), 0, width - 1)
            yi = np.clip(np.round(y_pix).astype(int, copy=False), 0, height - 1)
            mask_ok = np.zeros(len(galaxy_clusters), dtype=bool)
            edge_distance = np.full(len(galaxy_clusters), np.nan, dtype=float)
            mask_ok[inside] = prepared["valid_mask"][yi[inside], xi[inside]]
            edge_distance[inside] = prepared["edge_distance_map"][yi[inside], xi[inside]]
            measurable = inside & mask_ok & np.isfinite(edge_distance) & (edge_distance >= 4.0)
            metrics = pd.DataFrame()
            if measurable.any():
                positions = np.column_stack([x_pix[measurable], y_pix[measurable]])
                metrics = _measure_residual_positions(
                    prepared["ref_sub"],
                    prepared["cmp_scaled"],
                    prepared["diff_image"],
                    prepared["valid_mask"],
                    positions,
                    radius=3.5,
                )
                metrics.index = np.flatnonzero(measurable)
            for idx, cluster in galaxy_clusters.iterrows():
                status = "MEASURED" if measurable[idx] else ("MASKED_OR_EDGE" if inside[idx] else "OUTSIDE_FOOTPRINT")
                row = {
                    "cluster_id": str(cluster["cluster_id"]),
                    "galaxy_name": str(cluster["galaxy_name"]),
                    "pair_id": str(pair["pair_id"]),
                    "pre_obs_id": str(pair["obs_id_1"]),
                    "post_obs_id": str(pair["obs_id_2"]),
                    "pre_filter": str(pair["filter_1"]),
                    "post_filter": str(pair["filter_2"]),
                    "pre_instrument": str(pair["instrument_1"]),
                    "post_instrument": str(pair["instrument_2"]),
                    "band_group": _band_group(str(pair["filter_1"]), str(pair["instrument_1"])),
                    "baseline_days": float(pair["baseline_days"]),
                    "pre_mjd": float("nan"),
                    "post_mjd": float("nan"),
                    "coverage_status": status,
                    "x_pix": float(x_pix[idx]) if np.isfinite(x_pix[idx]) else float("nan"),
                    "y_pix": float(y_pix[idx]) if np.isfinite(y_pix[idx]) else float("nan"),
                    "edge_distance_px": float(edge_distance[idx]) if np.isfinite(edge_distance[idx]) else float("nan"),
                    "registration_residual_px": float(pair["registration_residual_px"]) if pd.notna(pair["registration_residual_px"]) else float("nan"),
                    "registration_n_matches": float(pair["registration_n_matches"]) if pd.notna(pair["registration_n_matches"]) else float("nan"),
                    "pair_stable_fraction": float(pair["pair_stable_fraction"]) if pd.notna(pair["pair_stable_fraction"]) else float("nan"),
                    "pair_extreme_fade_fraction": float(pair["pair_extreme_fade_fraction"]) if pd.notna(pair["pair_extreme_fade_fraction"]) else float("nan"),
                    "pair_failure_reason": None,
                }
                if measurable[idx]:
                    metric = metrics.loc[idx]
                    row.update(
                        {
                            "diff_flux": float(metric["flux_diff"]),
                            "diff_err": float(metric["err_diff"]),
                            "diff_sigma": float(metric["diff_sigma"]),
                            "flux_pre": float(metric["flux_pre"]),
                            "flux_post": float(metric["flux_post"]),
                            "err_pre": float(metric["err_pre"]),
                            "err_post": float(metric["err_post"]),
                            "post_to_pre_ratio": float(metric["post_to_pre_ratio"]) if math.isfinite(metric["post_to_pre_ratio"]) else float("nan"),
                        }
                    )
                else:
                    row.update(
                        {
                            "diff_flux": float("nan"),
                            "diff_err": float("nan"),
                            "diff_sigma": float("nan"),
                            "flux_pre": float("nan"),
                            "flux_post": float("nan"),
                            "err_pre": float("nan"),
                            "err_post": float("nan"),
                            "post_to_pre_ratio": float("nan"),
                        }
                    )
                rows.append(row)
        except Exception as exc:
            _log(output_dir, f"forced pair failed pair_id={pair['pair_id']} error={exc}")
            for _, cluster in galaxy_clusters.iterrows():
                rows.append(
                    {
                        "cluster_id": str(cluster["cluster_id"]),
                        "galaxy_name": str(cluster["galaxy_name"]),
                        "pair_id": str(pair["pair_id"]),
                        "pre_obs_id": str(pair["obs_id_1"]),
                        "post_obs_id": str(pair["obs_id_2"]),
                        "pre_filter": str(pair["filter_1"]),
                        "post_filter": str(pair["filter_2"]),
                        "pre_instrument": str(pair["instrument_1"]),
                        "post_instrument": str(pair["instrument_2"]),
                        "band_group": _band_group(str(pair["filter_1"]), str(pair["instrument_1"])),
                        "baseline_days": float(pair["baseline_days"]),
                        "pre_mjd": float("nan"),
                        "post_mjd": float("nan"),
                        "coverage_status": "PAIR_PREP_FAILED",
                        "x_pix": float("nan"),
                        "y_pix": float("nan"),
                        "edge_distance_px": float("nan"),
                        "registration_residual_px": float("nan"),
                        "registration_n_matches": float("nan"),
                        "pair_stable_fraction": float("nan"),
                        "pair_extreme_fade_fraction": float("nan"),
                        "pair_failure_reason": str(exc),
                        "diff_flux": float("nan"),
                        "diff_err": float("nan"),
                        "diff_sigma": float("nan"),
                        "flux_pre": float("nan"),
                        "flux_post": float("nan"),
                        "err_pre": float("nan"),
                        "err_post": float("nan"),
                        "post_to_pre_ratio": float("nan"),
                    }
                )
        write_progress(
            progress_path,
            "stage01_forced_measurements",
            "running",
            {
                "n_clusters": int(len(clusters)),
                "n_scanned_pairs": int(len(scanned_pairs)),
                "n_pairs_processed": int(pair_index),
                "n_forced_rows_written": int(len(rows)),
                "active_pair_id": str(pair["pair_id"]),
            },
        )

    forced = pd.DataFrame(rows).sort_values(["cluster_id", "pair_id"]).reset_index(drop=True)
    forced.to_parquet(output_dir / "forced_measurements.parquet", index=False)
    forced.to_csv(output_dir / "forced_measurements.csv", index=False)
    write_progress(
        progress_path,
        "stage01_forced_measurements",
        "completed",
        {
            "n_clusters": int(len(clusters)),
            "n_scanned_pairs": int(len(scanned_pairs)),
            "n_forced_rows_written": int(len(forced)),
        },
    )
    return forced


def _per_band_support(measured: pd.DataFrame, band_group: str) -> tuple[float, float]:
    band = measured[measured["band_group"].astype(str) == band_group]
    if band.empty:
        return float("nan"), 0.0
    low_ratio = float((pd.to_numeric(band["post_to_pre_ratio"], errors="coerce") <= 0.35).mean())
    positive = float((pd.to_numeric(band["diff_sigma"], errors="coerce") >= 2.0).mean())
    return float(np.mean([low_ratio, positive])), float(len(band))


def _late_stage_fractions(measured: pd.DataFrame) -> tuple[float, float]:
    if measured.empty:
        return 0.0, 0.0
    ordered = measured.sort_values("baseline_days")
    tail = ordered.iloc[max(int(len(ordered) * 2 / 3), 0) :]
    if tail.empty:
        tail = ordered
    fade = (
        (pd.to_numeric(tail["diff_sigma"], errors="coerce") >= 2.0)
        & (pd.to_numeric(tail["post_to_pre_ratio"], errors="coerce") <= 0.50)
    ).mean()
    ret = (
        (pd.to_numeric(tail["diff_sigma"], errors="coerce") <= -2.0)
        | (pd.to_numeric(tail["post_to_pre_ratio"], errors="coerce") > 0.60)
    ).mean()
    return float(fade), float(ret)


def _reason_codes(row: pd.Series, *, benchmark_ok: bool) -> list[str]:
    codes: list[str] = []
    if benchmark_ok:
        codes.append("BENCHMARK_GUARDRAIL_OK")
    if int(row["n_support_pairs"]) >= 2:
        codes.append("MULTI_PAIR_LOCALIZATION")
    if float(row["detection_baseline_span_days"]) >= 365.0:
        codes.append("LONG_BASELINE_SUPPORT")
    if float(row["forced_photometry_coherence"]) >= 0.60:
        codes.append("FORCED_COHERENCE_STRONG")
    if float(row["fade_persistence"]) >= 0.60:
        codes.append("PERSISTENT_FADE_SUPPORT")
    if float(row["partial_branch_score"]) >= 0.55:
        codes.append("PARTIAL_SUPPRESSION_PATTERN")
    if float(row["mir_fade_support"]) >= 0.45:
        codes.append("MIR_FADE_SUPPORT")
    if float(row["dust_survivor_score"]) >= 0.55:
        codes.append("DUST_COMPETITOR_PRESENT")
    if float(row["rebrightening_penalty"]) >= 0.40:
        codes.append("REBRIGHTENING_OR_RETURN")
    if float(row["artifact_penalty"]) >= 0.50:
        codes.append("ARTIFACT_RISK_HIGH")
    if float(row["coverage_completeness"]) < 0.35:
        codes.append("COVERAGE_LIMITED")
    if float(row["registration_confidence"]) < 0.45:
        codes.append("REGISTRATION_WEAK")
    if float(row["detector_confidence"]) >= 0.75:
        codes.append("DETECTOR_SIGNAL_STRONG")
    return codes


def _classify_cluster(row: pd.Series) -> tuple[str, str]:
    export_score = float(row["export_failure_score"])
    intermediate_score = float(row["partial_branch_score"])
    dust_score = float(row["dust_survivor_score"])
    systematic_score = float(row["systematic_risk"])
    unresolved_score = float(row["unresolved_variability_score"])

    score_map = [
        (CLASS_STRONG_EXPORT, export_score),
        (CLASS_INTERMEDIATE, intermediate_score),
        (CLASS_DUST, dust_score),
        (CLASS_SYSTEMATIC, systematic_score),
        (CLASS_UNRESOLVED, unresolved_score),
    ]
    ordered = sorted(score_map, key=lambda item: item[1], reverse=True)
    best_label, best_score = ordered[0]
    second_score = ordered[1][1] if len(ordered) > 1 else 0.0
    margin = best_score - second_score
    coverage = float(row.get("coverage_completeness", 0.0) or 0.0)
    forced_pairs_measured = int(row.get("forced_n_pairs_measured", 0) or 0)
    n_support_pairs = int(row.get("n_support_pairs", 0) or 0)
    n_pass_members = int(row.get("n_pass_members", 0) or 0)
    support_ok = coverage >= 0.45 and forced_pairs_measured >= 2
    strong_support = n_support_pairs >= 3 or n_pass_members >= 1
    intermediate_support = n_support_pairs >= 2 or n_pass_members >= 1

    if systematic_score >= 0.62 and systematic_score >= export_score + 0.05:
        label = CLASS_SYSTEMATIC
    elif dust_score >= 0.58 and dust_score >= export_score + 0.05:
        label = CLASS_DUST
    elif (
        export_score >= 0.72
        and float(row["coverage_completeness"]) >= 0.50
        and float(row["forced_photometry_coherence"]) >= 0.55
        and float(row["registration_confidence"]) >= 0.55
        and float(row["fade_persistence"]) >= 0.55
        and float(row["artifact_penalty"]) <= 0.40
        and strong_support
    ):
        label = CLASS_STRONG_EXPORT
    elif (
        intermediate_score >= 0.55
        and float(row["coverage_completeness"]) >= 0.45
        and float(row["forced_photometry_coherence"]) >= 0.45
        and float(row["artifact_penalty"]) <= 0.45
        and float(row["dust_survivor_score"]) <= 0.55
        and intermediate_support
    ):
        label = CLASS_INTERMEDIATE
    else:
        label = CLASS_UNRESOLVED

    if margin >= 0.16 and support_ok and intermediate_support:
        confidence = "HIGH"
    elif margin >= 0.08 and float(row["coverage_completeness"]) >= 0.30:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"
    return label, confidence


def _build_cluster_scores(
    *,
    clusters: pd.DataFrame,
    members: pd.DataFrame,
    forced: pd.DataFrame,
    pair_catalog: pd.DataFrame,
    benchmark_summary: dict[str, Any],
) -> pd.DataFrame:
    galaxy_pair_stats = (
        pair_catalog.groupby("galaxy_name", dropna=False)["status"]
        .value_counts()
        .unstack(fill_value=0)
        .rename(columns=str)
        .reset_index()
    )
    galaxy_pair_stats["n_pairs_selected_galaxy"] = galaxy_pair_stats.drop(columns=["galaxy_name"]).sum(axis=1)
    galaxy_pair_stats["pair_failure_fraction_galaxy"] = galaxy_pair_stats.apply(
        lambda row: _frac(float(row.get("FAILED", 0.0)), float(row["n_pairs_selected_galaxy"])),
        axis=1,
    )
    galaxy_pair_stats["pair_skip_fraction_galaxy"] = galaxy_pair_stats.apply(
        lambda row: _frac(float(row.get("SKIPPED", 0.0)), float(row["n_pairs_selected_galaxy"])),
        axis=1,
    )

    forced_rows: list[dict[str, Any]] = []
    for cluster_id, cluster_forced in forced.groupby("cluster_id", dropna=False):
        measured = cluster_forced[cluster_forced["coverage_status"].astype(str) == "MEASURED"].copy()
        n_considered = int(len(cluster_forced))
        n_measured = int(len(measured))
        positive = measured[
            (pd.to_numeric(measured["diff_sigma"], errors="coerce") >= 2.0)
            & (pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce") <= 0.50)
        ]
        severe = measured[
            (pd.to_numeric(measured["diff_sigma"], errors="coerce") >= 3.0)
            & (pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce") <= 0.05)
        ]
        partial = measured[
            (pd.to_numeric(measured["diff_sigma"], errors="coerce") >= 2.0)
            & (pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce") > 0.05)
            & (pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce") <= 0.50)
        ]
        returns = measured[
            (pd.to_numeric(measured["diff_sigma"], errors="coerce") <= -2.0)
            | (pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce") > 0.60)
        ]
        late_fade_fraction, late_return_fraction = _late_stage_fractions(measured)
        optical_support, optical_n = _per_band_support(measured, "OPTICAL")
        nir_support, nir_n = _per_band_support(measured, "NIR")
        mir_support, mir_n = _per_band_support(measured, "MIR")
        support_rows = {
            "cluster_id": str(cluster_id),
            "forced_n_pairs_considered": n_considered,
            "forced_n_pairs_measured": n_measured,
            "forced_n_no_coverage": int((cluster_forced["coverage_status"].astype(str) == "OUTSIDE_FOOTPRINT").sum()),
            "forced_n_masked_or_edge": int((cluster_forced["coverage_status"].astype(str) == "MASKED_OR_EDGE").sum()),
            "forced_n_pair_prep_failed": int((cluster_forced["coverage_status"].astype(str) == "PAIR_PREP_FAILED").sum()),
            "forced_positive_fraction": _frac(float(len(positive)), float(n_measured)),
            "forced_severe_fraction": _frac(float(len(severe)), float(n_measured)),
            "forced_partial_fraction": _frac(float(len(partial)), float(n_measured)),
            "forced_return_fraction": _frac(float(len(returns)), float(n_measured)),
            "forced_ratio_median": float(pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce").median()) if n_measured else float("nan"),
            "forced_ratio_iqr": float(
                pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce").quantile(0.75)
                - pd.to_numeric(measured["post_to_pre_ratio"], errors="coerce").quantile(0.25)
            )
            if n_measured
            else float("nan"),
            "forced_diff_sigma_median": float(pd.to_numeric(measured["diff_sigma"], errors="coerce").median()) if n_measured else float("nan"),
            "forced_diff_sigma_max": float(pd.to_numeric(measured["diff_sigma"], errors="coerce").max()) if n_measured else float("nan"),
            "forced_registration_residual_median": float(
                pd.to_numeric(measured["registration_residual_px"], errors="coerce").median()
            )
            if n_measured
            else float("nan"),
            "forced_registration_matches_median": float(
                pd.to_numeric(measured["registration_n_matches"], errors="coerce").median()
            )
            if n_measured
            else float("nan"),
            "forced_pair_stable_fraction_median": float(
                pd.to_numeric(measured["pair_stable_fraction"], errors="coerce").median()
            )
            if n_measured
            else float("nan"),
            "forced_pair_extreme_fade_fraction_median": float(
                pd.to_numeric(measured["pair_extreme_fade_fraction"], errors="coerce").median()
            )
            if n_measured
            else float("nan"),
            "forced_baseline_span_days": float(pd.to_numeric(measured["baseline_days"], errors="coerce").max()) if n_measured else float("nan"),
            "forced_late_fade_fraction": late_fade_fraction,
            "forced_late_return_fraction": late_return_fraction,
            "optical_fade_support": optical_support,
            "nir_fade_support": nir_support,
            "mir_fade_support": mir_support,
            "optical_measured_pairs": optical_n,
            "nir_measured_pairs": nir_n,
            "mir_measured_pairs": mir_n,
        }
        forced_rows.append(support_rows)
    forced_summary = pd.DataFrame(forced_rows)

    scored = clusters.merge(forced_summary, on="cluster_id", how="left").merge(
        galaxy_pair_stats[
            [
                "galaxy_name",
                "n_pairs_selected_galaxy",
                "pair_failure_fraction_galaxy",
                "pair_skip_fraction_galaxy",
            ]
        ],
        on="galaxy_name",
        how="left",
    )

    scored["support_pair_score"] = _norm_series(scored["n_support_pairs"], default=0.0)
    scored["member_signal_score"] = _norm_series(scored["member_diff_sigma_minmax"], log=True, default=0.0)
    scored["member_ratio_loss_score"] = _clip01(1.0 - pd.to_numeric(scored["member_ratio_median"], errors="coerce").fillna(1.0))
    scored["pass_fraction"] = scored.apply(
        lambda row: _frac(float(row["n_pass_members"]), float(row["n_members"])),
        axis=1,
    )
    scored["member_pre_snr_score"] = _norm_series(scored["member_pre_snr_median"], log=True, default=0.0)
    scored["forced_measurement_score"] = _norm_series(scored["forced_n_pairs_measured"], log=True, default=0.0)
    scored["measured_fraction_score"] = scored.apply(
        lambda row: _frac(float(row["forced_n_pairs_measured"]), float(row["forced_n_pairs_considered"])),
        axis=1,
    )
    scored["forced_ratio_loss_score"] = _clip01(1.0 - pd.to_numeric(scored["forced_ratio_median"], errors="coerce").fillna(1.0))
    scored["ratio_consistency_score"] = 1.0 - normalize_01(
        pd.to_numeric(scored["forced_ratio_iqr"], errors="coerce").clip(lower=0.0, upper=1.5)
    ).fillna(0.0)
    scored["member_reg_score"] = 1.0 - normalize_01(
        pd.to_numeric(scored["member_registration_residual_median"], errors="coerce").clip(lower=0.0, upper=1.2)
    ).fillna(0.0)
    scored["forced_reg_score"] = 1.0 - normalize_01(
        pd.to_numeric(scored["forced_registration_residual_median"], errors="coerce").clip(lower=0.0, upper=1.2)
    ).fillna(0.0)
    scored["anchor_score"] = _norm_series(scored["forced_registration_matches_median"].combine_first(scored["member_registration_matches_median"]), log=True, default=0.0)
    scored["cluster_scatter_penalty"] = normalize_01(
        pd.to_numeric(scored["member_centroid_scatter_arcsec"], errors="coerce").clip(lower=0.0, upper=1.0)
    ).fillna(0.0)
    scored["pair_stability_score"] = (
        0.55 * pd.to_numeric(scored["member_pair_stable_fraction_median"], errors="coerce").fillna(0.0)
        + 0.45 * pd.to_numeric(scored["forced_pair_stable_fraction_median"], errors="coerce").fillna(0.0)
    ).clip(0.0, 1.0)
    scored["pair_extreme_penalty"] = (
        0.50 * pd.to_numeric(scored["member_pair_extreme_fade_fraction_median"], errors="coerce").fillna(0.0)
        + 0.50 * pd.to_numeric(scored["forced_pair_extreme_fade_fraction_median"], errors="coerce").fillna(0.0)
    ).clip(0.0, 1.0)
    scored["warning_penalty"] = (
        0.35 * scored["member_warning_flags_json"].map(lambda value: "BLEND_RISK" in _safe_json_list(value)).astype(float)
        + 0.25 * scored["member_warning_flags_json"].map(lambda value: "CROWDED" in _safe_json_list(value)).astype(float)
        + 0.20 * scored["member_warning_flags_json"].map(lambda value: "EDGE_NEAR_MASK" in _safe_json_list(value)).astype(float)
        + 0.20 * scored["member_warning_flags_json"].map(lambda value: "ASTROMETRY_RESIDUAL" in _safe_json_list(value)).astype(float)
    ).clip(0.0, 1.0)
    scored["blend_penalty"] = normalize_01(
        pd.to_numeric(scored["member_blend_risk_median"], errors="coerce").clip(lower=0.0, upper=1.0)
    ).fillna(0.0)
    scored["crowding_penalty"] = normalize_01(
        pd.to_numeric(scored["member_crowding_median"], errors="coerce").clip(lower=0.0, upper=8.0)
    ).fillna(0.0)
    scored["edge_penalty"] = 1.0 - normalize_01(
        pd.to_numeric(scored["member_edge_distance_median"], errors="coerce").clip(lower=0.0, upper=40.0)
    ).fillna(0.0)
    scored["filter_diversity_score"] = _norm_series(scored["n_support_filters"], default=0.0)
    scored["baseline_span_score"] = _norm_series(
        pd.to_numeric(scored["detection_baseline_span_days"], errors="coerce").combine_first(
            pd.to_numeric(scored["forced_baseline_span_days"], errors="coerce")
        ),
        log=True,
        default=0.0,
    )
    scored["optical_fade_support"] = pd.to_numeric(scored["optical_fade_support"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    scored["nir_fade_support"] = pd.to_numeric(scored["nir_fade_support"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    scored["mir_fade_support"] = pd.to_numeric(scored["mir_fade_support"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    scored["has_ir_coverage"] = (
        (pd.to_numeric(scored["nir_measured_pairs"], errors="coerce").fillna(0.0) > 0)
        | (pd.to_numeric(scored["mir_measured_pairs"], errors="coerce").fillna(0.0) > 0)
    )
    scored["multi_band_export_support"] = np.maximum(scored["nir_fade_support"], scored["mir_fade_support"])
    scored["dust_bias_score"] = np.where(
        scored["has_ir_coverage"],
        _clip01(scored["optical_fade_support"] - np.maximum(scored["nir_fade_support"], scored["mir_fade_support"])),
        0.0,
    )

    scored["detector_confidence"] = (
        0.28 * scored["member_signal_score"]
        + 0.20 * scored["support_pair_score"]
        + 0.12 * scored["pass_fraction"]
        + 0.12 * scored["member_ratio_loss_score"]
        + 0.14 * scored["forced_measurement_score"]
        + 0.14 * scored["member_pre_snr_score"]
    ).clip(0.0, 1.0)
    scored["registration_confidence"] = (
        0.42 * scored["member_reg_score"]
        + 0.28 * scored["forced_reg_score"]
        + 0.18 * scored["anchor_score"]
        + 0.12 * (1.0 - scored["cluster_scatter_penalty"])
    ).clip(0.0, 1.0)
    scored["subtraction_cleanliness"] = (
        0.30 * scored["pair_stability_score"]
        + 0.20 * (1.0 - scored["pair_extreme_penalty"])
        + 0.20 * scored["ratio_consistency_score"]
        + 0.15 * (1.0 - scored["warning_penalty"])
        + 0.15 * (1.0 - scored["edge_penalty"])
    ).clip(0.0, 1.0)
    scored["coverage_completeness"] = (
        0.50 * scored["measured_fraction_score"]
        + 0.20 * scored["filter_diversity_score"]
        + 0.15 * scored["baseline_span_score"]
        + 0.15 * (1.0 - pd.to_numeric(scored["pair_failure_fraction_galaxy"], errors="coerce").fillna(1.0).clip(0.0, 1.0))
    ).clip(0.0, 1.0)
    scored["forced_photometry_coherence"] = (
        0.28 * pd.to_numeric(scored["forced_positive_fraction"], errors="coerce").fillna(0.0)
        + 0.24 * scored["forced_ratio_loss_score"]
        + 0.20 * scored["ratio_consistency_score"]
        + 0.14 * scored["forced_measurement_score"]
        + 0.14 * (1.0 - pd.to_numeric(scored["forced_return_fraction"], errors="coerce").fillna(0.0))
    ).clip(0.0, 1.0)
    scored["fade_persistence"] = (
        0.30 * pd.to_numeric(scored["forced_severe_fraction"], errors="coerce").fillna(0.0)
        + 0.20 * pd.to_numeric(scored["forced_positive_fraction"], errors="coerce").fillna(0.0)
        + 0.20 * pd.to_numeric(scored["forced_late_fade_fraction"], errors="coerce").fillna(0.0)
        + 0.15 * scored["baseline_span_score"]
        + 0.15 * (1.0 - pd.to_numeric(scored["forced_late_return_fraction"], errors="coerce").fillna(0.0))
    ).clip(0.0, 1.0)
    scored["rebrightening_penalty"] = (
        0.60 * pd.to_numeric(scored["forced_return_fraction"], errors="coerce").fillna(0.0)
        + 0.40 * pd.to_numeric(scored["forced_late_return_fraction"], errors="coerce").fillna(0.0)
    ).clip(0.0, 1.0)
    scored["partial_branch_score"] = (
        0.40 * pd.to_numeric(scored["forced_partial_fraction"], errors="coerce").fillna(0.0)
        + 0.20 * scored["forced_photometry_coherence"]
        + 0.15 * scored["coverage_completeness"]
        + 0.10 * scored["registration_confidence"]
        + 0.10 * (1.0 - scored["rebrightening_penalty"])
        + 0.05 * scored["baseline_span_score"]
        - 0.10 * scored["warning_penalty"]
    ).clip(0.0, 1.0)
    scored["artifact_penalty"] = (
        0.24 * scored["warning_penalty"]
        + 0.16 * scored["blend_penalty"]
        + 0.10 * scored["crowding_penalty"]
        + 0.15 * (1.0 - scored["registration_confidence"])
        + 0.10 * scored["edge_penalty"]
        + 0.10 * pd.to_numeric(scored["pair_failure_fraction_galaxy"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
        + 0.15 * scored["cluster_scatter_penalty"]
    ).clip(0.0, 1.0)
    scored["dust_survivor_score"] = (
        0.45 * scored["dust_bias_score"]
        + 0.20 * np.where(scored["has_ir_coverage"], 1.0 - scored["multi_band_export_support"], 0.0)
        + 0.20 * scored["rebrightening_penalty"]
        + 0.15 * np.where(scored["has_ir_coverage"], 0.0, 0.35 * scored["optical_fade_support"])
    ).clip(0.0, 1.0)
    scored["systematic_risk"] = (
        0.45 * scored["artifact_penalty"]
        + 0.20 * (1.0 - scored["registration_confidence"])
        + 0.20 * (1.0 - scored["subtraction_cleanliness"])
        + 0.15 * (1.0 - scored["detector_confidence"])
    ).clip(0.0, 1.0)
    scored["unresolved_variability_score"] = (
        0.35 * scored["rebrightening_penalty"]
        + 0.25 * (1.0 - scored["coverage_completeness"])
        + 0.20 * (1.0 - scored["forced_photometry_coherence"])
        + 0.20 * (1.0 - scored["fade_persistence"])
    ).clip(0.0, 1.0)
    scored["export_failure_score"] = (
        0.22 * scored["detector_confidence"]
        + 0.16 * scored["registration_confidence"]
        + 0.10 * scored["subtraction_cleanliness"]
        + 0.17 * scored["forced_photometry_coherence"]
        + 0.17 * scored["fade_persistence"]
        + 0.08 * scored["coverage_completeness"]
        + 0.05 * scored["multi_band_export_support"]
        + 0.05 * scored["partial_branch_score"]
        + 0.05 * scored["baseline_span_score"]
        - 0.10 * scored["rebrightening_penalty"]
        - 0.13 * scored["artifact_penalty"]
        - 0.12 * scored["dust_survivor_score"]
    ).clip(0.0, 1.0)

    scored["support_depth_score"] = (
        0.55 * scored["support_pair_score"]
        + 0.45 * scored["pass_fraction"]
    ).clip(0.0, 1.0)
    scored["single_pair_penalty"] = (pd.to_numeric(scored["n_support_pairs"], errors="coerce").fillna(0.0) <= 1.0).astype(float)
    scored["branch_rank_raw"] = (
        scored["export_failure_score"]
        + 0.35 * scored["partial_branch_score"]
        + 0.20 * scored["coverage_completeness"]
        + 0.15 * scored["support_depth_score"]
        - 0.45 * scored["systematic_risk"]
        - 0.35 * scored["dust_survivor_score"]
        - 0.20 * scored["single_pair_penalty"]
    )
    scored["branch_rank_score"] = normalize_01(scored["branch_rank_raw"]).fillna(0.0)

    benchmark_ok = float(benchmark_summary.get("recovery_fraction", 0.0) or 0.0) >= 0.70
    scored["benchmark_recovery_fraction"] = float(benchmark_summary.get("recovery_fraction", 0.0) or 0.0)
    scored["benchmark_median_recovered_sep_arcsec"] = float(
        benchmark_summary.get("median_recovered_sep_arcsec", float("nan")) or float("nan")
    )
    scored["benchmark_guardrail_pass"] = benchmark_ok

    classifications: list[str] = []
    confidences: list[str] = []
    reasons: list[str] = []
    for _, row in scored.iterrows():
        label, confidence = _classify_cluster(row)
        classifications.append(label)
        confidences.append(confidence)
        reasons.append(_json_list(_reason_codes(row, benchmark_ok=benchmark_ok)))
    scored["branch_class"] = classifications
    scored["branch_confidence"] = confidences
    scored["reason_codes_json"] = reasons
    scored = scored.sort_values(
        ["branch_rank_score", "export_failure_score", "n_support_pairs", "member_diff_sigma_minmax"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    scored["branch_rank"] = np.arange(1, len(scored) + 1, dtype=int)
    return scored


def _explanation_markdown(row: pd.Series) -> str:
    return "\n".join(
        [
            f"# {row['cluster_id']}",
            "",
            f"- Branch class: `{row['branch_class']}`",
            f"- Confidence: `{row['branch_confidence']}`",
            f"- Branch rank score: `{row['branch_rank_score']:.3f}`",
            f"- Export-failure score: `{row['export_failure_score']:.3f}`",
            f"- Intermediate-branch score: `{row['partial_branch_score']:.3f}`",
            f"- Dust-survivor score: `{row['dust_survivor_score']:.3f}`",
            f"- Systematic risk: `{row['systematic_risk']:.3f}`",
            f"- Unresolved/variability score: `{row['unresolved_variability_score']:.3f}`",
            "",
            "## Score Decomposition",
            f"- Detector confidence: `{row['detector_confidence']:.3f}`",
            f"- Registration confidence: `{row['registration_confidence']:.3f}`",
            f"- Subtraction cleanliness: `{row['subtraction_cleanliness']:.3f}`",
            f"- Coverage completeness: `{row['coverage_completeness']:.3f}`",
            f"- Forced photometry coherence: `{row['forced_photometry_coherence']:.3f}`",
            f"- Fade persistence: `{row['fade_persistence']:.3f}`",
            f"- Rebrightening penalty: `{row['rebrightening_penalty']:.3f}`",
            f"- Artifact penalty: `{row['artifact_penalty']:.3f}`",
            "",
            "## Raw Support",
            f"- Members: `{int(row['n_members'])}` total, `{int(row['n_pass_members'])}` PASS, `{int(row['n_review_members'])}` REVIEW",
            f"- Support pairs: `{int(row['n_support_pairs'])}`",
            f"- Support filters: `{int(row['n_support_filters'])}`",
            f"- Detection baseline span days: `{float(row['detection_baseline_span_days']):.1f}`",
            f"- Forced measured pairs: `{int(row['forced_n_pairs_measured'])}` / `{int(row['forced_n_pairs_considered'])}`",
            f"- Forced severe fraction: `{float(row['forced_severe_fraction']):.3f}`",
            f"- Forced partial fraction: `{float(row['forced_partial_fraction']):.3f}`",
            f"- Forced return fraction: `{float(row['forced_return_fraction']):.3f}`",
            f"- Optical fade support: `{float(row['optical_fade_support']):.3f}`",
            f"- NIR fade support: `{float(row['nir_fade_support']):.3f}`",
            f"- MIR fade support: `{float(row['mir_fade_support']):.3f}`",
            "",
            "## Reason Codes",
            *[f"- `{code}`" for code in _safe_json_list(row["reason_codes_json"])],
            "",
        ]
    )


def run_branch_scoring(
    *,
    root_dir: Path,
    strict_output_dir: Path,
    output_dir: Path | None = None,
    benchmark_dir: Path | None = None,
    cluster_radius_arcsec: float = 0.25,
    precomputed_dir: Path | None = None,
) -> dict[str, Path]:
    out_dir = ensure_dir(output_dir or (strict_output_dir / "branch_aware"))
    packets_dir = ensure_dir(out_dir / "packets")
    predictaroni_dir = ensure_dir(root_dir / "predictaroni")
    progress_path = out_dir / "progress.json"
    _write_run_metadata(
        out_dir,
        root_dir,
        {
            "created_utc": utc_stamp(),
            "mode": "branch_scoring",
            "strict_output_dir": str(strict_output_dir),
            "benchmark_dir": str(benchmark_dir) if benchmark_dir is not None else None,
            "cluster_radius_arcsec": float(cluster_radius_arcsec),
            "precomputed_dir": str(precomputed_dir) if precomputed_dir is not None else None,
        },
    )
    _log(out_dir, f"run start mode=branch_scoring strict_output_dir={strict_output_dir}")

    detections = pd.read_parquet(strict_output_dir / "fade_candidates.parquet").copy()
    detections["reason_codes"] = detections["reason_codes_json"].map(_safe_json_list)
    detections["warning_flags"] = detections["warning_flags_json"].map(_safe_json_list)

    if precomputed_dir is None:
        write_progress(progress_path, "stage00_cluster_detections", "running", {"n_detections": int(len(detections))})
        clusters, members = _build_cluster_tables(detections, radius_arcsec=cluster_radius_arcsec)
        clusters.to_parquet(out_dir / "cluster_index.parquet", index=False)
        clusters.to_csv(out_dir / "cluster_index.csv", index=False)
        members.to_parquet(out_dir / "cluster_members.parquet", index=False)
        members.to_csv(out_dir / "cluster_members.csv", index=False)
        write_progress(
            progress_path,
            "stage00_cluster_detections",
            "completed",
            {
                "n_detections": int(len(detections)),
                "n_clusters": int(len(clusters)),
            },
        )
        _log(out_dir, f"clustered detections n_detections={len(detections)} n_clusters={len(clusters)}")
        pair_catalog = _load_pair_catalog(strict_output_dir)
        forced = _force_measure_cluster_centers(
            root_dir=root_dir,
            strict_output_dir=strict_output_dir,
            output_dir=out_dir,
            clusters=clusters,
            pair_catalog=pair_catalog,
        )
    else:
        write_progress(progress_path, "stage00_reuse_precomputed", "running", {"source_dir": str(precomputed_dir)})
        clusters = pd.read_parquet(precomputed_dir / "cluster_index.parquet")
        members = pd.read_parquet(precomputed_dir / "cluster_members.parquet")
        forced = pd.read_parquet(precomputed_dir / "forced_measurements.parquet")
        clusters.to_parquet(out_dir / "cluster_index.parquet", index=False)
        clusters.to_csv(out_dir / "cluster_index.csv", index=False)
        members.to_parquet(out_dir / "cluster_members.parquet", index=False)
        members.to_csv(out_dir / "cluster_members.csv", index=False)
        forced.to_parquet(out_dir / "forced_measurements.parquet", index=False)
        forced.to_csv(out_dir / "forced_measurements.csv", index=False)
        pair_catalog = _load_pair_catalog(strict_output_dir)
        write_progress(
            progress_path,
            "stage00_reuse_precomputed",
            "completed",
            {
                "source_dir": str(precomputed_dir),
                "n_clusters": int(len(clusters)),
                "n_forced_rows_written": int(len(forced)),
            },
        )
        _log(out_dir, f"reused precomputed branch inputs source_dir={precomputed_dir} n_clusters={len(clusters)} forced_rows={len(forced)}")

    write_progress(progress_path, "stage02_score_clusters", "running", {"n_clusters": int(len(clusters))})
    benchmark_summary = {
        "recovery_fraction": 0.0,
        "median_recovered_sep_arcsec": float("nan"),
    }
    if benchmark_dir is not None and (benchmark_dir / "benchmark_summary.json").exists():
        benchmark_summary = json.loads((benchmark_dir / "benchmark_summary.json").read_text())
    scored = _build_cluster_scores(
        clusters=clusters,
        members=members,
        forced=forced,
        pair_catalog=pair_catalog,
        benchmark_summary=benchmark_summary,
    )
    scored.to_parquet(out_dir / "branch_cluster_scores.parquet", index=False)
    scored.to_csv(out_dir / "branch_cluster_scores.csv", index=False)

    detection_cluster = members[["candidate_id", "cluster_id"]].drop_duplicates("candidate_id")
    detection_scored = detections.merge(detection_cluster, on="candidate_id", how="left").merge(
        scored[
            [
                "cluster_id",
                "branch_rank",
                "branch_rank_score",
                "branch_class",
                "branch_confidence",
                "export_failure_score",
                "partial_branch_score",
                "dust_survivor_score",
                "systematic_risk",
                "unresolved_variability_score",
            ]
        ],
        on="cluster_id",
        how="left",
    )
    detection_scored.to_parquet(out_dir / "branch_detection_scores.parquet", index=False)
    detection_scored.to_csv(out_dir / "branch_detection_scores.csv", index=False)

    for _, row in scored.iterrows():
        packet_dir = ensure_dir(packets_dir / str(row["cluster_id"]))
        packet_summary = {
            "cluster_id": str(row["cluster_id"]),
            "galaxy_name": str(row["galaxy_name"]),
            "ra_deg": float(row["ra_deg"]),
            "dec_deg": float(row["dec_deg"]),
            "branch_rank": int(row["branch_rank"]),
            "branch_rank_score": float(row["branch_rank_score"]),
            "branch_class": str(row["branch_class"]),
            "branch_confidence": str(row["branch_confidence"]),
            "export_failure_score": float(row["export_failure_score"]),
            "partial_branch_score": float(row["partial_branch_score"]),
            "dust_survivor_score": float(row["dust_survivor_score"]),
            "systematic_risk": float(row["systematic_risk"]),
            "unresolved_variability_score": float(row["unresolved_variability_score"]),
            "reason_codes": _safe_json_list(row["reason_codes_json"]),
        }
        write_json(packet_dir / "summary.json", packet_summary)
        write_text(packet_dir / "decision.md", _explanation_markdown(row))

    class_counts = scored["branch_class"].value_counts(dropna=False).to_dict()
    summary = {
        "created_utc": utc_stamp(),
        "n_clusters": int(len(scored)),
        "n_detections": int(len(detections)),
        "class_counts": {str(key): int(value) for key, value in class_counts.items()},
        "benchmark_recovery_fraction": float(benchmark_summary.get("recovery_fraction", 0.0) or 0.0),
        "benchmark_median_recovered_sep_arcsec": float(
            benchmark_summary.get("median_recovered_sep_arcsec", float("nan")) or float("nan")
        ),
        "benchmark_guardrail_pass": bool(float(benchmark_summary.get("recovery_fraction", 0.0) or 0.0) >= 0.70),
        "top_strong_export": scored.loc[
            scored["branch_class"].astype(str).eq(CLASS_STRONG_EXPORT),
            ["cluster_id", "galaxy_name", "branch_rank_score", "export_failure_score"],
        ]
        .head(10)
        .to_dict(orient="records"),
        "top_intermediate": scored.loc[
            scored["branch_class"].astype(str).eq(CLASS_INTERMEDIATE),
            ["cluster_id", "galaxy_name", "branch_rank_score", "partial_branch_score"],
        ]
        .head(10)
        .to_dict(orient="records"),
        "top_dust": scored.loc[
            scored["branch_class"].astype(str).eq(CLASS_DUST),
            ["cluster_id", "galaxy_name", "branch_rank_score", "dust_survivor_score"],
        ]
        .head(10)
        .to_dict(orient="records"),
        "top_systematic": scored.loc[
            scored["branch_class"].astype(str).eq(CLASS_SYSTEMATIC),
            ["cluster_id", "galaxy_name", "branch_rank_score", "systematic_risk"],
        ]
        .head(10)
        .to_dict(orient="records"),
    }
    write_json(out_dir / "branch_summary.json", summary)

    report_lines = [
        "# Branch-Aware Export-Failure Report",
        "",
        "## Summary",
        f"- Strict input run: `{strict_output_dir}`",
        f"- Clusters scored: `{len(scored)}`",
        f"- Input detections: `{len(detections)}`",
        f"- Benchmark recovery fraction carried into this run: `{summary['benchmark_recovery_fraction']:.2f}`",
        f"- Benchmark guardrail pass: `{summary['benchmark_guardrail_pass']}`",
        "",
        "## Class Counts",
    ]
    for label in [CLASS_STRONG_EXPORT, CLASS_INTERMEDIATE, CLASS_DUST, CLASS_SYSTEMATIC, CLASS_UNRESOLVED]:
        report_lines.append(f"- `{label}`: `{summary['class_counts'].get(label, 0)}`")
    report_lines.extend(["", "## Top Branch-Ranked Clusters"])
    for _, row in scored.head(15).iterrows():
        report_lines.append(
            f"- `{row['cluster_id']}` in `{row['galaxy_name']}`: class `{row['branch_class']}`, "
            f"rank `{row['branch_rank_score']:.3f}`, export `{row['export_failure_score']:.3f}`, "
            f"partial `{row['partial_branch_score']:.3f}`, dust `{row['dust_survivor_score']:.3f}`, systematic `{row['systematic_risk']:.3f}`"
        )
    write_text(out_dir / "branch_report.md", "\n".join(report_lines) + "\n")

    manifest = {
        "branch_cluster_scores": str(out_dir / "branch_cluster_scores.parquet"),
        "branch_detection_scores": str(out_dir / "branch_detection_scores.parquet"),
        "cluster_index": str(out_dir / "cluster_index.parquet"),
        "cluster_members": str(out_dir / "cluster_members.parquet"),
        "forced_measurements": str(out_dir / "forced_measurements.parquet"),
        "branch_summary": str(out_dir / "branch_summary.json"),
        "branch_report": str(out_dir / "branch_report.md"),
        "packets_dir": str(packets_dir),
        "predictaroni_dir": str(predictaroni_dir),
    }
    write_json(out_dir / "branch_manifest.json", manifest)
    write_progress(
        progress_path,
        "stage02_score_clusters",
        "completed",
        {
            "n_clusters": int(len(scored)),
            "n_strong_export": int((scored["branch_class"] == CLASS_STRONG_EXPORT).sum()),
            "n_intermediate": int((scored["branch_class"] == CLASS_INTERMEDIATE).sum()),
            "n_dust": int((scored["branch_class"] == CLASS_DUST).sum()),
            "n_systematic": int((scored["branch_class"] == CLASS_SYSTEMATIC).sum()),
            "n_unresolved": int((scored["branch_class"] == CLASS_UNRESOLVED).sum()),
        },
    )
    _log(out_dir, f"run completed n_clusters={len(scored)}")
    return {
        "branch_cluster_scores": out_dir / "branch_cluster_scores.parquet",
        "branch_cluster_scores_csv": out_dir / "branch_cluster_scores.csv",
        "branch_detection_scores": out_dir / "branch_detection_scores.parquet",
        "cluster_index": out_dir / "cluster_index.parquet",
        "cluster_members": out_dir / "cluster_members.parquet",
        "forced_measurements": out_dir / "forced_measurements.parquet",
        "branch_summary": out_dir / "branch_summary.json",
        "branch_report": out_dir / "branch_report.md",
        "branch_manifest": out_dir / "branch_manifest.json",
        "packets_dir": packets_dir,
    }
