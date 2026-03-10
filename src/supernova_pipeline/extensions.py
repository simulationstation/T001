from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa
from reproject import reproject_interp
from scipy.ndimage import gaussian_filter

from .archive_matrix import infer_filter_central_um
from .candidate_followup import _candidate_cutout
from .pixel_search import (
    _aperture_fluxes,
    _background_stats,
    _detect_sources,
    _fwhm_guess,
    _load_science_image,
)
from .utils import ensure_dir, normalize_01, robust_log10, write_json


Irsa.ROW_LIMIT = 5


def _infer_obs_collection_from_path(path: Path) -> str:
    parts = {part.upper() for part in path.parts}
    if "JWST" in parts:
        return "JWST"
    return "HST"


def _infer_instrument_name(obs_id: str) -> str:
    text = str(obs_id).lower()
    if "wfc3_ir" in text:
        return "WFC3/IR"
    if "wfc3_uvis" in text:
        return "WFC3/UVIS"
    if "acs" in text:
        return "ACS/WFC"
    if "wfpc2" in text:
        return "WFPC2"
    if "nircam" in text:
        return "NIRCam"
    if "miri" in text:
        return "MIRI"
    return "UNKNOWN"


def _decode_provenance(candidate: pd.Series) -> dict[str, Any]:
    value = candidate.get("provenance_json")
    if isinstance(value, str) and value.strip():
        return json.loads(value)
    return {}


def _candidate_world(candidate: pd.Series) -> SkyCoord:
    return SkyCoord(float(candidate["ra_deg"]) * u.deg, float(candidate["dec_deg"]) * u.deg, frame="icrs")


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _ensure_columns(df: pd.DataFrame | None, columns: list[str]) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame(columns=columns)
    out = df.copy()
    for column in columns:
        if column not in out.columns:
            out[column] = np.nan
    return out


def _weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not np.any(good):
        return float("nan"), float("nan")
    weights = 1.0 / np.square(errors[good])
    mean = float(np.sum(values[good] * weights) / np.sum(weights))
    err = float(np.sqrt(1.0 / np.sum(weights)))
    return mean, err


def _robust_corrcoef(x: np.ndarray, y: np.ndarray) -> float:
    good = np.isfinite(x) & np.isfinite(y)
    if np.sum(good) < 3:
        return float("nan")
    x_rank = pd.Series(x[good]).rank(method="average").to_numpy(dtype=float)
    y_rank = pd.Series(y[good]).rank(method="average").to_numpy(dtype=float)
    corr = np.corrcoef(x_rank, y_rank)[0, 1]
    return float(corr) if math.isfinite(corr) else float("nan")


def _extract_scene_metrics(candidate: pd.Series, *, cutout_arcsec: float = 90.0) -> tuple[dict[str, Any], pd.DataFrame]:
    provenance = _decode_provenance(candidate)
    product_info = provenance.get("product_info", {})
    path_1 = Path(str(product_info.get("product_path_1", "")))
    path_2 = Path(str(product_info.get("product_path_2", "")))
    if not path_1.exists() or not path_2.exists():
        raise RuntimeError("Pair products are not present in the local cache.")

    obs_id_1 = str(candidate.get("pre_obs_id", product_info.get("obs_id_1", "pre_obs")))
    obs_id_2 = str(candidate.get("post_obs_id", product_info.get("obs_id_2", "post_obs")))
    image_1 = _load_science_image(
        path_1,
        obs_collection=_infer_obs_collection_from_path(path_1),
        obs_id=obs_id_1,
        filter_name=str(candidate.get("pre_filter", "")),
        instrument_name=_infer_instrument_name(obs_id_1),
    )
    image_2 = _load_science_image(
        path_2,
        obs_collection=_infer_obs_collection_from_path(path_2),
        obs_id=obs_id_2,
        filter_name=str(candidate.get("post_filter", "")),
        instrument_name=_infer_instrument_name(obs_id_2),
    )

    world = _candidate_world(candidate)
    image_1 = _candidate_cutout(image_1, world, size_arcsec=cutout_arcsec, min_pixels=120)
    image_2 = _candidate_cutout(image_2, world, size_arcsec=cutout_arcsec, min_pixels=120)
    post_reproj, footprint = reproject_interp((image_2.data, image_2.wcs), image_1.wcs, shape_out=image_1.data.shape)
    valid_mask = np.isfinite(image_1.data) & np.isfinite(post_reproj) & (footprint > 0)
    if valid_mask.sum() < 1_000:
        raise RuntimeError("Insufficient overlap in candidate scene cutout.")

    pre_bg, _ = _background_stats(image_1.data, valid_mask)
    post_bg, _ = _background_stats(post_reproj, valid_mask)
    pre_sub = image_1.data - pre_bg
    post_sub = post_reproj - post_bg
    x_cand_arr, y_cand_arr = image_1.wcs.world_to_pixel(world)
    x_cand = float(np.asarray(x_cand_arr).reshape(-1)[0])
    y_cand = float(np.asarray(y_cand_arr).reshape(-1)[0])
    if not (math.isfinite(x_cand) and math.isfinite(y_cand)):
        raise RuntimeError("Candidate projection failed in the scene model.")

    sources = _detect_sources(pre_sub, valid_mask, pixel_scale_arcsec=image_1.pixel_scale_arcsec, max_sources=400)
    if sources.empty:
        raise RuntimeError("No comparison stars detected in the scene model cutout.")

    radius = max(2.5, 1.2 * _fwhm_guess(image_1.pixel_scale_arcsec))
    positions = sources[["x", "y"]].to_numpy(dtype=float)
    flux_pre, err_pre, _ = _aperture_fluxes(pre_sub, valid_mask, positions, radius=radius)
    flux_post, err_post, _ = _aperture_fluxes(post_sub, valid_mask, positions, radius=radius)

    cand_flux_pre, cand_err_pre, _ = _aperture_fluxes(pre_sub, valid_mask, np.array([[x_cand, y_cand]], dtype=float), radius=radius)
    cand_flux_post, cand_err_post, _ = _aperture_fluxes(post_sub, valid_mask, np.array([[x_cand, y_cand]], dtype=float), radius=radius)

    source_table = sources.copy()
    source_table["flux_pre"] = flux_pre
    source_table["flux_post"] = flux_post
    source_table["ratio"] = np.divide(
        flux_post,
        flux_pre,
        out=np.full_like(flux_post, np.nan, dtype=float),
        where=np.isfinite(flux_pre) & (np.abs(flux_pre) > 0),
    )
    source_table["fade_fraction"] = 1.0 - source_table["ratio"]
    source_table["distance_pix"] = np.hypot(source_table["x"] - float(x_cand), source_table["y"] - float(y_cand))
    candidate_flux_pre = float(cand_flux_pre[0])
    brightness_ratio = np.divide(
        np.abs(source_table["flux_pre"].to_numpy(dtype=float)),
        max(abs(candidate_flux_pre), 1e-9),
    )
    comp_mask = (
        (source_table["distance_pix"] <= 40.0)
        & (source_table["distance_pix"] >= 4.0)
        & np.isfinite(source_table["fade_fraction"])
        & np.isfinite(source_table["flux_pre"])
        & (brightness_ratio >= 0.25)
        & (brightness_ratio <= 4.0)
    )
    comparators = source_table.loc[comp_mask].copy()
    neighbor_median_fade = float(comparators["fade_fraction"].median()) if not comparators.empty else float("nan")
    neighbor_scatter = (
        float(np.nanmedian(np.abs(comparators["fade_fraction"].to_numpy(dtype=float) - neighbor_median_fade)))
        if not comparators.empty and math.isfinite(neighbor_median_fade)
        else float("nan")
    )
    if not math.isfinite(neighbor_scatter) or neighbor_scatter <= 0:
        neighbor_scatter = float(comparators["fade_fraction"].std(ddof=0)) if len(comparators) > 1 else float("nan")

    candidate_ratio = float(cand_flux_post[0] / cand_flux_pre[0]) if math.isfinite(cand_flux_pre[0]) and abs(cand_flux_pre[0]) > 0 else float("nan")
    candidate_fade = 1.0 - candidate_ratio if math.isfinite(candidate_ratio) else float("nan")
    uniqueness_sigma = float(abs(candidate_fade - neighbor_median_fade) / max(neighbor_scatter, 1e-3)) if math.isfinite(candidate_fade) and math.isfinite(neighbor_median_fade) else float("nan")
    candidate_sigma = float((cand_flux_pre[0] - cand_flux_post[0]) / math.sqrt(max(cand_err_pre[0] ** 2 + cand_err_post[0] ** 2, 1e-9)))
    residual_std = float(np.nanstd((pre_sub - post_sub)[valid_mask]))
    candidate_residual = float(pre_sub[int(round(y_cand)), int(round(x_cand))] - post_sub[int(round(y_cand)), int(round(x_cand))])

    followup_path = Path("sed") / str(candidate["candidate_id"]) / "all_photometry.parquet"
    scene_mode = "pair_proxy"
    n_epochs_measured = 2
    monotonic_score = float("nan")
    pre_stage_flux = candidate_flux_pre
    post_stage_flux = cand_flux_post[0]
    if followup_path.exists():
        phot = pd.read_parquet(followup_path)
        measured = phot[phot["status"].astype(str) == "MEASURED"].copy()
        measured = measured[
            np.isfinite(measured["flux_jy"])
            & np.isfinite(measured["err_jy"])
            & (measured["err_jy"] > 0)
            & np.isfinite(measured["t_min"])
        ].copy()
        if not measured.empty:
            scene_mode = "all_epoch"
            n_epochs_measured = int(len(measured))
            target_filter = str(candidate.get("pre_filter", ""))
            if target_filter:
                filter_measured = measured[measured["filter_name"].astype(str) == target_filter].copy()
                if not filter_measured.empty:
                    measured = filter_measured
            if "stage" in measured.columns:
                pre_mean, _ = _weighted_mean(
                    measured.loc[measured["stage"] == "pre", "flux_jy"].to_numpy(dtype=float),
                    measured.loc[measured["stage"] == "pre", "err_jy"].to_numpy(dtype=float),
                )
                post_mean, _ = _weighted_mean(
                    measured.loc[measured["stage"] == "post", "flux_jy"].to_numpy(dtype=float),
                    measured.loc[measured["stage"] == "post", "err_jy"].to_numpy(dtype=float),
                )
                if math.isfinite(pre_mean):
                    pre_stage_flux = pre_mean
                if math.isfinite(post_mean):
                    post_stage_flux = post_mean
            monotonic_score = _robust_corrcoef(
                measured["t_min"].to_numpy(dtype=float),
                -measured["flux_jy"].to_numpy(dtype=float),
            )

    scene_metrics = {
        "candidate_id": str(candidate["candidate_id"]),
        "galaxy_name": str(candidate["galaxy_name"]),
        "scene_mode": scene_mode,
        "n_scene_epochs": int(n_epochs_measured),
        "scene_pre_flux": float(pre_stage_flux) if math.isfinite(pre_stage_flux) else float("nan"),
        "scene_post_flux": float(post_stage_flux) if math.isfinite(post_stage_flux) else float("nan"),
        "scene_pair_pre_flux_native": candidate_flux_pre,
        "scene_pair_post_flux_native": float(cand_flux_post[0]),
        "scene_fade_fraction": float(1.0 - post_stage_flux / pre_stage_flux) if math.isfinite(pre_stage_flux) and abs(pre_stage_flux) > 0 and math.isfinite(post_stage_flux) else candidate_fade,
        "scene_fade_sigma": candidate_sigma,
        "scene_neighbor_median_fade": neighbor_median_fade,
        "scene_neighbor_scatter": neighbor_scatter,
        "scene_unique_fade_sigma": uniqueness_sigma,
        "scene_local_comparators": int(len(comparators)),
        "scene_monotonic_score": monotonic_score,
        "scene_residual_center": candidate_residual,
        "scene_residual_std": residual_std,
    }
    if not comparators.empty:
        comparators = comparators.assign(candidate_id=str(candidate["candidate_id"]))
    return scene_metrics, comparators


def _run_scene_model_stage(root_dir: Path, ledger: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Path]]:
    out_dir = ensure_dir(root_dir / "extensions" / "scene_model")
    rows: list[dict[str, Any]] = []
    neighbor_frames: list[pd.DataFrame] = []
    failures: list[dict[str, str]] = []
    for _, candidate in ledger.iterrows():
        try:
            metrics, comparators = _extract_scene_metrics(candidate)
            rows.append(metrics)
            if not comparators.empty:
                neighbor_frames.append(comparators)
        except Exception as exc:
            failures.append({"candidate_id": str(candidate["candidate_id"]), "error": str(exc)})
    scene_columns = [
        "candidate_id",
        "galaxy_name",
        "scene_mode",
        "n_scene_epochs",
        "scene_pre_flux",
        "scene_post_flux",
        "scene_pair_pre_flux_native",
        "scene_pair_post_flux_native",
        "scene_fade_fraction",
        "scene_fade_sigma",
        "scene_neighbor_median_fade",
        "scene_neighbor_scatter",
        "scene_unique_fade_sigma",
        "scene_local_comparators",
        "scene_monotonic_score",
        "scene_residual_center",
        "scene_residual_std",
    ]
    scene_df = pd.DataFrame(rows, columns=scene_columns)
    if not scene_df.empty:
        scene_df = scene_df.sort_values(["scene_unique_fade_sigma", "scene_fade_sigma"], ascending=[False, False]).reset_index(drop=True)
    neighbors_df = pd.concat(neighbor_frames, ignore_index=True) if neighbor_frames else pd.DataFrame()
    scene_path = out_dir / "scene_model_metrics.csv"
    neighbors_path = out_dir / "scene_model_neighbors.parquet"
    failures_path = out_dir / "scene_model_failures.csv"
    scene_df.to_csv(scene_path, index=False)
    if not neighbors_df.empty:
        neighbors_df.to_parquet(neighbors_path, index=False)
    pd.DataFrame(failures).to_csv(failures_path, index=False)
    write_json(
        out_dir / "summary.json",
        {
            "n_candidates_attempted": int(len(ledger)),
            "n_candidates_completed": int(len(scene_df)),
            "n_failures": int(len(failures)),
        },
    )
    return scene_df, neighbors_df, {
        "scene_metrics": scene_path,
        "scene_neighbors": neighbors_path,
        "scene_failures": failures_path,
        "summary": out_dir / "summary.json",
    }


def _deep_reference_images_for_galaxy(ledger: pd.DataFrame, galaxy_name: str) -> list[dict[str, Any]]:
    work = ledger[ledger["galaxy_name"].astype(str) == str(galaxy_name)].copy()
    if work.empty:
        return []
    work["filter_key"] = np.where(
        work["pre_filter"].astype(str) == work["post_filter"].astype(str),
        work["pre_filter"].astype(str),
        work["pre_filter"].astype(str) + "__" + work["post_filter"].astype(str),
    )
    filter_key = work.groupby("filter_key").size().sort_values(ascending=False).index[0]
    use = work[work["filter_key"].astype(str) == str(filter_key)].copy()
    images: dict[str, dict[str, Any]] = {}
    for _, row in use.iterrows():
        provenance = _decode_provenance(row)
        product_info = provenance.get("product_info", {})
        for suffix, obs_id_col, mjd_col, filter_col in [
            ("1", "pre_obs_id", "pre_mjd", "pre_filter"),
            ("2", "post_obs_id", "post_mjd", "post_filter"),
        ]:
            path = Path(str(product_info.get(f"product_path_{suffix}", "")))
            if not path.exists():
                continue
            obs_id = str(row[obs_id_col])
            images.setdefault(
                obs_id,
                {
                    "path": path,
                    "obs_id": obs_id,
                    "filter_name": str(row[filter_col]),
                    "epoch_mjd": float(row[mjd_col]),
                },
            )
    return sorted(images.values(), key=lambda item: item["epoch_mjd"])


def _run_deep_reference_stage(root_dir: Path, ledger: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "deep_reference")
    host_priority = ["NGC 5861", "MESSIER 101", "MESSIER 077"]
    hosts = [host for host in host_priority if host in set(ledger["galaxy_name"].astype(str))]
    summary_rows: list[dict[str, Any]] = []
    for galaxy_name in hosts:
        galaxy_dir = ensure_dir(out_dir / galaxy_name.replace(" ", "_"))
        image_rows = _deep_reference_images_for_galaxy(ledger, galaxy_name)
        if len(image_rows) < 2:
            continue
        galaxy_candidates = ledger[ledger["galaxy_name"].astype(str) == galaxy_name].copy()
        center = SkyCoord(
            float(galaxy_candidates["ra_deg"].median()) * u.deg,
            float(galaxy_candidates["dec_deg"].median()) * u.deg,
            frame="icrs",
        )
        loaded = []
        for item in image_rows:
            image = _load_science_image(
                item["path"],
                obs_collection=_infer_obs_collection_from_path(item["path"]),
                obs_id=item["obs_id"],
                filter_name=item["filter_name"],
                instrument_name=_infer_instrument_name(item["obs_id"]),
            )
            image = _candidate_cutout(image, center, size_arcsec=120.0, min_pixels=160)
            loaded.append((item, image))
        ref_item, ref_image = loaded[0]
        aligned_stack = [np.asarray(ref_image.data, dtype=float)]
        image_flux_tables: list[pd.DataFrame] = []
        epoch_rows: list[dict[str, Any]] = []
        for item, image in loaded[1:]:
            reproj, footprint = reproject_interp((image.data, image.wcs), ref_image.wcs, shape_out=ref_image.data.shape)
            reproj[footprint <= 0] = np.nan
            aligned_stack.append(np.asarray(reproj, dtype=float))
        deep_ref = np.nanmedian(np.stack(aligned_stack), axis=0)
        valid_ref = np.isfinite(deep_ref)
        sources = _detect_sources(deep_ref, valid_ref, pixel_scale_arcsec=ref_image.pixel_scale_arcsec, max_sources=300)
        if sources.empty:
            continue
        positions = sources[["x", "y"]].to_numpy(dtype=float)
        radius = max(2.5, 1.2 * _fwhm_guess(ref_image.pixel_scale_arcsec))
        for item, image in loaded:
            if image is ref_image:
                aligned = np.asarray(image.data, dtype=float)
            else:
                aligned, footprint = reproject_interp((image.data, image.wcs), ref_image.wcs, shape_out=ref_image.data.shape)
                aligned[footprint <= 0] = np.nan
            valid = np.isfinite(aligned)
            bg, _ = _background_stats(aligned, valid)
            flux, err, _ = _aperture_fluxes(aligned - bg, valid, positions, radius=radius)
            frame = pd.DataFrame(
                {
                    "source_id": np.arange(len(sources)),
                    "x": sources["x"],
                    "y": sources["y"],
                    "obs_id": item["obs_id"],
                    "epoch_mjd": item["epoch_mjd"],
                    "filter_name": item["filter_name"],
                    "flux_native": flux,
                    "err_native": err,
                }
            )
            image_flux_tables.append(frame)
            epoch_rows.append(
                {
                    "obs_id": item["obs_id"],
                    "epoch_mjd": item["epoch_mjd"],
                    "filter_name": item["filter_name"],
                }
            )
        forced = pd.concat(image_flux_tables, ignore_index=True)
        first_epoch = forced["epoch_mjd"].min()
        last_epoch = forced["epoch_mjd"].max()
        first_flux = forced.loc[forced["epoch_mjd"] == first_epoch, ["source_id", "flux_native", "err_native"]].rename(columns={"flux_native": "flux_first", "err_native": "err_first"})
        last_flux = forced.loc[forced["epoch_mjd"] == last_epoch, ["source_id", "flux_native", "err_native"]].rename(columns={"flux_native": "flux_last", "err_native": "err_last"})
        ranking = sources.copy()
        ranking["source_id"] = np.arange(len(sources))
        ranking = ranking.merge(first_flux, on="source_id", how="left").merge(last_flux, on="source_id", how="left")
        ranking["fade_fraction"] = 1.0 - np.divide(
            ranking["flux_last"],
            ranking["flux_first"],
            out=np.full(len(ranking), np.nan, dtype=float),
            where=np.isfinite(ranking["flux_first"]) & (np.abs(ranking["flux_first"]) > 0),
        )
        ranking["fade_sigma"] = np.divide(
            ranking["flux_first"] - ranking["flux_last"],
            np.sqrt(np.square(ranking["err_first"]) + np.square(ranking["err_last"])),
            out=np.full(len(ranking), np.nan, dtype=float),
            where=np.isfinite(ranking["err_first"]) & np.isfinite(ranking["err_last"]),
        )
        matched_ids: list[str | None] = []
        for _, source in ranking.iterrows():
            world = ref_image.wcs.pixel_to_world(float(source["x"]), float(source["y"]))
            sep = SkyCoord(galaxy_candidates["ra_deg"].to_numpy(dtype=float) * u.deg, galaxy_candidates["dec_deg"].to_numpy(dtype=float) * u.deg).separation(world).arcsec
            if len(sep) and float(np.nanmin(sep)) <= 0.4:
                matched_ids.append(str(galaxy_candidates.iloc[int(np.nanargmin(sep))]["candidate_id"]))
            else:
                matched_ids.append(None)
        ranking["matched_candidate_id"] = matched_ids
        ranking = ranking.sort_values(["fade_sigma", "fade_fraction"], ascending=[False, False]).reset_index(drop=True)
        reference_path = galaxy_dir / "reference_catalog.parquet"
        forced_path = galaxy_dir / "forced_photometry.parquet"
        ranking_path = galaxy_dir / "fade_ranking.csv"
        pd.DataFrame(epoch_rows).to_csv(galaxy_dir / "epochs.csv", index=False)
        sources.assign(source_id=np.arange(len(sources))).to_parquet(reference_path, index=False)
        forced.to_parquet(forced_path, index=False)
        ranking.to_csv(ranking_path, index=False)
        write_json(
            galaxy_dir / "summary.json",
            {
                "galaxy_name": galaxy_name,
                "n_images": int(len(loaded)),
                "n_sources": int(len(sources)),
                "filter_name": str(image_rows[0]["filter_name"]),
            },
        )
        summary_rows.append(
            {
                "galaxy_name": galaxy_name,
                "n_images": int(len(loaded)),
                "n_sources": int(len(sources)),
                "ranking_path": str(ranking_path),
            }
        )
    summary_path = out_dir / "summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    return {"summary": summary_path}


def _mid_ir_query(candidate: pd.Series) -> dict[str, Any]:
    coord = _candidate_world(candidate)
    row: dict[str, Any] = {
        "candidate_id": str(candidate["candidate_id"]),
        "galaxy_name": str(candidate["galaxy_name"]),
        "irsa_catalog": "allwise_p3as_psd",
        "status": "NO_MATCH",
        "match_sep_arcsec": float("nan"),
        "w1mpro": float("nan"),
        "w2mpro": float("nan"),
        "w3mpro": float("nan"),
        "w4mpro": float("nan"),
        "w1_w2": float("nan"),
        "mid_ir_survivor_flag": False,
        "dusty_color_flag": False,
        "error": None,
    }
    try:
        table = Irsa.query_region(coord, catalog="allwise_p3as_psd", spatial="Cone", radius=2.5 * u.arcsec)
        if len(table) == 0:
            return row
        df = table.to_pandas()
        sep = coord.separation(SkyCoord(df["ra"].to_numpy(dtype=float) * u.deg, df["dec"].to_numpy(dtype=float) * u.deg)).arcsec
        idx = int(np.nanargmin(sep))
        match = df.iloc[idx]
        row.update(
            {
                "status": "MATCHED",
                "match_sep_arcsec": float(sep[idx]),
                "w1mpro": _safe_float(match.get("w1mpro")),
                "w2mpro": _safe_float(match.get("w2mpro")),
                "w3mpro": _safe_float(match.get("w3mpro")),
                "w4mpro": _safe_float(match.get("w4mpro")),
            }
        )
        if math.isfinite(row["w1mpro"]) and math.isfinite(row["w2mpro"]):
            row["w1_w2"] = float(row["w1mpro"] - row["w2mpro"])
        row["dusty_color_flag"] = bool(math.isfinite(row["w1_w2"]) and row["w1_w2"] >= 0.5)
        row["mid_ir_survivor_flag"] = bool(
            row["match_sep_arcsec"] <= 1.5
            and (
                (math.isfinite(row["w1mpro"]) and row["w1mpro"] <= 17.0)
                or (math.isfinite(row["w2mpro"]) and row["w2mpro"] <= 16.5)
            )
        )
    except Exception as exc:
        row["status"] = "FAILED"
        row["error"] = str(exc)
    return row


def _run_mid_ir_stage(root_dir: Path, ledger: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "mid_ir")
    rows = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(_mid_ir_query, row): str(row["candidate_id"]) for _, row in ledger.iterrows()}
        for future in as_completed(futures):
            rows.append(future.result())
    df = pd.DataFrame(rows).sort_values(["mid_ir_survivor_flag", "match_sep_arcsec"], ascending=[False, True]).reset_index(drop=True)
    output = out_dir / "mid_ir_crossmatch.csv"
    df.to_csv(output, index=False)
    write_json(
        out_dir / "summary.json",
        {
            "n_candidates": int(len(df)),
            "n_matches": int((df["status"] == "MATCHED").sum()),
            "n_survivor_flags": int(df["mid_ir_survivor_flag"].fillna(False).sum()),
            "n_dusty_flags": int(df["dusty_color_flag"].fillna(False).sum()),
        },
    )
    return {"crossmatch": output, "summary": out_dir / "summary.json"}


def _run_neighborhood_stage(root_dir: Path, scene_df: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "neighborhood")
    df = scene_df.copy()
    if df.empty:
        df = pd.DataFrame(columns=["candidate_id"])
    if "scene_neighbor_scatter" in df.columns:
        scatter = pd.to_numeric(df["scene_neighbor_scatter"], errors="coerce").fillna(np.nan)
        df["local_unique_sigma"] = pd.to_numeric(df["scene_unique_fade_sigma"], errors="coerce")
        df["local_systematic_risk"] = 1.0 / (1.0 + np.clip(df["local_unique_sigma"], 0.0, np.inf))
        df["local_coherent_trend"] = pd.to_numeric(df["scene_neighbor_median_fade"], errors="coerce").fillna(0.0)
        df["local_neighbor_count"] = pd.to_numeric(df["scene_local_comparators"], errors="coerce").fillna(0).astype(int)
        df["local_verdict"] = np.where(
            (df["local_unique_sigma"].fillna(0.0) >= 3.0) & (df["local_neighbor_count"] >= 2),
            "UNIQUE",
            np.where(df["local_coherent_trend"] >= 0.2, "SYSTEMATIC_TREND", "WEAK"),
        )
    output = out_dir / "neighborhood_scores.csv"
    keep_cols = [
        "candidate_id",
        "galaxy_name",
        "local_unique_sigma",
        "local_systematic_risk",
        "local_coherent_trend",
        "local_neighbor_count",
        "local_verdict",
    ]
    df.reindex(columns=keep_cols).to_csv(output, index=False)
    write_json(
        out_dir / "summary.json",
        {
            "n_candidates": int(len(df)),
            "n_unique": int((df.get("local_verdict", pd.Series(dtype=str)) == "UNIQUE").sum()) if not df.empty else 0,
            "n_systematic_trend": int((df.get("local_verdict", pd.Series(dtype=str)) == "SYSTEMATIC_TREND").sum()) if not df.empty else 0,
        },
    )
    return {"scores": output, "summary": out_dir / "summary.json"}


def _run_anomaly_stage(root_dir: Path, ledger: pd.DataFrame, scene_df: pd.DataFrame, neighborhood_df: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "anomaly")
    scene_df = _ensure_columns(
        scene_df,
        [
            "candidate_id",
            "galaxy_name",
            "scene_fade_sigma",
            "scene_fade_fraction",
            "scene_unique_fade_sigma",
            "scene_monotonic_score",
        ],
    )
    neighborhood_df = _ensure_columns(neighborhood_df, ["candidate_id", "local_unique_sigma", "local_systematic_risk"])
    merged = ledger.merge(scene_df, on=["candidate_id", "galaxy_name"], how="left").merge(
        neighborhood_df[["candidate_id", "local_unique_sigma", "local_systematic_risk"]],
        on="candidate_id",
        how="left",
    )
    feature_cols = [
        "priority_score",
        "fade_sigma",
        "fade_fraction",
        "crowding_index",
        "blend_risk",
        "host_bg_penalty",
        "depth_margin_post",
        "astrometric_residual",
        "cross_reduction_agreement",
        "scene_fade_sigma",
        "scene_fade_fraction",
        "scene_unique_fade_sigma",
        "scene_monotonic_score",
        "local_unique_sigma",
        "local_systematic_risk",
    ]
    X = merged.reindex(columns=feature_cols).apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(numeric_only=True)).fillna(0.0)
    x_arr = X.to_numpy(dtype=float)
    x_mean = x_arr.mean(axis=0, keepdims=True)
    x_std = x_arr.std(axis=0, keepdims=True)
    x_std[x_std <= 0] = 1.0
    z = (x_arr - x_mean) / x_std
    if len(z) >= 3 and z.shape[1] >= 2:
        u, s, vt = np.linalg.svd(z, full_matrices=False)
        k = int(max(1, min(3, min(z.shape) - 1)))
        z_hat = (u[:, :k] * s[:k]) @ vt[:k]
        recon = np.sqrt(np.sum(np.square(z - z_hat), axis=1))
        latent = u[:, :k] * s[:k]
        latent_std = latent.std(axis=0, keepdims=True)
        latent_std[latent_std <= 0] = 1.0
        dist = np.sqrt(np.sum(np.square(latent / latent_std), axis=1))
    else:
        recon = np.zeros(len(z), dtype=float)
        dist = np.sqrt(np.sum(np.square(z), axis=1))
    recon_s = pd.Series(recon)
    dist_s = pd.Series(dist)
    merged["anomaly_reconstruction"] = recon
    merged["anomaly_distance"] = dist
    merged["anomaly_score"] = normalize_01(recon_s) * 0.6 + normalize_01(dist_s) * 0.4
    feature_path = out_dir / "anomaly_features.parquet"
    score_path = out_dir / "anomaly_scores.csv"
    merged[["candidate_id", "galaxy_name", *feature_cols, "anomaly_reconstruction", "anomaly_distance", "anomaly_score"]].to_parquet(feature_path, index=False)
    merged[["candidate_id", "galaxy_name", "status", "priority_score", "anomaly_score"]].sort_values("anomaly_score", ascending=False).to_csv(score_path, index=False)
    write_json(out_dir / "summary.json", {"n_candidates": int(len(merged)), "max_anomaly_score": float(merged["anomaly_score"].max()) if len(merged) else 0.0})
    return {"features": feature_path, "scores": score_path, "summary": out_dir / "summary.json"}


def _run_population_prior_stage(root_dir: Path, ledger: pd.DataFrame, galaxy_master: pd.DataFrame, scene_df: pd.DataFrame, neighborhood_df: pd.DataFrame, anomaly_scores: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "population_prior")
    scene_df = _ensure_columns(scene_df, ["candidate_id", "scene_unique_fade_sigma", "scene_monotonic_score"])
    neighborhood_df = _ensure_columns(neighborhood_df, ["candidate_id", "local_unique_sigma", "local_systematic_risk"])
    anomaly_scores = _ensure_columns(anomaly_scores, ["candidate_id", "anomaly_score"])
    host_cols = ["galaxy_name", "distance_mpc", "mw_av_mag", "searchability_score", "sfr_proxy_msun_per_yr", "stellar_mass_msun"]
    host = galaxy_master.reindex(columns=host_cols).copy()
    merged = ledger.merge(host, on="galaxy_name", how="left")
    merged = merged.merge(scene_df[["candidate_id", "scene_unique_fade_sigma", "scene_monotonic_score"]], on="candidate_id", how="left")
    merged = merged.merge(neighborhood_df[["candidate_id", "local_unique_sigma", "local_systematic_risk"]], on="candidate_id", how="left")
    merged = merged.merge(anomaly_scores[["candidate_id", "anomaly_score"]], on="candidate_id", how="left")
    merged["baseline_days"] = pd.to_numeric(merged["post_mjd"], errors="coerce") - pd.to_numeric(merged["pre_mjd"], errors="coerce")
    merged["distance_prior"] = 1.0 - normalize_01(pd.to_numeric(merged["distance_mpc"], errors="coerce"))
    merged["host_prior"] = (
        0.45 * normalize_01(pd.to_numeric(merged["searchability_score"], errors="coerce"))
        + 0.25 * normalize_01(robust_log10(pd.to_numeric(merged["sfr_proxy_msun_per_yr"], errors="coerce").fillna(1e-6) + 1e-6))
        + 0.20 * merged["distance_prior"]
        + 0.10 * (1.0 - normalize_01(pd.to_numeric(merged["mw_av_mag"], errors="coerce")))
    )
    merged["candidate_prior"] = (
        0.30 * normalize_01(pd.to_numeric(merged["fade_fraction"], errors="coerce"))
        + 0.20 * normalize_01(pd.to_numeric(merged["fade_sigma"], errors="coerce"))
        + 0.15 * normalize_01(pd.to_numeric(merged["scene_unique_fade_sigma"], errors="coerce"))
        + 0.15 * normalize_01(pd.to_numeric(merged["local_unique_sigma"], errors="coerce"))
        + 0.10 * normalize_01(robust_log10(pd.to_numeric(merged["baseline_days"], errors="coerce").fillna(1.0)))
        + 0.10 * (1.0 - normalize_01(pd.to_numeric(merged["anomaly_score"], errors="coerce")))
    )
    merged["population_prior_score"] = 0.45 * merged["host_prior"] + 0.55 * merged["candidate_prior"]
    out_path = out_dir / "population_prior_scores.csv"
    merged[
        [
            "candidate_id",
            "galaxy_name",
            "status",
            "population_prior_score",
            "host_prior",
            "candidate_prior",
            "distance_mpc",
            "searchability_score",
            "sfr_proxy_msun_per_yr",
            "baseline_days",
        ]
    ].sort_values("population_prior_score", ascending=False).to_csv(out_path, index=False)
    write_json(out_dir / "summary.json", {"n_candidates": int(len(merged))})
    return {"scores": out_path, "summary": out_dir / "summary.json"}


def _gaussian_stamp(size: int, sigma_pix: float, amplitude: float) -> np.ndarray:
    y, x = np.mgrid[:size, :size]
    c = (size - 1) / 2.0
    return amplitude * np.exp(-((x - c) ** 2 + (y - c) ** 2) / (2.0 * sigma_pix**2))


def _inject_point_source(data: np.ndarray, x: float, y: float, amplitude: float, sigma_pix: float) -> np.ndarray:
    stamp_size = int(max(9, math.ceil(8 * sigma_pix)))
    if stamp_size % 2 == 0:
        stamp_size += 1
    stamp = _gaussian_stamp(stamp_size, sigma_pix, amplitude)
    out = np.array(data, dtype=float, copy=True)
    half = stamp_size // 2
    x0 = int(round(x)) - half
    y0 = int(round(y)) - half
    x1 = x0 + stamp_size
    y1 = y0 + stamp_size
    xs0 = max(x0, 0)
    ys0 = max(y0, 0)
    xs1 = min(x1, out.shape[1])
    ys1 = min(y1, out.shape[0])
    if xs0 >= xs1 or ys0 >= ys1:
        return out
    sx0 = xs0 - x0
    sy0 = ys0 - y0
    sx1 = sx0 + (xs1 - xs0)
    sy1 = sy0 + (ys1 - ys0)
    out[ys0:ys1, xs0:xs1] += stamp[sy0:sy1, sx0:sx1]
    return out


def _injection_trials_for_candidate(candidate: pd.Series, n_trials_per_ratio: int = 8) -> pd.DataFrame:
    provenance = _decode_provenance(candidate)
    product_info = provenance.get("product_info", {})
    path_1 = Path(str(product_info.get("product_path_1", "")))
    path_2 = Path(str(product_info.get("product_path_2", "")))
    if not path_1.exists() or not path_2.exists():
        return pd.DataFrame()
    obs_id_1 = str(candidate.get("pre_obs_id", product_info.get("obs_id_1", "pre_obs")))
    obs_id_2 = str(candidate.get("post_obs_id", product_info.get("obs_id_2", "post_obs")))
    image_1 = _load_science_image(
        path_1,
        obs_collection=_infer_obs_collection_from_path(path_1),
        obs_id=obs_id_1,
        filter_name=str(candidate.get("pre_filter", "")),
        instrument_name=_infer_instrument_name(obs_id_1),
    )
    image_2 = _load_science_image(
        path_2,
        obs_collection=_infer_obs_collection_from_path(path_2),
        obs_id=obs_id_2,
        filter_name=str(candidate.get("post_filter", "")),
        instrument_name=_infer_instrument_name(obs_id_2),
    )
    world = _candidate_world(candidate)
    image_1 = _candidate_cutout(image_1, world, size_arcsec=90.0, min_pixels=120)
    image_2 = _candidate_cutout(image_2, world, size_arcsec=90.0, min_pixels=120)
    post_reproj, footprint = reproject_interp((image_2.data, image_2.wcs), image_1.wcs, shape_out=image_1.data.shape)
    valid = np.isfinite(image_1.data) & np.isfinite(post_reproj) & (footprint > 0)
    if valid.sum() < 1_000:
        return pd.DataFrame()
    pre_bg, _ = _background_stats(image_1.data, valid)
    post_bg, _ = _background_stats(post_reproj, valid)
    pre_sub = image_1.data - pre_bg
    post_sub = post_reproj - post_bg
    sources = _detect_sources(pre_sub, valid, pixel_scale_arcsec=image_1.pixel_scale_arcsec, max_sources=250)
    if sources.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(17)
    isolated = sources.copy()
    isolated["edge"] = np.minimum.reduce([isolated["x"], isolated["y"], pre_sub.shape[1] - 1 - isolated["x"], pre_sub.shape[0] - 1 - isolated["y"]])
    isolated = isolated[isolated["edge"] >= 10].copy()
    if isolated.empty:
        return pd.DataFrame()
    sample = isolated.sample(n=min(12, len(isolated)), random_state=17)
    radius = max(2.5, 1.2 * _fwhm_guess(image_1.pixel_scale_arcsec))
    positions = sample[["x", "y"]].to_numpy(dtype=float)
    flux_pre, err_pre, _ = _aperture_fluxes(pre_sub, valid, positions, radius=radius)
    sigma_pix = max(_fwhm_guess(image_1.pixel_scale_arcsec) / 2.355, 1.0)
    ratios = [1.0, 0.2, 0.05, 0.0]
    rows = []
    for ratio in ratios:
        for trial_idx, (pos, base_flux) in enumerate(zip(positions, flux_pre, strict=False)):
            amplitude = float(max(abs(base_flux), np.nanmedian(abs(flux_pre)))) / max(2.0 * np.pi * sigma_pix**2, 1.0)
            amplitude *= float(rng.uniform(0.8, 1.2))
            pre_inj = _inject_point_source(pre_sub, pos[0], pos[1], amplitude, sigma_pix)
            post_inj = _inject_point_source(post_sub, pos[0], pos[1], amplitude * ratio, sigma_pix)
            inj_pre, inj_pre_err, _ = _aperture_fluxes(pre_inj, valid, np.array([pos]), radius=radius)
            inj_post, inj_post_err, _ = _aperture_fluxes(post_inj, valid, np.array([pos]), radius=radius)
            fade = 1.0 - (inj_post[0] / inj_pre[0]) if math.isfinite(inj_pre[0]) and abs(inj_pre[0]) > 0 else float("nan")
            sigma = (inj_pre[0] - inj_post[0]) / math.sqrt(max(inj_pre_err[0] ** 2 + inj_post_err[0] ** 2, 1e-9))
            recovered = bool(math.isfinite(fade) and math.isfinite(sigma) and fade >= (1.0 - ratio) * 0.6 and sigma >= 5.0)
            rows.append(
                {
                    "candidate_id": str(candidate["candidate_id"]),
                    "galaxy_name": str(candidate["galaxy_name"]),
                    "trial_index": int(trial_idx),
                    "injected_ratio": float(ratio),
                    "injected_fade_fraction": float(1.0 - ratio),
                    "recovered_fade_fraction": float(fade) if math.isfinite(fade) else float("nan"),
                    "recovered_sigma": float(sigma) if math.isfinite(sigma) else float("nan"),
                    "recovered": recovered,
                }
            )
            if trial_idx + 1 >= n_trials_per_ratio:
                break
    return pd.DataFrame(rows)


def _run_injection_stage(root_dir: Path, ledger: pd.DataFrame) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "injection_recovery")
    selected_ids = []
    for galaxy_name in ["NGC 5861", "MESSIER 101", "MESSIER 077"]:
        sub = ledger[ledger["galaxy_name"].astype(str) == galaxy_name].sort_values("priority_score", ascending=False).head(1)
        if not sub.empty:
            selected_ids.append(str(sub.iloc[0]["candidate_id"]))
    selected = ledger[ledger["candidate_id"].astype(str).isin(selected_ids)].copy()
    frames = [_injection_trials_for_candidate(row) for _, row in selected.iterrows()]
    trials = pd.concat([frame for frame in frames if not frame.empty], ignore_index=True) if any(not frame.empty for frame in frames) else pd.DataFrame()
    trials_path = out_dir / "injection_trials.parquet"
    summary_path = out_dir / "recovery_summary.csv"
    if not trials.empty:
        trials.to_parquet(trials_path, index=False)
        summary = (
            trials.groupby(["candidate_id", "injected_ratio"], dropna=False)
            .agg(n_trials=("recovered", "size"), recovery_rate=("recovered", "mean"), median_sigma=("recovered_sigma", "median"))
            .reset_index()
        )
        summary.to_csv(summary_path, index=False)
        recommended = summary.groupby("injected_ratio", dropna=False)["recovery_rate"].median().reset_index()
        detectable = recommended[recommended["recovery_rate"] >= 0.8].sort_values("injected_ratio")
        threshold_ratio = float(detectable["injected_ratio"].max()) if not detectable.empty else 0.05
    else:
        pd.DataFrame().to_csv(summary_path, index=False)
        threshold_ratio = 0.05
    write_json(
        out_dir / "summary.json",
        {
            "n_candidate_pairs": int(len(selected)),
            "n_trials": int(len(trials)),
            "recommended_max_post_flux_ratio": threshold_ratio,
        },
    )
    return {"trials": trials_path, "summary": summary_path, "json": out_dir / "summary.json"}


def _run_combined_stage(
    root_dir: Path,
    ledger: pd.DataFrame,
    scene_df: pd.DataFrame,
    neighborhood_df: pd.DataFrame,
    anomaly_scores: pd.DataFrame,
    population_scores: pd.DataFrame,
    mid_ir_df: pd.DataFrame,
) -> dict[str, Path]:
    out_dir = ensure_dir(root_dir / "extensions" / "combined")
    scene_df = _ensure_columns(scene_df, ["candidate_id", "scene_unique_fade_sigma", "scene_monotonic_score"])
    neighborhood_df = _ensure_columns(neighborhood_df, ["candidate_id", "local_unique_sigma", "local_systematic_risk", "local_verdict"])
    anomaly_scores = _ensure_columns(anomaly_scores, ["candidate_id", "anomaly_score"])
    population_scores = _ensure_columns(population_scores, ["candidate_id", "population_prior_score"])
    mid_ir_df = _ensure_columns(mid_ir_df, ["candidate_id", "mid_ir_survivor_flag", "dusty_color_flag"])
    combined = ledger.copy()
    combined = combined.merge(scene_df[["candidate_id", "scene_unique_fade_sigma", "scene_monotonic_score"]], on="candidate_id", how="left")
    combined = combined.merge(neighborhood_df[["candidate_id", "local_unique_sigma", "local_systematic_risk", "local_verdict"]], on="candidate_id", how="left")
    combined = combined.merge(anomaly_scores[["candidate_id", "anomaly_score"]], on="candidate_id", how="left")
    combined = combined.merge(population_scores[["candidate_id", "population_prior_score"]], on="candidate_id", how="left")
    combined = combined.merge(mid_ir_df[["candidate_id", "mid_ir_survivor_flag", "dusty_color_flag"]], on="candidate_id", how="left")
    combined["extended_rank_score"] = (
        0.28 * normalize_01(pd.to_numeric(combined["priority_score"], errors="coerce"))
        + 0.20 * normalize_01(pd.to_numeric(combined["population_prior_score"], errors="coerce"))
        + 0.17 * normalize_01(pd.to_numeric(combined["scene_unique_fade_sigma"], errors="coerce"))
        + 0.10 * normalize_01(pd.to_numeric(combined["scene_monotonic_score"], errors="coerce"))
        + 0.10 * normalize_01(pd.to_numeric(combined["local_unique_sigma"], errors="coerce"))
        + 0.10 * (1.0 - normalize_01(pd.to_numeric(combined["anomaly_score"], errors="coerce")))
        + 0.05 * (1.0 - pd.to_numeric(combined["mid_ir_survivor_flag"], errors="coerce").fillna(0.0))
    )
    combined["extension_verdict"] = np.where(
        combined["mid_ir_survivor_flag"].fillna(False),
        "MID_IR_SURVIVOR",
        np.where(
            combined["local_verdict"].astype(str).eq("UNIQUE") & (combined["extended_rank_score"] >= combined["extended_rank_score"].median()),
            "UPRANK",
            "KEEP",
        ),
    )
    output = out_dir / "candidate_extension_ranking.csv"
    combined.sort_values(["extended_rank_score", "priority_score"], ascending=[False, False]).to_csv(output, index=False)
    write_json(
        out_dir / "summary.json",
        {
            "n_candidates": int(len(combined)),
            "n_uprank": int((combined["extension_verdict"] == "UPRANK").sum()),
            "n_mid_ir_survivor": int((combined["extension_verdict"] == "MID_IR_SURVIVOR").sum()),
        },
    )
    return {"ranking": output, "summary": out_dir / "summary.json"}


def run_extensions(
    *,
    root_dir: Path,
    candidate_ledger_path: Path,
    observation_matrix_path: Path,
    galaxy_master_path: Path,
) -> dict[str, Path]:
    ensure_dir(root_dir / "extensions")
    ledger = pd.read_parquet(candidate_ledger_path).copy()
    observations = pd.read_parquet(observation_matrix_path)
    galaxy_master = pd.read_parquet(galaxy_master_path)
    del observations  # reserved for future extension stages

    scene_df, neighbors_df, scene_paths = _run_scene_model_stage(root_dir, ledger)
    deep_paths = _run_deep_reference_stage(root_dir, ledger)
    mid_paths = _run_mid_ir_stage(root_dir, ledger)
    neighborhood_paths = _run_neighborhood_stage(root_dir, scene_df)

    neighborhood_df = pd.read_csv(neighborhood_paths["scores"]) if Path(neighborhood_paths["scores"]).exists() else pd.DataFrame()
    anomaly_paths = _run_anomaly_stage(root_dir, ledger, scene_df, neighborhood_df)
    anomaly_df = pd.read_csv(anomaly_paths["scores"]) if Path(anomaly_paths["scores"]).exists() else pd.DataFrame()
    population_paths = _run_population_prior_stage(root_dir, ledger, galaxy_master, scene_df, neighborhood_df, anomaly_df)
    population_df = pd.read_csv(population_paths["scores"]) if Path(population_paths["scores"]).exists() else pd.DataFrame()
    injection_paths = _run_injection_stage(root_dir, ledger)
    mid_ir_df = pd.read_csv(mid_paths["crossmatch"]) if Path(mid_paths["crossmatch"]).exists() else pd.DataFrame()
    combined_paths = _run_combined_stage(root_dir, ledger, scene_df, neighborhood_df, anomaly_df, population_df, mid_ir_df)

    manifest = {
        "extended_plan": str((root_dir / "extended_plan.md").resolve()),
        "scene_metrics": str(Path(scene_paths["scene_metrics"]).resolve()),
        "deep_reference_summary": str(Path(deep_paths["summary"]).resolve()),
        "mid_ir_crossmatch": str(Path(mid_paths["crossmatch"]).resolve()),
        "anomaly_scores": str(Path(anomaly_paths["scores"]).resolve()),
        "population_prior_scores": str(Path(population_paths["scores"]).resolve()),
        "injection_summary": str(Path(injection_paths["summary"]).resolve()),
        "neighborhood_scores": str(Path(neighborhood_paths["scores"]).resolve()),
        "combined_ranking": str(Path(combined_paths["ranking"]).resolve()),
    }
    manifest_path = root_dir / "extensions" / "manifest.json"
    write_json(manifest_path, manifest)
    return {
        "manifest": manifest_path,
        "scene_metrics": Path(scene_paths["scene_metrics"]),
        "deep_reference_summary": Path(deep_paths["summary"]),
        "mid_ir_crossmatch": Path(mid_paths["crossmatch"]),
        "anomaly_scores": Path(anomaly_paths["scores"]),
        "population_prior_scores": Path(population_paths["scores"]),
        "injection_summary": Path(injection_paths["summary"]),
        "neighborhood_scores": Path(neighborhood_paths["scores"]),
        "combined_ranking": Path(combined_paths["ranking"]),
    }
