from __future__ import annotations

import concurrent.futures as cf
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests

from .utils import ensure_dir, utc_stamp, write_json


OBS_COLLECTIONS = {"HST", "JWST", "HLSP"}
MAST_INVOKE_URL = "https://mast.stsci.edu/api/v0/invoke"


def _mast_cone_query(ra_deg: float, dec_deg: float, radius_deg: float, *, page_size: int = 2000) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    page = 1
    pages_filtered = 1
    while page <= pages_filtered:
        payload = {
            "service": "Mast.Caom.Cone",
            "params": {
                "ra": ra_deg,
                "dec": dec_deg,
                "radius": radius_deg,
                "cachebreaker": utc_stamp(),
            },
            "format": "json",
            "pagesize": page_size,
            "page": page,
            "removenullcolumns": True,
            "timeout": 30,
        }
        resp = requests.post(MAST_INVOKE_URL, data={"request": json.dumps(payload)}, timeout=60)
        resp.raise_for_status()
        out = resp.json()
        rows.extend(out.get("data", []))
        paging = out.get("paging", {})
        pages_filtered = int(paging.get("pagesFiltered", 1))
        page += 1
    return pd.DataFrame(rows)


def infer_filter_central_um(filter_name: str | None, instrument_name: str | None, obs_collection: str | None) -> float | None:
    filter_name = str(filter_name or "")
    instrument_name = str(instrument_name or "").upper()
    obs_collection = str(obs_collection or "").upper()
    tokens = re.split(r"[;,/ ]+", filter_name)
    candidate = next((tok for tok in tokens if tok.upper().startswith("F") and any(ch.isdigit() for ch in tok)), "")
    if not candidate:
        return None
    match = re.search(r"F(\d{2,4})", candidate.upper())
    if not match:
        return None
    value = int(match.group(1))
    if obs_collection == "JWST" or "NIRCAM" in instrument_name or "NIRISS" in instrument_name or "NIRSPEC" in instrument_name or "MIRI" in instrument_name:
        return value / 100.0
    if "WFC3/IR" in instrument_name or "NICMOS" in instrument_name or value < 200:
        return value / 100.0
    return value / 1000.0


def filter_compatibility(
    filter_1: str,
    filter_2: str,
    instrument_1: str,
    instrument_2: str,
    collection_1: str,
    collection_2: str,
) -> tuple[str, float]:
    lam1 = infer_filter_central_um(filter_1, instrument_1, collection_1)
    lam2 = infer_filter_central_um(filter_2, instrument_2, collection_2)
    if filter_1 and filter_2 and filter_1 == filter_2:
        return "exact", 1.0
    if lam1 is None or lam2 is None:
        return "unknown", 0.2
    delta = abs(lam1 - lam2)
    if delta <= 0.10:
        return "very_similar", 0.85
    if delta <= 0.25:
        return "similar", 0.65
    if (lam1 < 0.9 and lam2 < 0.9) or (0.9 <= lam1 < 3.0 and 0.9 <= lam2 < 3.0) or (lam1 >= 3.0 and lam2 >= 3.0):
        return "same_band", 0.40
    return "incompatible", 0.0


def query_galaxy_observations(row: pd.Series, *, radius_deg: float) -> pd.DataFrame:
    coord = SkyCoord(float(row["ra"]) * u.deg, float(row["dec"]) * u.deg, frame="icrs")
    df = _mast_cone_query(float(row["ra"]), float(row["dec"]), radius_deg)
    if df.empty:
        return df

    df = df[df["intentType"].astype(str).str.lower() == "science"].copy()
    df = df[df["obs_collection"].astype(str).isin(OBS_COLLECTIONS)].copy()
    if "dataproduct_type" in df.columns:
        df = df[df["dataproduct_type"].astype(str).str.lower().isin({"image", "cube"})].copy()
    if df.empty:
        return df

    keep_cols = [
        c
        for c in [
            "obsid",
            "obs_id",
            "target_name",
            "obs_collection",
            "proposal_id",
            "filters",
            "instrument_name",
            "t_min",
            "t_max",
            "s_ra",
            "s_dec",
            "calib_level",
            "intentType",
            "dataRights",
            "dataproduct_type",
            "obs_title",
        ]
        if c in df.columns
    ]
    df = df[keep_cols].copy()
    sep = SkyCoord(df["s_ra"].astype(float).to_numpy() * u.deg, df["s_dec"].astype(float).to_numpy() * u.deg).separation(coord)
    df["separation_arcmin"] = sep.arcmin
    df["galaxy_name"] = str(row["galaxy_name"])
    df["galaxy_rank"] = int(row["rank"])
    df["galaxy_ra"] = float(row["ra"])
    df["galaxy_dec"] = float(row["dec"])
    df["search_radius_deg"] = float(radius_deg)
    df["filters"] = df["filters"].astype(str).str.strip()
    df = df[~df["filters"].str.lower().isin({"detection", "none", "", "nan"})].copy()
    df["filter_central_um"] = [
        infer_filter_central_um(f, i, c) for f, i, c in zip(df["filters"], df["instrument_name"], df["obs_collection"], strict=True)
    ]
    df["is_public"] = df["dataRights"].astype(str).str.upper().ne("PROPRIETARY")
    df["epoch_mjd"] = np.floor(pd.to_numeric(df["t_min"], errors="coerce")).astype("Int64")
    dedup_cols = [c for c in ["galaxy_name", "obs_collection", "filters", "instrument_name", "epoch_mjd"] if c in df.columns]
    if dedup_cols:
        df = df.sort_values(["t_min", "t_max"]).drop_duplicates(subset=dedup_cols, keep="first").reset_index(drop=True)
    return df


def build_epoch_pairs(obs: pd.DataFrame, *, min_baseline_days: float) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for galaxy_name, group in obs.groupby("galaxy_name"):
        group = group.sort_values("t_min").reset_index(drop=True)
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                left = group.iloc[i]
                right = group.iloc[j]
                baseline_days = abs(float(right["t_min"]) - float(left["t_min"]))
                if baseline_days < min_baseline_days:
                    continue
                compat_label, compat_score = filter_compatibility(
                    str(left["filters"]),
                    str(right["filters"]),
                    str(left["instrument_name"]),
                    str(right["instrument_name"]),
                    str(left["obs_collection"]),
                    str(right["obs_collection"]),
                )
                if compat_score <= 0.0:
                    continue
                pair_score = compat_score + 0.15 * min(math.log10(baseline_days + 1.0), 2.0)
                if str(left["obs_collection"]) == str(right["obs_collection"]):
                    pair_score += 0.1
                rows.append(
                    {
                        "pair_id": f"{galaxy_name}::{left['obs_id']}::{right['obs_id']}",
                        "galaxy_name": galaxy_name,
                        "galaxy_rank": int(left["galaxy_rank"]),
                        "obsid_1": int(left["obsid"]) if pd.notna(left.get("obsid")) else pd.NA,
                        "obsid_2": int(right["obsid"]) if pd.notna(right.get("obsid")) else pd.NA,
                        "obs_id_1": str(left["obs_id"]),
                        "obs_id_2": str(right["obs_id"]),
                        "obs_collection_1": str(left["obs_collection"]),
                        "obs_collection_2": str(right["obs_collection"]),
                        "filter_1": str(left["filters"]),
                        "filter_2": str(right["filters"]),
                        "instrument_1": str(left["instrument_name"]),
                        "instrument_2": str(right["instrument_name"]),
                        "baseline_days": baseline_days,
                        "compatibility": compat_label,
                        "compatibility_score": compat_score,
                        "pair_score": pair_score,
                        "public_both": bool(left["is_public"] and right["is_public"]),
                    }
                )
    pairs = pd.DataFrame(rows)
    if not pairs.empty:
        pairs = pairs.sort_values(["pair_score", "baseline_days"], ascending=[False, False]).reset_index(drop=True)
    return pairs


def build_archive_products(
    *,
    root_dir: Path,
    galaxy_master_path: Path,
    archive_top_n: int,
    radius_deg: float,
    min_baseline_days: float,
    max_workers: int = 4,
) -> dict[str, Path]:
    archive_dir = ensure_dir(root_dir / "archive")
    galaxies = pd.read_parquet(galaxy_master_path).head(archive_top_n).copy()

    frames: list[pd.DataFrame] = []
    with cf.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(query_galaxy_observations, row, radius_deg=radius_deg): row["galaxy_name"] for _, row in galaxies.iterrows()}
        for future in cf.as_completed(futures):
            df = future.result()
            if not df.empty:
                frames.append(df)

    observations = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not observations.empty:
        observations = observations.sort_values(["galaxy_rank", "t_min", "obs_collection", "filters"]).reset_index(drop=True)
    epoch_pairs = build_epoch_pairs(observations, min_baseline_days=min_baseline_days) if not observations.empty else pd.DataFrame()

    observation_matrix_path = archive_dir / "observation_matrix.parquet"
    epoch_pairs_path = archive_dir / "epoch_pairs.parquet"
    summary_path = archive_dir / "archive_summary.json"

    observations.to_parquet(observation_matrix_path, index=False)
    epoch_pairs.to_parquet(epoch_pairs_path, index=False)
    write_json(
        summary_path,
        {
            "created_utc": utc_stamp(),
            "archive_top_n": int(archive_top_n),
            "radius_deg": float(radius_deg),
            "min_baseline_days": float(min_baseline_days),
            "n_galaxies_queried": int(len(galaxies)),
            "n_observations": int(len(observations)),
            "n_epoch_pairs": int(len(epoch_pairs)),
            "collections": observations["obs_collection"].value_counts().to_dict() if not observations.empty else {},
            "top_pair_galaxies": epoch_pairs["galaxy_name"].head(20).tolist() if not epoch_pairs.empty else [],
        },
    )

    return {
        "observation_matrix": observation_matrix_path,
        "epoch_pairs": epoch_pairs_path,
        "summary": summary_path,
    }
