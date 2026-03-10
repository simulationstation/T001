from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u

from .utils import angular_sep_deg, decode_object_series, ensure_dir, normalize_01, robust_log10, utc_stamp, write_json


NED_LVS_URL = "https://ned.ipac.caltech.edu/NED::LVS/fits/AsPublished/"
DEFAULT_GLADEPLUS_INDEX_DIR = Path(
    "/home/primary/Dark_Siren_Ladder_Audit/data/processed/galaxies/gladeplus/index_nside128_wlumB_zmax0.3"
)
NED_LVS_COLUMNS = [
    "objname",
    "ra",
    "dec",
    "objtype",
    "z",
    "DistMpc",
    "DistMpc_unc",
    "ebv",
    "SFR_W4",
    "SFR_W4_unc",
    "SFR_hybrid",
    "SFR_hybrid_unc",
    "ET_flag",
    "Mstar",
    "Mstar_unc",
    "Lum_W1",
    "Lum_W1_unc",
]

EXPLICIT_EXCLUSION_PATTERNS = [
    "M31",
    "MESSIER 031",
    "ANDROMEDA",
    "LARGE MAGELLANIC",
    "SMALL MAGELLANIC",
    "LMC",
    "SMC",
]


@dataclass(slots=True)
class GalaxyBuildConfig:
    max_distance_mpc: float = 40.0
    max_av_mag: float = 0.5
    min_sfr: float = 0.0
    top_n: int = 130
    gladeplus_radius_deg: float = 0.25
    gladeplus_z_window: float = 0.003


class GladePlusLocalIndex:
    def __init__(self, index_dir: Path) -> None:
        self.index_dir = index_dir
        meta_path = index_dir / "meta.json"
        if not meta_path.exists():
            raise FileNotFoundError(f"Missing GLADE+ meta file: {meta_path}")
        self.meta = pd.read_json(meta_path, typ="series")
        self.nside = int(self.meta["nside"])
        self.nest = bool(self.meta["nest"])
        self.hp = HEALPix(nside=self.nside, order="nested" if self.nest else "ring", frame="icrs")
        self.offsets = np.load(index_dir / "hpix_offsets.npy", mmap_mode="r")
        self.ra = np.load(index_dir / "ra_deg.npy", mmap_mode="r")
        self.dec = np.load(index_dir / "dec_deg.npy", mmap_mode="r")
        self.z = np.load(index_dir / "z.npy", mmap_mode="r")
        self.w = np.load(index_dir / "w.npy", mmap_mode="r")

    def cone_stats(self, ra_deg: float, dec_deg: float, *, radius_deg: float, z_center: float | None, z_window: float) -> dict[str, float]:
        pixels = self.hp.cone_search_lonlat(ra_deg * u.deg, dec_deg * u.deg, radius_deg * u.deg)
        total = 0
        weight_sum = 0.0
        for pixel in np.asarray(pixels, dtype=int):
            start = int(self.offsets[pixel])
            stop = int(self.offsets[pixel + 1])
            if stop <= start:
                continue
            ra_chunk = self.ra[start:stop]
            dec_chunk = self.dec[start:stop]
            z_chunk = self.z[start:stop]
            w_chunk = self.w[start:stop]
            mask = angular_sep_deg(ra_deg, dec_deg, ra_chunk, dec_chunk) <= radius_deg
            if z_center is not None and np.isfinite(z_center) and z_center > 0:
                mask &= np.abs(z_chunk - z_center) <= z_window
            if not np.any(mask):
                continue
            total += int(mask.sum())
            weight_sum += float(np.sum(w_chunk[mask]))
        area_sq_deg = math.pi * radius_deg * radius_deg
        return {
            "gladeplus_neighbor_count": float(total),
            "gladeplus_weight_sum": float(weight_sum),
            "gladeplus_surface_density_sqdeg": float(total / area_sq_deg if area_sq_deg > 0 else 0.0),
        }


def download_ned_lvs(cache_dir: Path, *, url: str = NED_LVS_URL) -> Path:
    ensure_dir(cache_dir)
    out_path = cache_dir / "NEDLVS_20210922.fits"
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with out_path.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    fh.write(chunk)
    return out_path


def load_ned_lvs(path: Path) -> pd.DataFrame:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        payload: dict[str, Any] = {}
        for col in NED_LVS_COLUMNS:
            arr = np.array(data[col])
            if hasattr(arr.dtype, "byteorder") and arr.dtype.byteorder not in ("=", "|"):
                arr = arr.byteswap().view(arr.dtype.newbyteorder("="))
            payload[col] = arr
    df = pd.DataFrame(payload)
    for col in df.columns:
        df[col] = decode_object_series(df[col])
    return df


def _is_explicit_exclusion(name: str) -> bool:
    upper = str(name).upper()
    return any(pattern in upper for pattern in EXPLICIT_EXCLUSION_PATTERNS)


def build_galaxy_master(
    *,
    root_dir: Path,
    config: GalaxyBuildConfig,
    gladeplus_index_dir: Path | None = DEFAULT_GLADEPLUS_INDEX_DIR,
) -> dict[str, Path]:
    cache_dir = ensure_dir(root_dir / "data" / "cache" / "ned_lvs")
    catalog_dir = ensure_dir(root_dir / "catalogs")

    ned_path = download_ned_lvs(cache_dir)
    df = load_ned_lvs(ned_path)

    df = df.rename(
        columns={
            "objname": "galaxy_name",
            "DistMpc": "distance_mpc",
            "DistMpc_unc": "distance_mpc_unc",
            "Mstar": "stellar_mass_msun",
            "Mstar_unc": "stellar_mass_msun_unc",
            "SFR_hybrid": "sfr_hybrid_msun_per_yr",
            "SFR_hybrid_unc": "sfr_hybrid_msun_per_yr_unc",
            "SFR_W4": "sfr_w4_msun_per_yr",
            "SFR_W4_unc": "sfr_w4_msun_per_yr_unc",
            "Lum_W1": "lum_w1_lsun",
            "Lum_W1_unc": "lum_w1_lsun_unc",
        }
    )

    df["galaxy_name"] = df["galaxy_name"].astype(str).str.strip()
    df["mw_av_mag"] = 3.1 * pd.to_numeric(df["ebv"], errors="coerce")
    sfr_hybrid = pd.to_numeric(df["sfr_hybrid_msun_per_yr"], errors="coerce")
    sfr_w4 = pd.to_numeric(df["sfr_w4_msun_per_yr"], errors="coerce")
    df["sfr_proxy_msun_per_yr"] = sfr_hybrid.where(np.isfinite(sfr_hybrid) & (sfr_hybrid > 0), sfr_w4)
    df["distance_mpc"] = pd.to_numeric(df["distance_mpc"], errors="coerce")
    df["stellar_mass_msun"] = pd.to_numeric(df["stellar_mass_msun"], errors="coerce")
    df["lum_w1_lsun"] = pd.to_numeric(df["lum_w1_lsun"], errors="coerce")
    df["z"] = pd.to_numeric(df["z"], errors="coerce")
    df["ET_flag"] = df["ET_flag"].astype(bool)
    df["explicit_exclusion"] = df["galaxy_name"].map(_is_explicit_exclusion)

    df["mass_proxy"] = df["stellar_mass_msun"].where(np.isfinite(df["stellar_mass_msun"]), df["lum_w1_lsun"])
    df["score_mass"] = normalize_01(robust_log10(df["mass_proxy"], floor=1.0))
    df["score_sfr"] = normalize_01(robust_log10(df["sfr_proxy_msun_per_yr"], floor=1e-5))
    df["score_distance"] = 1.0 - normalize_01(df["distance_mpc"])
    df["score_extinction"] = 1.0 - normalize_01(df["mw_av_mag"])
    df["searchability_score"] = (
        0.40 * df["score_mass"]
        + 0.30 * df["score_sfr"]
        + 0.20 * df["score_distance"]
        + 0.10 * df["score_extinction"]
    )

    eligible_mask = (
        np.isfinite(df["distance_mpc"])
        & (df["distance_mpc"] > 0)
        & (df["distance_mpc"] <= config.max_distance_mpc)
        & np.isfinite(df["mw_av_mag"])
        & (df["mw_av_mag"] <= config.max_av_mag)
        & ~df["ET_flag"]
        & ~df["explicit_exclusion"]
    )
    eligible = df.loc[eligible_mask].copy()
    eligible = eligible.sort_values(["searchability_score", "mass_proxy", "sfr_proxy_msun_per_yr"], ascending=[False, False, False])
    eligible["rank"] = np.arange(1, len(eligible) + 1)

    glade_annotated = False
    if gladeplus_index_dir and gladeplus_index_dir.exists():
        index = GladePlusLocalIndex(gladeplus_index_dir)
        rows = []
        top_slice = eligible.head(max(config.top_n, 250)).copy()
        for row in top_slice.itertuples(index=False):
            stats = index.cone_stats(
                float(row.ra),
                float(row.dec),
                radius_deg=config.gladeplus_radius_deg,
                z_center=float(row.z) if np.isfinite(row.z) else None,
                z_window=config.gladeplus_z_window,
            )
            rows.append(stats)
        if rows:
            stats_df = pd.DataFrame(rows, index=top_slice.index)
            for col in stats_df.columns:
                eligible.loc[stats_df.index, col] = stats_df[col]
            eligible["score_environment"] = 1.0 - normalize_01(eligible["gladeplus_surface_density_sqdeg"])
            eligible["searchability_score"] = 0.92 * eligible["searchability_score"] + 0.08 * eligible["score_environment"].fillna(0.5)
            eligible = eligible.sort_values(["searchability_score", "mass_proxy", "sfr_proxy_msun_per_yr"], ascending=[False, False, False])
            eligible["rank"] = np.arange(1, len(eligible) + 1)
            glade_annotated = True

    exclusions = df.copy()
    exclusions["exclusion_reason"] = np.select(
        [
            ~np.isfinite(exclusions["distance_mpc"]) | (exclusions["distance_mpc"] <= 0),
            exclusions["distance_mpc"] > config.max_distance_mpc,
            exclusions["mw_av_mag"] > config.max_av_mag,
            exclusions["ET_flag"],
            exclusions["explicit_exclusion"],
        ],
        [
            "missing_or_bad_distance",
            "too_distant",
            "too_extincted",
            "early_type",
            "explicit_local_group_exclusion",
        ],
        default="eligible",
    )

    aliases = eligible[["galaxy_name", "ra", "dec", "z"]].copy()
    aliases["canonical_name"] = aliases["galaxy_name"]

    galaxy_master_path = catalog_dir / "galaxy_master.parquet"
    galaxy_scores_path = catalog_dir / "galaxy_scores.parquet"
    galaxy_aliases_path = catalog_dir / "galaxy_aliases.parquet"
    galaxy_exclusions_path = catalog_dir / "galaxy_exclusions.csv"
    summary_path = catalog_dir / "galaxy_selection_summary.json"

    eligible.to_parquet(galaxy_master_path, index=False)
    eligible[
        [
            "galaxy_name",
            "rank",
            "searchability_score",
            "distance_mpc",
            "mw_av_mag",
            "mass_proxy",
            "sfr_proxy_msun_per_yr",
            "score_mass",
            "score_sfr",
            "score_distance",
            "score_extinction",
        ]
        + [c for c in ["gladeplus_neighbor_count", "gladeplus_weight_sum", "gladeplus_surface_density_sqdeg", "score_environment"] if c in eligible.columns]
    ].to_parquet(galaxy_scores_path, index=False)
    aliases.to_parquet(galaxy_aliases_path, index=False)
    exclusions.loc[exclusions["exclusion_reason"] != "eligible", ["galaxy_name", "distance_mpc", "mw_av_mag", "ET_flag", "exclusion_reason"]].to_csv(
        galaxy_exclusions_path, index=False
    )
    write_json(
        summary_path,
        {
            "created_utc": utc_stamp(),
            "source_catalog": str(ned_path),
            "n_total": int(len(df)),
            "n_eligible": int(len(eligible)),
            "top_n_requested": int(config.top_n),
            "gladeplus_annotation": bool(glade_annotated),
            "gladeplus_index_dir": str(gladeplus_index_dir) if gladeplus_index_dir else None,
            "max_distance_mpc": float(config.max_distance_mpc),
            "max_av_mag": float(config.max_av_mag),
            "top_ranked": eligible.head(config.top_n)["galaxy_name"].tolist(),
        },
    )

    return {
        "galaxy_master": galaxy_master_path,
        "galaxy_scores": galaxy_scores_path,
        "galaxy_aliases": galaxy_aliases_path,
        "galaxy_exclusions": galaxy_exclusions_path,
        "summary": summary_path,
    }
