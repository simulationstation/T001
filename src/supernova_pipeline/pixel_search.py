from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import re
import threading
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote
from uuid import uuid4

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture
from photutils.detection import DAOStarFinder
from reproject import reproject_interp
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree

from .candidate_ledger import CANDIDATE_COLUMNS
from .utils import ensure_dir, json_list, utc_stamp, write_json


MAST_INVOKE_URL = "https://mast.stsci.edu/api/v0/invoke"
MAST_DOWNLOAD_URL = "https://mast.stsci.edu/api/v0.1/Download/file"
_FIGURE_LOCK = threading.Lock()

warnings.filterwarnings("ignore", category=FITSFixedWarning)


@dataclass(slots=True)
class ScienceImage:
    path: Path
    data: np.ndarray
    wcs: WCS
    pixel_scale_arcsec: float
    obs_collection: str
    obs_id: str
    filter_name: str
    instrument_name: str
    extname: str


def _safe_divide(num: np.ndarray, den: np.ndarray | float, *, default: float = np.nan) -> np.ndarray:
    num_arr = np.asarray(num, dtype=float)
    den_arr = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num_arr, den_arr).shape, default, dtype=float)
    ok = np.isfinite(num_arr) & np.isfinite(den_arr) & (np.abs(den_arr) > 0.0)
    out[ok] = num_arr[ok] / den_arr[ok]
    return out


def _slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return text.strip("_") or "item"


def _product_query(obsid: int) -> pd.DataFrame:
    payload = {
        "service": "Mast.Caom.Products",
        "params": {"obsid": str(int(obsid))},
        "format": "json",
        "pagesize": 2000,
        "page": 1,
        "removenullcolumns": True,
        "timeout": 30,
    }
    resp = requests.post(MAST_INVOKE_URL, data={"request": json.dumps(payload)}, timeout=60)
    resp.raise_for_status()
    return pd.DataFrame(resp.json().get("data", []))


def _product_priority(row: pd.Series, *, obs_collection: str) -> float:
    filename = str(row.get("productFilename", "")).lower()
    subgroup = str(row.get("productSubGroupDescription", "")).upper()
    datatype = str(row.get("dataproduct_type", "")).lower()
    description = str(row.get("description", "")).lower()
    rights = str(row.get("dataRights", "")).upper()
    if rights not in {"", "PUBLIC"}:
        return -1e9
    if not (filename.endswith(".fits") or filename.endswith(".fits.gz")):
        return -1e9
    if any(token in filename for token in ["_x1d", "_s2d", "_rateints", "_segm", "_cat", "_asn", "_trl", "_jpg", "_png"]):
        return -1e9
    score = 0.0
    score += 8.0 * float(row.get("calib_level", 0) or 0)
    if "minimum recommended" in str(row.get("productGroupDescription", "")).lower():
        score += 25.0
    if datatype == "image":
        score += 10.0
    if "fits" in description:
        score += 2.0
    if obs_collection == "HST":
        if subgroup in {"DRC", "DRZ"} or filename.endswith(("_drc.fits", "_drz.fits", "_drc.fits.gz", "_drz.fits.gz")):
            score += 120.0
        elif subgroup in {"FLC", "FLT"} or filename.endswith(("_flc.fits", "_flt.fits", "_flc.fits.gz", "_flt.fits.gz")):
            score += 80.0
        elif subgroup in {"C1M", "C0M", "CRJ", "IMA"}:
            score += 40.0
    elif obs_collection == "JWST":
        if subgroup == "I2D" or filename.endswith(("_i2d.fits", "_i2d.fits.gz")):
            score += 120.0
        elif subgroup == "CAL" or filename.endswith(("_cal.fits", "_cal.fits.gz")):
            score += 80.0
        elif subgroup == "RATE" or filename.endswith(("_rate.fits", "_rate.fits.gz")):
            score += 55.0
    return score


def _select_best_product(products: pd.DataFrame, *, obs_collection: str) -> pd.Series | None:
    if products.empty:
        return None
    work = products.copy()
    work["priority"] = [_product_priority(row, obs_collection=obs_collection) for _, row in work.iterrows()]
    work = work[work["priority"] > 0].copy()
    if work.empty:
        return None
    return work.sort_values(["priority", "size"], ascending=[False, False]).iloc[0]


def _download_product(product: pd.Series, *, root_dir: Path, obs_collection: str, obs_id: str) -> Path:
    cache_dir = ensure_dir(root_dir / "data" / "cache" / "mast_products" / obs_collection / _slugify(obs_id))
    filename = str(product["productFilename"])
    destination = cache_dir / filename
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    uri = str(product["dataURI"])
    tmp_destination = destination.with_name(f".{destination.name}.{uuid4().hex}.part")
    resp = requests.get(f"{MAST_DOWNLOAD_URL}?uri={quote(uri, safe=':/')}", stream=True, timeout=(30, 300))
    resp.raise_for_status()
    with tmp_destination.open("wb") as handle:
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    if destination.exists() and destination.stat().st_size > 0:
        tmp_destination.unlink(missing_ok=True)
        return destination
    tmp_destination.replace(destination)
    return destination


def _load_science_image(path: Path, *, obs_collection: str, obs_id: str, filter_name: str, instrument_name: str) -> ScienceImage:
    with fits.open(path, memmap=False) as hdul:
        best_score = -1.0
        best_data: np.ndarray | None = None
        best_wcs: WCS | None = None
        best_extname = "PRIMARY"
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if not isinstance(data, np.ndarray) or data.ndim != 2:
                continue
            if not np.isfinite(data).any():
                continue
            try:
                wcs = WCS(hdu.header).celestial
            except Exception:
                continue
            if not wcs.has_celestial:
                continue
            extname = (getattr(hdu, "name", "") or "PRIMARY").upper()
            score = float(np.isfinite(data).sum())
            if extname == "SCI":
                score += 1e9
            elif extname == "PRIMARY":
                score += 5e8
            if score > best_score:
                best_score = score
                best_data = np.asarray(data, dtype=float)
                best_wcs = wcs
                best_extname = extname

    if best_data is None or best_wcs is None:
        raise RuntimeError(f"No celestial 2D science image found in {path}")

    pix_scale = float(np.nanmedian(proj_plane_pixel_scales(best_wcs)) * 3600.0)
    return ScienceImage(
        path=path,
        data=best_data,
        wcs=best_wcs,
        pixel_scale_arcsec=pix_scale if math.isfinite(pix_scale) and pix_scale > 0 else 0.05,
        obs_collection=obs_collection,
        obs_id=obs_id,
        filter_name=filter_name,
        instrument_name=instrument_name,
        extname=best_extname,
    )


def _cutout_if_needed(image: ScienceImage, *, ra_deg: float, dec_deg: float, max_pixels: int = 1800, size_arcsec: float = 240.0) -> ScienceImage:
    if image.data.shape[0] <= max_pixels and image.data.shape[1] <= max_pixels:
        return image
    size_pix = max(256, min(max_pixels, int(size_arcsec / max(image.pixel_scale_arcsec, 1e-3))))
    position = SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg, frame="icrs")
    cutout = Cutout2D(image.data, position, size=(size_pix, size_pix), wcs=image.wcs, mode="trim", fill_value=np.nan)
    return ScienceImage(
        path=image.path,
        data=np.asarray(cutout.data, dtype=float),
        wcs=cutout.wcs,
        pixel_scale_arcsec=image.pixel_scale_arcsec,
        obs_collection=image.obs_collection,
        obs_id=image.obs_id,
        filter_name=image.filter_name,
        instrument_name=image.instrument_name,
        extname=image.extname,
    )


def _finite_mask(data: np.ndarray) -> np.ndarray:
    return np.isfinite(data) & np.isreal(data)


def _bbox_from_mask(mask: np.ndarray, *, padding: int = 24) -> tuple[slice, slice]:
    ys, xs = np.where(mask)
    if len(xs) == 0 or len(ys) == 0:
        return slice(0, mask.shape[0]), slice(0, mask.shape[1])
    y0 = max(int(ys.min()) - padding, 0)
    y1 = min(int(ys.max()) + padding + 1, mask.shape[0])
    x0 = max(int(xs.min()) - padding, 0)
    x1 = min(int(xs.max()) + padding + 1, mask.shape[1])
    return slice(y0, y1), slice(x0, x1)


def _background_stats(data: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    good = mask & np.isfinite(data)
    if good.sum() < 100:
        return 0.0, float(np.nanstd(np.where(good, data, np.nan)))
    _, median, std = sigma_clipped_stats(data[good], sigma=3.0, maxiters=5)
    return float(median), float(std if math.isfinite(std) and std > 0 else np.nanstd(data[good]))


def _fwhm_guess(pixel_scale_arcsec: float) -> float:
    return float(np.clip(0.18 / max(pixel_scale_arcsec, 0.02), 2.0, 4.5))


def _detect_sources(data: np.ndarray, valid_mask: np.ndarray, *, pixel_scale_arcsec: float, max_sources: int = 2500) -> pd.DataFrame:
    median, std = _background_stats(data, valid_mask)
    if not math.isfinite(std) or std <= 0:
        return pd.DataFrame(columns=["x", "y", "flux", "peak"])

    search_image = np.where(valid_mask & np.isfinite(data), data - median, 0.0)
    fwhm = _fwhm_guess(pixel_scale_arcsec)
    attempts: list[tuple[np.ndarray, float]] = [(search_image, 6.0 * std)]

    filtered = search_image - gaussian_filter(search_image, sigma=max(fwhm, 2.0))
    good = valid_mask & np.isfinite(filtered)
    if good.sum() > 0:
        _, filtered_median, filtered_std = sigma_clipped_stats(filtered[good], sigma=3.0, maxiters=5)
        filtered_std = float(filtered_std if math.isfinite(filtered_std) and filtered_std > 0 else np.nanstd(filtered[good]))
        if math.isfinite(filtered_std) and filtered_std > 0:
            attempts.append((np.where(good, filtered - filtered_median, 0.0), 4.5 * filtered_std))

    for candidate_image, threshold in attempts:
        finder = DAOStarFinder(
            fwhm=fwhm,
            threshold=threshold,
            exclude_border=True,
            sharplo=0.02,
            roundlo=-2.5,
            roundhi=2.5,
        )
        table = finder(candidate_image)
        if table is None or len(table) == 0:
            continue
        df = table.to_pandas().rename(columns={"xcentroid": "x", "ycentroid": "y"})
        df = df[np.isfinite(df["x"]) & np.isfinite(df["y"]) & np.isfinite(df["flux"])].copy()
        if df.empty:
            continue
        xi = np.clip(np.round(df["x"]).astype(int), 0, valid_mask.shape[1] - 1)
        yi = np.clip(np.round(df["y"]).astype(int), 0, valid_mask.shape[0] - 1)
        df = df[valid_mask[yi, xi]].copy()
        if len(df) > max_sources:
            df = df.sort_values("flux", ascending=False).head(max_sources).copy()
        return df.reset_index(drop=True)

    return pd.DataFrame(columns=["x", "y", "flux", "peak"])


def _aperture_fluxes(
    data: np.ndarray,
    valid_mask: np.ndarray,
    positions: np.ndarray,
    *,
    radius: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aperture = CircularAperture(positions, r=radius)
    annulus = CircularAnnulus(positions, r_in=radius * 2.0, r_out=radius * 3.0)
    ap_stats = ApertureStats(data, aperture, mask=~valid_mask)
    ann_stats = ApertureStats(data, annulus, mask=~valid_mask)
    area = np.asarray(getattr(ap_stats.sum_aper_area, "value", ap_stats.sum_aper_area), dtype=float)
    bg = np.asarray(ann_stats.median, dtype=float)
    bg_std = np.asarray(ann_stats.std, dtype=float)
    flux = np.asarray(ap_stats.sum, dtype=float) - bg * area
    err = np.abs(bg_std) * np.sqrt(np.clip(area, 1.0, None))
    return flux, err, bg


def _estimate_scale_factor(pre_data: np.ndarray, post_data: np.ndarray, valid_mask: np.ndarray, positions: np.ndarray) -> tuple[float, int]:
    if len(positions) == 0:
        return 1.0, 0
    flux_pre, err_pre, _ = _aperture_fluxes(pre_data, valid_mask, positions, radius=3.0)
    flux_post, err_post, _ = _aperture_fluxes(post_data, valid_mask, positions, radius=3.0)
    snr_pre = _safe_divide(np.abs(flux_pre), err_pre, default=0.0)
    snr_post = _safe_divide(np.abs(flux_post), err_post, default=0.0)
    good = (
        np.isfinite(flux_pre)
        & np.isfinite(flux_post)
        & (flux_pre > 0)
        & (flux_post > 0)
        & (snr_pre >= 8.0)
        & (snr_post >= 5.0)
    )
    if good.sum() < 8:
        return 1.0, int(good.sum())
    ratios = flux_post[good] / flux_pre[good]
    lo, hi = np.nanpercentile(ratios, [20, 80])
    keep = ratios[(ratios >= lo) & (ratios <= hi)]
    if len(keep) < 5:
        keep = ratios
    scale = float(np.nanmedian(keep))
    if not math.isfinite(scale) or scale <= 0:
        scale = 1.0
    return float(np.clip(scale, 0.05, 20.0)), int(good.sum())


def _astrometric_residual(pre_data: np.ndarray, post_data: np.ndarray, valid_mask: np.ndarray, *, pixel_scale_arcsec: float) -> float:
    pre_src = _detect_sources(pre_data, valid_mask, pixel_scale_arcsec=pixel_scale_arcsec, max_sources=300)
    post_src = _detect_sources(post_data, valid_mask, pixel_scale_arcsec=pixel_scale_arcsec, max_sources=300)
    if pre_src.empty or post_src.empty:
        return float("nan")
    tree = cKDTree(post_src[["x", "y"]].to_numpy())
    dist, _ = tree.query(pre_src[["x", "y"]].to_numpy(), distance_upper_bound=2.5)
    good = np.isfinite(dist) & (dist < 2.5)
    if good.sum() < 5:
        return float("nan")
    return float(np.nanmedian(dist[good]))


def _priority_score(
    *,
    s_min_loo: float,
    fade_fraction: float,
    depth_margin_post: float,
    crowding_index: float,
    blend_risk: float,
    host_bg_penalty: float,
    agreement: float,
) -> float:
    penalty = 1.0 + 0.2 * crowding_index + 1.5 * blend_risk + 8.0 * max(host_bg_penalty, 0.0)
    depth_term = float(np.clip(depth_margin_post / 3.0, 0.2, 2.0))
    agreement_term = float(np.clip(0.5 + 0.5 * agreement, 0.2, 1.5))
    return float(max(s_min_loo, 0.0) * max(fade_fraction, 0.0) * depth_term * agreement_term / penalty)


def _pair_stability_metrics(source_catalog: pd.DataFrame) -> tuple[float, float, bool]:
    usable = source_catalog[
        np.isfinite(source_catalog["flux_ratio"])
        & np.isfinite(source_catalog["flux_pre"])
        & (source_catalog["flux_pre"] > 0)
    ].copy()
    if "edge_distance_px" in usable.columns:
        usable = usable[usable["edge_distance_px"].fillna(0.0) >= 6.0].copy()
    if usable.empty:
        return float("nan"), float("nan"), False
    stable_fraction = float(((usable["flux_ratio"] >= 0.8) & (usable["flux_ratio"] <= 1.2)).mean())
    extreme_fade_fraction = float((usable["flux_ratio"] <= 0.05).mean())
    pair_is_stable = stable_fraction >= 0.60 and extreme_fade_fraction <= 0.02
    return stable_fraction, extreme_fade_fraction, pair_is_stable


def _stretch_limits(data: np.ndarray) -> tuple[float, float]:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return -1.0, 1.0
    lo, hi = np.nanpercentile(finite, [2, 98])
    if math.isclose(lo, hi):
        hi = lo + 1.0
    return float(lo), float(hi)


def _save_candidate_figure(
    packet_dir: Path,
    candidate_id: str,
    pre_cutout: np.ndarray,
    post_cutout: np.ndarray,
    diff_cutout: np.ndarray,
    *,
    title: str,
) -> Path:
    output = packet_dir / f"{candidate_id}_cutout.png"
    with _FIGURE_LOCK:
        fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), constrained_layout=True)
        for ax, image, label in zip(axes, [pre_cutout, post_cutout, diff_cutout], ["Pre", "Post", "Pre - Post"], strict=True):
            lo, hi = _stretch_limits(image)
            ax.imshow(image, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            ax.set_title(label)
            ax.set_xticks([])
            ax.set_yticks([])
            cy = image.shape[0] / 2.0
            cx = image.shape[1] / 2.0
            ax.plot(cx, cy, marker="+", color="tab:red", ms=10, mew=1.5)
        fig.suptitle(title, fontsize=10)
        fig.savefig(output, dpi=180)
        plt.close(fig)
    return output


def _build_candidate_packet(
    *,
    root_dir: Path,
    candidate_id: str,
    pair: pd.Series,
    world: SkyCoord,
    pre_crop: np.ndarray,
    post_crop: np.ndarray,
    diff_crop: np.ndarray,
    x_pix: float,
    y_pix: float,
    cutout_half_size: int = 18,
    payload: dict[str, Any],
) -> Path:
    packet_dir = ensure_dir(root_dir / "candidate_packets" / candidate_id)
    x0 = max(int(round(x_pix)) - cutout_half_size, 0)
    x1 = min(int(round(x_pix)) + cutout_half_size + 1, pre_crop.shape[1])
    y0 = max(int(round(y_pix)) - cutout_half_size, 0)
    y1 = min(int(round(y_pix)) + cutout_half_size + 1, pre_crop.shape[0])
    title = f"{pair['galaxy_name']}  {pair['filter_1']}->{pair['filter_2']}"
    figure_path = _save_candidate_figure(
        packet_dir,
        candidate_id,
        pre_crop[y0:y1, x0:x1],
        post_crop[y0:y1, x0:x1],
        diff_crop[y0:y1, x0:x1],
        title=title,
    )
    summary_path = packet_dir / "summary.json"
    write_json(
        summary_path,
        {
            **payload,
            "candidate_id": candidate_id,
            "galaxy_name": str(pair["galaxy_name"]),
            "ra_deg": float(world.ra.deg),
            "dec_deg": float(world.dec.deg),
            "pre_obs_id": str(pair["obs_id_1"]),
            "post_obs_id": str(pair["obs_id_2"]),
            "pre_filter": str(pair["filter_1"]),
            "post_filter": str(pair["filter_2"]),
            "baseline_days": float(pair["baseline_days"]),
            "figure_path": str(figure_path),
            "created_utc": utc_stamp(),
        },
    )
    return figure_path


def _status_for_candidate(
    *,
    flux_ratio: float,
    s_raw: float,
    s_min_loo: float,
    depth_margin_post: float,
    blend_risk: float,
    host_bg_penalty: float,
    edge_distance_px: float,
) -> str | None:
    if (
        math.isfinite(flux_ratio)
        and flux_ratio <= 0.05
        and flux_ratio > -0.10
        and s_min_loo >= 5.0
        and depth_margin_post >= 2.0
        and blend_risk < 0.6
        and host_bg_penalty < 0.35
        and edge_distance_px >= 10.0
    ):
        return "PASS"
    if (
        math.isfinite(flux_ratio)
        and flux_ratio <= 0.50
        and s_raw >= 5.0
        and depth_margin_post >= 1.0
        and edge_distance_px >= 6.0
    ):
        return "REVIEW"
    return None


def _resolve_pair_products(pair: pd.Series, *, root_dir: Path) -> tuple[dict[str, Any], ScienceImage, ScienceImage]:
    products_1 = _product_query(int(pair["obsid_1"]))
    products_2 = _product_query(int(pair["obsid_2"]))
    product_1 = _select_best_product(products_1, obs_collection=str(pair["obs_collection_1"]))
    product_2 = _select_best_product(products_2, obs_collection=str(pair["obs_collection_2"]))
    if product_1 is None or product_2 is None:
        raise RuntimeError("No suitable calibrated FITS product found for one or both observations.")

    path_1 = _download_product(product_1, root_dir=root_dir, obs_collection=str(pair["obs_collection_1"]), obs_id=str(pair["obs_id_1"]))
    path_2 = _download_product(product_2, root_dir=root_dir, obs_collection=str(pair["obs_collection_2"]), obs_id=str(pair["obs_id_2"]))
    image_1 = _load_science_image(
        path_1,
        obs_collection=str(pair["obs_collection_1"]),
        obs_id=str(pair["obs_id_1"]),
        filter_name=str(pair["filter_1"]),
        instrument_name=str(pair["instrument_1"]),
    )
    image_2 = _load_science_image(
        path_2,
        obs_collection=str(pair["obs_collection_2"]),
        obs_id=str(pair["obs_id_2"]),
        filter_name=str(pair["filter_2"]),
        instrument_name=str(pair["instrument_2"]),
    )
    product_info = {
        "obsid_1": int(pair["obsid_1"]),
        "obsid_2": int(pair["obsid_2"]),
        "obs_id_1": str(pair["obs_id_1"]),
        "obs_id_2": str(pair["obs_id_2"]),
        "product_filename_1": str(product_1["productFilename"]),
        "product_filename_2": str(product_2["productFilename"]),
        "product_uri_1": str(product_1["dataURI"]),
        "product_uri_2": str(product_2["dataURI"]),
        "product_path_1": str(path_1),
        "product_path_2": str(path_2),
        "product_priority_1": float(product_1["priority"]),
        "product_priority_2": float(product_2["priority"]),
        "image_ext_1": image_1.extname,
        "image_ext_2": image_2.extname,
    }
    return product_info, image_1, image_2


def _scan_pair(pair: pd.Series, obs_meta: dict[str, dict[str, Any]], *, root_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]], pd.DataFrame]:
    meta_1 = obs_meta[str(pair["obs_id_1"])]
    meta_2 = obs_meta[str(pair["obs_id_2"])]
    product_info, image_1, image_2 = _resolve_pair_products(pair, root_dir=root_dir)

    galaxy_ra = float(meta_1["galaxy_ra"])
    galaxy_dec = float(meta_1["galaxy_dec"])
    search_radius_deg = max(float(meta_1.get("search_radius_deg", 0.0) or 0.0), float(meta_2.get("search_radius_deg", 0.0) or 0.0))
    cutout_arcsec = max(240.0, 2.2 * search_radius_deg * 3600.0)
    image_1 = _cutout_if_needed(
        image_1,
        ra_deg=galaxy_ra,
        dec_deg=galaxy_dec,
        size_arcsec=cutout_arcsec,
        max_pixels=6000,
    )
    image_2 = _cutout_if_needed(
        image_2,
        ra_deg=galaxy_ra,
        dec_deg=galaxy_dec,
        size_arcsec=cutout_arcsec,
        max_pixels=6000,
    )

    post_reproj, footprint = reproject_interp((image_2.data, image_2.wcs), image_1.wcs, shape_out=image_1.data.shape)
    overlap_mask = _finite_mask(image_1.data) & _finite_mask(post_reproj) & (footprint > 0)
    if overlap_mask.sum() < 10_000:
        raise RuntimeError("Insufficient image overlap after reprojection.")

    ys, xs = _bbox_from_mask(overlap_mask)
    pre_crop = np.asarray(image_1.data[ys, xs], dtype=float)
    post_crop = np.asarray(post_reproj[ys, xs], dtype=float)
    valid_crop = np.asarray(overlap_mask[ys, xs], dtype=bool)

    pre_bg, pre_std = _background_stats(pre_crop, valid_crop)
    post_bg, post_std = _background_stats(post_crop, valid_crop)
    pre_sub = pre_crop - pre_bg
    post_sub = post_crop - post_bg
    edge_distance_map = distance_transform_edt(valid_crop)

    seed_sources = _detect_sources(pre_sub, valid_crop, pixel_scale_arcsec=image_1.pixel_scale_arcsec, max_sources=1500)
    if seed_sources.empty:
        raise RuntimeError("No point-like sources detected in the pre image.")

    positions = seed_sources[["x", "y"]].to_numpy()
    scale_factor, n_scale_sources = _estimate_scale_factor(pre_sub, post_sub, valid_crop, positions)
    post_scaled = post_sub / scale_factor
    diff_image = pre_sub - post_scaled

    astrom_resid = _astrometric_residual(pre_sub, post_scaled, valid_crop, pixel_scale_arcsec=image_1.pixel_scale_arcsec)

    radii = [2.5, 3.5, 4.5]
    measures: list[pd.DataFrame] = []
    for radius in radii:
        flux_pre, err_pre, bg_pre = _aperture_fluxes(pre_sub, valid_crop, positions, radius=radius)
        flux_post, err_post, bg_post = _aperture_fluxes(post_scaled, valid_crop, positions, radius=radius)
        diff = flux_pre - flux_post
        diff_err = np.sqrt(np.square(err_pre) + np.square(err_post))
        sigma = _safe_divide(diff, diff_err, default=np.nan)
        ratio = _safe_divide(flux_post, flux_pre, default=np.nan)
        fade_fraction = 1.0 - ratio
        frame = pd.DataFrame(
            {
                "radius": radius,
                "flux_pre": flux_pre,
                "flux_post": flux_post,
                "err_pre": err_pre,
                "err_post": err_post,
                "bg_pre": bg_pre,
                "bg_post": bg_post,
                "sigma": sigma,
                "ratio": ratio,
                "fade_fraction": fade_fraction,
            }
        )
        measures.append(frame)

    metrics = measures[1].copy()
    metrics["sigma_min"] = np.nanmin(np.vstack([m["sigma"].to_numpy() for m in measures]), axis=0)
    metrics["sigma_max"] = np.nanmax(np.vstack([m["sigma"].to_numpy() for m in measures]), axis=0)
    metrics["fade_median"] = np.nanmedian(np.vstack([m["fade_fraction"].to_numpy() for m in measures]), axis=0)
    metrics["fade_std"] = np.nanstd(np.vstack([m["fade_fraction"].to_numpy() for m in measures]), axis=0)
    metrics["ratio_default"] = metrics["ratio"]

    src_positions = seed_sources[["x", "y"]].to_numpy()
    tree = cKDTree(src_positions)
    dist, _ = tree.query(src_positions, k=2)
    nearest = dist[:, 1]
    crowding_index = np.array([len(tree.query_ball_point(p, r=10.0)) - 1 for p in src_positions], dtype=float)
    blend_risk = np.clip(1.0 - nearest / 10.0, 0.0, 1.0)
    agreement = np.clip(1.0 - _safe_divide(metrics["fade_std"].to_numpy(), np.abs(metrics["fade_median"].to_numpy()), default=1.0), 0.0, 1.0)
    depth_margin_post = _safe_divide(metrics["flux_pre"].to_numpy(), 5.0 * np.clip(metrics["err_post"].to_numpy(), 1e-12, None), default=0.0)
    host_bg_penalty = _safe_divide(np.abs(metrics["bg_post"].to_numpy()), np.abs(metrics["flux_pre"].to_numpy()), default=1.0)

    source_catalog = seed_sources.copy()
    source_catalog["flux_pre"] = metrics["flux_pre"].to_numpy()
    source_catalog["flux_post"] = metrics["flux_post"].to_numpy()
    source_catalog["fade_sigma"] = metrics["sigma"].to_numpy()
    source_catalog["fade_sigma_min"] = metrics["sigma_min"].to_numpy()
    source_catalog["fade_sigma_max"] = metrics["sigma_max"].to_numpy()
    source_catalog["fade_fraction"] = metrics["fade_median"].to_numpy()
    source_catalog["flux_ratio"] = metrics["ratio_default"].to_numpy()
    source_catalog["crowding_index"] = crowding_index
    source_catalog["blend_risk"] = blend_risk
    source_catalog["depth_margin_post"] = depth_margin_post
    source_catalog["host_bg_penalty"] = host_bg_penalty
    source_catalog["agreement"] = agreement
    xi = np.clip(np.round(source_catalog["x"]).astype(int), 0, edge_distance_map.shape[1] - 1)
    yi = np.clip(np.round(source_catalog["y"]).astype(int), 0, edge_distance_map.shape[0] - 1)
    source_catalog["edge_distance_px"] = edge_distance_map[yi, xi]
    source_catalog["pair_id"] = str(pair["pair_id"])

    stable_fraction, extreme_fade_fraction, pair_is_stable = _pair_stability_metrics(source_catalog)

    candidates: list[dict[str, Any]] = []
    for idx, source in source_catalog.iterrows():
        if not pair_is_stable:
            continue
        status = _status_for_candidate(
            flux_ratio=float(source["flux_ratio"]) if math.isfinite(source["flux_ratio"]) else float("nan"),
            s_raw=float(source["fade_sigma"]),
            s_min_loo=float(source["fade_sigma_min"]),
            depth_margin_post=float(source["depth_margin_post"]),
            blend_risk=float(source["blend_risk"]),
            host_bg_penalty=float(source["host_bg_penalty"]),
            edge_distance_px=float(source["edge_distance_px"]),
        )
        if status is None:
            continue

        world = image_1.wcs.pixel_to_world(float(source["x"] + xs.start), float(source["y"] + ys.start))
        fade_fraction = float(source["fade_fraction"])
        flux_ratio = float(source["flux_ratio"]) if math.isfinite(source["flux_ratio"]) else float("nan")
        s_raw = float(source["fade_sigma"])
        s_min_loo = float(source["fade_sigma_min"])
        agreement_val = float(source["agreement"])
        crowding_val = float(source["crowding_index"])
        blend_val = float(source["blend_risk"])
        host_bg_val = float(source["host_bg_penalty"])
        depth_val = float(source["depth_margin_post"])
        edge_distance_val = float(source["edge_distance_px"])
        priority = _priority_score(
            s_min_loo=s_min_loo,
            fade_fraction=fade_fraction,
            depth_margin_post=depth_val,
            crowding_index=crowding_val,
            blend_risk=blend_val,
            host_bg_penalty=host_bg_val,
            agreement=agreement_val,
        )

        reason_codes: list[str] = []
        warning_flags: list[str] = []
        if math.isfinite(flux_ratio) and flux_ratio <= 0.05:
            reason_codes.append("POST_LT_5PCT")
        if s_min_loo >= 5.0:
            reason_codes.append("ROBUST_FADE_SIG")
        if depth_val >= 2.0:
            reason_codes.append("POST_DEPTH_OK")
        if crowding_val >= 4.0:
            warning_flags.append("CROWDED")
        if blend_val >= 0.6:
            warning_flags.append("BLEND_RISK")
        if host_bg_val >= 0.35:
            warning_flags.append("HOST_BACKGROUND")
        if edge_distance_val < 10.0:
            warning_flags.append("EDGE_NEAR_MASK")
        if math.isfinite(astrom_resid) and astrom_resid >= 1.5:
            warning_flags.append("ASTROMETRY_RESIDUAL")
        if n_scale_sources < 15:
            warning_flags.append("LOW_SCALE_ANCHORS")

        candidate_id = _slugify(f"{pair['galaxy_name']}_{pair['obs_id_1']}_{pair['obs_id_2']}_{idx:04d}")
        figure_path = _build_candidate_packet(
            root_dir=root_dir,
            candidate_id=candidate_id,
            pair=pair,
            world=world,
            pre_crop=pre_sub,
            post_crop=post_scaled,
            diff_crop=diff_image,
            x_pix=float(source["x"]),
            y_pix=float(source["y"]),
            payload={
                "status": status,
                "priority_score": priority,
                "fade_sigma": s_raw,
                "fade_fraction": fade_fraction,
                "flux_ratio": flux_ratio,
                "S_min_LOO": s_min_loo,
                "crowding_index": crowding_val,
                "blend_risk": blend_val,
                "host_bg_penalty": host_bg_val,
                "edge_distance_px": edge_distance_val,
                "depth_margin_post": depth_val,
                "astrometric_residual": astrom_resid,
                "cross_reduction_agreement": agreement_val,
                "product_info": product_info,
                "warning_flags": warning_flags,
                "reason_codes": reason_codes,
                "pre_mjd": float(meta_1["t_min"]),
                "post_mjd": float(meta_2["t_min"]),
            },
        )

        candidates.append(
            {
                "candidate_id": candidate_id,
                "pair_id": str(pair["pair_id"]),
                "galaxy_name": str(pair["galaxy_name"]),
                "ra_deg": float(world.ra.deg),
                "dec_deg": float(world.dec.deg),
                "status": status,
                "priority_score": priority,
                "fade_sigma": s_raw,
                "fade_fraction": fade_fraction,
                "S_raw": s_raw,
                "S_min_LOO": s_min_loo,
                "S_robust": s_min_loo,
                "d_max": float(source_catalog.loc[idx, "fade_sigma_max"] - source_catalog.loc[idx, "fade_sigma_min"]),
                "depth_margin_post": depth_val,
                "crowding_index": crowding_val,
                "blend_risk": blend_val,
                "host_bg_penalty": host_bg_val,
                "astrometric_residual": float(astrom_resid) if math.isfinite(astrom_resid) else np.nan,
                "cross_reduction_agreement": agreement_val,
                "pre_obsid": int(pair["obsid_1"]),
                "post_obsid": int(pair["obsid_2"]),
                "pre_obs_id": str(pair["obs_id_1"]),
                "post_obs_id": str(pair["obs_id_2"]),
                "pre_filter": str(pair["filter_1"]),
                "post_filter": str(pair["filter_2"]),
                "pre_mjd": float(meta_1["t_min"]),
                "post_mjd": float(meta_2["t_min"]),
                "cutout_path": str(figure_path),
                "reason_codes_json": json_list(reason_codes),
                "warning_flags_json": json_list(warning_flags),
                "provenance_json": json.dumps(
                    {
                        "pair_score": float(pair["pair_score"]),
                        "compatibility": str(pair["compatibility"]),
                        "baseline_days": float(pair["baseline_days"]),
                        "product_info": product_info,
                        "scale_factor_post_to_pre": float(scale_factor),
                        "n_scale_sources": int(n_scale_sources),
                        "edge_distance_px": edge_distance_val,
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
        "scale_factor_post_to_pre": float(scale_factor),
        "n_scale_sources": int(n_scale_sources),
        "astrometric_residual_px": float(astrom_resid) if math.isfinite(astrom_resid) else np.nan,
        "n_sources_detected": int(len(seed_sources)),
        "pair_stable_fraction": stable_fraction,
        "pair_extreme_fade_fraction": extreme_fade_fraction,
        "pair_is_stable": pair_is_stable,
        "n_candidates": int(len(candidates)),
        "n_pass": int(sum(c["status"] == "PASS" for c in candidates)),
        "n_review": int(sum(c["status"] == "REVIEW" for c in candidates)),
        **product_info,
    }
    return pair_summary, candidates, source_catalog


def _select_pixel_pairs(
    pairs: pd.DataFrame,
    *,
    max_pairs: int,
    per_galaxy: int,
    compatibilities: set[str],
    include_cross_collection: bool,
) -> pd.DataFrame:
    work = pairs.copy()
    work = work[work["public_both"].fillna(False)].copy()
    work = work[work["compatibility"].astype(str).isin(compatibilities)].copy()
    work = work[~work["obs_collection_1"].astype(str).eq("HLSP") & ~work["obs_collection_2"].astype(str).eq("HLSP")].copy()
    work = work[work["obsid_1"].notna() & work["obsid_2"].notna()].copy()
    same_collection = work["obs_collection_1"].astype(str).eq(work["obs_collection_2"].astype(str))
    if not include_cross_collection:
        work = work[same_collection].copy()
    reducible_1 = work["obs_collection_1"].astype(str).eq("JWST") | work["obs_id_1"].astype(str).str.contains("skycell|wfc3|acs|jw", case=False, na=False)
    reducible_2 = work["obs_collection_2"].astype(str).eq("JWST") | work["obs_id_2"].astype(str).str.contains("skycell|wfc3|acs|jw", case=False, na=False)
    imaging_like_1 = ~work["instrument_1"].astype(str).str.contains("NIRSPEC", case=False, na=False) & ~work["filter_1"].astype(str).str.contains(";", regex=False, na=False)
    imaging_like_2 = ~work["instrument_2"].astype(str).str.contains("NIRSPEC", case=False, na=False) & ~work["filter_2"].astype(str).str.contains(";", regex=False, na=False)
    work = work[reducible_1 & reducible_2 & imaging_like_1 & imaging_like_2].copy()
    same_instrument = work["instrument_1"].astype(str).eq(work["instrument_2"].astype(str))
    work["pixel_priority"] = (
        work["pair_score"].astype(float)
        + 0.25 * same_collection.astype(float)
        + 0.10 * work["compatibility"].astype(str).eq("exact").astype(float)
        + 0.35 * same_instrument.astype(float)
        + 0.20 * work["obs_collection_1"].astype(str).eq("JWST").astype(float)
        + 0.20 * work["obs_collection_2"].astype(str).eq("JWST").astype(float)
        + 0.05 * np.clip(np.log10(work["baseline_days"].astype(float) + 1.0), 0.0, 5.0)
    )
    work = work.sort_values(["pixel_priority", "baseline_days"], ascending=[False, False]).copy()
    if per_galaxy > 0:
        work = work.groupby("galaxy_name", group_keys=False).head(per_galaxy).copy()
    if max_pairs > 0:
        work = work.head(max_pairs).copy()
    return work.reset_index(drop=True)


def run_pixel_search(
    *,
    root_dir: Path,
    epoch_pairs_path: Path,
    observation_matrix_path: Path,
    max_pairs: int = 12,
    per_galaxy: int = 2,
    compatibilities: tuple[str, ...] = ("exact", "very_similar"),
    include_cross_collection: bool = True,
    max_workers: int = 1,
) -> dict[str, Path]:
    archive_dir = ensure_dir(root_dir / "archive")
    candidate_dir = ensure_dir(root_dir / "candidates")

    pairs = pd.read_parquet(epoch_pairs_path)
    observations = pd.read_parquet(observation_matrix_path)
    obs_meta = observations.drop_duplicates("obs_id").set_index("obs_id").to_dict(orient="index")

    selected = _select_pixel_pairs(
        pairs,
        max_pairs=max_pairs,
        per_galaxy=per_galaxy,
        compatibilities=set(compatibilities),
        include_cross_collection=include_cross_collection,
    )

    queue_path = archive_dir / "pixel_pair_queue.parquet"
    selected.to_parquet(queue_path, index=False)

    pair_rows: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    source_frames: list[pd.DataFrame] = []

    def _failed_row(pair_row: pd.Series, exc: Exception) -> dict[str, Any]:
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
            "n_candidates": 0,
            "n_pass": 0,
            "n_review": 0,
        }

    workers = max(int(max_workers), 1)
    if workers == 1:
        iterator = selected.iterrows()
        for _, pair in iterator:
            try:
                pair_summary, pair_candidates, source_catalog = _scan_pair(pair, obs_meta, root_dir=root_dir)
                pair_summary["status"] = "SCANNED"
                pair_summary["error"] = None
                pair_rows.append(pair_summary)
                candidates.extend(pair_candidates)
                source_frames.append(source_catalog)
            except Exception as exc:
                pair_rows.append(_failed_row(pair, exc))
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_scan_pair, pair.copy(), obs_meta, root_dir=root_dir): pair.copy()
                for _, pair in selected.iterrows()
            }
            for future in as_completed(futures):
                pair = futures[future]
                try:
                    pair_summary, pair_candidates, source_catalog = future.result()
                    pair_summary["status"] = "SCANNED"
                    pair_summary["error"] = None
                    pair_rows.append(pair_summary)
                    candidates.extend(pair_candidates)
                    source_frames.append(source_catalog)
                except Exception as exc:
                    pair_rows.append(_failed_row(pair, exc))

    pair_summary_df = pd.DataFrame(pair_rows).sort_values(["status", "n_candidates"], ascending=[True, False])
    pair_summary_path = candidate_dir / "pixel_pair_summary.parquet"
    pair_summary_csv = candidate_dir / "pixel_pair_summary.csv"
    pair_summary_df.to_parquet(pair_summary_path, index=False)
    pair_summary_df.to_csv(pair_summary_csv, index=False)

    source_catalog_df = pd.concat(source_frames, ignore_index=True) if source_frames else pd.DataFrame()
    source_path = candidate_dir / "source_measurements.parquet"
    source_catalog_df.to_parquet(source_path, index=False)

    ledger = pd.DataFrame(candidates)
    if ledger.empty:
        ledger = pd.DataFrame(columns=CANDIDATE_COLUMNS)
    else:
        ledger = ledger.reindex(columns=CANDIDATE_COLUMNS)
        ledger = ledger.sort_values(["status", "priority_score"], ascending=[True, False]).reset_index(drop=True)

    ledger_path = candidate_dir / "candidate_ledger.parquet"
    ledger_csv = candidate_dir / "candidate_ledger.csv"
    ledger.to_parquet(ledger_path, index=False)
    ledger.to_csv(ledger_csv, index=False)

    review = ledger[ledger["status"].astype(str).eq("REVIEW")].copy() if not ledger.empty else pd.DataFrame(columns=ledger.columns)
    review_queue = review[["candidate_id", "galaxy_name", "priority_score", "status", "cutout_path"]].copy() if not review.empty else pd.DataFrame(columns=["candidate_id", "galaxy_name", "priority_score", "status", "cutout_path"])
    review_queue_path = candidate_dir / "review_queue.csv"
    review_queue.to_csv(review_queue_path, index=False)

    summary = {
        "created_utc": utc_stamp(),
        "n_pairs_selected": int(len(selected)),
        "n_pairs_scanned": int((pair_summary_df["status"] == "SCANNED").sum()) if not pair_summary_df.empty else 0,
        "n_pairs_failed": int((pair_summary_df["status"] == "FAILED").sum()) if not pair_summary_df.empty else 0,
        "n_sources_measured": int(len(source_catalog_df)),
        "n_candidates_total": int(len(ledger)),
        "n_pass": int((ledger["status"] == "PASS").sum()) if not ledger.empty else 0,
        "n_review": int((ledger["status"] == "REVIEW").sum()) if not ledger.empty else 0,
        "compatibilities": list(compatibilities),
        "include_cross_collection": bool(include_cross_collection),
        "max_workers": int(workers),
    }
    summary_path = candidate_dir / "candidate_summary.json"
    write_json(summary_path, summary)

    return {
        "pixel_pair_queue": queue_path,
        "pair_summary": pair_summary_path,
        "pair_summary_csv": pair_summary_csv,
        "source_measurements": source_path,
        "candidate_ledger": ledger_path,
        "candidate_ledger_csv": ledger_csv,
        "review_queue": review_queue_path,
        "summary": summary_path,
    }
