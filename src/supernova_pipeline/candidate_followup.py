from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import threading
import warnings
from pathlib import Path
from typing import Any

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import log as astropy_log
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import FITSFixedWarning
from scipy.optimize import minimize

from .archive_matrix import infer_filter_central_um
from .pixel_search import (
    ScienceImage,
    _aperture_fluxes,
    _background_stats,
    _download_product,
    _fwhm_guess,
    _load_science_image,
    _product_query,
    _select_best_product,
)
from .utils import ensure_dir, utc_stamp, write_json


HST_C_ANG_PER_S = 2.99792458e18
_PLOT_LOCK = threading.Lock()
warnings.filterwarnings("ignore", category=FITSFixedWarning)
astropy_log.setLevel("ERROR")
LIGHTCURVE_COLS = [
    "candidate_id",
    "galaxy_name",
    "obsid",
    "obs_id",
    "obs_collection",
    "instrument_name",
    "filter_name",
    "filter_central_um",
    "t_min",
    "t_max",
    "stage",
    "native_flux",
    "native_err",
    "flux_jy",
    "err_jy",
    "snr",
    "pixel_scale_arcsec",
    "radius_pix",
    "x_pix",
    "y_pix",
    "product_path",
    "status",
    "error",
]


def extinction_curve_ccm89(wavelength_um: float, rv: float = 3.1) -> float:
    x = 1.0 / max(float(wavelength_um), 1e-6)
    if x < 0.3:
        a = 0.0
        b = 0.0
    elif x < 1.1:
        a = 0.574 * x**1.61
        b = -0.527 * x**1.61
    elif x < 3.3:
        y = x - 1.82
        a = (
            1
            + 0.17699 * y
            - 0.50447 * y**2
            - 0.02427 * y**3
            + 0.72085 * y**4
            + 0.01979 * y**5
            - 0.77530 * y**6
            + 0.32999 * y**7
        )
        b = (
            1.41338 * y
            + 2.28305 * y**2
            + 1.07233 * y**3
            - 5.38434 * y**4
            - 0.62251 * y**5
            + 5.30260 * y**6
            - 2.09002 * y**7
        )
    elif x < 8.0:
        if x < 5.9:
            fa = 0.0
            fb = 0.0
        else:
            fa = -0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3
            fb = 0.2130 * (x - 5.9) ** 2 + 0.1207 * (x - 5.9) ** 3
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) + fa
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + fb
    else:
        x2 = x - 8.0
        a = -1.073 - 0.628 * x2 + 0.137 * x2**2 - 0.070 * x2**3
        b = 13.670 + 4.257 * x2 - 0.420 * x2**2 + 0.374 * x2**3
    return float(a + b / rv)


def _planck_fnu_jy(wavelength_um: np.ndarray, temperature_k: float) -> np.ndarray:
    h = 6.62607015e-34
    c = 2.99792458e8
    k = 1.380649e-23
    lam_m = np.asarray(wavelength_um, dtype=float) * 1e-6
    expo = h * c / np.clip(lam_m * k * temperature_k, 1e-30, None)
    expo = np.clip(expo, None, 700.0)
    b_lambda = (2.0 * h * c**2) / np.clip(lam_m**5, 1e-99, None) / np.expm1(expo)
    f_lambda = np.pi * b_lambda
    f_nu = f_lambda * lam_m**2 / c
    return f_nu / 1e-26


def _read_calibration_header(path: Path, extname: str) -> fits.Header:
    with fits.open(path, memmap=False) as hdul:
        header = fits.Header()
        header.update(hdul[0].header)
        for hdu in hdul:
            if (getattr(hdu, "name", "") or "PRIMARY").upper() == str(extname).upper():
                header.update(hdu.header)
                break
    return header


def _candidate_cutout(image: ScienceImage, world: SkyCoord, *, size_arcsec: float = 60.0, min_pixels: int = 96) -> ScienceImage:
    size_pix = max(min_pixels, int(size_arcsec / max(image.pixel_scale_arcsec, 1e-3)))
    cutout = Cutout2D(image.data, world, size=(size_pix, size_pix), wcs=image.wcs, mode="partial", fill_value=np.nan)
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


def _jy_from_native(
    *,
    native_flux: float,
    native_err: float,
    image: ScienceImage,
    header: fits.Header,
) -> tuple[float, float]:
    bunit = str(header.get("BUNIT", "")).strip().upper()
    if image.obs_collection == "HST":
        photflam = header.get("PHOTFLAM")
        photplam = header.get("PHOTPLAM")
        if photflam and photplam and float(photflam) > 0 and float(photplam) > 0:
            flam = native_flux * float(photflam)
            flam_err = abs(native_err) * float(photflam)
            lam_ang = float(photplam)
            flux_jy = flam * lam_ang**2 / HST_C_ANG_PER_S * 1e23
            err_jy = flam_err * lam_ang**2 / HST_C_ANG_PER_S * 1e23
            return float(flux_jy), float(err_jy)
    if "MJY/SR" in bunit or image.obs_collection == "JWST":
        pixel_area_sr = header.get("PIXAR_SR")
        if pixel_area_sr is None or not math.isfinite(float(pixel_area_sr)) or float(pixel_area_sr) <= 0:
            pixel_area_sr = (image.pixel_scale_arcsec / 206265.0) ** 2
        flux_jy = native_flux * float(pixel_area_sr) * 1e6
        err_jy = abs(native_err) * float(pixel_area_sr) * 1e6
        return float(flux_jy), float(err_jy)
    return float("nan"), float("nan")


def _measure_forced_photometry(
    *,
    image: ScienceImage,
    world: SkyCoord,
) -> dict[str, float]:
    image = _candidate_cutout(image, world)
    x_pix, y_pix = image.wcs.world_to_pixel(world)
    if not (np.isfinite(x_pix) and np.isfinite(y_pix)):
        raise RuntimeError("Candidate position could not be projected into image coordinates.")
    if x_pix < 0 or y_pix < 0 or x_pix >= image.data.shape[1] or y_pix >= image.data.shape[0]:
        raise RuntimeError("Candidate position lies outside image footprint.")
    valid_mask = np.isfinite(image.data)
    radius_pix = max(2.5, 1.2 * _fwhm_guess(image.pixel_scale_arcsec))
    flux_native, err_native, _ = _aperture_fluxes(
        image.data,
        valid_mask,
        np.array([[x_pix, y_pix]], dtype=float),
        radius=radius_pix,
    )
    return {
        "native_flux": float(flux_native[0]),
        "native_err": float(err_native[0]),
        "pixel_scale_arcsec": float(image.pixel_scale_arcsec),
        "radius_pix": float(radius_pix),
        "x_pix": float(x_pix),
        "y_pix": float(y_pix),
    }


def _candidate_stage(epoch_mjd: float, *, pre_mjd: float, post_mjd: float) -> str:
    if not math.isfinite(epoch_mjd):
        return "unknown"
    if epoch_mjd <= pre_mjd + 0.5:
        return "pre"
    if epoch_mjd >= post_mjd - 0.5:
        return "post"
    return "window"


def _fit_stage_sed(stage_df: pd.DataFrame, *, mw_av_mag: float) -> dict[str, Any]:
    use = stage_df[
        np.isfinite(stage_df["wavelength_um"])
        & np.isfinite(stage_df["flux_jy"])
        & np.isfinite(stage_df["err_jy"])
        & (stage_df["err_jy"] > 0)
    ].copy()
    if len(use) < 2:
        return {"fit_ok": False, "n_points": int(len(use))}

    wavelength_um = use["wavelength_um"].to_numpy(dtype=float)
    flux_jy = use["flux_jy"].to_numpy(dtype=float)
    err_jy = use["err_jy"].to_numpy(dtype=float)
    init_scale = max(np.nanmax(np.abs(flux_jy)), 1e-9)

    def objective(theta: np.ndarray) -> float:
        log_t, log_scale, av_int = theta
        t = 10 ** log_t
        scale = 10 ** log_scale
        ext = np.array(
            [
                10 ** (-0.4 * (mw_av_mag + av_int) * extinction_curve_ccm89(lam))
                for lam in wavelength_um
            ],
            dtype=float,
        )
        model = scale * _planck_fnu_jy(wavelength_um, t) * ext
        resid = (model - flux_jy) / err_jy
        return float(np.sum(np.square(resid)))

    result = minimize(
        objective,
        x0=np.array([math.log10(6000.0), math.log10(init_scale) - 15.0, 0.2], dtype=float),
        bounds=[(3.2, 5.0), (-35.0, 10.0), (0.0, 6.0)],
        method="L-BFGS-B",
    )
    log_t, log_scale, av_int = result.x
    best_t = 10 ** log_t
    best_scale = 10 ** log_scale
    ext = np.array(
        [10 ** (-0.4 * (mw_av_mag + av_int) * extinction_curve_ccm89(lam)) for lam in wavelength_um],
        dtype=float,
    )
    model_flux = best_scale * _planck_fnu_jy(wavelength_um, best_t) * ext
    chi2 = float(np.sum(np.square((model_flux - flux_jy) / err_jy)))
    return {
        "fit_ok": bool(result.success),
        "n_points": int(len(use)),
        "temperature_k": float(best_t),
        "scale": float(best_scale),
        "internal_av_mag": float(av_int),
        "mw_av_mag": float(mw_av_mag),
        "chi2": chi2,
        "dof": int(max(len(use) - 3, 0)),
        "fit_message": str(result.message),
    }


def _plot_lightcurve(df: pd.DataFrame, output: Path, *, title: str) -> None:
    with _PLOT_LOCK:
        fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
        for filter_name, group in df.groupby("filter_name"):
            ax.errorbar(
                group["t_min"],
                group["flux_jy"],
                yerr=group["err_jy"],
                fmt="o",
                ms=4,
                capsize=2,
                label=str(filter_name),
                alpha=0.85,
            )
        ax.set_xlabel("MJD")
        ax.set_ylabel("Flux density (Jy)")
        ax.set_title(title)
        if df["filter_name"].nunique() <= 12:
            ax.legend(fontsize=8, ncol=2)
        ax.grid(alpha=0.25)
        fig.savefig(output, dpi=180)
        plt.close(fig)


def _plot_sed(stage_table: pd.DataFrame, fits_by_stage: dict[str, dict[str, Any]], output: Path, *, title: str) -> None:
    with _PLOT_LOCK:
        fig, ax = plt.subplots(figsize=(8.0, 5.4), constrained_layout=True)
        colors = {"pre": "tab:blue", "post": "tab:red"}
        for stage in ["pre", "post"]:
            group = stage_table[stage_table["stage"] == stage].copy()
            if group.empty:
                continue
            ax.errorbar(
                group["wavelength_um"],
                group["flux_jy"],
                yerr=group["err_jy"],
                fmt="o",
                ms=6,
                capsize=2,
                label=f"{stage.title()} photometry",
                color=colors[stage],
            )
            fit_payload = fits_by_stage.get(stage, {})
            if fit_payload.get("fit_ok"):
                lam = np.logspace(np.log10(max(group["wavelength_um"].min() * 0.8, 0.2)), np.log10(max(group["wavelength_um"].max() * 1.2, 30.0)), 200)
                total_av = float(fit_payload["mw_av_mag"]) + float(fit_payload["internal_av_mag"])
                ext = np.array([10 ** (-0.4 * total_av * extinction_curve_ccm89(x)) for x in lam], dtype=float)
                model = float(fit_payload["scale"]) * _planck_fnu_jy(lam, float(fit_payload["temperature_k"])) * ext
                ax.plot(lam, model, color=colors[stage], alpha=0.8)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-8)
        ax.set_xlabel("Wavelength (micron)")
        ax.set_ylabel("Flux density (Jy)")
        ax.set_title(title)
        ax.grid(alpha=0.25, which="both")
        ax.legend()
        fig.savefig(output, dpi=180)
        plt.close(fig)


def _aggregate_stage_photometry(measurements: pd.DataFrame) -> pd.DataFrame:
    usable = measurements[
        np.isfinite(measurements["flux_jy"])
        & np.isfinite(measurements["err_jy"])
        & (measurements["err_jy"] > 0)
        & measurements["stage"].isin(["pre", "post"])
    ].copy()
    rows: list[dict[str, Any]] = []
    for (stage, filter_name), group in usable.groupby(["stage", "filter_name"]):
        weight = 1.0 / np.square(group["err_jy"].to_numpy(dtype=float))
        flux = float(np.sum(group["flux_jy"] * weight) / np.sum(weight))
        err = float(np.sqrt(1.0 / np.sum(weight)))
        rows.append(
            {
                "stage": stage,
                "filter_name": filter_name,
                "wavelength_um": float(np.nanmedian(group["filter_central_um"])),
                "flux_jy": flux,
                "err_jy": err,
                "n_epochs": int(len(group)),
            }
        )
    return pd.DataFrame(rows).sort_values(["stage", "wavelength_um", "filter_name"]).reset_index(drop=True) if rows else pd.DataFrame(columns=["stage", "filter_name", "wavelength_um", "flux_jy", "err_jy", "n_epochs"])


def _follow_candidate(
    candidate: pd.Series,
    *,
    root_dir: Path,
    observations: pd.DataFrame,
    galaxy_master: pd.DataFrame,
) -> dict[str, Any]:
    candidate_id = str(candidate["candidate_id"])
    output_dir = ensure_dir(root_dir / "sed" / candidate_id)
    world = SkyCoord(float(candidate["ra_deg"]) * u.deg, float(candidate["dec_deg"]) * u.deg, frame="icrs")

    galaxy_name = str(candidate["galaxy_name"])
    obs_rows = observations[observations["galaxy_name"].astype(str) == galaxy_name].copy()
    obs_rows = obs_rows.sort_values(["t_min", "obs_collection", "filters"]).reset_index(drop=True)

    measurements: list[dict[str, Any]] = []
    for obs in obs_rows.itertuples(index=False):
        row = {
            "candidate_id": candidate_id,
            "galaxy_name": galaxy_name,
            "obsid": int(obs.obsid) if pd.notna(getattr(obs, "obsid", pd.NA)) else pd.NA,
            "obs_id": str(obs.obs_id),
            "obs_collection": str(obs.obs_collection),
            "instrument_name": str(obs.instrument_name),
            "filter_name": str(obs.filters),
            "filter_central_um": infer_filter_central_um(str(obs.filters), str(obs.instrument_name), str(obs.obs_collection)),
            "t_min": float(obs.t_min) if pd.notna(obs.t_min) else np.nan,
            "t_max": float(obs.t_max) if pd.notna(obs.t_max) else np.nan,
            "stage": _candidate_stage(float(obs.t_min) if pd.notna(obs.t_min) else np.nan, pre_mjd=float(candidate["pre_mjd"]), post_mjd=float(candidate["post_mjd"])),
            "status": "FAILED",
            "error": None,
        }
        try:
            products = _product_query(int(obs.obsid))
            product = _select_best_product(products, obs_collection=str(obs.obs_collection))
            if product is None:
                raise RuntimeError("No suitable calibrated FITS product found.")
            path = _download_product(product, root_dir=root_dir, obs_collection=str(obs.obs_collection), obs_id=str(obs.obs_id))
            image = _load_science_image(
                path,
                obs_collection=str(obs.obs_collection),
                obs_id=str(obs.obs_id),
                filter_name=str(obs.filters),
                instrument_name=str(obs.instrument_name),
            )
            forced = _measure_forced_photometry(image=image, world=world)
            header = _read_calibration_header(path, image.extname)
            flux_jy, err_jy = _jy_from_native(
                native_flux=float(forced["native_flux"]),
                native_err=float(forced["native_err"]),
                image=image,
                header=header,
            )
            row.update(
                {
                    **forced,
                    "flux_jy": flux_jy,
                    "err_jy": err_jy,
                    "snr": float(flux_jy / err_jy) if math.isfinite(flux_jy) and math.isfinite(err_jy) and err_jy > 0 else np.nan,
                    "product_path": str(path),
                    "status": "MEASURED",
                }
            )
        except Exception as exc:
            row.update(
                {
                    "native_flux": np.nan,
                    "native_err": np.nan,
                    "flux_jy": np.nan,
                    "err_jy": np.nan,
                    "snr": np.nan,
                    "pixel_scale_arcsec": np.nan,
                    "radius_pix": np.nan,
                    "x_pix": np.nan,
                    "y_pix": np.nan,
                    "product_path": None,
                    "error": str(exc),
                }
            )
        measurements.append(row)

    measurements_df = pd.DataFrame(measurements).reindex(columns=LIGHTCURVE_COLS)
    photometry_path = output_dir / "all_photometry.parquet"
    photometry_csv = output_dir / "all_photometry.csv"
    measurements_df.to_parquet(photometry_path, index=False)
    measurements_df.to_csv(photometry_csv, index=False)

    stage_table = _aggregate_stage_photometry(measurements_df)
    stage_table_path = output_dir / "stage_sed.csv"
    stage_table.to_csv(stage_table_path, index=False)

    galaxy_row = galaxy_master[galaxy_master["galaxy_name"].astype(str) == galaxy_name].head(1)
    mw_av_mag = float(galaxy_row["mw_av_mag"].iloc[0]) if not galaxy_row.empty and pd.notna(galaxy_row["mw_av_mag"].iloc[0]) else 0.0
    fits_by_stage = {stage: _fit_stage_sed(stage_table[stage_table["stage"] == stage], mw_av_mag=mw_av_mag) for stage in ["pre", "post"]}

    measured = measurements_df[measurements_df["status"] == "MEASURED"].copy()
    lightcurve_path = output_dir / "lightcurve.png"
    if not measured.empty and measured["flux_jy"].notna().any():
        _plot_lightcurve(measured, lightcurve_path, title=f"{galaxy_name} | {candidate_id} light curve")

    sed_plot_path = output_dir / "sed_fit.png"
    if not stage_table.empty:
        _plot_sed(stage_table, fits_by_stage, sed_plot_path, title=f"{galaxy_name} | {candidate_id} pre/post SED")

    summary = {
        "candidate_id": candidate_id,
        "galaxy_name": galaxy_name,
        "created_utc": utc_stamp(),
        "fade_window_mjd": [float(candidate["pre_mjd"]), float(candidate["post_mjd"])],
        "mw_av_mag": mw_av_mag,
        "n_observations_considered": int(len(obs_rows)),
        "n_measurements": int((measurements_df["status"] == "MEASURED").sum()),
        "n_stage_points": int(len(stage_table)),
        "fits": fits_by_stage,
        "lightcurve_path": str(lightcurve_path) if lightcurve_path.exists() else None,
        "sed_plot_path": str(sed_plot_path) if sed_plot_path.exists() else None,
        "photometry_path": str(photometry_path),
        "stage_sed_path": str(stage_table_path),
    }
    summary_path = output_dir / "followup_summary.json"
    write_json(summary_path, summary)
    return {
        "candidate_id": candidate_id,
        "galaxy_name": galaxy_name,
        "status": str(candidate["status"]),
        "n_measurements": int((measurements_df["status"] == "MEASURED").sum()),
        "n_stage_points": int(len(stage_table)),
        "followup_summary_path": str(summary_path),
        "lightcurve_path": str(lightcurve_path) if lightcurve_path.exists() else None,
        "sed_plot_path": str(sed_plot_path) if sed_plot_path.exists() else None,
    }


def _follow_candidate_worker(candidate_payload: dict[str, Any], *, root_dir: str, observations: pd.DataFrame, galaxy_master: pd.DataFrame) -> dict[str, Any]:
    return _follow_candidate(pd.Series(candidate_payload), root_dir=Path(root_dir), observations=observations, galaxy_master=galaxy_master)


def run_candidate_followup(
    *,
    root_dir: Path,
    candidate_ledger_path: Path,
    observation_matrix_path: Path,
    galaxy_master_path: Path,
    statuses: tuple[str, ...] = ("PASS", "REVIEW"),
    max_candidates: int = 0,
    max_workers: int = 1,
) -> dict[str, Path]:
    sed_dir = ensure_dir(root_dir / "sed")
    ledger = pd.read_parquet(candidate_ledger_path)
    observations = pd.read_parquet(observation_matrix_path)
    galaxy_master = pd.read_parquet(galaxy_master_path)

    work = ledger[ledger["status"].astype(str).isin([str(x) for x in statuses])].copy()
    work = work.sort_values(["status", "priority_score"], ascending=[True, False]).reset_index(drop=True)
    if max_candidates > 0:
        work = work.head(max_candidates).copy()

    workers = max(int(max_workers), 1)
    if workers == 1 or len(work) <= 1:
        rows = [_follow_candidate(candidate, root_dir=root_dir, observations=observations, galaxy_master=galaxy_master) for _, candidate in work.iterrows()]
    else:
        rows = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    _follow_candidate_worker,
                    candidate.to_dict(),
                    root_dir=str(root_dir),
                    observations=observations,
                    galaxy_master=galaxy_master,
                ): str(candidate["candidate_id"])
                for _, candidate in work.iterrows()
            }
            for future in as_completed(futures):
                rows.append(future.result())
    index_df = pd.DataFrame(rows)
    index_path = sed_dir / "followup_index.csv"
    index_df.to_csv(index_path, index=False)

    summary_path = sed_dir / "followup_summary.json"
    write_json(
        summary_path,
        {
            "created_utc": utc_stamp(),
            "n_candidates_requested": int(len(work)),
            "statuses": list(statuses),
            "max_candidates": int(max_candidates),
            "max_workers": int(workers),
            "index_path": str(index_path),
        },
    )
    return {
        "followup_index": index_path,
        "summary": summary_path,
    }
