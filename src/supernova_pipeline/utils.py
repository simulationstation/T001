from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SUTC")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, payload: dict[str, Any]) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def decode_object_series(series: pd.Series) -> pd.Series:
    if series.dtype != object:
        return series

    def _decode(value: Any) -> Any:
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="ignore").strip("\x00 ").strip()
        return value

    return series.map(_decode)


def robust_log10(series: pd.Series, *, floor: float = 1e-12) -> pd.Series:
    return np.log10(np.clip(pd.to_numeric(series, errors="coerce"), floor, np.inf))


def normalize_01(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite = np.isfinite(values)
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if not np.any(finite):
        return out.fillna(0.0)
    lo = float(values[finite].min())
    hi = float(values[finite].max())
    if math.isclose(lo, hi):
        return out.fillna(1.0)
    out.loc[finite] = (values.loc[finite] - lo) / (hi - lo)
    return out.fillna(0.0)


def angular_sep_deg(ra1_deg: float, dec1_deg: float, ra2_deg: np.ndarray, dec2_deg: np.ndarray) -> np.ndarray:
    ra1 = np.deg2rad(ra1_deg)
    dec1 = np.deg2rad(dec1_deg)
    ra2 = np.deg2rad(ra2_deg)
    dec2 = np.deg2rad(dec2_deg)
    cos_sep = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.rad2deg(np.arccos(cos_sep))


def json_list(values: list[str]) -> str:
    return json.dumps(sorted(values))

