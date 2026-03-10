from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir, utc_stamp, write_json


CANDIDATE_COLUMNS = [
    "candidate_id",
    "pair_id",
    "galaxy_name",
    "ra_deg",
    "dec_deg",
    "status",
    "priority_score",
    "fade_sigma",
    "fade_fraction",
    "S_raw",
    "S_min_LOO",
    "S_robust",
    "d_max",
    "depth_margin_post",
    "crowding_index",
    "blend_risk",
    "host_bg_penalty",
    "astrometric_residual",
    "cross_reduction_agreement",
    "pre_obsid",
    "post_obsid",
    "pre_obs_id",
    "post_obs_id",
    "pre_filter",
    "post_filter",
    "pre_mjd",
    "post_mjd",
    "cutout_path",
    "reason_codes_json",
    "warning_flags_json",
    "provenance_json",
]


def build_detection_queue(epoch_pairs_path: Path, *, root_dir: Path, top_per_galaxy: int = 5) -> dict[str, Path]:
    candidate_dir = ensure_dir(root_dir / "candidates")
    pairs = pd.read_parquet(epoch_pairs_path)
    if pairs.empty:
        queue = pd.DataFrame(
            columns=[
                "queue_rank",
                "pair_id",
                "galaxy_name",
                "galaxy_rank",
                "obs_id_1",
                "obs_id_2",
                "pair_score",
                "compatibility",
                "baseline_days",
                "scan_status",
            ]
        )
    else:
        queue = pairs.groupby("galaxy_name", group_keys=False).head(top_per_galaxy).copy()
        queue = queue.sort_values(["pair_score", "baseline_days"], ascending=[False, False]).reset_index(drop=True)
        queue.insert(0, "queue_rank", np.arange(1, len(queue) + 1))
        queue["scan_status"] = "UNSCANNED"

    queue_path = candidate_dir / "detection_queue.parquet"
    queue_csv_path = candidate_dir / "detection_queue.csv"
    queue.to_parquet(queue_path, index=False)
    queue.to_csv(queue_csv_path, index=False)

    ledger = pd.DataFrame(columns=CANDIDATE_COLUMNS)
    ledger_path = candidate_dir / "candidate_ledger.parquet"
    ledger.to_parquet(ledger_path, index=False)

    review_queue = pd.DataFrame(columns=["candidate_id", "galaxy_name", "priority_score", "status"])
    review_queue_path = candidate_dir / "review_queue.csv"
    review_queue.to_csv(review_queue_path, index=False)

    summary_path = candidate_dir / "candidate_summary.json"
    write_json(
        summary_path,
        {
            "created_utc": utc_stamp(),
            "n_detection_queue_rows": int(len(queue)),
            "n_candidate_rows": int(len(ledger)),
            "top_per_galaxy": int(top_per_galaxy),
            "note": "Detection queue is populated from epoch-pair opportunities. Candidate ledger remains empty until pixel-level or catalog-level fading discovery begins.",
        },
    )

    return {
        "detection_queue": queue_path,
        "detection_queue_csv": queue_csv_path,
        "candidate_ledger": ledger_path,
        "review_queue": review_queue_path,
        "summary": summary_path,
    }
