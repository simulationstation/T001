from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
from PIL import Image

from supernova_pipeline.uber_nova_site import _build_gif_and_panel, _highest_value, _merge_benchmark_recovered


class UberNovaSiteTests(unittest.TestCase):
    def test_build_gif_and_panel_creates_expected_assets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cutout = tmp_path / "triptych.png"
            image = Image.new("RGB", (900, 300))
            for x in range(300):
                for y in range(300):
                    image.putpixel((x, y), (255, 0, 0))
                    image.putpixel((x + 300, y), (0, 255, 0))
                    image.putpixel((x + 600, y), (0, 0, 255))
            image.save(cutout)

            gif_path = tmp_path / "blink.gif"
            panel_path = tmp_path / "panel.png"
            diff_path = tmp_path / "diff.png"
            _build_gif_and_panel(
                cutout_path=cutout,
                gif_path=gif_path,
                panel_path=panel_path,
                diff_path=diff_path,
            )

            self.assertTrue(gif_path.exists())
            self.assertTrue(panel_path.exists())
            self.assertTrue(diff_path.exists())

    def test_merge_benchmark_recovered_keeps_only_recovered_truth_groups(self) -> None:
        evaluation = pd.DataFrame(
            [
                {"galaxy_name": "MESSIER 066", "sn_name": "SN-A", "recovered": True, "best_sep_arcsec": 0.42, "matched_pair_id": "pair-a", "matched_event_sign": "fade", "matched_diff_sigma": 11.0},
                {"galaxy_name": "MESSIER 077", "sn_name": "SN-B", "recovered": False, "best_sep_arcsec": 2.1, "matched_pair_id": "pair-b", "matched_event_sign": "brighten", "matched_diff_sigma": -8.0},
            ]
        )
        pairs = pd.DataFrame(
            [
                {"pair_id": "pair-a", "galaxy_name": "MESSIER 066", "sn_name": "SN-A", "truth_ra_deg": 10.0, "truth_dec_deg": 20.0, "filter_1": "F555W", "filter_2": "F555W"},
                {"pair_id": "pair-b", "galaxy_name": "MESSIER 077", "sn_name": "SN-B", "truth_ra_deg": 30.0, "truth_dec_deg": 40.0, "filter_1": "F814W", "filter_2": "F814W"},
            ]
        )
        truth = pd.DataFrame(
            [
                {"galaxy_name": "MESSIER 066", "sn_name": "SN-A", "truth_ra_deg": 10.0, "truth_dec_deg": 20.0},
                {"galaxy_name": "MESSIER 077", "sn_name": "SN-B", "truth_ra_deg": 30.0, "truth_dec_deg": 40.0},
            ]
        )
        pair_summary = pd.DataFrame(
            [
                {"pair_id": "pair-a", "registration_residual_px": 0.12, "n_residual_detections": 50, "n_survivors": 2},
                {"pair_id": "pair-b", "registration_residual_px": 0.34, "n_residual_detections": 12, "n_survivors": 0},
            ]
        )
        merged = _merge_benchmark_recovered(evaluation, pairs, truth, pair_summary)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged.iloc[0]["pair_id"], "pair-a")
        self.assertEqual(merged.iloc[0]["sn_name"], "SN-A")
        self.assertAlmostEqual(float(merged.iloc[0]["registration_residual_px"]), 0.12, places=6)

    def test_highest_value_requires_real_support(self) -> None:
        cluster = pd.Series(
            {
                "branch_class": "STRONG_EXPORT_FAILURE_LIKE",
                "branch_rank_score": 0.91,
                "export_failure_score": 0.81,
                "n_support_pairs": 4,
                "systematic_risk": 0.08,
                "forced_late_return_fraction": 0.0,
                "branch_confidence": "MEDIUM",
            }
        )
        self.assertTrue(_highest_value(cluster))
        cluster["n_support_pairs"] = 1
        self.assertFalse(_highest_value(cluster))


if __name__ == "__main__":
    unittest.main()
