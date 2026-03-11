from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
from PIL import Image

from supernova_pipeline.uber_nova_site import _blind_match, _build_gif_and_panel, _highest_value


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

    def test_blind_match_uses_radius_threshold(self) -> None:
        truth = pd.DataFrame(
            [
                {"galaxy_name": "MESSIER 061", "sn_name": "SN-A", "truth_ra_deg": 10.0, "truth_dec_deg": 20.0},
            ]
        )
        cluster = pd.Series({"galaxy_name": "MESSIER 061", "ra_deg": 10.0, "dec_deg": 20.0})
        payload = _blind_match(cluster, truth)
        self.assertTrue(payload["matched_to_blind"])
        self.assertEqual(payload["blind_match_name"], "SN-A")

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
