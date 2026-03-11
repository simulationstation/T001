from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from PIL import Image

from supernova_pipeline.uber_nova_site import _build_gif_and_panel


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


if __name__ == "__main__":
    unittest.main()
