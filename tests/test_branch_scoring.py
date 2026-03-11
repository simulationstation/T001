from __future__ import annotations

import unittest

import pandas as pd

from supernova_pipeline.branch_scoring import (
    CLASS_STRONG_EXPORT,
    CLASS_SYSTEMATIC,
    _band_group,
    _classify_cluster,
    _cluster_single_galaxy,
)


class BranchScoringTests(unittest.TestCase):
    def test_band_group_prefers_instrument_family(self) -> None:
        self.assertEqual(_band_group("F1000W", "MIRI/IMAGE"), "MIR")
        self.assertEqual(_band_group("F200W", "NIRCAM/IMAGE"), "NIR")
        self.assertEqual(_band_group("F555W", "WFC3/UVIS"), "OPTICAL")

    def test_cluster_single_galaxy_groups_close_positions(self) -> None:
        frame = pd.DataFrame(
            [
                {"galaxy_name": "TEST", "ra_deg": 10.0, "dec_deg": 20.0, "diff_sigma": 8.0, "diff_sigma_min": 7.0},
                {"galaxy_name": "TEST", "ra_deg": 10.0 + 0.00001, "dec_deg": 20.0 + 0.00001, "diff_sigma": 7.5, "diff_sigma_min": 6.5},
                {"galaxy_name": "TEST", "ra_deg": 10.01, "dec_deg": 20.01, "diff_sigma": 6.0, "diff_sigma_min": 5.5},
            ]
        )
        clusters = _cluster_single_galaxy(frame, radius_arcsec=0.25)
        self.assertEqual(len(clusters), 2)
        sizes = sorted(len(cluster.member_indices) for cluster in clusters)
        self.assertEqual(sizes, [1, 2])

    def test_classify_cluster_prefers_strong_export_when_clean(self) -> None:
        row = pd.Series(
            {
                "export_failure_score": 0.82,
                "partial_branch_score": 0.41,
                "dust_survivor_score": 0.12,
                "systematic_risk": 0.18,
                "unresolved_variability_score": 0.20,
                "forced_photometry_coherence": 0.63,
                "registration_confidence": 0.71,
                "artifact_penalty": 0.22,
                "coverage_completeness": 0.65,
                "forced_n_pairs_measured": 4,
                "fade_persistence": 0.69,
                "n_support_pairs": 3,
                "n_pass_members": 1,
            }
        )
        label, confidence = _classify_cluster(row)
        self.assertEqual(label, CLASS_STRONG_EXPORT)
        self.assertIn(confidence, {"HIGH", "MEDIUM"})

    def test_classify_cluster_prefers_systematic_when_penalty_dominates(self) -> None:
        row = pd.Series(
            {
                "export_failure_score": 0.38,
                "partial_branch_score": 0.22,
                "dust_survivor_score": 0.18,
                "systematic_risk": 0.74,
                "unresolved_variability_score": 0.35,
                "forced_photometry_coherence": 0.19,
                "registration_confidence": 0.28,
                "artifact_penalty": 0.71,
                "coverage_completeness": 0.46,
                "forced_n_pairs_measured": 3,
            }
        )
        label, _ = _classify_cluster(row)
        self.assertEqual(label, CLASS_SYSTEMATIC)


if __name__ == "__main__":
    unittest.main()
