# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Unit tests for scripts/visualize_alignments.py."""

from __future__ import annotations

import json
from pathlib import Path

import visualize_alignments


def test_query_coverage_formula_and_edge_cases() -> None:
    cov_pct, cov_str = visualize_alignments.compute_query_coverage(
        {"qStartPos": 1, "qEndPos": 60, "qLen": 60}
    )
    assert cov_pct == 100.0
    assert cov_str == "100.0%"

    # Partial span: residues 11..40 on 60-residue query = 30/60 = 50.0%
    cov_pct2, cov_str2 = visualize_alignments.compute_query_coverage(
        {"qStartPos": 11, "qEndPos": 40, "qLen": 60}
    )
    assert cov_pct2 == 50.0
    assert cov_str2 == "50.0%"

    # Cap over-100% anomalies cleanly at 100.0%
    cov_pct3, cov_str3 = visualize_alignments.compute_query_coverage(
        {"qStartPos": 1, "qEndPos": 75, "qLen": 60}
    )
    assert cov_pct3 == 100.0
    assert cov_str3 == "100.0%"

    # Invalid / missing qLen or inverted coordinates -> N/A
    assert visualize_alignments.compute_query_coverage({"qStartPos": 10, "qEndPos": 10, "qLen": 60}) == (
        0.0,
        "N/A",
    )
    assert visualize_alignments.compute_query_coverage({"qStartPos": 1, "qEndPos": 20, "qLen": 0}) == (
        0.0,
        "N/A",
    )


def test_seq_identity_normalization() -> None:
    assert visualize_alignments.normalize_seq_identity({"seqId": 0.214}) == (21.4, "21.4%")
    assert visualize_alignments.normalize_seq_identity({"seqIdentity": 46.7}) == (46.7, "46.7%")
    assert visualize_alignments.normalize_seq_identity({"fident": 1.0, "missmatches": 0, "alnLength": 60}) == (
        100.0,
        "100.0%",
    )
    assert visualize_alignments.normalize_seq_identity({"seqId": 1.0, "missmatches": 55, "alnLength": 60}) == (
        1.0,
        "1.0%",
    )


def test_homology_regime_classification() -> None:
    assert (
        visualize_alignments.classify_homology_regime(0.95, 52.0, 95.0)
        == visualize_alignments.REGIME_HIGH_IDENTITY
    )
    assert (
        visualize_alignments.classify_homology_regime(0.94, 17.5, 92.0)
        == visualize_alignments.REGIME_TWILIGHT_ZONE
    )
    assert (
        visualize_alignments.classify_homology_regime(0.68, 25.0, 40.0)
        == visualize_alignments.REGIME_PARTIAL_MOTIF
    )
    assert (
        visualize_alignments.classify_homology_regime(0.32, 19.0, 30.0)
        == visualize_alignments.REGIME_LOW_CONFIDENCE
    )


def test_sorting_and_300_hit_cap() -> None:
    synthetic_alignments = [
        {
            "target": f"SYN_{i:04d}_A Synthetic test family protein [Synthetic]",
            "prob": 0.85,
            "eval": float(400 - i) * 1e-5,
            "seqId": 0.25,
            "qStartPos": 1,
            "qEndPos": 50,
            "qLen": 60,
            "alnLength": 50,
        }
        for i in range(350)
    ]
    raw_payload = {"results": [{"db": "pdb100", "alignments": [synthetic_alignments]}]}
    hits = visualize_alignments.extract_and_normalize_hits(raw_payload)
    assert len(hits) == 300
    assert hits[0]["rank"] == 1
    assert hits[-1]["rank"] == 300
    evalues = [h["evalue"] for h in hits]
    assert evalues == sorted(evalues)


def test_sample_fixture_analytics_and_clusters(sample_results_path: Path) -> None:
    raw = json.loads(sample_results_path.read_text(encoding="utf-8"))
    analytics = visualize_alignments.build_alignment_analytics(raw, fallback_query_len=60)

    hits = analytics["hits"]
    assert len(hits) == 20
    # Verify all four regimes are represented in synthetic fixture
    regimes = {h["homology_regime"] for h in hits}
    assert visualize_alignments.REGIME_HIGH_IDENTITY in regimes
    assert visualize_alignments.REGIME_TWILIGHT_ZONE in regimes
    assert visualize_alignments.REGIME_PARTIAL_MOTIF in regimes
    assert visualize_alignments.REGIME_LOW_CONFIDENCE in regimes

    profile = analytics["coverage_profile"]
    assert len(profile) == 60
    assert max(profile) >= 15

    clusters = analytics["functional_clusters"]
    assert len(clusters) >= 2
    assert all(c["count"] >= 1 and 0.0 <= c["mean_prob"] <= 1.0 for c in clusters)
