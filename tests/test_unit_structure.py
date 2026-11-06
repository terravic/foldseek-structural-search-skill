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

"""Unit tests for scripts/visualize_structure.py."""

from __future__ import annotations

from pathlib import Path
import pytest

import visualize_structure


def test_parse_synthetic_pdb_metadata(sample_pdb_path: Path) -> None:
    meta = visualize_structure.parse_structure_metadata(sample_pdb_path)
    assert meta["filename"] == "sample_query.pdb"
    assert meta["extension"] == ".pdb"
    assert meta["viewer_format"] == "pdb"
    assert meta["chains"] == ["A"]
    assert meta["residue_count"] == 60
    assert meta["atom_count"] == 475
    assert meta["min_res_seq"] == 1
    assert meta["max_res_seq"] == 60
    assert 55.0 <= meta["mean_bfactor"] <= 98.0
    assert len(meta["ca_trace"]) == 60
    assert len(meta["atoms"]) == 475
    assert {pt["ss"] for pt in meta["ca_trace"]} == {"H", "E", "C"}
    assert "ATOM" in meta["raw_coordinates"]


def test_parse_synthetic_cif_metadata(sample_cif_path: Path) -> None:
    meta = visualize_structure.parse_structure_metadata(sample_cif_path)
    assert meta["filename"] == "sample_query.cif"
    assert meta["extension"] == ".cif"
    assert meta["viewer_format"] == "cif"
    assert meta["chains"] == ["A"]
    assert meta["residue_count"] == 60
    assert meta["atom_count"] == 475
    assert len(meta["ca_trace"]) == 60
    assert len(meta["atoms"]) == 475
    assert {pt["ss"] for pt in meta["ca_trace"]} == {"H", "E", "C"}


def test_parse_second_synthetic_dimer_pdb_and_cif(project_root: Path) -> None:
    pdb2 = project_root / "sample_data" / "sample_query_2.pdb"
    cif2 = project_root / "sample_data" / "sample_query_2.cif"
    meta_pdb = visualize_structure.parse_structure_metadata(pdb2)
    meta_cif = visualize_structure.parse_structure_metadata(cif2)

    for meta, expected_fmt in ((meta_pdb, "pdb"), (meta_cif, "cif")):
        assert meta["viewer_format"] == expected_fmt
        assert meta["chains"] == ["A", "B"]
        assert meta["residue_count"] == 96
        assert meta["atom_count"] == 748
        assert len(meta["ca_trace"]) == 96
        assert len(meta["atoms"]) == 748
        assert 18.0 <= meta["mean_bfactor"] <= 45.0
        assert {pt["ss"] for pt in meta["ca_trace"]} == {"H", "E", "C"}


def test_mmcif_extension_supported(sample_cif_path: Path, tmp_path: Path) -> None:
    mmcif_copy = tmp_path / "synthetic_model.mmcif"
    mmcif_copy.write_text(sample_cif_path.read_text(encoding="utf-8"), encoding="utf-8")
    meta = visualize_structure.parse_structure_metadata(mmcif_copy)
    assert meta["extension"] == ".mmcif"
    assert meta["viewer_format"] == "cif"
    assert meta["residue_count"] == 60


def test_nmr_multimodel_pdb_counts_first_model_only(tmp_path: Path) -> None:
    nmr_pdb = tmp_path / "ensemble.pdb"
    nmr_pdb.write_text(
        "\n".join(
            [
                "MODEL        1",
                "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.00           C",
                "ATOM      2  CA  GLY A   2       4.800   2.100   3.200  1.00 87.00           C",
                "ENDMDL",
                "MODEL        2",
                "ATOM      3  CA  ALA A   1       1.100   2.100   3.100  1.00 84.00           C",
                "ATOM      4  CA  GLY A   2       4.900   2.200   3.300  1.00 86.00           C",
                "ENDMDL",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    meta = visualize_structure.parse_structure_metadata(nmr_pdb)
    assert meta["residue_count"] == 2
    assert meta["atom_count"] == 2
    assert len(meta["ca_trace"]) == 2


@pytest.mark.parametrize(
    "bad_input,match_substr",
    [
        ("P69905", "cannot search by protein name"),
        ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYN", "Raw amino acid sequence"),
        ("", "Empty input path"),
        ("nonexistent_sequence.fasta", "Unsupported file extension"),
    ],
)
def test_reject_non_coordinate_inputs(bad_input: str, match_substr: str) -> None:
    with pytest.raises(ValueError, match=match_substr):
        visualize_structure.validate_structure_file(bad_input)


def test_reject_missing_pdb_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing_coords.pdb"
    with pytest.raises(FileNotFoundError, match="not found"):
        visualize_structure.validate_structure_file(missing)


def test_reject_empty_and_atomless_pdb(tmp_path: Path) -> None:
    empty_pdb = tmp_path / "empty.pdb"
    empty_pdb.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty"):
        visualize_structure.parse_structure_metadata(empty_pdb)

    header_only = tmp_path / "header_only.pdb"
    header_only.write_text("HEADER    SYNTHETIC EMPTY STRUCTURE\nEND\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No valid macromolecular ATOM"):
        visualize_structure.parse_structure_metadata(header_only)


def test_infer_target_structure_urls() -> None:
    pdb_url, pdb_fmt = visualize_structure.infer_target_structure_url("1UBQ_A", "pdb100")
    assert pdb_url == "https://files.rcsb.org/download/1UBQ.cif"
    assert pdb_fmt == "cif"

    af_url, af_fmt = visualize_structure.infer_target_structure_url("AF-P69905-F1", "afdb50")
    assert af_url == "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.cif"
    assert af_fmt == "cif"

    cath_url, cath_fmt = visualize_structure.infer_target_structure_url(
        "cath|4_3_0|1ubqA00/1-76",
        "cath50",
        "cath|4_3_0|1ubqA00/1-76 Ubiquitin-like superfamily",
    )
    assert cath_url == "https://files.rcsb.org/download/1UBQ.cif"
    assert cath_fmt == "cif"

    none_url, none_fmt = visualize_structure.infer_target_structure_url("MGYP001234567890", "mgnify_esm30")
    assert none_url is None
    assert none_fmt is None
