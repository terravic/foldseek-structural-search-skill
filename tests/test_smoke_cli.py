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

"""CLI smoke and negative validation tests for foldseek-structural-search."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


def _run_cli(args: list[str], cwd: Path, ws_root: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["FOLDSEEK_WORKSPACE_ROOT"] = str(ws_root)
    return subprocess.run(
        [sys.executable, *args],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_backward_compatible_search_stdout_redirection(
    project_root: Path,
    sample_pdb_path: Path,
    isolated_workspace: Path,
    tmp_path: Path,
) -> None:
    out_json = tmp_path / "foldseek_compat.json"
    proc = _run_cli(
        [
            str(project_root / "scripts" / "search.py"),
            str(sample_pdb_path),
            "-o",
            str(out_json),
            "--databases",
            "pdb100,afdb50",
            "--mock",
        ],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc.returncode == 0, proc.stderr
    assert out_json.is_file()
    saved = json.loads(out_json.read_text(encoding="utf-8"))
    assert "results" in saved

    assert "### All Structural Matches" in proc.stdout
    assert "| Target ID | Q-Cov | Prob | E-value | Seq Identity | Aln Length |" in proc.stdout
    assert "van Kempen" in proc.stdout


def test_monorepo_symlink_cli_compatibility(
    project_root: Path,
     isolated_workspace: Path,
    tmp_path: Path,
) -> None:
    symlink_script = project_root / "skills" / "foldseek_structural_search" / "scripts" / "search.py"
    symlink_sample = (
        project_root / "skills" / "foldseek_structural_search" / "sample_data" / "sample_query.pdb"
    )
    assert symlink_script.is_file()
    assert symlink_sample.is_file()

    out_json = tmp_path / "symlink_out.json"
    proc = _run_cli(
        [str(symlink_script), str(symlink_sample), "-o", str(out_json), "--mock"],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc.returncode == 0, proc.stderr
    assert "### All Structural Matches" in proc.stdout


def test_invalid_database_exits_with_code_1_and_lists_allowed(
    project_root: Path,
    sample_pdb_path: Path,
    isolated_workspace: Path,
    tmp_path: Path,
) -> None:
    proc = _run_cli(
        [
            str(project_root / "scripts" / "search.py"),
            str(sample_pdb_path),
            "-o",
            str(tmp_path / "unused.json"),
            "--databases",
            "invalid_db",
            "--mock",
        ],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc.returncode == 1
    for allowed_db in (
        "afdb50",
        "afdb-swissprot",
        "pdb100",
        "BFVD",
        "mgnify_esm30",
        "cath50",
        "gmgcl_id",
        "bfmd",
        "afdb-proteome",
    ):
        assert allowed_db in proc.stderr


def test_unsupported_fasta_and_missing_file_exit_code_1(
    project_root: Path,
    isolated_workspace: Path,
    tmp_path: Path,
) -> None:
    fasta_file = tmp_path / "sequence.fasta"
    fasta_file.write_text(">SYN1\nMQIFVKTLTGKTITLEVEPS\n", encoding="utf-8")

    proc_fasta = _run_cli(
        [str(project_root / "scripts" / "search.py"), str(fasta_file), "--mock"],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc_fasta.returncode == 1
    assert ".pdb" in proc_fasta.stderr

    proc_missing = _run_cli(
        [str(project_root / "scripts" / "search.py"), str(tmp_path / "missing.pdb"), "--mock"],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc_missing.returncode == 1

    proc_uniprot = _run_cli(
        [str(project_root / "scripts" / "run_pipeline.py"), "P69905", "-o", str(tmp_path / "out")],
        cwd=project_root,
        ws_root=isolated_workspace,
    )
    assert proc_uniprot.returncode == 1
