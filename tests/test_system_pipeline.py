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

"""End-to-end system, mocked HTTP client, dashboard schema, and repository validation tests."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any
import pytest

import generate_dashboard
import run_pipeline
import search


class FakeFoldseekHttpClient:
    """Deterministic fake HttpClient simulating Foldseek ticket submission and polling."""

    def __init__(self, sample_payload: dict[str, Any], fail_ticket: bool = False) -> None:
        self.sample_payload = sample_payload
        self.fail_ticket = fail_ticket
        self.calls: list[dict[str, Any]] = []
        self._poll_count = 0

    def fetch_json(
        self,
        url: str,
        *,
        method: str = "GET",
        headers: dict[str, str] | None = None,
        data: bytes | None = None,
        timeout: float = 30.0,
    ) -> Any:
        self.calls.append(
            {"url": url, "method": method, "headers": headers, "data": data, "timeout": timeout}
        )
        if url.endswith("/api/ticket") and method == "POST":
            assert data is not None
            assert b'name="mode"\r\n\r\n3diaa' in data
            assert b'name="database[]"' in data
            assert b'name="q"; filename="' in data
            return {"id": "synth-ticket-001", "status": "PENDING"}

        if "/api/ticket/synth-ticket-001" in url:
            self._poll_count += 1
            if self.fail_ticket:
                return {"id": "synth-ticket-001", "status": "ERROR"}
            if self._poll_count < 2:
                return {"id": "synth-ticket-001", "status": "RUNNING"}
            return {"id": "synth-ticket-001", "status": "COMPLETE"}

        if "/api/result/synth-ticket-001/0" in url:
            return self.sample_payload

        raise AssertionError(f"Unexpected URL requested: {url}")


def test_live_ticket_polling_workflow_with_fake_client(
    sample_pdb_path: Path,
    sample_results_path: Path,
) -> None:
    sample_payload = json.loads(sample_results_path.read_text(encoding="utf-8"))
    fake_client = FakeFoldseekHttpClient(sample_payload)

    result = search.submit_and_poll_foldseek(
        sample_pdb_path,
        ["pdb100", "afdb50"],
        client=fake_client,
        poll_interval_secs=0.0,
        finalization_wait_secs=0.0,
    )
    assert result == sample_payload
    assert len(fake_client.calls) == 4  # 1 POST + 2 status GETs + 1 result GET


def test_ticket_error_status_raises_runtime_error(
    sample_pdb_path: Path,
    sample_results_path: Path,
) -> None:
    sample_payload = json.loads(sample_results_path.read_text(encoding="utf-8"))
    failing_client = FakeFoldseekHttpClient(sample_payload, fail_ticket=True)

    with pytest.raises(RuntimeError, match="ERROR"):
        search.submit_and_poll_foldseek(
            sample_pdb_path,
            ["pdb100"],
            client=failing_client,
            poll_interval_secs=0.0,
            finalization_wait_secs=0.0,
        )


def test_end_to_end_pipeline_and_dashboard_verification(
    project_root: Path,
    sample_pdb_path: Path,
    sample_cif_path: Path,
    isolated_workspace: Path,
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "pipeline_pdb_out"
    artifacts = run_pipeline.run_pipeline(
        input_file=sample_pdb_path,
        output_dir=out_dir,
        databases="pdb100,afdb50,afdb-swissprot,cath50",
        mock=True,
        workspace_root=isolated_workspace,
    )

    # 1. Verify license notification guard created timestamped license file
    license_file = isolated_workspace / ".licenses" / "foldseek_structural_search_LICENSE.txt"
    assert license_file.is_file()
    lic_text = license_file.read_text(encoding="utf-8")
    assert "https://search.foldseek.com/search" in lic_text
    assert "https://github.com/steineggerlab/foldseek" in lic_text
    assert "Timestamp:" in lic_text

    # 2. Verify generated outputs exist and are non-empty
    assert artifacts["json"].name == "sample_query_foldseek_results.json"
    assert artifacts["markdown"].name == "sample_query_foldseek_results.md"
    assert artifacts["dashboard"].name == "dashboard.html"
    for key in ("json", "markdown", "dashboard"):
        assert artifacts[key].is_file()
        assert artifacts[key].stat().st_size > 0

    # 3. Verify Markdown synthesis sections
    md_text = artifacts["markdown"].read_text(encoding="utf-8")
    assert "## 1. Executive Summary" in md_text
    assert "## 2. Top Structural Homologues" in md_text
    assert "## 3. Remote Twilight-Zone Structural Homologues" in md_text
    assert "### All Structural Matches" in md_text
    assert "van Kempen" in md_text

    # 4. Verify standalone dashboard HTML structure, theme tokens, and embedded payload
    html_text = artifacts["dashboard"].read_text(encoding="utf-8")
    for required_token in (
        "3Dmol",
        "highlightResidueRange",
        "clearResidueHighlight",
        'id="btn-reset-profile"',
        'id="software-viewport"',
        'id="theme-btn-light"',
        'id="theme-btn-dark"',
        'id="theme-btn-system"',
        'value="cartoon_sticks"',
        'value="cartoon"',
        'value="surface"',
        'value="stick"',
        'value="sphere"',
        'value="spectrum"',
        'value="ss"',
        'value="bfactor"',
        'value="conservation"',
        'value="chain"',
        "removeAllSurfaces",
        "#FAFAF9",
        "#0D0E11",
        "#16171D",
        "#0D9488",
        "#D97706",
        "#4F46E5",
        "#059669",
        "#E11D48",
        "van Kempen",
        "https://search.foldseek.com/search",
        "https://github.com/steineggerlab/foldseek",
    ):
        assert required_token in html_text, f"Missing token in dashboard.html: {required_token}"

    match = re.search(
        r'<script id="foldseek-dashboard-payload" type="application/json">(.*?)</script>',
        html_text,
        flags=re.DOTALL,
    )
    assert match is not None
    embedded_payload = json.loads(match.group(1))
    schema = json.loads(
        (project_root / "assets" / "dashboard_payload_schema.json").read_text(encoding="utf-8")
    )
    for req_key in schema["required"]:
        assert req_key in embedded_payload
    assert embedded_payload["query_structure"]["residue_count"] == 60
    assert embedded_payload["query_structure"]["atom_count"] == 475
    assert len(embedded_payload["query_structure"]["atoms"]) == 475
    assert {pt["ss"] for pt in embedded_payload["query_structure"]["ca_trace"]} == {"H", "E", "C"}
    assert "ATOM" in embedded_payload["query_structure"]["raw_coordinates"]
    assert len(embedded_payload["hits"]) == 20

    # 5. Also verify mmCIF (.cif) input through full pipeline
    cif_out = tmp_path / "pipeline_cif_out"
    cif_artifacts = run_pipeline.run_pipeline(
        input_file=sample_cif_path,
        output_dir=cif_out,
        databases="pdb100,afdb50",
        mock=True,
        workspace_root=isolated_workspace,
    )
    assert cif_artifacts["dashboard"].is_file()


def test_repository_license_and_encoding_integrity(project_root: Path) -> None:
    skill_licenses = (project_root / "SKILL_LICENSES.md").read_text(encoding="utf-8")
    assert "foldseek-structural-search" in skill_licenses
    assert "https://search.foldseek.com/search" in skill_licenses
    assert "https://github.com/steineggerlab/foldseek" in skill_licenses

    for py_script in sorted((project_root / "scripts").glob("*.py")):
        txt = py_script.read_text(encoding="utf-8")
        assert "Copyright 2026 Google LLC" in txt, f"Missing Copyright header in {py_script.name}"
        assert "http://www.apache.org/licenses/LICENSE-2.0" in txt

    checked_files = (
        list(project_root.glob("*.md"))
        + list(project_root.glob("*.html"))
        + [project_root / "LICENSE", project_root / "requirements.txt"]
        + list((project_root / "documents").glob("*.md"))
        + list((project_root / "references").glob("*"))
        + list((project_root / "assets").glob("*"))
        + list((project_root / "sample_data").glob("*"))
        + list((project_root / "scripts").glob("*.py"))
        + list((project_root / "tests").glob("*.py"))
    )
    assert len(checked_files) >= 20
    binary_image_suffixes = {".png", ".jpg", ".jpeg", ".webp", ".gif", ".ico"}
    local_root_str = str(project_root)

    for fpath in checked_files:
        if not fpath.is_file() or fpath.suffix.lower() in binary_image_suffixes:
            continue
        rel_name = fpath.relative_to(project_root).as_posix()
        content = fpath.read_text(encoding="utf-8")
        assert content.isascii(), f"Non-ASCII byte sequence found in {rel_name}"
        assert local_root_str not in content, f"Local absolute workspace path found in {rel_name}"


def test_structure_only_dashboard_and_on_demand_pdb_loader(
    sample_pdb_path: Path,
    tmp_path: Path,
) -> None:
    out_html = tmp_path / "pdb_only_dashboard.html"
    compiled = generate_dashboard.build_dashboard(
        input_structure=sample_pdb_path,
        results_json=None,
        output_html=out_html,
        title="On-Demand PDB Structure Inspector",
    )
    assert compiled.is_file()
    html_text = compiled.read_text(encoding="utf-8")
    for required_hook in (
        'id="input-pdb-file"',
        'id="btn-load-pdb"',
        'id="btn-card-load-pdb"',
        'id="btn-run-foldseek"',
        'id="struct-details-pane"',
        'id="struct-chains-tbody"',
        'id="struct-sse-tbody"',
        'id="struct-seq-box"',
        "window.loadStructureFromText",
        "window.loadHitsJsonFromText",
    ):
        assert required_hook in html_text, f"Missing on-demand PDB loader hook: {required_hook}"

    match = re.search(
        r'<script id="foldseek-dashboard-payload" type="application/json">(.*?)</script>',
        html_text,
        flags=re.DOTALL,
    )
    assert match is not None
    payload = json.loads(match.group(1))
    assert payload["hits"] == []
    assert payload["summary_kpis"]["total_hits"] == 0
    assert payload["query_structure"]["residue_count"] == 60
