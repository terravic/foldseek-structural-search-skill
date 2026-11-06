#!/usr/bin/env python3
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
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///

"""Standalone HTML5 and WebGL Foldseek dashboard compiler."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

_SCRIPTS_DIR = Path(__file__).resolve().parent
_SKILL_ROOT = _SCRIPTS_DIR.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import visualize_alignments  # noqa: E402
import visualize_structure  # noqa: E402

LOGGER = logging.getLogger(__name__)

FOLDSEEK_CITATION_TEXT = (
    "van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., "
    "Gilchrist, C.L.M., Soding, J., & Steinegger, M. (2023). Fast and accurate "
    "protein structure search with Foldseek. Nature Biotechnology, 41(12), "
    "1673-1682. doi:10.1038/s41587-023-01773-0."
)
FOLDSEEK_SERVER_URL = "https://search.foldseek.com/search"
FOLDSEEK_GITHUB_URL = "https://github.com/steineggerlab/foldseek"
RESEARCH_DISCLAIMER = (
    "Provided AS-IS for scientific research and structural bioinformatics analysis. "
    "Verify structural homology hypotheses experimentally before clinical or therapeutic use."
)


def detect_host_system_theme() -> str:
    """Detect native host operating system light or dark appearance preference."""
    env_theme = os.environ.get("GTK_THEME", "").lower()
    if "dark" in env_theme:
        return "dark"

    colorfgbg = os.environ.get("COLORFGBG", "")
    if ";" in colorfgbg:
        bg_token = colorfgbg.rsplit(";", 1)[-1].strip()
        if bg_token.isdigit() and int(bg_token) in {0, 8}:
            return "dark"

    try:
        proc = subprocess.run(
            ["gsettings", "get", "org.gnome.desktop.interface", "color-scheme"],
            capture_output=True,
            text=True,
            timeout=1.5,
            check=False,
        )
        if proc.returncode == 0 and "dark" in proc.stdout.lower():
            return "dark"
    except (OSError, subprocess.SubprocessError):
        pass

    kde_cfg = Path.home() / ".config" / "kdeglobals"
    if kde_cfg.is_file():
        try:
            if "dark" in kde_cfg.read_text(encoding="utf-8", errors="ignore").lower():
                return "dark"
        except OSError:
            pass

    return "light"


def build_dashboard_payload(
    input_structure: str | Path,
    results_data: dict[str, Any] | str | Path | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Assemble and validate the normalized dashboard payload dictionary.

    Args:
        input_structure: Path to the query `.pdb`, `.cif`, or `.mmcif` coordinate file.
        results_data: Optional parsed Foldseek JSON dictionary or path to `.json` file.
        title: Optional custom dashboard title.

    Returns:
        Normalized payload dictionary conforming to `assets/dashboard_payload_schema.json`.
    """
    query_meta = visualize_structure.parse_structure_metadata(input_structure)

    if results_data is None:
        raw_json: dict[str, Any] = {"results": []}
    elif isinstance(results_data, (str, Path)):
        res_path = Path(results_data)
        if not res_path.is_file():
            raise FileNotFoundError(
                f"Foldseek JSON results file not found at '{results_data}'."
            )
        raw_json = json.loads(res_path.read_text(encoding="utf-8"))
    elif isinstance(results_data, dict):
        raw_json = results_data
    else:
        raise ValueError("results_data must be a dict, path to a JSON file, or None.")

    analytics = visualize_alignments.build_alignment_analytics(
        raw_json=raw_json,
        fallback_query_len=query_meta["residue_count"],
    )

    resolved_title = (
        title.strip()
        if title and title.strip()
        else f"Foldseek Structural Homology Analysis - {query_meta['filename']}"
    )

    payload: dict[str, Any] = {
        "title": resolved_title,
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "system_theme": detect_host_system_theme(),
        "query_structure": query_meta,
        "summary_kpis": analytics["summary_kpis"],
        "hits": analytics["hits"],
        "coverage_profile": analytics["coverage_profile"],
        "database_breakdown": analytics["database_breakdown"],
        "functional_clusters": analytics["functional_clusters"],
        "provenance": {
            "citation_text": FOLDSEEK_CITATION_TEXT,
            "foldseek_url": FOLDSEEK_SERVER_URL,
            "github_url": FOLDSEEK_GITHUB_URL,
            "disclaimer": RESEARCH_DISCLAIMER,
        },
    }
    return payload


def _safe_json_for_script_tag(data: dict[str, Any]) -> str:
    """Serialize JSON safely for embedding inside an HTML `<script type="application/json">` tag."""
    raw = json.dumps(data, ensure_ascii=True, separators=(",", ":"))
    return raw.replace("</", "<\\/").replace("<!--", "<\\!--")


def build_dashboard(
    input_structure: str | Path,
    results_json: dict[str, Any] | str | Path | None = None,
    output_html: str | Path = "dashboard.html",
    title: str | None = None,
    template_path: Path | None = None,
) -> Path:
    """Compile query coordinates and Foldseek alignments into a standalone HTML5 dashboard.

    Args:
        input_structure: Path to query `.pdb`, `.cif`, or `.mmcif` coordinate file.
        results_json: Optional parsed Foldseek JSON dict or path to saved Foldseek JSON file.
        output_html: Destination `.html` file path.
        title: Optional custom heading title.
        template_path: Optional override path to `dashboard_template.html`.

    Returns:
        Path to the written standalone `.html` dashboard file.
    """
    tpl_file = template_path or (_SKILL_ROOT / "assets" / "dashboard_template.html")
    if not tpl_file.is_file():
        raise FileNotFoundError(f"Dashboard HTML template missing at '{tpl_file}'.")

    payload = build_dashboard_payload(
        input_structure=input_structure,
        results_data=results_json,
        title=title,
    )

    template_html = tpl_file.read_text(encoding="utf-8")
    escaped_title = html.escape(payload["title"], quote=True)
    safe_json = _safe_json_for_script_tag(payload)

    rendered_html = template_html.replace("__DASHBOARD_TITLE__", escaped_title).replace(
        "__DASHBOARD_PAYLOAD_JSON__", safe_json
    )

    out_path = Path(output_html)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered_html, encoding="utf-8")
    LOGGER.info("Generated standalone Foldseek dashboard at %s", out_path)
    return out_path


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for standalone dashboard compilation."""
    parser = argparse.ArgumentParser(
        description="Compile query 3D coordinates and optional Foldseek JSON results into a standalone HTML5 dashboard."
    )
    parser.add_argument(
        "--input_structure",
        required=True,
        help="Path to the query .pdb, .cif, or .mmcif coordinate file.",
    )
    parser.add_argument(
        "--results_json",
        required=False,
        default=None,
        help="Optional path to the Foldseek JSON results file.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional custom dashboard title.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Destination path for standalone dashboard HTML file.",
    )
    args = parser.parse_args(argv)

    try:
        build_dashboard(
            input_structure=args.input_structure,
            results_json=args.results_json,
            output_html=args.output,
            title=args.title,
        )
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
