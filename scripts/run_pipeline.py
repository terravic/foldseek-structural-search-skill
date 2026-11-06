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
# dependencies = [
#   "polite-http",
# ]
# ///

"""Unified end-to-end Foldseek structural search, Markdown synthesis, and dashboard pipeline."""

from __future__ import annotations

import argparse
from functools import partial
import http.server
import json
from pathlib import Path
import socketserver
import sys
from typing import Any

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import generate_dashboard  # noqa: E402
import search  # noqa: E402
import visualize_alignments  # noqa: E402
import visualize_structure  # noqa: E402


def compose_synthesis_markdown(
    query_meta: dict[str, Any],
    analytics: dict[str, Any],
    raw_json: dict[str, Any],
) -> str:
    """Compose a comprehensive Markdown structural homology report and full alignment table."""
    kpis = analytics["summary_kpis"]
    hits: list[dict[str, Any]] = analytics["hits"]
    clusters: list[dict[str, Any]] = analytics["functional_clusters"]

    lines: list[str] = [
        f"# Foldseek 3D Structural Search Report: `{query_meta['filename']}`",
        "",
        "## 1. Executive Summary",
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Query Coordinate File | `{query_meta['filename']}` (`{query_meta['extension']}`) |",
        f"| Detected Chains / Length | {', '.join(query_meta['chains']) or 'A'} ({kpis['query_len']} residues, {query_meta['atom_count']} atoms) |",
        f"| Mean Query B-factor / pLDDT | {query_meta['mean_bfactor']:.2f} (range {query_meta['min_bfactor']:.2f} - {query_meta['max_bfactor']:.2f}) |",
        f"| Searched Databases | {', '.join(kpis['databases_searched']) or 'N/A'} |",
        f"| Total Structural Hits | {kpis['total_hits']} |",
        f"| High-Confidence Homologs (`Prob >= 0.90`) | {kpis['high_confidence_hits']} |",
        f"| Remote Twilight-Zone Homologs (`Prob >= 0.80`, `Seq ID < 30%`) | {kpis['twilight_zone_hits']} |",
        f"| Best Hit (`Rank 1`) | `{kpis['best_hit_accession']}` (E-value `{kpis['best_hit_evalue']}`) |",
        "",
        "## 2. Top Structural Homologues",
        "",
    ]

    top_hits = hits[:5]
    if not top_hits:
        lines.append("No significant structural alignments were returned for this query.")
    else:
        for h in top_hits:
            lines.append(
                f"- **Rank {h['rank']}: `{h['accession']}`** (`{h['database']}`) - "
                f"{h['description']}. "
                f"`Prob = {h['prob_display']}`, `Q-Cov = {h['q_cov_display']}` "
                f"(residues {h['q_start']}-{h['q_end']}), "
                f"`E-value = {h['evalue_display']}`, `Seq Identity = {h['seq_identity_display']}`, "
                f"Regime: *{h['homology_regime']}*."
            )

    twilight_hits = [
        h
        for h in hits
        if h["homology_regime"] == visualize_alignments.REGIME_TWILIGHT_ZONE
    ][:5]
    lines.extend(["", "## 3. Remote Twilight-Zone Structural Homologues", ""])
    if not twilight_hits:
        lines.append(
            "No remote twilight-zone homologues (`Prob >= 0.80` with `Seq Identity < 30%`) detected."
        )
    else:
        for h in twilight_hits:
            lines.append(
                f"- **`{h['accession']}`** (`{h['database']}`): {h['description']} - "
                f"`Prob = {h['prob_display']}`, `Seq Identity = {h['seq_identity_display']}`, "
                f"`Q-Cov = {h['q_cov_display']}`, `E-value = {h['evalue_display']}`."
            )

    lines.extend(["", "## 4. Functional Family & Domain Cluster Breakdown", ""])
    if not clusters:
        lines.append("Insufficient target descriptions to derive recurring functional clusters.")
    else:
        lines.extend(
            [
                "| Functional Theme / Protein Family | Hit Count | Mean Prob | Representative Target |",
                "|---|---|---|---|",
            ]
        )
        for cl in clusters:
            lines.append(
                f"| {cl['label']} | {cl['count']} | {cl['mean_prob']:.3f} | `{cl.get('top_accession', '')}` |"
            )

    lines.extend(
        [
            "",
            search.format_markdown_matches_table(
                raw_json,
                max_hits=search.MAX_ALIGNMENT_HITS,
                include_citation=True,
            ).rstrip("\n"),
            "",
        ]
    )
    return "\n".join(lines)


def run_pipeline(
    input_file: str | Path,
    output_dir: str | Path,
    databases: str = search.DEFAULT_DATABASES,
    mock: bool = False,
    title: str | None = None,
    *,
    workspace_root: Path | None = None,
) -> dict[str, Path]:
    """Execute full Foldseek search, Markdown synthesis, and dashboard generation pipeline."""
    validated_structure = visualize_structure.validate_structure_file(input_file)
    query_meta = visualize_structure.parse_structure_metadata(validated_structure)
    selected_dbs = search.validate_databases(databases)

    search.ensure_license_notification(workspace_root=workspace_root)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = validated_structure.stem

    json_path = out_dir / f"{stem}_foldseek_results.json"
    md_path = out_dir / f"{stem}_foldseek_results.md"
    dashboard_path = out_dir / "dashboard.html"

    if mock:
        raw_results = search.load_mock_results(selected_dbs)
    else:
        raw_results = search.submit_and_poll_foldseek(validated_structure, selected_dbs)

    json_path.write_text(json.dumps(raw_results, indent=2) + "\n", encoding="utf-8")

    analytics = visualize_alignments.build_alignment_analytics(
        raw_json=raw_results,
        fallback_query_len=query_meta["residue_count"],
    )
    md_content = compose_synthesis_markdown(query_meta, analytics, raw_results)
    md_path.write_text(md_content, encoding="utf-8")

    generate_dashboard.build_dashboard(
        input_structure=validated_structure,
        results_json=raw_results,
        output_html=dashboard_path,
        title=title,
    )

    return {
        "json": json_path,
        "markdown": md_path,
        "dashboard": dashboard_path,
    }


def _serve_directory(directory: Path, port: int) -> None:
    """Serve the output directory over local HTTP for browser inspection."""
    handler = partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
    with socketserver.TCPServer(("127.0.0.1", port), handler) as httpd:
        print(f"Serving dashboard at http://127.0.0.1:{port}/dashboard.html (Ctrl+C to stop)")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


def build_arg_parser() -> argparse.ArgumentParser:
    """Construct CLI parser for `run_pipeline.py`."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the end-to-end Foldseek 3D structural search pipeline, emitting JSON, "
            "Markdown functional synthesis, and standalone HTML5 dashboard artifacts."
        )
    )
    parser.add_argument(
        "input_file",
        help="Path to query 3D protein coordinate file (.pdb, .cif, or .mmcif).",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        "--output-dir",
        dest="output_dir",
        required=True,
        help="Output directory where JSON, Markdown, and dashboard.html will be saved.",
    )
    parser.add_argument(
        "--databases",
        default=search.DEFAULT_DATABASES,
        help=(
            f"Comma-separated Foldseek databases (default: {search.DEFAULT_DATABASES}). "
            f"Allowed: {', '.join(search.ALLOWED_DATABASES)}"
        ),
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Execute deterministically offline using synthetic Foldseek results.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional custom title for the generated HTML5 dashboard.",
    )
    parser.add_argument(
        "--serve",
        action="store_true",
        help="Start a local HTTP server in --output_dir after generating artifacts.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8765,
        help="Local HTTP port when --serve is enabled (default: 8765).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for unified Foldseek pipeline."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        outputs = run_pipeline(
            input_file=args.input_file,
            output_dir=args.output_dir,
            databases=args.databases,
            mock=args.mock,
            title=args.title,
        )
    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        TimeoutError,
        search.HttpError,
        OSError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print("Foldseek structural search pipeline completed.")
    print(f"  JSON payload:      {outputs['json']}")
    print(f"  Markdown report:   {outputs['markdown']}")
    print(f"  HTML5 dashboard:   {outputs['dashboard']}")
    print(f"  Citation:          {generate_dashboard.FOLDSEEK_CITATION_TEXT}")

    if args.serve:
        _serve_directory(Path(args.output_dir).resolve(), args.port)
    return 0


if __name__ == "__main__":
    sys.exit(main())
