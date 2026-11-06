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
#
# Modifications Notice (Apache License 2.0, Section 4(b)):
# Extended from the baseline foldseek_structural_search search.py script to add
# workspace license notification enforcement, deterministic --mock offline mode,
# --dashboard standalone HTML5 visualization generation, and modular alignment
# normalization while preserving backward-compatible CLI arguments and stdout
# Markdown formatting.
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "polite-http",
# ]
# ///

"""Core Foldseek 3Di+AA REST API client, license guard, and Markdown/JSON emitter."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import sys
import time
from typing import Any
import urllib.error
import urllib.request
import uuid

_SCRIPTS_DIR = Path(__file__).resolve().parent
_SKILL_ROOT = _SCRIPTS_DIR.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import generate_dashboard  # noqa: E402
import visualize_alignments  # noqa: E402
import visualize_structure  # noqa: E402

try:
    from polite_http import HttpError  # type: ignore[import-not-found]
    from polite_http import http_client  # type: ignore[import-not-found]
except ModuleNotFoundError:  # pragma: no cover - fallback for bare system Python without uv

    class HttpError(Exception):
        """Fallback HTTP error matching polite_http.HttpError."""

        def __init__(
            self,
            message: str,
            *,
            status_code: int | None = None,
            body: bytes | None = None,
            url: str | None = None,
        ) -> None:
            super().__init__(message)
            self.status_code = status_code
            self.body = body
            self.url = url

    class _FallbackHttpClient:
        """Rate-limited stdlib urllib fallback matching polite_http.HttpClient."""

        def __init__(self, base_url: str, qps: float = 1.0) -> None:
            self.base_url = base_url
            self.min_interval = 1.0 / max(qps, 0.01)
            self._last_call = 0.0

        def fetch_json(
            self,
            url: str,
            *,
            method: str = "GET",
            headers: dict[str, str] | None = None,
            data: bytes | None = None,
            timeout: float = 30.0,
        ) -> Any:
            elapsed = time.monotonic() - self._last_call
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self._last_call = time.monotonic()
            req_headers = {"User-Agent": "foldseek-structural-search/1.0"}
            if headers:
                req_headers.update(headers)
            req = urllib.request.Request(url, data=data, headers=req_headers, method=method)
            try:
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    raw = resp.read().decode("utf-8")
                    return json.loads(raw)
            except urllib.error.HTTPError as exc:
                err_body = exc.read()
                raise HttpError(
                    f"HTTP Error {exc.code} while fetching {url}",
                    status_code=exc.code,
                    body=err_body,
                    url=url,
                ) from exc
            except (urllib.error.URLError, OSError, json.JSONDecodeError) as exc:
                raise HttpError(f"Network error fetching {url}: {exc}", url=url) from exc

    class _FallbackModule:
        HttpClient = _FallbackHttpClient

    http_client = _FallbackModule()

LOGGER = logging.getLogger(__name__)

FOLDSEEK_API_BASE = "https://search.foldseek.com"
ALLOWED_DATABASES: list[str] = [
    "afdb50",
    "afdb-swissprot",
    "pdb100",
    "BFVD",
    "mgnify_esm30",
    "cath50",
    "gmgcl_id",
    "bfmd",
    "afdb-proteome",
]
ALLOWED_EXTENSIONS: set[str] = {".pdb", ".cif", ".mmcif"}
MAX_ALIGNMENT_HITS: int = 300
DEFAULT_EVALUE: float = 1000.0
DEFAULT_DATABASES: str = "pdb100,afdb50"

CLIENT = http_client.HttpClient(FOLDSEEK_API_BASE, qps=1)

LICENSE_REL_PATH = Path(".licenses") / "foldseek_structural_search_LICENSE.txt"
LICENSE_NOTICE_LINE = (
    "Notification: Please review the Foldseek terms of use and licensing at "
    "https://search.foldseek.com/search and https://github.com/steineggerlab/foldseek "
    "before using structural search results."
)


def resolve_workspace_root() -> Path:
    """Resolve workspace root directory for `.licenses/` storage without hardcoded paths."""
    env_override = os.environ.get("FOLDSEEK_WORKSPACE_ROOT", "").strip()
    if env_override:
        return Path(env_override).resolve()

    # Check if installed inside `<monorepo>/skills/foldseek_structural_search`
    monorepo_candidate = _SKILL_ROOT.parent.parent
    if (
        _SKILL_ROOT.parent.name == "skills"
        and (monorepo_candidate / "SKILL_LICENSES.md").is_file()
    ):
        return monorepo_candidate

    return _SKILL_ROOT


def ensure_license_notification(workspace_root: Path | None = None) -> Path:
    """Verify `.licenses/foldseek_structural_search_LICENSE.txt` exists; notify and create if absent."""
    root = workspace_root or resolve_workspace_root()
    license_file = root / LICENSE_REL_PATH
    if license_file.is_file():
        return license_file

    timestamp_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    content = f"{LICENSE_NOTICE_LINE}\nTimestamp: {timestamp_utc}\n"

    print(
        "\n[LICENSE NOTIFICATION] " + LICENSE_NOTICE_LINE + "\n",
        file=sys.stderr,
    )
    license_file.parent.mkdir(parents=True, exist_ok=True)
    license_file.write_text(content, encoding="utf-8")
    return license_file


def validate_databases(databases_arg: str) -> list[str]:
    """Validate comma-separated database names against ALLOWED_DATABASES."""
    selected = [db.strip() for db in (databases_arg or "").split(",") if db.strip()]
    if not selected:
        raise ValueError(
            "No target databases specified. Allowed databases: "
            + ", ".join(ALLOWED_DATABASES)
        )

    invalid = [db for db in selected if db not in ALLOWED_DATABASES]
    if invalid:
        raise ValueError(
            f"Unsupported Foldseek database(s): {', '.join(invalid)}. "
            f"Allowed databases are: {', '.join(ALLOWED_DATABASES)}"
        )

    # Deduplicate while preserving user order
    deduped: list[str] = []
    for db in selected:
        if db not in deduped:
            deduped.append(db)
    return deduped


def build_multipart_payload(
    fields: dict[str, str | list[str]],
    files: dict[str, str | Path],
) -> tuple[bytes, str]:
    """Construct an RFC 7578 `multipart/form-data` byte body and Content-Type header."""
    boundary = f"----FoldseekSkillBoundary{uuid.uuid4().hex}"
    chunks: list[bytes] = []

    for field_name, raw_val in fields.items():
        values = raw_val if isinstance(raw_val, list) else [raw_val]
        for val in values:
            chunks.append(f"--{boundary}\r\n".encode("utf-8"))
            chunks.append(
                f'Content-Disposition: form-data; name="{field_name}"\r\n\r\n'.encode(
                    "utf-8"
                )
            )
            chunks.append(f"{val}\r\n".encode("utf-8"))

    for file_field, fpath in files.items():
        path_obj = Path(fpath)
        file_bytes = path_obj.read_bytes()
        safe_filename = path_obj.name.replace('"', "_")
        chunks.append(f"--{boundary}\r\n".encode("utf-8"))
        chunks.append(
            (
                f'Content-Disposition: form-data; name="{file_field}"; '
                f'filename="{safe_filename}"\r\n'
                "Content-Type: application/octet-stream\r\n\r\n"
            ).encode("utf-8")
        )
        chunks.append(file_bytes)
        chunks.append(b"\r\n")

    chunks.append(f"--{boundary}--\r\n".encode("utf-8"))
    content_type = f"multipart/form-data; boundary={boundary}"
    return (b"".join(chunks), content_type)


def load_mock_results(selected_dbs: list[str]) -> dict[str, Any]:
    """Load deterministic offline synthetic Foldseek results for `--mock` execution."""
    mock_path = _SKILL_ROOT / "sample_data" / "sample_foldseek_results.json"
    if not mock_path.is_file():
        raise FileNotFoundError(
            f"Synthetic mock results fixture missing at '{mock_path}'."
        )
    raw: dict[str, Any] = json.loads(mock_path.read_text(encoding="utf-8"))
    all_groups = raw.get("results", [])
    if isinstance(all_groups, list) and selected_dbs:
        filtered_groups = [
            grp for grp in all_groups if isinstance(grp, dict) and grp.get("db") in selected_dbs
        ]
        if filtered_groups:
            return {"query": raw.get("query", {}), "results": filtered_groups}
    return raw


def submit_and_poll_foldseek(
    structure_path: Path,
    selected_dbs: list[str],
    *,
    client: Any = None,
    poll_interval_secs: float = 10.0,
    finalization_wait_secs: float = 3.0,
    max_poll_secs: float | None = None,
) -> dict[str, Any]:
    """Submit structure file to Foldseek API in `3diaa` mode, poll ticket, and return JSON."""
    http = client or CLIENT
    effective_max_poll = (
        max_poll_secs
        if max_poll_secs is not None
        else float(os.environ.get("FOLDSEEK_MAX_POLL_SECONDS", "600"))
    )

    body_bytes, content_type = build_multipart_payload(
        fields={"mode": "3diaa", "database[]": selected_dbs},
        files={"q": structure_path},
    )

    ticket_url = f"{FOLDSEEK_API_BASE}/api/ticket"
    ticket_resp = http.fetch_json(
        ticket_url,
        method="POST",
        headers={"Content-Type": content_type},
        data=body_bytes,
        timeout=30,
    )
    if not isinstance(ticket_resp, dict) or not ticket_resp.get("id"):
        raise RuntimeError(
            f"Unexpected Foldseek ticket response payload: {ticket_resp!r}"
        )

    ticket_id = str(ticket_resp["id"])
    status_url = f"{FOLDSEEK_API_BASE}/api/ticket/{ticket_id}"
    started = time.monotonic()

    while True:
        status_resp = http.fetch_json(status_url, timeout=20)
        status = (
            str(status_resp.get("status", "")).upper()
            if isinstance(status_resp, dict)
            else ""
        )
        if status == "COMPLETE":
            break
        if status == "ERROR":
            raise RuntimeError(
                f"Foldseek server returned ERROR status for ticket '{ticket_id}'."
            )
        if (time.monotonic() - started) >= effective_max_poll:
            raise TimeoutError(
                f"Timed out waiting for Foldseek ticket '{ticket_id}' after {effective_max_poll:.0f}s."
            )
        if poll_interval_secs > 0:
            time.sleep(poll_interval_secs)

    if finalization_wait_secs > 0:
        time.sleep(finalization_wait_secs)

    result_url = f"{FOLDSEEK_API_BASE}/api/result/{ticket_id}/0"
    result_payload = http.fetch_json(result_url, timeout=120)
    if not isinstance(result_payload, dict):
        raise RuntimeError("Foldseek result endpoint did not return a JSON object.")
    return result_payload


def format_markdown_matches_table(
    raw_json: dict[str, Any],
    max_hits: int = MAX_ALIGNMENT_HITS,
    include_citation: bool = True,
) -> str:
    """Format Foldseek JSON payload into the canonical Markdown alignment table."""
    hits = visualize_alignments.extract_and_normalize_hits(raw_json, max_hits=max_hits)
    lines = [
        "### All Structural Matches",
        "| Target ID | Q-Cov | Prob | E-value | Seq Identity | Aln Length |",
        "|---|---|---|---|---|---|",
    ]
    for hit in hits:
        safe_target = str(hit["target_id"]).replace("|", "\\|")
        lines.append(
            f"| {safe_target} | {hit['q_cov_display']} | {hit['prob_display']} | "
            f"{hit['evalue_display']} | {hit['seq_identity_display']} | {hit['aln_len']} |"
        )

    if include_citation:
        lines.extend(
            [
                "",
                "---",
                f"**Citation:** {generate_dashboard.FOLDSEEK_CITATION_TEXT}",
            ]
        )
    return "\n".join(lines) + "\n"


def execute_search(
    input_file: str | Path,
    output_json: str | Path = "foldseek_results.json",
    databases: str = DEFAULT_DATABASES,
    dashboard_path: str | Path | None = None,
    mock: bool = False,
    *,
    client: Any = None,
    workspace_root: Path | None = None,
) -> tuple[dict[str, Any], str]:
    """Validate inputs, run Foldseek search (or mock), write outputs, and return `(json_dict, markdown)`."""
    validated_structure = visualize_structure.validate_structure_file(input_file)
    # Also verify non-zero coordinate records before submitting network requests
    visualize_structure.parse_structure_metadata(validated_structure)
    selected_dbs = validate_databases(databases)

    ensure_license_notification(workspace_root=workspace_root)

    if mock:
        raw_results = load_mock_results(selected_dbs)
    else:
        raw_results = submit_and_poll_foldseek(
            validated_structure,
            selected_dbs,
            client=client,
        )

    out_json_path = Path(output_json)
    out_json_path.parent.mkdir(parents=True, exist_ok=True)
    out_json_path.write_text(json.dumps(raw_results, indent=2) + "\n", encoding="utf-8")

    if dashboard_path:
        generate_dashboard.build_dashboard(
            input_structure=validated_structure,
            results_json=raw_results,
            output_html=dashboard_path,
        )

    markdown_report = format_markdown_matches_table(
        raw_results,
        max_hits=MAX_ALIGNMENT_HITS,
        include_citation=True,
    )
    return (raw_results, markdown_report)


def build_arg_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser compatible with baseline `scripts/search.py`."""
    parser = argparse.ArgumentParser(
        description=(
            "Submit a 3D protein coordinate file (.pdb, .cif, or .mmcif) to the "
            "Foldseek REST API in 3diaa mode and output sorted structural alignments."
        )
    )
    parser.add_argument(
        "input_file",
        help="Path to physical 3D protein structure file (.pdb, .cif, or .mmcif).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="foldseek_results.json",
        help="Path to save the raw Foldseek JSON payload (default: foldseek_results.json).",
    )
    parser.add_argument(
        "--databases",
        default=DEFAULT_DATABASES,
        help=(
            f"Comma-separated Foldseek databases (default: {DEFAULT_DATABASES}). "
            f"Allowed: {', '.join(ALLOWED_DATABASES)}"
        ),
    )
    parser.add_argument(
        "--dashboard",
        default=None,
        help="Optional path to generate a standalone interactive HTML5 dashboard.",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Run offline using deterministic synthetic Foldseek alignment data.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for Foldseek structural search."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        _raw_json, markdown_table = execute_search(
            input_file=args.input_file,
            output_json=args.output,
            databases=args.databases,
            dashboard_path=args.dashboard,
            mock=args.mock,
        )
    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        TimeoutError,
        HttpError,
        OSError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    sys.stdout.write(markdown_table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
