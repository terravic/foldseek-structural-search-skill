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

"""Foldseek alignment normalization, coverage depth profiling, and functional clustering."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

# Support both package-relative and direct script execution
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import visualize_structure  # noqa: E402

MAX_ALIGNMENT_HITS: int = 300
DEFAULT_EVALUE: float = 1000.0

REGIME_HIGH_IDENTITY = "High-Identity Homolog"
REGIME_TWILIGHT_ZONE = "Remote Structural Homolog (Twilight Zone)"
REGIME_PARTIAL_MOTIF = "Local / Partial Motif Match"
REGIME_LOW_CONFIDENCE = "Low Confidence"

_UNIPROT_KV_RE = re.compile(r"\b(?:OS|OX|GN|PE|SV)=[^\s]+(?:\s+(?![A-Z]{2}=)[^\s]+)*")
_UNIPROT_OS_CAPTURE_RE = re.compile(r"\bOS=([^=]+?)(?:\s+(?:OX|GN|PE|SV)=|$)")
_BRACKET_ORG_RE = re.compile(r"\[([^\]]+)\]\s*$")
_PAREN_TRAIL_RE = re.compile(r"\(([^)]+)\)\s*$")
_CATH_CODE_RE = re.compile(r"\b\d+\.\d+\.\d+\.\d+\b")
_WORD_TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9\-]{2,}")

_ANNOTATION_STOPWORDS: frozenset[str] = frozenset(
    {
        "protein",
        "proteins",
        "domain",
        "domains",
        "subunit",
        "chain",
        "chains",
        "like",
        "putative",
        "uncharacterized",
        "hypothetical",
        "predicted",
        "synthetic",
        "construct",
        "model",
        "organism",
        "bacterium",
        "archaeon",
        "eukaryote",
        "metazoan",
        "fungus",
        "phage",
        "microbe",
        "reference",
        "design",
        "fragment",
        "segment",
        "module",
        "fold",
        "superfamily",
        "family",
        "related",
        "containing",
        "type",
        "class",
        "member",
        "isoform",
        "variant",
        "precursor",
        "probable",
        "similar",
        "the",
        "and",
        "for",
        "with",
        "from",
        "terminal",
        "core",
        "short",
    }
)


def extract_evalue(hit: dict[str, Any], default: float = DEFAULT_EVALUE) -> float:
    """Extract numeric E-value from a Foldseek hit dictionary."""
    for key in ("eval", "eValue", "evalue"):
        if key in hit and hit[key] is not None:
            try:
                val = float(hit[key])
                if math.isnan(val) or val < 0.0:
                    return default
                return val
            except (TypeError, ValueError):
                continue
    return default


def format_evalue(evalue: float) -> str:
    """Format E-value compactly for Markdown tables and UI badges."""
    if evalue == 0.0:
        return "0.0"
    if evalue < 1e-3 or evalue >= 1e4:
        return f"{evalue:.2e}"
    return f"{evalue:.4g}"


def compute_query_coverage(hit: dict[str, Any]) -> tuple[float, str]:
    """Compute bounded query coverage percentage and display string per Foldseek spec.

    Formula: `min(((qEndPos - qStartPos + 1) / qLen) * 100, 100.0)` when
    `qLen > 0` and `qEndPos > qStartPos`.
    """
    try:
        q_len = int(hit.get("qLen", hit.get("qlen", 0)) or 0)
        q_start = int(hit.get("qStartPos", 0) or 0)
        q_end = int(hit.get("qEndPos", 0) or 0)
    except (TypeError, ValueError):
        return (0.0, "N/A")

    if q_len > 0 and q_end > q_start:
        aligned_q_residues = (q_end - q_start) + 1
        cov_percentage = min((aligned_q_residues / q_len) * 100.0, 100.0)
        return (round(cov_percentage, 2), f"{cov_percentage:.1f}%")
    return (0.0, "N/A")


def normalize_seq_identity(hit: dict[str, Any]) -> tuple[float, str]:
    """Normalize fractional `[0, 1]` or percentage `[0, 100]` sequence identity."""
    raw_val: Any = None
    for key in ("seqId", "seqIdentity", "fident", "pident"):
        if key in hit and hit[key] is not None:
            raw_val = hit[key]
            break

    if raw_val is None:
        return (0.0, "0.0%")

    try:
        val = float(raw_val)
    except (TypeError, ValueError):
        return (0.0, "0.0%")

    if val < 0.0 or math.isnan(val):
        return (0.0, "0.0%")

    if val <= 1.0:
        # Distinguish exact 1.0 (100% identity) from 1.0% when mismatch metadata is present
        mismatches = hit.get("missmatches", hit.get("mismatches"))
        aln_len = hit.get("alnLength", hit.get("alnLen", hit.get("alnlen")))
        if val == 1.0 and mismatches is not None and aln_len:
            try:
                if int(mismatches) > 0 and int(aln_len) > 10:
                    pct = 1.0
                else:
                    pct = 100.0
            except (TypeError, ValueError):
                pct = 100.0
        else:
            pct = val * 100.0
    else:
        pct = min(val, 100.0)

    pct_rounded = round(pct, 2)
    return (pct_rounded, f"{pct_rounded:.1f}%")


def classify_homology_regime(
    prob: float,
    seq_identity_pct: float,
    q_cov_pct: float,
) -> str:
    """Classify a normalized structural match into one of four canonical homology regimes."""
    if prob < 0.5:
        return REGIME_LOW_CONFIDENCE
    if prob >= 0.8:
        if seq_identity_pct >= 30.0:
            return REGIME_HIGH_IDENTITY
        return REGIME_TWILIGHT_ZONE
    # 0.5 <= prob < 0.8 (including partial motif matches with q_cov_pct < 50.0)
    _ = q_cov_pct
    return REGIME_PARTIAL_MOTIF


def parse_target_annotation(target_id: str) -> tuple[str, str, str]:
    """Split a raw Foldseek `target` string into `(accession, description, organism)`."""
    raw = (target_id or "").strip()
    if not raw:
        return ("unknown", "Unannotated structural target", "")

    parts = raw.split(None, 1)
    raw_acc = parts[0]
    rest = parts[1].strip() if len(parts) > 1 else ""

    # Clean common model suffixes in AlphaFold accessions while preserving core ID
    accession = re.sub(r"-model_v\d+(?:\.(?:cif|pdb))?(?:\.gz)?$", "", raw_acc)
    accession = re.sub(r"\.(?:cif|pdb|ent)(?:\.gz)?$", "", accession)

    organism = ""
    os_match = _UNIPROT_OS_CAPTURE_RE.search(rest)
    if os_match:
        organism = os_match.group(1).strip()
        rest = _UNIPROT_KV_RE.sub("", rest).strip()
    else:
        br_match = _BRACKET_ORG_RE.search(rest)
        if br_match:
            organism = br_match.group(1).strip()
            rest = _BRACKET_ORG_RE.sub("", rest).strip()

    description = rest if rest else accession
    return (accession, description, organism)


def flatten_raw_alignments(raw_json: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    """Flatten Foldseek `results[].alignments[]` across both list and nested-list layouts."""
    collected: list[tuple[str, dict[str, Any]]] = []
    if not isinstance(raw_json, dict):
        return collected

    results = raw_json.get("results")
    if not isinstance(results, list):
        return collected

    for group in results:
        if not isinstance(group, dict):
            continue
        db_name = str(group.get("db") or "unknown").strip() or "unknown"
        alignments = group.get("alignments")
        if not alignments:
            continue
        if isinstance(alignments, dict):
            iterable = list(alignments.values())
        elif isinstance(alignments, list):
            iterable = alignments
        else:
            continue

        for item in iterable:
            if isinstance(item, dict):
                collected.append((db_name, item))
            elif isinstance(item, list):
                for sub_item in item:
                    if isinstance(sub_item, dict):
                        collected.append((db_name, sub_item))
    return collected


def extract_and_normalize_hits(
    raw_json: dict[str, Any],
    max_hits: int = MAX_ALIGNMENT_HITS,
) -> list[dict[str, Any]]:
    """Extract, sort, cap, and normalize Foldseek alignment hits.

    Args:
        raw_json: Parsed Foldseek JSON response dictionary.
        max_hits: Maximum hits to retain after ascending E-value sort (default 300).

    Returns:
        List of normalized hit dictionaries matching `dashboard_payload_schema.json`.
    """
    flattened = flatten_raw_alignments(raw_json)
    annotated_pairs = [(extract_evalue(hit), idx, db_name, hit) for idx, (db_name, hit) in enumerate(flattened)]
    annotated_pairs.sort(key=lambda item: (item[0], item[1]))
    capped = annotated_pairs[: max(0, max_hits)]

    normalized: list[dict[str, Any]] = []
    for rank_idx, (evalue, _orig_idx, group_db, hit) in enumerate(capped, start=1):
        target_id = str(hit.get("target") or hit.get("tId") or "unknown").strip()
        database = str(hit.get("db") or group_db or "unknown").strip()
        accession, description, organism = parse_target_annotation(target_id)

        try:
            prob_raw = float(hit.get("prob", 0.0) or 0.0)
        except (TypeError, ValueError):
            prob_raw = 0.0
        prob = min(max(prob_raw, 0.0), 1.0)

        q_cov_pct, q_cov_display = compute_query_coverage(hit)
        seq_id_pct, seq_id_display = normalize_seq_identity(hit)

        try:
            q_len = max(0, int(hit.get("qLen", hit.get("qlen", 0)) or 0))
        except (TypeError, ValueError):
            q_len = 0
        try:
            q_start = max(0, int(hit.get("qStartPos", 0) or 0))
        except (TypeError, ValueError):
            q_start = 0
        try:
            q_end = max(0, int(hit.get("qEndPos", 0) or 0))
        except (TypeError, ValueError):
            q_end = 0

        try:
            t_start = max(0, int(hit.get("dbStartPos", hit.get("tStartPos", 0)) or 0))
            t_end = max(0, int(hit.get("dbEndPos", hit.get("tEndPos", 0)) or 0))
            t_len = max(0, int(hit.get("dbLen", hit.get("tLen", 0)) or 0))
        except (TypeError, ValueError):
            t_start, t_end, t_len = 0, 0, 0

        try:
            aln_len = max(
                0,
                int(hit.get("alnLength", hit.get("alnLen", hit.get("alnlen", 0))) or 0),
            )
        except (TypeError, ValueError):
            aln_len = 0

        try:
            score = float(hit.get("score", hit.get("bits", 0.0)) or 0.0)
        except (TypeError, ValueError):
            score = 0.0

        neg_log10_eval = round(-math.log10(max(evalue, 1e-300)), 3)
        regime = classify_homology_regime(prob, seq_id_pct, q_cov_pct)
        target_url, target_fmt = visualize_structure.infer_target_structure_url(
            accession=accession,
            database=database,
            target_id=target_id,
        )

        normalized.append(
            {
                "rank": rank_idx,
                "target_id": target_id,
                "accession": accession,
                "description": description,
                "organism": organism,
                "database": database,
                "prob": round(prob, 4),
                "prob_display": f"{prob:.3f}",
                "evalue": evalue,
                "evalue_display": format_evalue(evalue),
                "neg_log10_evalue": neg_log10_eval,
                "seq_identity_pct": seq_id_pct,
                "seq_identity_display": seq_id_display,
                "q_start": q_start,
                "q_end": q_end,
                "q_len": q_len,
                "q_cov_pct": q_cov_pct,
                "q_cov_display": q_cov_display,
                "t_start": t_start,
                "t_end": t_end,
                "t_len": t_len,
                "aln_len": aln_len,
                "score": round(score, 2),
                "homology_regime": regime,
                "target_structure_url": target_url,
                "target_structure_format": target_fmt,
            }
        )

    return normalized


def compute_residue_coverage_profile(hits: list[dict[str, Any]], q_len: int) -> list[int]:
    """Compute per-residue structural hit coverage depth across query positions `1..q_len`."""
    if q_len <= 0:
        inferred = [int(h.get("q_len", 0) or 0) for h in hits]
        q_len = max(inferred, default=0)
    if q_len <= 0:
        return []

    depth = [0] * q_len
    for hit in hits:
        q_start = int(hit.get("q_start", 0) or 0)
        q_end = int(hit.get("q_end", 0) or 0)
        if q_start <= 0 or q_end < q_start:
            continue
        start_idx = max(0, q_start - 1)
        end_idx = min(q_len, q_end)
        for pos in range(start_idx, end_idx):
            depth[pos] += 1
    return depth


def summarize_database_breakdown(hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Aggregate hit counts and mean homology probabilities per database."""
    grouped: dict[str, list[float]] = defaultdict(list)
    for hit in hits:
        db = str(hit.get("database") or "unknown")
        grouped[db].append(float(hit.get("prob", 0.0)))

    summary = [
        {
            "database": db,
            "count": len(probs),
            "mean_prob": round(sum(probs) / len(probs), 3) if probs else 0.0,
        }
        for db, probs in grouped.items()
    ]
    summary.sort(key=lambda item: (-item["count"], -item["mean_prob"], item["database"]))
    return summary


def _extract_theme_phrase(description: str, accession: str) -> list[str]:
    """Extract informative functional keyword tokens/phrases from a cleaned target description."""
    if not description or description == accession:
        return []
    cleaned = _PAREN_TRAIL_RE.sub(" ", description)
    cleaned = _CATH_CODE_RE.sub(" ", cleaned)
    raw_tokens = _WORD_TOKEN_RE.findall(cleaned)

    meaningful: list[str] = []
    for tok in raw_tokens:
        lower = tok.lower().strip("-")
        if len(lower) < 3 or lower in _ANNOTATION_STOPWORDS:
            continue
        if lower.isdigit():
            continue
        meaningful.append(lower)

    candidates: list[str] = []
    for i in range(len(meaningful) - 1):
        candidates.append(f"{meaningful[i]} {meaningful[i + 1]}")
    candidates.extend(meaningful)
    return candidates


def cluster_functional_annotations(
    hits: list[dict[str, Any]],
    max_clusters: int = 8,
) -> list[dict[str, Any]]:
    """Group structural hits into recurring protein family/function themes."""
    if not hits:
        return []

    phrase_hits: dict[str, list[dict[str, Any]]] = defaultdict(list)
    unannotated: list[dict[str, Any]] = []

    for hit in hits:
        desc = str(hit.get("description") or "")
        acc = str(hit.get("accession") or "")
        phrases = _extract_theme_phrase(desc, acc)
        if not phrases:
            unannotated.append(hit)
            continue
        for phrase in set(phrases):
            phrase_hits[phrase].append(hit)

    # Prefer multi-word bigrams when they capture the same hits as component unigrams
    scored_phrases: list[tuple[float, str]] = []
    for phrase, member_list in phrase_hits.items():
        count = len(member_list)
        bigram_bonus = 1.25 if " " in phrase else 1.0
        mean_p = sum(float(h.get("prob", 0.0)) for h in member_list) / count
        scored_phrases.append((count * bigram_bonus + mean_p * 0.5, phrase))
    scored_phrases.sort(key=lambda item: (-item[0], item[1]))

    assigned_ranks: set[int] = set()
    clusters: list[dict[str, Any]] = []

    for _score, phrase in scored_phrases:
        if len(clusters) >= max_clusters:
            break
        unassigned_members = [h for h in phrase_hits[phrase] if h["rank"] not in assigned_ranks]
        if not unassigned_members:
            continue
        if len(unassigned_members) == 1 and len(hits) > 4 and len(clusters) >= 4:
            continue
        for member in unassigned_members:
            assigned_ranks.add(member["rank"])
        probs = [float(m.get("prob", 0.0)) for m in unassigned_members]
        label = " ".join(
            word.upper() if len(word) <= 3 and any(c.isdigit() for c in word) else word.capitalize()
            for word in phrase.split()
        )
        clusters.append(
            {
                "label": label,
                "count": len(unassigned_members),
                "mean_prob": round(sum(probs) / len(probs), 3),
                "top_accession": str(unassigned_members[0].get("accession") or ""),
            }
        )

    remaining = [h for h in hits if h["rank"] not in assigned_ranks]
    if remaining and len(clusters) < max_clusters:
        probs = [float(m.get("prob", 0.0)) for m in remaining]
        clusters.append(
            {
                "label": "Other / Divergent Annotations",
                "count": len(remaining),
                "mean_prob": round(sum(probs) / len(probs), 3),
                "top_accession": str(remaining[0].get("accession") or ""),
            }
        )

    clusters.sort(key=lambda c: (-c["count"], -c["mean_prob"], c["label"]))
    return clusters[:max_clusters]


def build_alignment_analytics(
    raw_json: dict[str, Any],
    fallback_query_len: int = 0,
    max_hits: int = MAX_ALIGNMENT_HITS,
) -> dict[str, Any]:
    """Build complete alignment analytics bundle from a Foldseek JSON response."""
    hits = extract_and_normalize_hits(raw_json, max_hits=max_hits)
    inferred_q_len = max(
        [fallback_query_len] + [int(h.get("q_len", 0) or 0) for h in hits]
    )
    coverage_profile = compute_residue_coverage_profile(hits, inferred_q_len)
    db_breakdown = summarize_database_breakdown(hits)
    clusters = cluster_functional_annotations(hits)

    # Also collect any searched DBs from raw_json["results"] that returned 0 hits
    searched_dbs: list[str] = []
    if isinstance(raw_json, dict) and isinstance(raw_json.get("results"), list):
        for group in raw_json["results"]:
            if isinstance(group, dict) and group.get("db"):
                db_str = str(group["db"]).strip()
                if db_str and db_str not in searched_dbs:
                    searched_dbs.append(db_str)
    for item in db_breakdown:
        if item["database"] not in searched_dbs:
            searched_dbs.append(item["database"])

    high_conf = sum(1 for h in hits if h["prob"] >= 0.90)
    twilight = sum(
        1 for h in hits if h["homology_regime"] == REGIME_TWILIGHT_ZONE
    )
    best_acc = hits[0]["accession"] if hits else "None"
    best_eval = hits[0]["evalue_display"] if hits else "N/A"

    return {
        "hits": hits,
        "coverage_profile": coverage_profile,
        "database_breakdown": db_breakdown,
        "functional_clusters": clusters,
        "summary_kpis": {
            "query_len": inferred_q_len,
            "total_hits": len(hits),
            "high_confidence_hits": high_conf,
            "twilight_zone_hits": twilight,
            "best_hit_accession": best_acc,
            "best_hit_evalue": best_eval,
            "databases_searched": searched_dbs,
        },
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint to convert a Foldseek JSON payload into normalized analytics JSON."""
    parser = argparse.ArgumentParser(
        description="Compute normalized Foldseek alignment analytics, coverage profiles, and clusters."
    )
    parser.add_argument("results_json", help="Path to Foldseek raw JSON result file.")
    parser.add_argument("-o", "--output", help="Optional output path for normalized analytics JSON.")
    args = parser.parse_args(argv)

    in_path = Path(args.results_json)
    if not in_path.is_file():
        print(f"Error: Foldseek results JSON not found at '{args.results_json}'.", file=sys.stderr)
        return 1

    try:
        raw_data = json.loads(in_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        print(f"Error: Invalid JSON in '{args.results_json}': {exc}", file=sys.stderr)
        return 1

    analytics = build_alignment_analytics(raw_data)
    serialized = json.dumps(analytics, indent=2)
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(serialized + "\n", encoding="utf-8")
    else:
        print(serialized)
    return 0


if __name__ == "__main__":
    sys.exit(main())
