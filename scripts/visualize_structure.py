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

"""3D macromolecular coordinate parser and structural payload builder."""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import re
import shlex
import sys
from typing import Any

LOGGER = logging.getLogger(__name__)

ALLOWED_EXTENSIONS: frozenset[str] = frozenset({".pdb", ".cif", ".mmcif"})
MAX_EMBEDDED_ATOMS: int = 3500

_UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F\d+)?",
    re.IGNORECASE,
)
_PDB_CHAIN_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})(?:[_\-:.][A-Za-z0-9]+)?$")
_CATH_PDB_RE = re.compile(r"cath\|[^|]+\|([0-9][A-Za-z0-9]{3})", re.IGNORECASE)
_RAW_SEQUENCE_HINT_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{12,}$", re.IGNORECASE)


def validate_structure_file(file_path: str | Path) -> Path:
    """Validate that input points to an existing .pdb, .cif, or .mmcif coordinate file."""
    raw = str(file_path).strip()
    if not raw:
        raise ValueError(
            "Empty input path provided. Foldseek strictly requires a local .pdb, .cif, or .mmcif file."
        )

    candidate = Path(raw)
    suffix = candidate.suffix.lower()

    if not candidate.exists():
        if not suffix:
            if _RAW_SEQUENCE_HINT_RE.match(raw):
                raise ValueError(
                    "Raw amino acid sequence detected instead of a file path. "
                    "Foldseek strictly requires a physical 3D coordinate file (.pdb, .cif, or .mmcif)."
                )
            raise ValueError(
                f"Input '{raw}' is not an existing file and has no coordinate extension. "
                "Foldseek cannot search by protein name, gene symbol, or UniProt/PDB accession ID; "
                "provide a local .pdb, .cif, or .mmcif file."
            )
        if suffix not in ALLOWED_EXTENSIONS:
            allowed_str = ", ".join(sorted(ALLOWED_EXTENSIONS))
            raise ValueError(
                f"Unsupported file extension '{suffix}' for '{raw}'. "
                f"Allowed 3D coordinate extensions are: {allowed_str}."
            )
        raise FileNotFoundError(
            f"3D coordinate file not found at '{raw}'. Verify the file path in your workspace."
        )

    if not candidate.is_file():
        raise ValueError(
            f"Input path '{raw}' is a directory or special device, not a regular coordinate file."
        )

    if suffix not in ALLOWED_EXTENSIONS:
        allowed_str = ", ".join(sorted(ALLOWED_EXTENSIONS))
        raise ValueError(
            f"Unsupported file extension '{suffix}' for '{candidate.name}'. "
            f"Foldseek requires a physical 3D structure file ({allowed_str}); "
            "FASTA sequences and text annotations are not supported."
        )

    if candidate.stat().st_size == 0:
        raise ValueError(f"Coordinate file '{raw}' is empty (0 bytes).")

    return candidate


def _dist(a: dict[str, Any], b: dict[str, Any]) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _infer_ca_secondary_structure(ca_trace: list[dict[str, Any]]) -> dict[tuple[str, int], str]:
    """Infer H (helix) / E (sheet) / C (coil) from C-alpha inter-residue distances when headers are absent."""
    n = len(ca_trace)
    raw_labels = ["C"] * n
    for i in range(n):
        if i + 3 < n and ca_trace[i]["chain"] == ca_trace[i + 3]["chain"]:
            d13 = _dist(ca_trace[i], ca_trace[i + 2])
            d14 = _dist(ca_trace[i], ca_trace[i + 3])
            if 4.7 <= d14 <= 6.2 and 5.0 <= d13 <= 6.0:
                for k in range(i, min(n, i + 4)):
                    raw_labels[k] = "H"
                continue
            if 6.1 <= d13 <= 7.4 and 9.0 <= d14 <= 11.4:
                for k in range(i, min(n, i + 3)):
                    if raw_labels[k] == "C":
                        raw_labels[k] = "E"

    # Smooth isolated single-residue assignments
    for i in range(1, n - 1):
        if raw_labels[i] != "C" and raw_labels[i - 1] == "C" and raw_labels[i + 1] == "C":
            raw_labels[i] = "C"

    return {(pt["chain"], pt["res"]): raw_labels[idx] for idx, pt in enumerate(ca_trace)}


def _parse_pdb_coordinates(text: str) -> dict[str, Any]:
    """Extract structural summary, secondary structure, C-alpha trace, and heavy atoms from PDB text."""
    chains: list[str] = []
    seen_chains: set[str] = set()
    residues: set[tuple[str, int, str]] = set()
    res_numbers: list[int] = []
    bfactors: list[float] = []
    ss_map: dict[tuple[str, int], str] = {}
    ca_by_res: dict[tuple[str, int], dict[str, Any]] = {}
    o_by_res: dict[tuple[str, int], tuple[float, float, float]] = {}
    fallback_by_res: dict[tuple[str, int], dict[str, Any]] = {}
    atoms_list: list[dict[str, Any]] = []
    atom_count = 0
    in_first_model = True
    saw_atoms_in_model = False

    for line in text.splitlines():
        record = line[:6].strip().upper()
        if record == "HELIX" and len(line) >= 37:
            h_chain = line[19:20].strip() or "A"
            try:
                s_res = int(line[21:25].strip())
                e_res = int(line[33:37].strip())
                for r in range(s_res, e_res + 1):
                    ss_map[(h_chain, r)] = "H"
            except ValueError:
                pass
            continue
        if record == "SHEET" and len(line) >= 37:
            s_chain = line[21:22].strip() or "A"
            try:
                s_res = int(line[22:26].strip())
                e_res = int(line[33:37].strip())
                for r in range(s_res, e_res + 1):
                    ss_map.setdefault((s_chain, r), "E")
            except ValueError:
                pass
            continue
        if record == "ENDMDL":
            if saw_atoms_in_model:
                in_first_model = False
            continue
        if not in_first_model or record != "ATOM" or len(line) < 54:
            continue

        alt_loc = line[16:17].strip()
        if alt_loc not in ("", "A", "1"):
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue

        atom_name = line[12:16].strip()
        elem = (line[76:78].strip() if len(line) >= 78 else "") or atom_name[:1]
        if elem.upper() == "H":
            continue

        resn = line[17:20].strip() or "UNK"
        chain = line[21:22].strip() or "A"
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue
        ins_code = line[26:27].strip()
        bfac = 0.0
        if len(line) >= 66:
            try:
                bfac = float(line[60:66].strip())
            except ValueError:
                bfac = 0.0

        saw_atoms_in_model = True
        atom_count += 1
        if chain not in seen_chains:
            seen_chains.add(chain)
            chains.append(chain)
        residues.add((chain, res_seq, ins_code))
        res_numbers.append(res_seq)
        bfactors.append(bfac)

        key = (chain, res_seq)
        rx, ry, rz, rb = round(x, 3), round(y, 3), round(z, 3), round(bfac, 2)
        if atom_name == "O" and key not in o_by_res:
            o_by_res[key] = (rx, ry, rz)

        atom_entry = {
            "res": res_seq,
            "chain": chain,
            "resn": resn,
            "x": rx,
            "y": ry,
            "z": rz,
            "b": rb,
        }
        if atom_name in ("CA", "P") and key not in ca_by_res:
            ca_by_res[key] = atom_entry
        elif key not in fallback_by_res:
            fallback_by_res[key] = atom_entry

        if len(atoms_list) < MAX_EMBEDDED_ATOMS or atom_name in ("N", "CA", "C", "O", "CB"):
            if len(atoms_list) < MAX_EMBEDDED_ATOMS:
                atoms_list.append(
                    {
                        "n": atom_name,
                        "rn": resn,
                        "c": chain,
                        "r": res_seq,
                        "x": rx,
                        "y": ry,
                        "z": rz,
                        "b": rb,
                        "e": elem.upper(),
                    }
                )

    ordered_keys = sorted(
        set(ca_by_res.keys()) | set(fallback_by_res.keys()), key=lambda k: (k[0], k[1])
    )
    ca_trace = [
        dict(ca_by_res[k] if k in ca_by_res else fallback_by_res[k]) for k in ordered_keys
    ]

    if not ss_map and ca_trace:
        ss_map = _infer_ca_secondary_structure(ca_trace)

    for pt in ca_trace:
        key = (pt["chain"], pt["res"])
        pt["ss"] = ss_map.get(key, "C")
        if key in o_by_res:
            pt["ox"], pt["oy"], pt["oz"] = o_by_res[key]

    for at in atoms_list:
        at["ss"] = ss_map.get((at["c"], at["r"]), "C")

    return {
        "chains": chains,
        "residue_count": len(residues),
        "atom_count": atom_count,
        "res_numbers": res_numbers,
        "bfactors": bfactors,
        "ca_trace": ca_trace,
        "atoms": atoms_list,
    }


def _tokenize_cif_row(line: str) -> list[str]:
    """Tokenize a single mmCIF data row while respecting single/double quotes."""
    try:
        return shlex.split(line, comments=False, posix=True)
    except ValueError:
        return line.split()


def _parse_cif_coordinates(text: str) -> dict[str, Any]:
    """Extract structural summary, secondary structure, C-alpha trace, and heavy atoms from mmCIF."""
    lines = text.splitlines()
    n_lines = len(lines)
    idx = 0

    chains: list[str] = []
    seen_chains: set[str] = set()
    residues: set[tuple[str, int]] = set()
    res_numbers: list[int] = []
    bfactors: list[float] = []
    ss_map: dict[tuple[str, int], str] = {}
    ca_by_res: dict[tuple[str, int], dict[str, Any]] = {}
    o_by_res: dict[tuple[str, int], tuple[float, float, float]] = {}
    fallback_by_res: dict[tuple[str, int], dict[str, Any]] = {}
    atoms_list: list[dict[str, Any]] = []
    atom_count = 0
    first_model_num: str | None = None

    while idx < n_lines:
        stripped = lines[idx].strip()
        if stripped == "loop_":
            idx += 1
            headers: list[str] = []
            while idx < n_lines:
                hline = lines[idx].strip()
                if not hline or hline.startswith("#"):
                    idx += 1
                    continue
                if hline.startswith("_"):
                    headers.append(hline.split()[0])
                    idx += 1
                else:
                    break

            if not headers:
                continue

            col_map = {h: pos for pos, h in enumerate(headers)}
            n_cols = len(headers)

            def _first_col(*keys: str) -> int | None:
                for k in keys:
                    if k in col_map:
                        return col_map[k]
                return None

            # Parse secondary structure loops (_struct_conf / _struct_sheet_range)
            if any(h.startswith("_struct_conf.") for h in headers) or any(
                h.startswith("_struct_sheet_range.") for h in headers
            ):
                is_helix_loop = any(h.startswith("_struct_conf.") for h in headers)
                type_col = col_map.get("_struct_conf.conf_type_id")
                c_col = _first_col(
                    "_struct_conf.beg_auth_asym_id",
                    "_struct_conf.beg_label_asym_id",
                    "_struct_sheet_range.beg_auth_asym_id",
                    "_struct_sheet_range.beg_label_asym_id",
                )
                s_col = _first_col(
                    "_struct_conf.beg_label_seq_id",
                    "_struct_conf.beg_auth_seq_id",
                    "_struct_sheet_range.beg_label_seq_id",
                    "_struct_sheet_range.beg_auth_seq_id",
                )
                e_col = _first_col(
                    "_struct_conf.end_label_seq_id",
                    "_struct_conf.end_auth_seq_id",
                    "_struct_sheet_range.end_label_seq_id",
                    "_struct_sheet_range.end_auth_seq_id",
                )
                row_tokens: list[str] = []
                while idx < n_lines:
                    rline = lines[idx].strip()
                    if not rline:
                        idx += 1
                        continue
                    if (
                        rline.startswith("#")
                        or rline == "loop_"
                        or rline.startswith("_")
                        or rline.startswith("data_")
                    ):
                        break
                    row_tokens.extend(_tokenize_cif_row(rline))
                    idx += 1
                    while len(row_tokens) >= n_cols:
                        toks = row_tokens[:n_cols]
                        row_tokens = row_tokens[n_cols:]
                        if is_helix_loop and type_col is not None:
                            if not toks[type_col].upper().startswith("HELX"):
                                continue
                        ch = toks[c_col] if c_col is not None and toks[c_col] not in (".", "?") else "A"
                        try:
                            sr = int(toks[s_col]) if s_col is not None else 0
                            er = int(toks[e_col]) if e_col is not None else 0
                            for r in range(sr, er + 1):
                                ss_map[(ch, r)] = "H" if is_helix_loop else "E"
                        except ValueError:
                            pass
                continue

            if any(h.startswith("_atom_site.") for h in headers):
                group_col = col_map.get("_atom_site.group_PDB")
                elem_col = col_map.get("_atom_site.type_symbol")
                atom_col = _first_col("_atom_site.label_atom_id", "_atom_site.auth_atom_id")
                comp_col = _first_col("_atom_site.label_comp_id", "_atom_site.auth_comp_id")
                chain_col = (
                    col_map.get("_atom_site.auth_asym_id")
                    if "_atom_site.auth_asym_id" in col_map
                    else col_map.get("_atom_site.label_asym_id")
                )
                seq_col = (
                    col_map.get("_atom_site.label_seq_id")
                    if "_atom_site.label_seq_id" in col_map
                    else col_map.get("_atom_site.auth_seq_id")
                )
                auth_seq_col = col_map.get("_atom_site.auth_seq_id")
                alt_col = col_map.get("_atom_site.label_alt_id")
                x_col = col_map.get("_atom_site.Cartn_x")
                y_col = col_map.get("_atom_site.Cartn_y")
                z_col = col_map.get("_atom_site.Cartn_z")
                b_col = col_map.get("_atom_site.B_iso_or_equiv")
                model_col = col_map.get("_atom_site.pdbx_PDB_model_num")

                row_tokens = []
                while idx < n_lines:
                    rline = lines[idx].strip()
                    if not rline:
                        idx += 1
                        continue
                    if (
                        rline.startswith("#")
                        or rline == "loop_"
                        or rline.startswith("_")
                        or rline.startswith("data_")
                    ):
                        break
                    if rline.startswith(";"):
                        idx += 1
                        while idx < n_lines and not lines[idx].startswith(";"):
                            idx += 1
                        idx += 1
                        continue
                    row_tokens.extend(_tokenize_cif_row(rline))
                    idx += 1
                    while len(row_tokens) >= n_cols:
                        tokens = row_tokens[:n_cols]
                        row_tokens = row_tokens[n_cols:]
                        if group_col is not None and tokens[group_col].upper() != "ATOM":
                            continue
                        if model_col is not None:
                            mval = tokens[model_col]
                            if first_model_num is None:
                                first_model_num = mval
                            elif mval != first_model_num:
                                continue
                        if alt_col is not None and tokens[alt_col] not in (".", "?", "A", "1"):
                            continue
                        if x_col is None or y_col is None or z_col is None:
                            continue
                        try:
                            x = float(tokens[x_col])
                            y = float(tokens[y_col])
                            z = float(tokens[z_col])
                        except ValueError:
                            continue

                        raw_seq = (
                            tokens[seq_col]
                            if seq_col is not None and tokens[seq_col] not in (".", "?")
                            else (tokens[auth_seq_col] if auth_seq_col is not None else "0")
                        )
                        try:
                            res_seq = int(raw_seq)
                        except ValueError:
                            continue

                        atom_name = (
                            tokens[atom_col]
                            if atom_col is not None and tokens[atom_col] not in (".", "?")
                            else "CA"
                        )
                        elem = (
                            tokens[elem_col]
                            if elem_col is not None and tokens[elem_col] not in (".", "?")
                            else atom_name[:1]
                        )
                        if elem.upper() == "H":
                            continue

                        chain = (
                            tokens[chain_col]
                            if chain_col is not None and tokens[chain_col] not in (".", "?")
                            else "A"
                        )
                        resn = (
                            tokens[comp_col]
                            if comp_col is not None and tokens[comp_col] not in (".", "?")
                            else "UNK"
                        )
                        bfac = 0.0
                        if b_col is not None and tokens[b_col] not in (".", "?"):
                            try:
                                bfac = float(tokens[b_col])
                            except ValueError:
                                bfac = 0.0

                        atom_count += 1
                        if chain not in seen_chains:
                            seen_chains.add(chain)
                            chains.append(chain)
                        residues.add((chain, res_seq))
                        res_numbers.append(res_seq)
                        bfactors.append(bfac)

                        key = (chain, res_seq)
                        rx, ry, rz, rb = round(x, 3), round(y, 3), round(z, 3), round(bfac, 2)
                        if atom_name == "O" and key not in o_by_res:
                            o_by_res[key] = (rx, ry, rz)

                        atom_entry = {
                            "res": res_seq,
                            "chain": chain,
                            "resn": resn,
                            "x": rx,
                            "y": ry,
                            "z": rz,
                            "b": rb,
                        }
                        if atom_name in ("CA", "P") and key not in ca_by_res:
                            ca_by_res[key] = atom_entry
                        elif key not in fallback_by_res:
                            fallback_by_res[key] = atom_entry

                        if len(atoms_list) < MAX_EMBEDDED_ATOMS:
                            atoms_list.append(
                                {
                                    "n": atom_name,
                                    "rn": resn,
                                    "c": chain,
                                    "r": res_seq,
                                    "x": rx,
                                    "y": ry,
                                    "z": rz,
                                    "b": rb,
                                    "e": elem.upper(),
                                }
                            )
                continue
        idx += 1

    ordered_keys = sorted(
        set(ca_by_res.keys()) | set(fallback_by_res.keys()), key=lambda k: (k[0], k[1])
    )
    ca_trace = [
        dict(ca_by_res[k] if k in ca_by_res else fallback_by_res[k]) for k in ordered_keys
    ]

    if not ss_map and ca_trace:
        ss_map = _infer_ca_secondary_structure(ca_trace)

    for pt in ca_trace:
        key = (pt["chain"], pt["res"])
        pt["ss"] = ss_map.get(key, "C")
        if key in o_by_res:
            pt["ox"], pt["oy"], pt["oz"] = o_by_res[key]

    for at in atoms_list:
        at["ss"] = ss_map.get((at["c"], at["r"]), "C")

    return {
        "chains": chains,
        "residue_count": len(residues),
        "atom_count": atom_count,
        "res_numbers": res_numbers,
        "bfactors": bfactors,
        "ca_trace": ca_trace,
        "atoms": atoms_list,
    }


def parse_structure_metadata(file_path: str | Path) -> dict[str, Any]:
    """Parse a .pdb, .cif, or .mmcif file into dashboard metadata, all-atom records, and raw text."""
    validated = validate_structure_file(file_path)
    raw_text = validated.read_text(encoding="utf-8", errors="replace")
    suffix = validated.suffix.lower()
    viewer_format = "pdb" if suffix == ".pdb" else "cif"

    parsed = (
        _parse_pdb_coordinates(raw_text)
        if viewer_format == "pdb"
        else _parse_cif_coordinates(raw_text)
    )

    if parsed["atom_count"] == 0 or parsed["residue_count"] == 0:
        raise ValueError(
            f"No valid macromolecular ATOM coordinate records found in '{validated.name}'. "
            "Ensure the file contains standard 3D atomic coordinates."
        )

    bfactors: list[float] = parsed["bfactors"]
    res_numbers: list[int] = parsed["res_numbers"]
    mean_b = round(sum(bfactors) / len(bfactors), 2) if bfactors else 0.0
    min_b = round(min(bfactors), 2) if bfactors else 0.0
    max_b = round(max(bfactors), 2) if bfactors else 0.0

    return {
        "filename": validated.name,
        "extension": suffix,
        "viewer_format": viewer_format,
        "chains": parsed["chains"],
        "residue_count": parsed["residue_count"],
        "atom_count": parsed["atom_count"],
        "min_res_seq": min(res_numbers) if res_numbers else 1,
        "max_res_seq": max(res_numbers) if res_numbers else parsed["residue_count"],
        "mean_bfactor": mean_b,
        "min_bfactor": min_b,
        "max_bfactor": max_b,
        "ca_trace": parsed["ca_trace"],
        "atoms": parsed["atoms"],
        "raw_coordinates": raw_text,
    }


def infer_target_structure_url(
    accession: str,
    database: str = "",
    target_id: str = "",
) -> tuple[str | None, str | None]:
    """Construct public RCSB PDB or AlphaFold DB coordinate fetch URL for a target hit."""
    clean_acc = (accession or "").strip()
    db_lower = (database or "").strip().lower()

    af_match = _UNIPROT_ACCESSION_RE.match(clean_acc)
    if af_match and (clean_acc.upper().startswith("AF-") or db_lower.startswith("afdb")):
        uniprot_id = af_match.group(1).upper()
        return (
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif",
            "cif",
        )

    cath_match = _CATH_PDB_RE.search(target_id or clean_acc)
    if cath_match:
        pdb_id = cath_match.group(1).upper()
        return (f"https://files.rcsb.org/download/{pdb_id}.cif", "cif")

    pdb_match = _PDB_CHAIN_RE.match(clean_acc)
    if pdb_match and db_lower in ("pdb100", "cath50", ""):
        pdb_id = pdb_match.group(1).upper()
        return (f"https://files.rcsb.org/download/{pdb_id}.cif", "cif")

    return (None, None)


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for inspecting parsed 3D structure metadata."""
    parser = argparse.ArgumentParser(
        description="Validate and extract structural metadata from a .pdb, .cif, or .mmcif file."
    )
    parser.add_argument(
        "input_file",
        help="Path to input 3D coordinate file (.pdb, .cif, or .mmcif).",
    )
    parser.add_argument(
        "--include-coordinates",
        action="store_true",
        help="Include full raw coordinate text in JSON stdout output.",
    )
    args = parser.parse_args(argv)

    try:
        metadata = parse_structure_metadata(args.input_file)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if not args.include_coordinates:
        metadata = {k: v for k, v in metadata.items() if k != "raw_coordinates"}

    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
