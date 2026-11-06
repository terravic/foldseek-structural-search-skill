---
name: foldseek-structural-search
description: >
  Performs 3D structural searches of proteins against various databases (PDB,
  AlphaFold, CATH, MGnify, etc.) using the Foldseek API and generates an
  interactive 3D structure and alignment visualization dashboard. Use ONLY when
  the user provides a physical 3D coordinate file (.cif, .mmcif, or .pdb) and
  wants to find structurally similar proteins. Do NOT use if the user only
  provides a protein sequence, gene name, or UniProt ID.
---

# Foldseek 3D Structural Search and Visualization

## Prerequisites

1. **`uv`**: Ensure `uv` (or Python 3.10+ with dependencies from `requirements.txt`) is installed and available on `PATH`.
2. **User Notification**: If `.licenses/foldseek_structural_search_LICENSE.txt` does not already exist in the workspace root directory then (1) prominently notify the user to check the terms at https://search.foldseek.com/search and https://github.com/steineggerlab/foldseek, then (2) create the file recording the notification text and timestamp.

## Goal

Submit a user-provided 3D protein structure file (`.cif`, `.mmcif`, or `.pdb`) to the Foldseek web server API to find structurally similar proteins. Report the top structural hits, interpret key alignment metrics, summarize the inferred protein functions, save the Markdown-formatted table to a `.md` file, save the full detailed results to a local `.json` file, and generate an interactive HTML5 visualization dashboard (`dashboard.html`) featuring a 3Dmol.js macromolecular viewer, on-demand `.pdb` / `.cif` structure inspection, query domain coverage tracks, and homology scatter plots.

## Core Rules

- **File Requirement**: This tool cannot search by sequence, name, or accession ID. It strictly requires a `.pdb`, `.cif`, or `.mmcif` file path.
- **Strict Validation**: Never bypass the input validation or the database allowlist check.
- **Do Not Parse the JSON Manually in Context**: Rely on the generated `.md` file for immediate text summaries and the generated `dashboard.html` for visual exploration. The JSON file is saved for downstream tool and dashboard consumption.
- **No Raw Coordinate Parsing in Context**: Do not read raw 3D atomic coordinates into LLM context; pass the relative coordinate file path directly to the CLI scripts.
- **Data Governance Notice**: In live execution mode, the coordinate file is submitted over HTTPS to the public Foldseek server (`https://search.foldseek.com`). Obtain user confirmation before submitting confidential or unpublished structures.
- **Citation**: When reporting results from this skill, include the Foldseek citation (van Kempen et al., *Nature Biotechnology* 2023).

## Instructions

1. **Strict Input Validation**: Verify that the user has provided a valid path to a `.cif`, `.mmcif`, or `.pdb` file in the workspace.
   - If the user provided a protein name, an amino acid sequence, or an accession ID (e.g., a UniProt ID) without a coordinate file, **halt immediately** and do not run the script.
   - Inform the user that Foldseek requires a physical 3D coordinate file, and suggest downloading the structure first from RCSB PDB or the AlphaFold Protein Structure Database.
2. **Database Validation**: Check whether the user requested specific target databases.
   - **Allowed List**: `afdb50`, `afdb-swissprot`, `pdb100`, `BFVD`, `mgnify_esm30`, `cath50`, `gmgcl_id`, `bfmd`, `afdb-proteome`.
   - If the user requests a database outside this list, **halt immediately**, inform the user that the database is unsupported, and list the allowed identifiers.
3. **Execute Search and Generate Outputs** (resolve `scripts/` relative to the directory containing this `SKILL.md`):
   - **Option A - Unified Pipeline (JSON + Markdown + Dashboard)**:
     ```bash
     uv run scripts/run_pipeline.py <path-to-file> -o <output-dir> [--databases <db1,db2,db3>]
     ```
   - **Option B - Direct Search Script (Backward-compatible CLI)**:
     - **Default databases (`pdb100,afdb50`)**:
       ```bash
       uv run scripts/search.py <path-to-file> -o <generated-filename.json> --dashboard <generated-dashboard.html> > <generated-filename.md>
       ```
     - **Custom databases**:
       ```bash
       uv run scripts/search.py <path-to-file> -o <generated-filename.json> --databases <db1,db2,db3> --dashboard <generated-dashboard.html> > <generated-filename.md>
       ```
   - **Option C - Standalone Dashboard Generation (With Optional Results JSON)**:
     ```bash
     uv run scripts/generate_dashboard.py --input_structure <path-to-file> [--results_json <generated-filename.json>] -o <generated-dashboard.html>
     ```
   - **Offline Verification Mode**: Append `--mock` to `run_pipeline.py` or `search.py` when running offline verification with synthetic alignment data from `sample_data/sample_foldseek_results.json`.
4. **Read the Results**: Open and read the generated `.md` file to review the Markdown table and functional synthesis sections.
5. **Interpret the Metrics**: Summarize the top 3 to 5 annotated structural matches using:
   - **Prob (Probability)**: Values approaching `1.0` (`100%`) indicate high confidence that the target shares homologous 3D fold topology.
   - **Q-Cov (Query Coverage)**: High percentages indicate alignment across the global query fold rather than an isolated local motif.
   - **E-value and Seq Identity**: Provide evolutionary context, noting remote twilight-zone homologues (`Prob >= 0.80` paired with `Seq Identity < 30%`).
6. **Perform Functional Analysis**: Analyze the descriptions embedded within the `Target ID` column for the reported matches.
   - Report the specific protein names and functions of the top structural homologues.
   - Summarize the range of recurring functions, domains, or protein families represented across the hit list.
7. **Report Generated Files**: List all newly created files (`.json`, `.md`, and `dashboard.html`) using relative workspace paths alongside the Foldseek academic citation:
   > van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Soding, J., & Steinegger, M. (2023). *Fast and accurate protein structure search with Foldseek.* Nature Biotechnology, 41(12), 1673-1682. doi:10.1038/s41587-023-01773-0.
   If the API returns an error or the input path is missing, report the diagnostic message clearly and ask the user to verify the file path.

## Additional Reference Resources

- Database descriptions, `3diaa` alignment mechanics, and metric equations: `references/databases_and_metrics.md`
- Dashboard analytical JSON contract: `references/dashboard_schema.md` and `assets/dashboard_payload_schema.json`
- BibTeX record: `references/citation.bib`
