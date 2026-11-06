# Product Requirement Document: `foldseek-structural-search` Agent Skill

## 1. Document Metadata

| Field | Value |
| :--- | :--- |
| **Product Name** | `foldseek-structural-search` |
| **Skill Identifier** | `foldseek-structural-search` |
| **License** | Apache License 2.0 (Clean-room REST integration with Foldseek API) |
| **Primary Runtime** | Python 3.10+ (`uv` PEP 723 script execution and standard virtual environments) |
| **Target Environments** | Standardized SKILL harnesses, automated bioinformatics pipelines, standalone CLI, offline web browsers |

---

## 2. User Stories & Scope

### 2.1 Target Personas

1. **Computational Structural Biologists and Protein Engineers (Technical Users)**:
   - Execute command-line batch or single-structure homology searches on experimental (`.pdb`, `.cif`, `.mmcif`) or predicted 3D protein models.
   - Require rapid separation of high-sequence-identity orthologues from remote "twilight-zone" structural homologues (high structural probability `Prob >= 0.80` paired with low sequence identity `< 30%`).
   - Need interactive 3D macromolecular inspection of conserved structural cores, sidechain packing, and secondary structure architecture mapped directly onto atomic coordinates.
2. **Translational and Life-Science Researchers (Non-Technical / Semi-Technical Users)**:
   - Provide a local protein structure file path to a conversational agent or load a `.pdb`/`.cif` file directly in `dashboard.html` to inspect macromolecular architecture, evolutionary context, and structural homologues without writing code.
   - Rely on automated input validation, remote database querying, structured Markdown synthesis without raw coordinate dumps in prompt context, and an interactive standalone HTML dashboard (`dashboard.html`).
3. **Platform & Validation Engineers**:
   - Validate skill behavior deterministically in air-gapped or continuous testing environments using `--mock` mode and 100% synthetic Healthcare and Life Sciences (HCLS) compliant fixtures without external network calls.

### 2.2 User Stories

| ID | Persona | User Story | Acceptance Criteria |
| :--- | :--- | :--- | :--- |
| **US-01** | Life-Science Researcher | As an agent user, I want to supply a `.pdb`, `.cif`, or `.mmcif` file and receive a functional summary of structurally homologous proteins across PDB and AlphaFold databases. | Pipeline runs coordinate validation and Foldseek search, writes JSON, Markdown, and `dashboard.html` outputs, and summarizes top hits, twilight-zone homologues, and functional clusters with full Foldseek citation. |
| **US-02** | Life-Science Researcher | As an agent user who supplies only a UniProt ID, gene symbol, or FASTA sequence, I want immediate feedback explaining why Foldseek cannot run. | Both `SKILL.md` instructions and `scripts/search.py` reject non-coordinate inputs with exit code `1` and instruct the user to supply a physical 3D coordinate file. |
| **US-03** | Structural Biologist | As a researcher analyzing macromolecular structures, I want to switch between multiple 3D representation styles and color schemes, click any alignment hit or secondary structure element, and see the corresponding query residue span highlighted in 3D. | `dashboard.html` renders all-atom and ribbon representations via WebGL (`3Dmol.js`) and an offline Blinn-Phong shaded software 3D projection engine, supporting 5 representation styles, 5 color modes, and synchronized `qStartPos..qEndPos` highlighting. |
| **US-04** | Life-Science Researcher | As a browser user, I want to load any local `.pdb` or `.cif` file directly into `dashboard.html` on demand, inspect its structural details immediately, and optionally run a live Foldseek search from the UI. | Client-side PDB/mmCIF parser extracts chains, secondary structure, ligands, radius of gyration, bounding box, and sequence, switches to the `Structure Details` view when hits are absent, and supports live browser Foldseek submission or `.json` hit loading. |
| **US-05** | Protein Engineer | As a CLI user with existing pipelines redirecting `scripts/search.py` stdout to Markdown files, I want backward-compatible CLI flags and Markdown output formatting. | `scripts/search.py <file> -o <out.json> --databases <dbs> > <out.md>` emits the exact `### All Structural Matches` Markdown table of up to 300 E-value-sorted hits. |
| **US-06** | Platform Engineer | As a compliance engineer, I want upstream Foldseek licensing notices recorded before execution and deterministic offline testing with synthetic data. | Workspace license guard checks/creates `.licenses/foldseek_structural_search_LICENSE.txt`; `--mock` mode and `pytest` execute offline with exit code `0` on synthetic fixtures. |

### 2.3 Explicit Scope Boundaries

#### In-Scope Capabilities
- Validation of local 3D macromolecular coordinate files (`.pdb`, `.cif`, `.mmcif`, plus `.ent` in client-side browser inspection) and structural `ATOM` / `_atom_site` record presence.
- Allowlist validation across the 9 supported Foldseek server databases: `afdb50`, `afdb-swissprot`, `pdb100`, `BFVD`, `mgnify_esm30`, `cath50`, `gmgcl_id`, `bfmd`, and `afdb-proteome`.
- Polite rate-limited (`qps=1`) HTTP multipart submission (`3diaa` mode), asynchronous ticket polling, and JSON result retrieval from `https://search.foldseek.com/api/`.
- E-value sorting, hit capping (`MAX_ALIGNMENT_HITS = 300`), query coverage computation, homology regime classification, per-residue coverage depth profiling, and functional annotation keyword clustering.
- All-atom heavy-atom extraction, backbone C-alpha + carbonyl oxygen (`O`) ribbon plane extraction, and explicit (`HELIX`/`SHEET`, `_struct_conf`/`_struct_sheet_range`) or geometric C-alpha secondary structure classification (`H`/`E`/`C`).
- Standalone HTML5 dashboard generation (`dashboard.html`) compliant with `DESIGN.md` (Light/Dark/System themes, WebGL `3Dmol.js` + Blinn-Phong shaded software 3D projection engine, 5 representation styles, 5 color modes, on-demand `.pdb`/`.cif` file loading and drag-and-drop, `Structure Details` inspector, live browser `Search Foldseek` execution, RCSB PDB / AlphaFold DB target overlay, alignment track map, homology scatter plot, database/cluster distribution bars, filterable/sortable table with CSV/Markdown export).
- Structure-only dashboard compilation (`scripts/generate_dashboard.py --input_structure <file> -o dashboard.html` without `--results_json`) for immediate structural inspection prior to homology search.
- Deterministic `--mock` execution and CLI-executable unit, smoke, and system tests on 100% synthetic protein structures and alignment payloads.

#### Out-of-Scope Capabilities
- Searching by raw amino acid sequence (FASTA), gene name, or accession identifier (e.g., UniProt or PDB ID without a local coordinate file).
- De novo 3D protein structure prediction (e.g., running AlphaFold or ESMFold locally).
- Vendoring, copying, or linking GPL-3.0 Foldseek C++ source code or local binary database indexing.
- Storing, processing, or transmitting Protected Health Information (PHI) or Personally Identifiable Information (PII).

---

## 3. Technical Architecture

### 3.1 System Sequence & Dataflow Diagram

```text
+-------------------+       +------------------------+       +-------------------------------+
| End User / CLI    |       | Skill Orchestrator     |       | Workspace File System         |
+---------+---------+       | (SKILL.md)             |       | (.licenses/, assets/, out/)   |
          |                 +-----------+------------+       +---------------+---------------+
          | 1. Provide .pdb/.cif path   |                                    |
          +---------------------------->| 2. Validate file extension & DBs   |
          |                             | 3. Invoke run_pipeline.py /        |
          |                             |    search.py                       |
          |                             +-------------------+                |
          |                                                 |                |
          |                                                 v                |
          |                             +---------------------------------+  |
          |                             | scripts/run_pipeline.py         |  |
          |                             | scripts/search.py               |  |
          |                             +---+-------------------------+---+  |
          |                                 | 4. Check/Write License  |      |
          |                                 +-------------------------+----->| .licenses/foldseek_...txt
          |                                 |                         |      |
          |        [Live Mode]              |                         | [Mock Mode: --mock]
          |   +-----------------------------+                         +------------------+
          |   | 5a. POST /api/ticket (multipart 3diaa, qps=1)                            |
          |   | 5b. Poll GET /api/ticket/{id} (10s interval)                             | 5m. Load synthetic
          |   | 5c. Sleep 3s -> GET /api/result/{id}/0                                   |     sample_foldseek_results.json
          |   v                                                                          v
          | +-------------------------------+                        +----------------------------------+
          | | Foldseek Public REST API      |                        | sample_data/                     |
          | | https://search.foldseek.com   |                        | sample_foldseek_results.json     |
          | +---------------+---------------+                        +-----------------+----------------+
          |                 |                                                          |
          |                 +----------------------------+-----------------------------+
          |                                              | 6. Raw JSON Payload
          |                                              v
          |                 +----------------------------------------------------------+
          |                 | Analytical & Structural Processing Layer                 |
          |                 |  - scripts/visualize_structure.py (PDB/CIF -> chains,    |
          |                 |    B-factor/pLDDT, SS H/E/C, CA+O trace, heavy atoms)    |
          |                 |  - scripts/visualize_alignments.py (sort E-val <= 300,   |
          |                 |    Q-Cov, homology regimes, 1D depth, function clusters) |
          |                 +----------------------------+-----------------------------+
          |                                              |
          |                                              v
          |                 +----------------------------------------------------------+
          |                 | scripts/generate_dashboard.py + assets/                  |
          |                 |  - Validates payload structure (or structure-only mode)  |
          |                 |  - Compiles standalone HTML5 (DESIGN.md tokens,          |
          |                 |    3Dmol.js + software 3D renderer, SVG plots, loader)   |
          |                 +----------------------------+-----------------------------+
          |                                              | 7. Write Outputs
          |                                              +--------------------------------------------->|
          |                                              |   - <stem>_foldseek_results.json             |
          |                                              |   - <stem>_foldseek_results.md               |
          |                                              |   - dashboard.html                           |
          |                 +----------------------------+                                              |
          |                 | 8. Read .md summary only (no raw coordinates/JSON in prompt context)      |
          |                 v                                                                           |
          |<----------------+ 9. Return homologue synthesis, twilight-zone analysis & relative paths    |
```

### 3.2 Module Responsibilities

| File Path | Architectural Role |
| :--- | :--- |
| `SKILL.md` | Agent instruction specification, input/database guardrails, execution options, interpretation workflow, and citation enforcement. |
| `DESIGN.md` | Visual architecture, color tokens, typography scale, Bento workspace grid layout, and interactive component specifications for `dashboard.html`. |
| `scripts/search.py` | Core Foldseek REST client (`polite-http`), multipart form builder, license file guard, `--mock` loader, JSON writer, backward-compatible Markdown stdout emitter, and `--dashboard` hook. |
| `scripts/visualize_structure.py` | Parser for `.pdb`, `.cif`, and `.mmcif` coordinate files; extracts chains, residue bounds, atom counts, B-factor/pLDDT statistics, secondary structure (`H`/`E`/`C`), C-alpha + carbonyl `O` trace, heavy-atom coordinates, and target RCSB/AlphaFold fetch URLs. |
| `scripts/visualize_alignments.py` | Alignment normalizer computing `Q-Cov`, `-log10(E-value)`, normalized sequence identity percentages, homology regimes, per-residue coverage depth arrays, database distributions, and functional annotation clusters. |
| `scripts/generate_dashboard.py` | Compiles query structural data and optional normalized alignment analytics into a self-contained `dashboard.html` using `assets/dashboard_template.html` and `DESIGN.md` specifications. |
| `scripts/run_pipeline.py` | Unified end-to-end CLI entrypoint orchestrating search, Markdown functional synthesis report generation, dashboard compilation, and optional local HTTP serving. |
| `assets/dashboard_payload_schema.json` | JSON Schema defining the analytical data bundle contract injected into `dashboard.html`. |
| `assets/dashboard_template.html` | Parameterized HTML5, vanilla CSS custom properties, and client-side JavaScript workspace template including on-demand `.pdb`/`.cif` parsing and 3D rendering. |
| `assets/foldseek_skill_overview.png` | Architectural workflow diagram illustrating inputs, processing stages, interactive dashboard panels, and structured outputs for `README.md`. |
| `skills/foldseek_structural_search/` | Standardized skill package directory containing relative symlinks to root `SKILL.md`, `scripts/`, `references/`, `assets/`, `documents/`, and `sample_data/`. |

---

## 4. Functional Requirements

### 4.1 Input & Parameter Validation (`FR-01` to `FR-04`)

- **FR-01 (Coordinate File Enforcement)**: All CLI entrypoints (`search.py`, `visualize_structure.py`, `generate_dashboard.py`, `run_pipeline.py`) must verify that the input path points to an existing regular file with an extension in `{".pdb", ".cif", ".mmcif"}` (case-insensitive) containing parseable `ATOM` or `_atom_site` coordinate records. Passing a raw sequence string, gene symbol, accession ID, missing path, or unsupported file extension must terminate immediately with exit code `1`.
- **FR-02 (Database Allowlist Enforcement)**: Requested databases (`--databases`, comma-separated, default `"pdb100,afdb50"`) must be validated against the exact allowlist:
  `["afdb50", "afdb-swissprot", "pdb100", "BFVD", "mgnify_esm30", "cath50", "gmgcl_id", "bfmd", "afdb-proteome"]`. Any unrecognized database identifier must terminate execution immediately with exit code `1` prior to network activity.
- **FR-03 (Workspace License Notification Guard)**: Before executing a Foldseek search or pipeline run, the tool must check for `<workspace_root>/.licenses/foldseek_structural_search_LICENSE.txt` (overridable via `FOLDSEEK_WORKSPACE_ROOT` for isolated testing). If absent, it must emit a notice to `sys.stderr` referencing `https://search.foldseek.com/search` and `https://github.com/steineggerlab/foldseek` and write the file with an ISO 8601 UTC timestamp.
- **FR-04 (Offline Mock Execution)**: Passing `--mock` to `search.py` or `run_pipeline.py` must bypass all outbound HTTP calls and load deterministic synthetic alignment records from `sample_data/sample_foldseek_results.json`.

### 4.2 Foldseek API Interaction & Alignment Math (`FR-05` to `FR-09`)

- **FR-05 (Multipart Submission)**: `build_multipart_payload(fields, files)` must construct RFC 7578 compliant `multipart/form-data` bytes encoding `mode="3diaa"`, one `database[]` part per selected database, and file upload field `q` containing the coordinate file bytes.
- **FR-06 (Ticket Polling Lifecycle)**:
  1. Submit `POST https://search.foldseek.com/api/ticket` (`timeout=30` seconds).
  2. Poll `GET https://search.foldseek.com/api/ticket/{ticket_id}` every `10` seconds (`timeout=20` seconds) until `status == "COMPLETE"`. Halt with exit code `1` if `status == "ERROR"` or if maximum polling duration (`600` seconds default, configurable via `FOLDSEEK_MAX_POLL_SECONDS`) is exceeded.
  3. Pause `3` seconds after `"COMPLETE"`, then fetch `GET https://search.foldseek.com/api/result/{ticket_id}/0` (`timeout=120` seconds).
- **FR-07 (Sorting, Capping & Query Coverage Formula)**:
  - Flatten nested alignments across `results[*].alignments[*]`, attaching `result_group.get("db", "unknown")` to each hit.
  - Sort ascending by numeric E-value (`eval` / `eValue` / `evalue`, defaulting to `1000.0`), slicing at `MAX_ALIGNMENT_HITS = 300`.
  - Compute Query Coverage (`Q-Cov`) using:
    `cov_percentage = min(((qEndPos - qStartPos + 1) / qLen) * 100, 100.0)` when `qLen > 0` and `qEndPos > qStartPos` (formatting as `f"{cov_percentage:.1f}%"`, else `"N/A"`).
- **FR-08 (Homology Regime Classification)**:
  Each normalized hit is assigned to one of four mutually exclusive regimes in `classify_homology_regime()`:
  1. `"High-Identity Homolog"`: `prob >= 0.8` and `seq_identity_pct >= 30.0`
  2. `"Remote Structural Homolog (Twilight Zone)"`: `prob >= 0.8` and `seq_identity_pct < 30.0`
  3. `"Local / Partial Motif Match"`: `0.5 <= prob < 0.8` (capturing partial domain and local motif alignments)
  4. `"Low Confidence"`: `prob < 0.5`
- **FR-09 (Conservation Depth & Functional Clustering)**:
  - Generate a 1-indexed per-residue coverage depth profile array of length `qLen` where element `i` counts hits spanning query residue `i + 1`.
  - Tokenize target descriptions (excluding accession tokens, bracketed/parenthesized taxonomy strings, `OS=`/`OX=` UniProt tags, and English/structural stopwords) into ranked protein family/function clusters with hit counts and mean `Prob`.

### 4.3 Structural Coordinate Extraction, Dashboard & Reporting (`FR-10` to `FR-14`)

- **FR-10 (Backward-Compatible Markdown Table)**: `scripts/search.py` stdout must emit `### All Structural Matches` followed by `| Target ID | Q-Cov | Prob | E-value | Seq Identity | Aln Length |` and the Foldseek citation block.
- **FR-11 (Pipeline Synthesis Report)**: `scripts/run_pipeline.py` must write `<stem>_foldseek_results.json`, `<stem>_foldseek_results.md` (containing summary metrics, top homologue interpretations, twilight-zone candidates, functional family clusters, the full Markdown table, and academic citation), and `dashboard.html`.
- **FR-12 (All-Atom & Secondary Structure Extraction)**: `scripts/visualize_structure.py` must extract both `ca_trace` (per-residue `x, y, z`, optional peptide carbonyl oxygen `ox, oy, oz`, B-factor/pLDDT `b`, and secondary structure `ss` in `{"H", "E", "C"}` from `HELIX`/`SHEET` or `_struct_conf`/`_struct_sheet_range` records with C-alpha distance-geometry fallback) and `atoms` (up to `3,500` non-hydrogen heavy atoms for CLI embedding).
- **FR-13 (Interactive Standalone Dashboard & 3D Rendering)**: `dashboard.html` must operate directly in a web browser from local disk or HTTP without build steps and implement:
  - Light, Dark, and System theme switching (`#theme-btn-light`, `#theme-btn-dark`, `#theme-btn-system`) with exact `DESIGN.md` surface tokens and OS theme tracking.
  - KPI summary header cards reflecting either search statistics + query metadata (when hits are loaded) or macromolecular biophysical properties (Radius of Gyration, Bounding Box, Alpha-Helix %, Beta-Sheet %, Bound Ligands when inspecting a structure without hits).
  - Dual-engine 3D macromolecular viewport supporting WebGL (`3Dmol.js`) and an automatic Blinn-Phong shaded software 3D projection engine (`#software-viewport`) with Catmull-Rom spline ribbons, helical ribbons, pleated beta-strand arrows, covalent sidechain bond inference, 5 representation styles (`Cartoon + Sidechains`, `Cartoon Ribbon`, `Surface + Cartoon`, `Ball & Stick`, `Spacefill (VDW)`), and 5 color modes (`N-to-C Spectrum`, `Secondary Structure`, `B-Factor / pLDDT`, `Alignment Depth`, `Chain / Heteroatom`).
  - Dual-mode right-hand inspector toggling between `Alignment Tracks` (with 1D `Alignment Depth` or `B-Factor / pLDDT + SS` profile) and `Structure Details` (8-metric biophysical summary, clickable Chain Inventory, clickable Secondary Structure Elements table, and clickable 1-letter Primary Sequence viewer).
  - Homology Landscape SVG scatter plot (`Seq Identity %` vs. `-log10(E-value)` / `Prob`), Database / Functional Cluster bar charts, and real-time filterable/sortable Structural Hits Explorer table with CSV and Markdown export controls.
  - Compliance footer displaying the `vanKempen2023` citation, Foldseek service/repository links, database attributions, and research disclaimer.
- **FR-14 (On-Demand Browser Coordinate Loader & Structure-Only CLI Mode)**:
  - `dashboard.html` must allow users to load `.pdb`, `.ent`, `.cif`, or `.mmcif` files on demand via file picker (`Load PDB / CIF`, `Open .PDB`), drag-and-drop onto the 3D viewport, or programmatic `window.loadStructureFromText(filename, rawText)` (parsing up to `8,000` heavy atoms client-side), optionally load a Foldseek `.json` payload (`Load Hits .JSON` / `window.loadHitsJsonFromText`), trigger a live browser `Search Foldseek` request, or click `Restore Default` to return to the compiled payload.
  - `scripts/generate_dashboard.py` must support `--input_structure <file> -o dashboard.html` with `--results_json` omitted, compiling a structure-inspection dashboard pre-focused on the `Structure Details` view.

---

## 5. Non-Functional Requirements

### 5.1 Performance & Rate-Limiting Latency
- **Polite API Throttling**: Outbound requests to `search.foldseek.com` are governed by `polite-http` at `qps=1` with a `10`-second polling interval and automatic exponential backoff on HTTP `429` or `5xx` responses.
- **Local Processing Latency**: Offline `--mock` pipeline execution, coordinate metadata extraction (up to 10,000 atoms), alignment normalization (300 hits), and `dashboard.html` compilation must complete in under `2.0` seconds on standard commodity hardware.
- **Prompt Context Protection**: Neither raw 3D coordinate records nor full Foldseek JSON alignment arrays are loaded into an agent's prompt context; only concise `.md` summaries are read during synthesis.

### 5.2 Security, HCLS Safety & Data Privacy Compliance
- **HIPAA & HCLS Privacy**: All bundled sample fixtures (`sample_data/sample_query.pdb`, `sample_data/sample_query.cif`, `sample_data/sample_foldseek_results.json`) and test cases are 100% synthetic de novo constructs containing zero Protected Health Information (PHI), patient identifiers, clinical metadata, or Personally Identifiable Information (PII).
- **Public Server Data Governance**: When executed in live mode (via CLI or the browser `Search Foldseek` button), the query coordinate file is transmitted over HTTPS (`TLS 1.2+`) to the public Foldseek server (`https://search.foldseek.com`). Users working with proprietary or pre-publication structures are informed in `README.md` and `SKILL.md` regarding external API transmission; local `.pdb`/`.cif` inspection in `dashboard.html` without clicking `Search Foldseek` runs 100% locally in browser memory.
- **XSS & Injection Prevention**: All target identifiers, protein annotations, filenames, and titles embedded into `dashboard.html` undergo strict JSON serialization (`</script>` escaping and ASCII-safe encoding) and DOM text-node/HTML entity escaping prior to rendering.
- **Path-Relative Portability**: Zero hardcoded local machine or user environment paths exist in source code or documentation; all asset, reference, sample, and license lookups resolve dynamically relative to the project root or explicit CLI/environment parameters.

### 5.3 Licensing & Extensibility
- **Clean-Room GPL-3.0 Separation**: All scripts in this repository are original Apache License 2.0 implementations communicating strictly with the public HTTP REST API; zero Foldseek GPL-3.0 source code is copied or linked.
- **Modular Extensibility**: `visualize_structure.py`, `visualize_alignments.py`, and `generate_dashboard.py` expose typed Python functions (`parse_structure_metadata`, `build_alignment_analytics`, `build_dashboard`) allowing downstream structural bioinformatics workflows to import them directly.

---

## 6. Error Handling & Edge Cases Matrix

| Failure Scenario / Edge Case | Detection Layer | System Behavior & Exit Code | User-Facing Diagnostic Message |
| :--- | :--- | :--- | :--- |
| **User passes UniProt ID, gene name, or FASTA sequence string (e.g., `P69905` or `MKWVTF...`)** | `SKILL.md` pre-check & `search.py` / `visualize_structure.py` path check | Halts immediately before API call; returns exit code `1`. | Explains that Foldseek requires a physical `.pdb`, `.cif`, or `.mmcif` coordinate file path and cannot search by sequence or accession ID. |
| **Unsupported file extension (e.g., `query.fasta`, `query.txt`, `query.json`)** | `validate_structure_file()` extension check | Halts immediately; returns exit code `1`. | Lists allowed coordinate extensions (`.pdb`, `.cif`, `.mmcif`) and rejects the supplied file. |
| **Input file does not exist or is a directory** | `Path.is_file()` validation | Halts immediately; returns exit code `1`. | Reports exact missing/non-file path and instructs user to verify workspace path. |
| **Empty or structurally invalid `.pdb`/`.cif` file (zero `ATOM`/`_atom_site` records)** | `parse_structure_metadata()` or client `parseStructureClientSide()` | CLI raises `ValueError` and exits with code `1`; browser displays inline error banner without losing current session. | Reports that zero valid macromolecular coordinate records were found in the file. |
| **Unrecognized database name in `--databases` (e.g., `nr`, `uniprot`, `invalid_db`)** | `validate_databases()` allowlist check | Halts immediately before API call; returns exit code `1`. | Identifies the invalid database token(s) and prints all 9 allowed Foldseek database names. |
| **Missing `.licenses/foldseek_structural_search_LICENSE.txt` in workspace** | `ensure_license_notification()` guard | Creates `.licenses/` directory and timestamped file; prints notice to `sys.stderr`; continues (`0`). | Displays Foldseek Terms of Use and repository URLs on `stderr` without corrupting `stdout` Markdown redirection. |
| **Foldseek ticket returns `status == "ERROR"`** | `search.py` or browser polling loop | Terminates polling immediately; CLI returns exit code `1`; browser displays status alert. | Reports remote Foldseek server ticket failure and ticket ID for troubleshooting. |
| **Network timeout, DNS failure, or HTTP `429`/`5xx` exhaustion** | `polite_http.HttpError` exception handler | Retries transient `429`/`5xx` per `polite-http` backoff; exits with code `1` if exhausted. | Logs HTTP status/URL failure cleanly to `stderr` without unhandled Python tracebacks. |
| **JSON result contains zero alignments (`alignments: null` or `[]`) or `--results_json` omitted** | `extract_and_normalize_hits()` / `build_dashboard_payload()` | Succeeds with exit code `0`; writes valid JSON/Markdown and valid `dashboard.html`. | Dashboard KPI cards display biophysical structure metrics and inspector opens `Structure Details` cleanly. |
| **Hit record missing `qLen` (`0`), inverted `qEndPos <= qStartPos`, or `seqId` in `[0, 1]` vs `[0, 100]`** | `extract_and_normalize_hits()` defensive math | Renders `Q-Cov` as `"N/A"` (`0.0` numeric), caps valid coverage at `100.0%`, and auto-normalizes `seqId <= 1.0` (`* 100.0`) unless exact integer `> 1`. | Prevents `ZeroDivisionError` and scale distortion in coverage tracks and scatter plots. |
| **Browser lacks hardware WebGL or runs offline without `3Dmol.js` CDN access** | `dashboard.html` client-side initialization check | Catches missing `$3Dmol` or WebGL context failure and activates interactive Blinn-Phong shaded software 3D projection engine (`#software-viewport`). | Renders interactive drag-rotatable query ribbon, sidechains, surfaces, spheres, and residue highlighting offline. |
