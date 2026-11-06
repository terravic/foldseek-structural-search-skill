# Dashboard Analytical Payload Schema Reference

`scripts/generate_dashboard.py` compiles a normalized analytical JSON bundle (`DASHBOARD_DATA`) and embeds it into the standalone `dashboard.html` output. The formal machine-readable JSON Schema is stored at `assets/dashboard_payload_schema.json`.

---

## 1. Top-Level Bundle Structure

| Property | Type | Required | Description |
| :--- | :--- | :--- | :--- |
| `title` | `string` | Yes | Human-readable dashboard title displayed in the header and `<title>` tag. |
| `generated_at_utc` | `string` | Yes | ISO 8601 UTC timestamp recording when the dashboard bundle was compiled. |
| `system_theme` | `string` | Yes | Host OS fallback theme (`"light"` or `"dark"`) resolved during CLI build per `DESIGN.md`. |
| `query_structure` | `object` | Yes | Structural metadata, C-alpha trace coordinates, heavy-atom records, and raw coordinate file payload for `3Dmol.js` and software 3D projection rendering. |
| `summary_kpis` | `object` | Yes | Top-level scalar metrics rendered in the header KPI strip. |
| `hits` | `array<object>` | Yes | Up to 300 E-value-sorted normalized alignment hit records. |
| `coverage_profile` | `array<integer>` | Yes | 1D array of length `q_len` containing hit coverage depth for residues `1..q_len`. |
| `database_breakdown` | `array<object>` | Yes | Per-database hit counts and mean homology probabilities. |
| `functional_clusters` | `array<object>` | Yes | Tokenized protein family/function annotation clusters with counts and mean `Prob`. |
| `provenance` | `object` | Yes | Citation text, BibTeX metadata, service links, and database attribution records. |

---

## 2. `query_structure` Object

| Field | Type | Description |
| :--- | :--- | :--- |
| `filename` | `string` | Basename of the input structure file (e.g., `sample_query.pdb`). |
| `extension` | `string` | Normalized lowercase file extension (`.pdb`, `.cif`, `.mmcif`). |
| `viewer_format` | `string` | `3Dmol.js` parser identifier (`"pdb"` or `"cif"`). |
| `chains` | `array<string>` | Ordered list of unique chain identifiers detected in `ATOM` / `_atom_site` records. |
| `residue_count` | `integer` | Total unique macromolecular residues parsed from coordinates. |
| `atom_count` | `integer` | Total macromolecular `ATOM` records parsed. |
| `min_res_seq` | `integer` | Minimum residue sequence index in structure. |
| `max_res_seq` | `integer` | Maximum residue sequence index in structure. |
| `mean_bfactor` | `number` | Arithmetic mean of isotropic B-factor or AlphaFold pLDDT column across parsed atoms. |
| `min_bfactor` | `number` | Minimum B-factor / pLDDT value. |
| `max_bfactor` | `number` | Maximum B-factor / pLDDT value. |
| `ca_trace` | `array<object>` | Array of `{res, chain, resn, x, y, z, b, ss, ox?, oy?, oz?}` records with secondary structure (`"H"`, `"E"`, `"C"`) and peptide-plane carbonyl coordinates used for 3D ribbon spline interpolation. |
| `atoms` | `array<object>` | Array of `{n, rn, c, r, x, y, z, b, e, ss}` heavy-atom records (up to 3,500 CLI-embedded or 8,000 client-parsed heavy atoms) used for covalent sidechain sticks, ball-and-stick, and VDW spacefill/surface rendering. |
| `raw_coordinates` | `string` | Complete text of the query `.pdb`/`.cif` file for client-side `3Dmol.js` and structure-details inspection. |

---

## 3. Normalized Hit Record (`hits[]`)

| Field | Type | Description |
| :--- | :--- | :--- |
| `rank` | `integer` | 1-based rank after ascending E-value sort. |
| `target_id` | `string` | Raw Foldseek `target` string including accession and description. |
| `accession` | `string` | Extracted primary accession (e.g., `1UBQ_A`, `AF-P69905-F1`). |
| `description` | `string` | Cleaned human-readable protein annotation parsed from `target_id`. |
| `organism` | `string` | Taxonomic organism name extracted from `OS=...` or bracketed/parenthesized suffix when present. |
| `database` | `string` | Source Foldseek database identifier (e.g., `pdb100`, `afdb50`). |
| `prob` | `number` | Homology probability in `[0.0, 1.0]`. |
| `prob_display` | `string` | Probability formatted to 3 decimal places (`"0.994"`). |
| `evalue` | `number` | Numeric expectation value (`eval` / `eValue` / `evalue`). |
| `evalue_display` | `string` | Scientific or compact string representation of E-value. |
| `neg_log10_evalue` | `number` | `-log10(max(evalue, 1e-300))` for scatter plot Y-axis scaling. |
| `seq_identity_pct` | `number` | Sequence identity normalized to `[0.0, 100.0]`. |
| `seq_identity_display` | `string` | Formatted sequence identity percentage (`"24.5%"`). |
| `q_start` | `integer` | 1-based start residue index on query (`qStartPos`). |
| `q_end` | `integer` | 1-based inclusive end residue index on query (`qEndPos`). |
| `q_len` | `integer` | Query sequence length reported in alignment (`qLen`). |
| `q_cov_pct` | `number` | Numeric query coverage percentage in `[0.0, 100.0]`. |
| `q_cov_display` | `string` | Formatted query coverage (`"95.0%"` or `"N/A"`). |
| `t_start` | `integer` | 1-based start residue index on target (`dbStartPos` / `tStartPos`). |
| `t_end` | `integer` | 1-based end residue index on target (`dbEndPos` / `tEndPos`). |
| `t_len` | `integer` | Target sequence length (`dbLen` / `tLen`). |
| `aln_len` | `integer` | Alignment column count (`alnLength` / `alnLen`). |
| `score` | `number` | Alignment bit score (`score` / `bits`). |
| `homology_regime` | `string` | One of `"High-Identity Homolog"`, `"Remote Structural Homolog (Twilight Zone)"`, `"Local / Partial Motif Match"`, or `"Low Confidence"`. |
| `target_structure_url` | `string \| null` | Optional HTTPS URL to RCSB PDB (`.cif`) or AlphaFold DB (`.cif`) coordinate file for on-demand 3D browser overlay. |
| `target_structure_format` | `string \| null` | Format of `target_structure_url` (`"cif"` or `"pdb"`). |
