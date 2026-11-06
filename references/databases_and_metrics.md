# Foldseek Structural Databases and Alignment Metrics Reference

## 1. Alignment Mode: `3diaa`

Foldseek encodes tertiary protein interactions by discretizing local backbone and virtual C-beta geometric conformations around each residue into a 20-state structural alphabet termed **3Di** (three-dimensional interactions). In `3diaa` alignment mode, Foldseek combines the 3Di tertiary interaction substitution score with standard 20-letter amino acid (**AA**) sequence substitution scoring. This hybrid representation preserves sensitivity across deep evolutionary timescales ("twilight zone" below 25-30% sequence identity) while maintaining specificity among closely related orthologues.

---

## 2. Supported Foldseek Target Databases

Requests submitted by `scripts/search.py` and `scripts/run_pipeline.py` are strictly validated against the following 9 target structural databases hosted by `https://search.foldseek.com`:

| Database Identifier | Full Name / Source | Scope & Coverage | Upstream Data License | Recommended Use Case |
| :--- | :--- | :--- | :--- | :--- |
| `pdb100` | Protein Data Bank (100% clustered) | Experimentally determined X-ray, cryo-EM, and NMR structures from wwPDB / RCSB PDB. | CC0 1.0 Public Domain | High-accuracy experimental homologues, ligand-bound states, and curated macromolecular complexes. Default target. |
| `afdb50` | AlphaFold Protein Structure Database (50% clustered) | Sequence-clustered (50% identity) representatives across 200M+ UniProt AlphaFold2 models. | CC-BY 4.0 | Broad evolutionary structural search across kingdoms without redundant sequence isoforms. Default target. |
| `afdb-swissprot` | AlphaFold DB - UniProtKB/Swiss-Prot | Manually curated, expert-annotated Swiss-Prot protein models. | CC-BY 4.0 | High-confidence functional inference from experimentally characterized and reviewed protein families. |
| `afdb-proteome` | AlphaFold DB - Reference Proteomes | Complete reference proteomes across model organisms and human pathogens. | CC-BY 4.0 | Organism-wide paralogue discovery and comparative structural genomics. |
| `cath50` | CATH-Gene3D Structural Domains (50% non-redundant) | Curated Class, Architecture, Topology, and Homologous superfamily domain definitions. | CC-BY 4.0 | Domain-boundary dissection and fold/superfamily classification of multi-domain queries. |
| `mgnify_esm30` | ESM Metagenomic Atlas (30% clustered) | Metagenomic protein structure predictions clustered at 30% sequence identity. | CC-BY 4.0 | Discovery of uncultured environmental and microbiome enzyme variants. |
| `BFVD` | Big Fantastic Virus Database | Predicted viral protein structures derived from global virome sequences. | CC-BY 4.0 | Viral capsid, polymerase, and host-interaction fold homology analysis. |
| `bfmd` | Big Fantastic Microbiome / Metagenomic Database | Specialized microbial and metagenomic structural assemblies. | CC-BY 4.0 | Environmental and prokaryotic fold diversity exploration. |
| `gmgcl_id` | Global Microbial Gene Catalog (GMGC) Structures | Non-redundant structural models from global microbial habitat catalogs. | CC-BY 4.0 | Microbiome-specific functional structural annotation. |

---

## 3. Mathematical Definitions of Alignment Metrics

### 3.1 Homology Probability (`Prob`)
- **Range**: `[0.000, 1.000]`
- **Definition**: Posterior probability estimated by Foldseek's empirical score distribution model that the target alignment represents a true structural homologue sharing homologous fold topology (comparable to TM-score >= 0.5 / shared SCOP/CATH superfamily) rather than a random geometric match.
- **Interpretation**:
  - `Prob >= 0.90`: High-confidence homologous fold.
  - `0.80 <= Prob < 0.90`: Strong structural similarity; inspect domain coverage and catalytic/core conservation.
  - `0.50 <= Prob < 0.80`: Moderate or partial local motif similarity; frequently arises in shared supersecondary motifs (e.g., isolated helix-turn-helix or beta-hairpin segments).
  - `Prob < 0.50`: Low confidence; likely non-homologous background geometry.

### 3.2 Expectation Value (`E-value`)
- **Range**: `[0.0, +inf)` (sorted ascending; fallback default `1000.0`)
- **Definition**: Expected number of chance alignments achieving an equivalent or higher 3Di+AA bit score in a database of the searched size.
- **Dashboard Transformation**: Plotted on logarithmic scale as:
  `neg_log10_evalue = -log10(max(evalue, 1e-300))`
- **Interpretation**:
  - `E-value < 1e-10` (`-log10(E) > 10`): Statistically decisive structural match.
  - `1e-10 <= E-value <= 1e-3`: Statistically significant match; cross-check `Prob` and `Q-Cov`.
  - `E-value > 0.01`: Weak statistical support unless corroborated by specific local motif geometry.

### 3.3 Query Coverage (`Q-Cov`)
- **Range**: `[0.0%, 100.0%]` (or `"N/A"` if `qLen == 0` or `qEndPos <= qStartPos`)
- **Formula**:
  ```text
  aligned_q_residues = (qEndPos - qStartPos) + 1
  Q-Cov (%) = min((aligned_q_residues / qLen) * 100.0, 100.0)
  ```
- **Interpretation**: Quantifies the contiguous span of the query polypeptide (`qStartPos` to `qEndPos`, 1-indexed) participating in the structural alignment relative to full query length `qLen`. Distinguishes global fold matches (`Q-Cov >= 70%`) from domain-specific or local motif alignments (`Q-Cov < 50%`).

### 3.4 Sequence Identity (`Seq Identity`)
- **Range**: `[0.0%, 100.0%]`
- **Normalization**: Foldseek API deployments emit `seqId` / `fident` either as a fraction in `[0.0, 1.0]` or as a percentage in `[0.0, 100.0]`. Values `<= 1.0` are normalized by multiplying by `100.0` (with `0.0` preserved as `0.0%`).
- **Interpretation**: Percentage of aligned residue pairs sharing identical amino acids. When `Prob >= 0.80` and `Seq Identity < 30.0%`, the target lies in the **Twilight Zone** of remote structural homology where primary sequence searches (e.g., BLASTp) routinely fail to detect evolutionary relatedness.

### 3.5 Alignment Length (`Aln Length`)
- **Definition**: Total number of aligned columns (`alnLength` / `alnLen`) across the local 3Di+AA alignment, including gap columns where applicable.

---

## 4. Homology Regime Classification Rules

| Homology Regime | Quantitative Criterion | Biological Interpretation |
| :--- | :--- | :--- |
| `High-Identity Homolog` | `Prob >= 0.80` and `Seq Identity >= 30.0%` | Close orthologue or paralogue detectable by both sequence and structure. |
| `Remote Structural Homolog (Twilight Zone)` | `Prob >= 0.80` and `Seq Identity < 30.0%` | Divergent evolutionary relative conserving 3D fold topology despite extensive sequence drift. |
| `Local / Partial Motif Match` | `Prob >= 0.50` and (`Q-Cov < 50.0%` or `Prob < 0.80`) | Shared sub-domain, repeat module, or local structural motif without high-confidence global fold equivalence. |
| `Low Confidence` | `Prob < 0.50` | Background geometric similarity; treat as putative false positive without orthogonal evidence. |

---

## 5. Primary Citation

van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Soding, J., & Steinegger, M. (2023). Fast and accurate protein structure search with Foldseek. *Nature Biotechnology*, 41(12), 1673-1682. https://doi.org/10.1038/s41587-023-01773-0
