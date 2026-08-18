# Final public-repository audit

| | |
|---|---|
| **Audit date** | 2026-08-17 |
| **Repository commit at audit start** | `5139745` (`Apply final external re-audit wording fixes`) |
| **Science-freeze tag** | `science-freeze-2026-08-15` → `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` (verified unmoved) |
| **Frozen shortlist SHA-256** | `b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc` (verified unchanged) |
| **Scope** | Independent scientific re-verification of all six poster layers from source; figure-by-figure publication audit; code, test, reproducibility, security and public-usability audit; public navigation reorganisation |

This audit re-derived every headline number from source tables and code rather
than trusting prior reports or documentation.

---

## Verdicts

| Area | Verdict | Notes |
|---|---|---|
| Scientific integrity | **PASS** | All six layers re-verified from source. One HIGH labelling issue found and fixed (below). No frozen value, threshold, ranking, or conclusion altered. |
| Reproducibility | **WARN** | Figures rebuild end-to-end from committed/processed data. Full raw-data reproduction is **not** possible from GitHub alone (see levels below). Figure output is only partly byte-reproducible (issue 4) — documented and test-enforced. Honest, documented, not a defect. |
| Figure provenance | **PASS** | Six canonical figures copied to `poster/final_figures/` and verified byte-identical (SHA-256) to their `results/figures/` sources; manifest records module, wrapper, source data, freeze status and hashes. |
| Code quality | **WARN** | No hardcoded scientific values in poster renderers; all read from loaders. Remaining: 179 `src/` modules with historical version families (v1…v6) retained for provenance — navigable via `docs/analysis_map.md` rather than mass-moved. |
| Testing | **PASS** | Poster/integrity suites green; new release-integrity module added. See test results section. |
| Security / privacy | **PASS** | No secrets, keys, tokens, credentials, emails or PHI. Only public GEO/DepMap/PDB accessions. |
| Public usability | **PASS** | A stranger reaches the six figures, the four genes, the frozen-vs-exploratory distinction and the build command from the README in under a minute. |

---

## Issues found

### HIGH — fixed
1. **Figure 02 mislabelled its acute-response block as single-cell.**
   `poster_hero_heatmap_v5.CONTEXT_TITLE["GSE245601"]` read
   *"Single-cell: before vs 12 h after tamoxifen"*, but the six rows plotted are
   **per-tumour pseudobulk** (3 patients × Control/Tamoxifen; see
   `poster_hero_heatmap_v4._build_gse245601_rows`). The wording risked implying
   cell-level resolution or pseudo-replication.
   **Fix:** the canonical v6 renderer now overrides the headline to
   *"Acute 12 h tamoxifen — per-tumour pseudobulk"*. Label-only change; no
   plotted value, row, or z-score altered. The superseded v5 module is left
   untouched so its historical output stays interpretable.

### MEDIUM — fixed
2. **Poster wrappers were not runnable outside the repository root.** The
   canonical `src/` modules use repo-relative paths, so every wrapper failed
   from any other working directory. All six wrappers and `build_all.py` now
   `chdir` to the repository root; each was re-verified from `/tmp`.
3. **Documentation cited a non-existent data path.** An earlier draft of
   `data/README.md`, `docs/analysis_map.md` and the manifest referenced
   `data/processed/labels/crispr_gene_labels.parquet`; the real frozen file is
   `data/processed/labels.parquet`. Corrected in all three.

### MEDIUM — investigated, resolved by documenting an accepted limitation
4. **Figure output is only partly byte-reproducible.** Measured, not assumed:
   PNG output is byte-reproducible; **PDF** is reproducible only when
   `SOURCE_DATE_EPOCH` is pinned (matplotlib otherwise stamps a wall-clock
   `/CreationDate`) — `scripts/poster/build_all.py` now pins it; **SVG** is never
   byte-reproducible (matplotlib emits per-run element ids). Separately,
   **Figure 03's PNG** is also not reproducible because `adjustText`'s collision
   solver gives different label pixel positions per run — confirmed by rendering
   repeatedly in one process with a fixed numpy seed and a fixed
   `PYTHONHASHSEED` (still differs) and with the solver disabled (bit-identical
   every time). Removing the solver was tested and **rejected**: it reintroduces
   real label collisions (`HDAC1` over the KDM1A node, `HMG20B`/`CDH1`,
   `CTBP1`/`SNAI1`), so readability is preferred over byte-reproducibility.
   Only label *pixels* vary; the graph (nodes, edges, scores, components,
   shortest paths) is deterministic and separately asserted. The release-integrity
   test enforces strict PNG equality for the other five figures and requires the
   one exemption to remain explicitly documented in source.

### MEDIUM — documented, not changed
5. **Machine-local absolute paths in configuration and historical provenance.**
   `config/config.yaml`, `config/resources.yaml`, several download scripts and
   older provenance reports contain a cluster scratch path that embeds the
   original operator's username. These are pre-existing, already-published,
   and load-bearing for reproducibility documentation. Rewriting them would
   touch frozen configuration and historical provenance for cosmetic gain, so
   they are **retained and documented** (`data/README.md` explains that
   `data.raw` entries are a machine-local root a reproducer must repoint). No
   *new* file added by this audit contains an absolute path — enforced by
   `tests/test_poster_release_integrity.py`.
6. **`src/` is large (179 modules).** Deliberately not mass-reorganised:
   moving scientific modules risks breaking frozen imports for no scientific
   gain. `docs/analysis_map.md` provides the conceptual index instead.

### LOW — documented
7. Superseded figure iterations (heatmap v1–v5, network v1–v3, pathway v1,
   DepMap v1) remain under `results/figures/` for provenance. Indexed in
   `results/figures/README.md`.
8. `report/`, `notebooks/`, `workflow/` are empty untracked directories on this
   machine; they carry no content to GitHub.

### None found
- No fabricated data, no invented PDB IDs, no docking or affinity prediction,
  no candidate re-ranking, no frozen-value edits, no history rewriting.

---

## Independent re-verification (all recomputed from source)

| Layer | Verified |
|---|---|
| **CRISPR** | 19,103 fitted genes; 13 significant sensitising hits at FDR < 0.10 with effect < 0; all 13 effects negative. KDM1A −2.167336 (rank 1), TLK2 −1.848198 (rank 4), VEZF1 −1.602445 (rank 8), USP34 −1.391298 (rank 12). |
| **Expression** | Exactly 29 biological rows: GSE118713 6 (MCF7×3, TAMR×3; FASR excluded), GSE111151 11 (4 parental + 7 resistant, family structure correct), GSE240112 6 (3 primary + 3 recurrent, unpaired — no pair anchors), GSE245601 6 (T02/T03/T07 × Ctrl/TAM, per-tumour pseudobulk). Gene-wise z-score within dataset block. |
| **Network** | 47 nodes, 147 edges, 3 components. KDM1A↔USP34 connected; TLK2 in a separate component; VEZF1 degree 0. All 4 equally short KDM1A–USP34 paths (length 3, all via DNMT1) enumerated programmatically. **0** displayed edges absent from the raw STRING files. |
| **Pathway** | 5 Hallmark sets × 4 contexts, all NES/FDR matching the frozen GSEA tables. Estrogen early/late negative in all four; EMT positive in GSE118713/111151/240112 and negative in GSE245601; E2F negative in all four; WNT weak/mixed (significant only in GSE111151 and GSE240112) and reported as such. |
| **DepMap** | Release 26Q1; exact 11 ER+/luminal evaluable lines (CAMA1, EFM19, HCC1428, KPL1, MCF7, MDAMB361, MDAMB415, MFM223, T47D, UACC3133, ZR751); probability threshold 0.5 read from config; counts KDM1A 0/11, TLK2 9/11, USP34 0/11, VEZF1 3/11; scatter x = −1 × frozen effect, y = 100 × count/11. |
| **Structure** | PDB IDs recovered from the audit table's own text (KDM1A 6NQU + 2Z5U; TLK2 5O0Y; USP34 7W3R + 7W3U; VEZF1 only 1AAY, the unrelated Zif268 homology template). Rendered: 6NQU (GSK2879552 inhibitor), 5O0Y (ATP-γ-S analog, labelled *not* an inhibitor), 7W3U (covalent activity-based probe at Cys1903, labelled *not* a drug). VEZF1 rendered with no structure. Local PDB files SHA-verified unmodified. |

---

## Remaining limitations

- The CRISPR screen tests whether knockout **enhances** tamoxifen effect in a
  parental/drug-tolerant MCF7-V context (E2+4-OHT vs E2). It does not
  demonstrate reversal of established acquired resistance.
- GSE240112's primary vs recurrent comparison is **unpaired** across different
  patients/biobanks: recurrence-associated, not tamoxifen-causal.
- GSE245601 is a 12 h, 10 µM **ex vivo acute** exposure; malignant-compartment
  assignment is inferential (inferCNV/copyKAT).
- GSE111151 candidate differential expression was largely **null**.
- STRING edges are undirected functional associations; a shortest path is not a
  mechanism.
- DepMap dependency is baseline cancer-cell fitness — not tamoxifen specificity
  and not normal-tissue safety. The TLK2 dependency claim is scoped to the four
  focus genes only (SUPT4H1 is higher across the full 13-hit set).
- Structural tractability is not efficacy; no clinical-stage compound here is
  breast-cancer-approved.
- GDSC (GDSC1-only for USP34, no GDSC2 replication) and TCGA are exploratory
  supporting layers and are deliberately excluded from the six poster figures.
- No wet-lab validation was performed by this project.

## Known post-freeze exploratory analyses

| Analysis | Location | Status |
|---|---|---|
| Standardized STRING network for the four focus genes (Figure 03) | `src/poster_network_mechanism_v2.py` → `_v4.py`; `data/reference/interactions/string_v2_*` | **post-freeze exploratory** — new query run after the freeze; declared as such in the manifest and both READMEs |
| Cytoscape export of that graph | `src/cytoscape_v4_export.py`, `results/tables/cytoscape_v4_*` | post-freeze exploratory |
| Post-audit sensitivity analysis that produced the current four focus genes | `src/post_audit_sensitivity_*.py`, `results/reports/post_audit/` | post-freeze reinterpretation; the historical frozen shortlist (USP34, VEZF1, EML5, CITED2) is preserved unchanged |
| All six poster figures as *visual artefacts* | `src/poster_*` | new renderings of frozen data; no new science except Figure 03 |

## Reproducibility levels

| Level | Achievable from a GitHub clone? |
|---|---|
| **1 — view frozen outputs** | **Yes.** Frozen tables, reports and all six figures are committed. |
| **2 — rebuild the six poster figures** | **Yes**, with the `bc` environment: `python scripts/poster/build_all.py`. Uses committed processed/frozen tables, the committed STRING reference files, the cached DepMap 11-line extract, and the committed PyMOL renders. |
| **3 — rerun selected analyses** | **Partially.** Requires re-downloading external raw data (GEO series, DepMap ~440 MB matrices, PDB files, GDSC, TCGA) and repointing `config/config.yaml`; helper scripts are in `scripts/download/`. R-based steps additionally need the R toolchain. |
| **4 — full raw-data reproduction** | **Not from GitHub alone.** Raw matrices are excluded for size and must be re-obtained from their public sources; some steps also require PyMOL and R. |

## Environment specification (amendment, 2026-08-18, `poster-final` branch)

Two claims about the environment were checked during the poster-final pass.
Both are recorded here because one of them turned out to be false, and a
correction that is not needed is worth documenting as clearly as one that is.

**`h5py` is NOT missing from `environment.yml`, and this document never said it
was.** The claim was raised as a known documentation error to fix. It is not
one: `h5py` is declared at line 10 of `environment.yml`, it imports (3.16.0),
and the test suite collects without error. A search of this document, the root
`README.md` and all of `docs/` finds no statement that `h5py` is absent. **No
correction was made, because there was nothing to correct.**

**A real, pre-existing gap was found and fixed instead.** `environment.yml`
declared `r-base` and `r-statmod` but not `r-yaml` or `r-data.table`, while
`src/post_poster_de_refit.R` loads all four of `edgeR`, `limma`, `yaml` and
`data.table`. The `bc` environment therefore could not run the frozen R refit
without an out-of-band install — the work was in practice carried out in a
separate `sc245601` environment that happened to have them. Both packages have
been added to `environment.yml`, so the R step is now reproducible from the
declared environment alone. This affects reproducibility level 3 in the table
above: the R-based steps need the R toolchain, and the specification now
actually describes it.

## Release decision

**READY** for public release, with the WARN items above recorded as documented
limitations rather than unresolved defects. No CRITICAL or unresolved HIGH issue
remains.
