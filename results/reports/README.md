# results/reports

Narrative `.md` reports and per-figure `NOTE.md` files, one subdirectory per
analysis phase. Start with the reports below.

## Most important reports

| Report | Why it matters |
|---|---|
| [`../../docs/THERAPEUTIC_SHORTLIST_FREEZE.md`](../../docs/THERAPEUTIC_SHORTLIST_FREEZE.md) | Defines the **historical frozen shortlist** (USP34, VEZF1, EML5, CITED2) and the gate that produced it |
| [`post_audit/`](post_audit/) | The post-audit sensitivity analysis that challenged the RNA-corroboration gate and produced the **current four poster focus genes** |
| [`../../docs/FINAL_PUBLIC_REPO_AUDIT.md`](../../docs/FINAL_PUBLIC_REPO_AUDIT.md) | The final independent public-release audit (scientific integrity, reproducibility, security, usability) |
| [`../../docs/DATA_PROVENANCE.md`](../../docs/DATA_PROVENANCE.md) | Where every dataset came from, with checksums |
| [`independent_validation/`](independent_validation/) | TCGA + DepMap orthogonal checks, including the DepMap 26Q1 access/verification status |

## Per-figure notes for the canonical poster figures

Each canonical figure has a short `NOTE.md` recording its exact sources,
derivations, conventions and caveats:

- `poster_crispr_discovery_v1/NOTE.md`
- `poster_hero_heatmap_v6/NOTE.md`
- `poster_network_mechanism_v4/NOTE.md`
- `poster_pathway_v2/NOTE.md`
- `poster_depmap_v2/NOTE.md`
- `poster_druggability_v1/NOTE.md`

Notes for superseded iterations (`..._v1`/`_v2`/`_v3`/`_v5`) are retained so the
reasoning behind each revision stays inspectable.
