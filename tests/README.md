# Tests

Every `src/` module has a pytest module that exercises its **logic**, not merely
that it runs (a project rule, see [`CLAUDE.md`](../CLAUDE.md)). Tests are the
main guard on scientific invariants: frozen values, frozen thresholds, data
provenance, and figure-to-source fidelity.

## Commands

```bash
# Full suite
python -m pytest -q

# The six canonical poster figures
python -m pytest tests/test_poster_crispr_discovery_v1.py \
                 tests/test_poster_hero_heatmap_v6.py \
                 tests/test_poster_network_mechanism_v4.py \
                 tests/test_poster_pathway_v2.py \
                 tests/test_poster_depmap_v2.py \
                 tests/test_poster_druggability_v1.py -q

# Freeze / scientific-integrity checks only
python -m pytest -q -k "frozen or freeze or sha256 or shortlist"

# Verify the published figures and their hashes without re-rendering
python scripts/poster/build_all.py --check
```

## What the test suite protects

| Category | Examples |
|---|---|
| **Freeze integrity** | frozen shortlist SHA-256 unchanged; `science-freeze-2026-08-15` tag resolution; frozen tables byte-stable |
| **CRISPR invariants** | 19,103 fitted genes; 13 significant sensitising hits; negative effect = sensitising; focus-gene ranks 1/4/8/12 |
| **Expression heatmap** | exact 29 biological rows; dataset block membership; within-dataset gene-wise z-scoring; FASR excluded |
| **Network provenance** | every displayed edge exists in the raw STRING files; 47 nodes / 147 edges / 3 components; VEZF1 degree 0; TLK2 in a separate component; shortest paths computed, not hardcoded |
| **Pathway fidelity** | every NES/FDR matches the frozen GSEA tables; no hand-typed values; contexts never mislabelled as resistance |
| **DepMap** | exact 11-line ER+/luminal subset; probability threshold 0.5 from config; counts 0/9/0/3; scatter x = −1 × frozen effect |
| **Structure** | PDB IDs traceable to the audited evidence table; **no docking/affinity code**; ATP analog and activity-based probe never described as drugs; VEZF1 given no fabricated structure |
| **Claim hygiene** | no evaluative wording drawn in figures; no absolute "druggable/undruggable"; no normal-tissue-safety implication |

## Notes

- Tests read real project data, so they must be run from the repository root
  with the `bc` environment active.
- **The figure tests re-render figures into `results/figures/`.** Because PDF and
  SVG output is not byte-reproducible (and figure 03's PNG is not either — see
  [`poster/README.md`](../poster/README.md)), running the suite leaves those
  committed artifacts showing as modified even though no figure changed. Discard
  the churn with:

  ```bash
  git checkout -- results/figures/
  ```

  The five reproducible PNGs are unaffected and are checked byte-for-byte by
  `tests/test_poster_release_integrity.py`.
- DepMap tests read a cached 11-line extract; if
  `results/tables/poster_depmap_v1/` is deleted they rebuild it from the external
  DepMap matrices, which requires those files locally (see
  [`data/README.md`](../data/README.md)).
