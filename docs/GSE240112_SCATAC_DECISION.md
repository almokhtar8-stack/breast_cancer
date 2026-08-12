# GSE240112 Phase 16: scATAC decision

**Decision: DEFERRED.**

**Reason:** GSE240112 includes 8 matched scATAC-seq libraries
(GSM7681693-700: 2 NT, 3 PT, 3 RT), but the public GEO supplementary files
are raw Cell Ranger ATAC outputs only (`*_fragments.tsv.gz`,
`*_peaks.bed.gz`, `*_singlecell.csv.gz` per sample) -- there is no
processed, cell-annotated ATAC object in GEO. A processed object
(`TTs.cancer.atac.motif.0615.rds`, motif-level, cancer cells only) is
linked from the author GitHub README on Dropbox, but was not downloaded
or inspected in this run: even in the best case it would still require
(1) confirming its barcode/sample labels map cleanly onto this analysis's
PT1-3/RT1-3 identity (not verified given the RNA object's own
raw-vs-processed sample-labeling quirks documented in
`docs/GSE240112_DATA_AUDIT.md` section 3), and (2) either using the
authors' own motif/gene-activity summarization as-is or performing
fresh peak calling / gene-activity scoring from the raw fragment files if
it does not cleanly support a USP34-locus or 13-candidate-loci query --
none of which is "straightforward within the existing environment" as
required by the Phase 16 gating criteria.

Per the user's explicit instruction, "a completed rigorous scRNA analysis
is much more important than partial scATAC" -- the scRNA primary and
secondary analyses (Phases 1-15) were prioritized to full rigor within
this run's time budget, and scATAC is deferred to a future, dedicated run
rather than attempted partially.
