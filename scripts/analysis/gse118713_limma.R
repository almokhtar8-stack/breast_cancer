# CLI entry point for GSE118713 Phase 2B limma differential expression.
#
# Usage:
#   Rscript scripts/analysis/gse118713_limma.R <expression_tsv_gz> <metadata_tsv> <output_tsv_gz> <blinded_gene_ids_csv> <redaction_record_tsv>
#
# <expression_tsv_gz>: filtered gene-level TPM matrix (gene_id,
#   gene_symbol, one column per sample_id), as written by
#   src/gse118713_expression_filter.py.
# <metadata_tsv>: sample metadata (sample_id, group, ...), as written by
#   src/gse118713_prep.py.
# <output_tsv_gz>: destination for the per-gene, per-contrast
#   differential-expression table.
# <blinded_gene_ids_csv>: comma-separated Ensembl gene IDs to withhold from
#   <output_tsv_gz> AFTER fitting and BH correction on the complete gene
#   set (see redact_blinded_genes() in gse118713_limma_lib.R -- the fit and
#   every other gene's statistics are unaffected). Pass an empty string to
#   explicitly redact nothing. This argument is mandatory (not optional)
#   precisely so no caller can produce a reportable result without making
#   an explicit blinding decision -- see PREANALYSIS.md S5 / CLAUDE.md.
# <redaction_record_tsv>: destination for the genes-fitted /
#   genes-withheld / genes-reported counts. Never contains gene identities.
#
# Any failure (bad arguments, invalid input, fitting error, or a
# blinded-gene-id that fails to match exactly one fitted gene) is reported
# to stderr with a nonzero exit status rather than a silent partial result.

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "gse118713_limma_lib.R"))

main <- function(argv) {
  if (length(argv) != 5) {
    stop(paste(
      "usage: Rscript gse118713_limma.R <expression_tsv_gz> <metadata_tsv> <output_tsv_gz>",
      "<blinded_gene_ids_csv> <redaction_record_tsv>",
      "(pass an empty string for blinded_gene_ids_csv to explicitly redact nothing)"
    ))
  }
  expression_path <- argv[1]
  metadata_path <- argv[2]
  output_path <- argv[3]
  blinded_gene_ids_csv <- argv[4]
  redaction_record_path <- argv[5]

  blinded_gene_ids <- if (nzchar(trimws(blinded_gene_ids_csv))) {
    trimws(strsplit(blinded_gene_ids_csv, ",")[[1]])
  } else {
    character(0)
  }
  if (any(!nzchar(blinded_gene_ids))) {
    stop("blinded_gene_ids_csv contains a blank entry (e.g. a stray comma) -- refusing to guess intent")
  }

  meta <- read.delim(metadata_path, stringsAsFactors = FALSE)
  expr <- read.delim(expression_path, stringsAsFactors = FALSE)

  de_table <- run_limma_de(expr, meta)
  redaction <- redact_blinded_genes(de_table, blinded_gene_ids)
  write_redaction_record(redaction, redaction_record_path)
  write_de_table(redaction$table, output_path)

  cat(sprintf(
    "gse118713_limma.R: fit %d genes, withheld %d preregistered blind gene(s), wrote %d rows (%d reported genes x %d contrasts) to %s\n",
    redaction$genes_fitted, redaction$genes_withheld, nrow(redaction$table), redaction$genes_reported,
    length(CONTRAST_NAMES), output_path
  ))
}

result <- tryCatch(
  {
    main(commandArgs(trailingOnly = TRUE))
    TRUE
  },
  error = function(e) {
    message("ERROR: gse118713_limma.R failed: ", conditionMessage(e))
    FALSE
  }
)

if (!isTRUE(result)) {
  quit(status = 1, save = "no")
}
