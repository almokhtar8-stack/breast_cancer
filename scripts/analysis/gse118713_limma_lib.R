# Limma differential-expression functions for GSE118713 Phase 2B.
#
# Source: filtered gene-level TPM matrix and sample metadata written by
# the Python Phase 2B pipeline (src/gse118713_expression_filter.py,
# src/gse118713_prep.py). Implements PREANALYSIS.md's 2026-08-05 Phase 2B
# statistical-plan amendment: log2(TPM + 1), unpaired three-group design
# (MCF7/TAMR/FASR, 3 replicates each, replicate numbers not matched
# blocks), limma lmFit/contrasts.fit/eBayes, three named contrasts, BH
# correction applied separately per contrast.
#
# This file defines functions only and has no top-level side effects
# beyond loading limma/statmod -- it is sourced by gse118713_limma.R and
# by tests, never executed directly.

suppressMessages(library(limma))
suppressMessages(library(statmod))

EXPECTED_GROUPS <- c("MCF7", "TAMR", "FASR")
EXPECTED_REPLICATES_PER_GROUP <- 3
CONTRAST_NAMES <- c("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR")

validate_metadata <- function(meta) {
  required_cols <- c("sample_id", "group")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(paste("metadata missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  expected_n <- length(EXPECTED_GROUPS) * EXPECTED_REPLICATES_PER_GROUP
  if (nrow(meta) != expected_n) {
    stop(paste("expected exactly", expected_n, "samples in metadata, got", nrow(meta)))
  }
  if (any(duplicated(meta$sample_id))) {
    stop(paste(sum(duplicated(meta$sample_id)), "duplicated sample_id values in metadata"))
  }

  group_counts <- table(meta$group)
  missing_groups <- setdiff(EXPECTED_GROUPS, names(group_counts))
  if (length(missing_groups) > 0) {
    stop(paste("metadata missing expected groups:", paste(missing_groups, collapse = ", ")))
  }
  extra_groups <- setdiff(names(group_counts), EXPECTED_GROUPS)
  if (length(extra_groups) > 0) {
    stop(paste("metadata contains unexpected groups:", paste(extra_groups, collapse = ", ")))
  }
  wrong <- group_counts[EXPECTED_GROUPS] != EXPECTED_REPLICATES_PER_GROUP
  if (any(wrong)) {
    stop(paste(
      "expected exactly", EXPECTED_REPLICATES_PER_GROUP, "replicates per group, got:",
      paste(names(group_counts)[wrong], "=", as.integer(group_counts[wrong]), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

validate_expression <- function(expr, meta) {
  if (!"gene_id" %in% colnames(expr)) {
    stop("expression table missing gene_id column")
  }
  if (any(duplicated(expr$gene_id))) {
    stop(paste(sum(duplicated(expr$gene_id)), "duplicated gene_id values in expression table"))
  }
  sample_cols <- meta$sample_id
  missing_cols <- setdiff(sample_cols, colnames(expr))
  if (length(missing_cols) > 0) {
    stop(paste("expression table missing sample columns:", paste(missing_cols, collapse = ", ")))
  }
  mat <- as.matrix(expr[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) {
    stop("non-finite TPM value in expression table")
  }
  if (any(mat < 0)) {
    stop("negative TPM value in expression table")
  }
  invisible(TRUE)
}

# Unpaired three-group design: group is the only term, replicate number is
# never used as a blocking factor.
build_design <- function(meta) {
  group <- factor(meta$group, levels = EXPECTED_GROUPS)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  rownames(design) <- meta$sample_id
  design
}

build_log2_matrix <- function(expr, meta) {
  sample_cols <- meta$sample_id
  mat <- as.matrix(expr[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- expr$gene_id
  log2(mat + 1)
}

fit_limma <- function(log2_mat, design) {
  fit <- lmFit(log2_mat, design)
  contrast_matrix <- makeContrasts(
    TAMR_vs_MCF7 = TAMR - MCF7,
    FASR_vs_MCF7 = FASR - MCF7,
    TAMR_vs_FASR = TAMR - FASR,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast_matrix)
  eBayes(fit2, trend = FALSE, robust = FALSE)
}

# Per-contrast table with BH correction computed independently for this
# contrast alone (topTable's adj.P.Val is computed within a single call,
# i.e. within a single coefficient/contrast).
extract_contrast <- function(fit2, contrast_name, expr) {
  tt <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "none")
  se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, contrast_name]
  direction <- ifelse(tt$logFC > 0, "up", ifelse(tt$logFC < 0, "down", "no_change"))
  data.frame(
    gene_id = expr$gene_id,
    gene_symbol = expr$gene_symbol,
    log2fc = tt$logFC,
    se = se,
    moderated_t = tt$t,
    p_value = tt$P.Value,
    fdr = tt$adj.P.Val,
    ave_expr = tt$AveExpr,
    contrast = contrast_name,
    direction = direction,
    stringsAsFactors = FALSE
  )
}

# Sort by: contrast (in the preregistered order), FDR ascending, absolute
# effect descending, gene_id ascending.
sort_de_table <- function(combined) {
  combined$contrast <- factor(combined$contrast, levels = CONTRAST_NAMES)
  ord <- order(combined$contrast, combined$fdr, -abs(combined$log2fc), combined$gene_id)
  combined <- combined[ord, ]
  combined$contrast <- as.character(combined$contrast)
  rownames(combined) <- NULL
  combined
}

run_limma_de <- function(expr, meta) {
  validate_metadata(meta)
  validate_expression(expr, meta)

  log2_mat <- build_log2_matrix(expr, meta)
  design <- build_design(meta)
  fit2 <- fit_limma(log2_mat, design)

  combined <- do.call(rbind, lapply(CONTRAST_NAMES, function(cn) extract_contrast(fit2, cn, expr)))
  sort_de_table(combined)
}

# Post-hoc reporting-stage redaction of preregistered blind-control genes.
# Must only ever be called AFTER run_limma_de() has produced a table fit and
# BH-corrected on the complete gene set -- removing rows here has no effect
# on any other gene's logFC, moderated t, p-value, or FDR, because those are
# already finalized. This function never receives or returns gene identity
# information beyond the count of genes withheld; callers must not log or
# print the removed rows.
#
# Fails loudly (count-only error message, no identity) rather than silently
# under-redacting: every requested ID must be a duplicate-free entry that
# matches exactly one gene actually present in de_table. A typo'd or
# duplicated configured ID is a hard error, not a quietly smaller redaction.
redact_blinded_genes <- function(de_table, blinded_gene_ids) {
  genes_fitted <- length(unique(de_table$gene_id))
  if (length(blinded_gene_ids) == 0) {
    return(list(table = de_table, genes_fitted = genes_fitted, genes_withheld = 0L, genes_reported = genes_fitted))
  }
  if (any(duplicated(blinded_gene_ids))) {
    stop(paste(sum(duplicated(blinded_gene_ids)), "duplicate entries in blinded gene id list"))
  }
  keep <- !(de_table$gene_id %in% blinded_gene_ids)
  reported <- de_table[keep, , drop = FALSE]
  rownames(reported) <- NULL
  genes_withheld <- genes_fitted - length(unique(reported$gene_id))
  if (genes_withheld != length(blinded_gene_ids)) {
    stop(sprintf(
      "blinded gene redaction count mismatch: expected to withhold %d configured gene id(s) but matched %d in the fitted table -- refusing to write a partially- or un-redacted result",
      length(blinded_gene_ids), genes_withheld
    ))
  }
  list(
    table = reported,
    genes_fitted = genes_fitted,
    genes_withheld = genes_withheld,
    genes_reported = length(unique(reported$gene_id))
  )
}

# Writes via a same-directory temp file + atomic rename so a killed or
# crashed process can never leave a truncated or partially-written file at
# output_path. file.rename() can return FALSE without raising (e.g. a
# cross-device destination or a permissions issue), which would otherwise
# leave a stale previous file silently in place while the caller believes
# the write succeeded -- so its result is checked explicitly.
.atomic_write <- function(output_path, write_fn) {
  tmp_path <- paste0(output_path, ".tmp", Sys.getpid())
  write_fn(tmp_path)
  renamed <- file.rename(tmp_path, output_path)
  if (!isTRUE(renamed)) {
    file.remove(tmp_path)
    stop(sprintf("failed to atomically write %s (temp file left removed)", output_path))
  }
}

write_de_table <- function(de_table, output_path) {
  .atomic_write(output_path, function(tmp_path) {
    con <- gzfile(tmp_path, "wt")
    on.exit(close(con))
    write.table(de_table, con, sep = "\t", row.names = FALSE, quote = FALSE)
  })
}

# Counts only -- never gene identities.
write_redaction_record <- function(redaction, output_path) {
  record <- data.frame(
    genes_fitted = redaction$genes_fitted,
    genes_withheld = redaction$genes_withheld,
    genes_reported = redaction$genes_reported
  )
  .atomic_write(output_path, function(tmp_path) {
    write.table(record, tmp_path, sep = "\t", row.names = FALSE, quote = FALSE)
  })
}
