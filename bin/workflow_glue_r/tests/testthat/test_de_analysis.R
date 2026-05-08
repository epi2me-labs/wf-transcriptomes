###
# Design parsing and validation
#
# Tests covering the experiment-design checks owned by supeRglue de_analysis:
# covariate parsing, required sample-sheet columns, alias matching, reference
# level selection, and the transcript metadata needed for DEXSeq grouping.

# Covariates CLI arg is comma-separated string with possible whitespace/empty values.
testthat::test_that("covariates parsed and trimmed", {
    testthat::expect_equal(
        de_parse_covariates(" batch, sex ,, site "),
        c("batch", "sex", "site")
    )
})

# Sample sheet must have 'alias' column matching SE colnames, condition column, and all covariates.
# Aliases must be unique. Missing columns should fail fast before passing to DESeq2/DRIMSeq.
testthat::test_that("sample sheet structure validated", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = NULL
    )

    missing_alias <- data.frame(condition = rep(c("control", "treated"), each = 3))
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, missing_alias, argv),
        "Sample sheet must contain an 'alias' column"
    )

    sample_df <- data.frame(
        alias = colnames(tx_se),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "Sample sheet must contain the 'condition' column"
    )

    sample_df$condition <- rep(c("control", "treated"), each = 3)
    argv$covariates <- "batch,site"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "Missing covariate columns: site"
    )

    duplicate_aliases <- sample_df
    duplicate_aliases$alias[2] <- duplicate_aliases$alias[1]
    argv$covariates <- "batch"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, duplicate_aliases, argv),
        "Sample sheet aliases must be unique"
    )
})

# R formulas break on spaces/hyphens in column names (e.g., ~ `time point` produces errors).
# Validate condition_column and covariate names are safe (alphanumeric + underscore).
testthat::test_that("formula-unsafe column names rejected", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "time point",
        covariates = "batch",
        reference_level = NULL
    )

    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "Design column names must be safe for R formulas"
    )

    argv$condition_column <- "condition"
    argv$covariates <- "batch-id"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "Invalid names: batch-id"
    )
})

# Sample aliases CAN have spaces/hyphens (they're not used in formulas, just for matching).
# Sample sheet rows can be in different order than SE columns - should reorder automatically.
testthat::test_that("non-syntactic aliases allowed, sheets reordered", {
    sample_names <- c(
        "1control",
        "control-2",
        "control rep 3",
        "4treated",
        "treated-5",
        "treated rep 6"
    )
    tx_se <- make_test_tx_se(sample_names = sample_names)
    gene_se <- make_test_gene_se(sample_names = sample_names)
    sample_df <- data.frame(
        barcode = sprintf("barcode%02d", seq_along(sample_names)),
        sample = paste0("sample_", seq_along(sample_names)),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        alias = sample_names,
        stringsAsFactors = FALSE
    )
    sample_df <- sample_df[c(6, 3, 5, 1, 4, 2), , drop = FALSE]
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control"
    )

    validated <- de_validate_inputs(tx_se, gene_se, sample_df, argv)

    testthat::expect_equal(validated$sample_df$alias, sample_names)
    testthat::expect_equal(
        as.character(validated$sample_df$condition),
        rep(c("control", "treated"), each = 3)
    )
})

# If --reference_level not specified, default to "control" if it exists in condition column.
# Sample aliases must exactly match SE column names (no missing, no extra).
testthat::test_that("reference level defaults to control", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = NULL
    )
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )

    validated <- de_validate_inputs(tx_se, gene_se, sample_df, argv)
    testthat::expect_equal(validated$reference_level, "control")

    bad_alias_df <- sample_df
    bad_alias_df$alias[1] <- "missing_alias"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, bad_alias_df, argv),
        "Sample sheet aliases do not match"
    )
})

# If no "control" level exists, user MUST provide --reference_level explicitly.
# Specified reference_level must exist in the condition column values.
testthat::test_that("explicit reference level required without control", {
    tx_se <- make_test_tx_se(
        sample_names = c(
            "baseline_rep1", "baseline_rep2", "baseline_rep3",
            "treated_rep1", "treated_rep2", "treated_rep3"
        )
    )
    gene_se <- make_test_gene_se(
        sample_names = c(
            "baseline_rep1", "baseline_rep2", "baseline_rep3",
            "treated_rep1", "treated_rep2", "treated_rep3"
        )
    )
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("baseline", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = NULL
    )

    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "Provide --reference_level"
    )

    argv$reference_level <- "does_not_exist"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, sample_df, argv),
        "requested reference level is not present"
    )
})

# Count matrices must be numeric, non-NA, non-zero totals per sample.
# Transcript and gene SE sample names must match exactly (same order, same values).
testthat::test_that("unusable count data rejected", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control"
    )
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )

    duplicated_tx_se <- tx_se
    colnames(duplicated_tx_se)[2] <- colnames(duplicated_tx_se)[1]
    testthat::expect_error(
        de_validate_inputs(duplicated_tx_se, gene_se, sample_df, argv),
        "Transcript RDS sample names must be unique"
    )

    mismatched_gene_se <- gene_se
    colnames(mismatched_gene_se)[1] <- "missing_from_transcripts"
    testthat::expect_error(
        de_validate_inputs(tx_se, mismatched_gene_se, sample_df, argv),
        "Transcript and gene RDS sample names must match"
    )

    zero_count_tx_se <- tx_se
    SummarizedExperiment::assay(zero_count_tx_se, "counts")[, "control_rep1"] <- 0
    testthat::expect_error(
        de_validate_inputs(zero_count_tx_se, gene_se, sample_df, argv),
        "zero total counts"
    )

    character_count_tx_se <- tx_se
    tx_counts_char <- SummarizedExperiment::assay(character_count_tx_se, "counts")
    storage.mode(tx_counts_char) <- "character"
    SummarizedExperiment::assay(character_count_tx_se, "counts") <- tx_counts_char
    testthat::expect_error(
        de_validate_inputs(character_count_tx_se, gene_se, sample_df, argv),
        "Count matrices must be numeric"
    )

    na_count_gene_se <- gene_se
    SummarizedExperiment::assay(na_count_gene_se, "counts")[1, 1] <- NA
    testthat::expect_error(
        de_validate_inputs(tx_se, na_count_gene_se, sample_df, argv),
        "must not contain NA values"
    )
})

# DRIMSeq (DTU analysis) requires GENEID metadata column to group transcripts by gene.
# Fail fast if bambu output is missing this column.
testthat::test_that("transcript SE without GENEID rejected", {
    fixture_dir <- tempfile("de-no-geneid-")
    dir.create(fixture_dir)

    tx_se <- make_test_tx_se(include_geneid = FALSE)
    gene_se <- make_test_gene_se()
    tx_path <- file.path(fixture_dir, "tx.rds")
    gene_path <- file.path(fixture_dir, "gene.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = colnames(tx_se),
            condition = rep(c("control", "treated"), each = 3),
            batch = rep(c("b1", "b2", "b1"), 2),
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    argv <- list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet,
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control",
        out_dir = file.path(fixture_dir, "out")
    )

    testthat::expect_error(
        main_run_de_analysis(
            argv,
            deseq_runner = function(...) stop("runner should not be called"),
            dexseq_runner = function(...) stop("runner should not be called")
        ),
        "Transcript rowData must contain GENEID"
    )
})

# DE analysis requires ≥2 replicates per condition level and ≥2 condition levels.
# RDS files must be valid R objects (not corrupt/malformed).
testthat::test_that("underspecified designs rejected", {
    fixture_dir <- tempfile("de-bad-data-")
    dir.create(fixture_dir)

    underspecified_samples <- c("control_rep1", "treated_rep1", "treated_rep2")
    tx_se <- make_test_tx_se(sample_names = underspecified_samples)
    gene_se <- make_test_gene_se(sample_names = underspecified_samples)
    tx_path <- file.path(fixture_dir, "tx.rds")
    gene_path <- file.path(fixture_dir, "gene.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = underspecified_samples,
            condition = c("control", "treated", "treated"),
            batch = c("b1", "b1", "b2"),
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    argv <- list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet,
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control",
        out_dir = file.path(fixture_dir, "out")
    )

    testthat::expect_error(
        suppressWarnings(main_run_de_analysis(
            argv,
            deseq_runner = function(...) stop("runner should not be called"),
            dexseq_runner = function(...) stop("runner should not be called")
        )),
        "fewer than 2 replicates"
    )

    single_condition_sheet <- file.path(fixture_dir, "single_condition.csv")
    utils::write.csv(
        data.frame(
            alias = colnames(make_test_tx_se()),
            condition = rep("control", 6),
            batch = rep(c("b1", "b2", "b1"), 2),
            stringsAsFactors = FALSE
        ),
        single_condition_sheet,
        row.names = FALSE,
        quote = FALSE
    )
    valid_tx_path <- file.path(fixture_dir, "valid_tx.rds")
    valid_gene_path <- file.path(fixture_dir, "valid_gene.rds")
    saveRDS(make_test_tx_se(), valid_tx_path)
    saveRDS(make_test_gene_se(), valid_gene_path)
    argv$transcript_rds <- valid_tx_path
    argv$gene_rds <- valid_gene_path
    argv$sample_sheet <- single_condition_sheet
    argv$out_dir <- file.path(fixture_dir, "single-out")
    testthat::expect_error(
        main_run_de_analysis(
            argv,
            deseq_runner = function(...) stop("runner should not be called"),
            dexseq_runner = function(...) stop("runner should not be called")
        ),
        "requires at least two condition levels"
    )

    malformed_rds <- file.path(fixture_dir, "not-an-rds.rds")
    writeLines("not an RDS file", malformed_rds)
    argv$transcript_rds <- malformed_rds
    argv$sample_sheet <- sample_sheet
    argv$out_dir <- file.path(fixture_dir, "malformed-out")
    testthat::expect_error(
        main_run_de_analysis(
            argv,
            deseq_runner = function(...) stop("runner should not be called"),
            dexseq_runner = function(...) stop("runner should not be called")
        )
    )
})

###
# Contrast planning and output writing
#
# Mock DESeq2/DRIMSeq to avoid slow runtime and test workflow logic:
# - One contrast created per non-reference condition level (treated vs control, treated2 vs control)
# - Each contrast subsets to only reference + target samples (not all samples)
# - Output directories created with correct naming

# Multi-level design expands to multiple pairwise contrasts (all vs reference).
# Each contrast subsets samples to just reference + target level.
testthat::test_that("contrasts expanded and samples subsetted", {
    fixture_dir <- tempfile("de-multi-")
    dir.create(fixture_dir)

    levels <- c("control", "treated", "treated2")
    bundle <- write_de_fixture_bundle(fixture_dir, levels = levels)
    calls <- new.env(parent = emptyenv())
    calls$deseq <- list()

    fake_deseq <- function(
        count_mat,
        coldata,
        target_level,
        reference_level,
        condition_column,
        covariates,
        out_dir,
        contrast_name
    ) {
        calls$deseq[[target_level]] <- coldata$alias
        result <- data.frame(
            baseMean = seq_len(nrow(count_mat)),
            log2FoldChange = rep(1, nrow(count_mat)),
            lfcSE = rep(0.1, nrow(count_mat)),
            stat = rep(1, nrow(count_mat)),
            pvalue = rep(0.05, nrow(count_mat)),
            padj = rep(0.05, nrow(count_mat)),
            row.names = rownames(count_mat)
        )
        list(dds = structure(list(), class = "fake_dds"), result = result)
    }

    fake_dexseq <- function(tx_counts, tx_meta, coldata, condition_column, covariates) {
        dxr <- data.frame(
            featureID = tx_meta$TXNAME,
            groupID = tx_meta$GENEID,
            log2fold = rep(0.5, nrow(tx_meta)),
            pvalue = rep(0.5, nrow(tx_meta)),
            padj = rep(0.5, nrow(tx_meta)),
            exonBaseMean = rep(10, nrow(tx_meta)),
            row.names = tx_meta$TXNAME
        )
        list(dxd = structure(list(), class = "fake_dxd"), dxr = dxr)
    }

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )

    main_run_de_analysis(
        argv,
        deseq_runner = fake_deseq,
        dexseq_runner = fake_dexseq,
        pdf_fn = function(path) file.create(path),
        dev_off_fn = function() NULL,
        plot_ma_fn = function(...) NULL,
        plot_disp_fn = function(...) NULL,
        per_gene_q_fn = function(dxr) stats::setNames(c(0.2, 0.3), c("gene1", "gene2")),
        placeholder_pdf_fn = function(path, label) file.create(path)
    )

    treated_dir <- file.path(argv$out_dir, "condition_treated_vs_control")
    treated2_dir <- file.path(argv$out_dir, "condition_treated2_vs_control")

    testthat::expect_true(dir.exists(treated_dir))
    testthat::expect_true(dir.exists(treated2_dir))
    testthat::expect_equal(
        sort(calls$deseq$treated),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated_rep1", "treated_rep2", "treated_rep3"))
    )
    testthat::expect_equal(
        sort(calls$deseq$treated2),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated2_rep1", "treated2_rep2", "treated2_rep3"))
    )
})

# Transcript IDs may contain '|' (e.g., Ensembl IDs like ENST00000123.4|ENSG00000456.7).
# TSV reading defaults to using '|' as separator - verify workflow preserves these IDs.
testthat::test_that("pipe characters in transcript IDs preserved", {
    fixture_dir <- tempfile("de-pipe-ids-")
    dir.create(fixture_dir)
    bundle <- write_de_fixture_bundle(fixture_dir, levels = c("control", "treated"))

    tx_se <- readRDS(bundle$transcript_rds)
    pipe_ids <- c(
        "ENSMUST000001|ENSMUSG000001|chr1:1-50",
        "ENSMUST000002|ENSMUSG000001|chr1:101-150",
        "ENSMUST000003|ENSMUSG000002|chr1:201-250",
        "ENSMUST000004|ENSMUSG000002|chr1:301-350"
    )
    rownames(tx_se) <- pipe_ids
    S4Vectors::mcols(SummarizedExperiment::rowRanges(tx_se))$TXNAME <- pipe_ids
    saveRDS(tx_se, bundle$transcript_rds)

    captured <- new.env(parent = emptyenv())
    fake_deseq <- function(
        count_mat,
        coldata,
        target_level,
        reference_level,
        condition_column,
        covariates,
        out_dir,
        contrast_name
    ) {
        result <- data.frame(
            baseMean = seq_len(nrow(count_mat)),
            log2FoldChange = rep(1, nrow(count_mat)),
            lfcSE = rep(0.1, nrow(count_mat)),
            stat = rep(1, nrow(count_mat)),
            pvalue = rep(0.05, nrow(count_mat)),
            padj = rep(0.05, nrow(count_mat)),
            row.names = rownames(count_mat)
        )
        list(dds = structure(list(), class = "fake_dds"), result = result)
    }
    fake_dexseq <- function(tx_counts, tx_meta, coldata, condition_column, covariates) {
        captured$feature_ids <- tx_meta$TXNAME
        dxr <- data.frame(
            featureID = tx_meta$TXNAME,
            groupID = tx_meta$GENEID,
            log2fold = rep(0.5, nrow(tx_meta)),
            pvalue = rep(0.5, nrow(tx_meta)),
            padj = rep(0.5, nrow(tx_meta)),
            exonBaseMean = rep(10, nrow(tx_meta)),
            row.names = tx_meta$TXNAME
        )
        list(dxd = structure(list(), class = "fake_dxd"), dxr = dxr)
    }

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )

    main_run_de_analysis(
        argv,
        deseq_runner = fake_deseq,
        dexseq_runner = fake_dexseq,
        pdf_fn = function(path) file.create(path),
        dev_off_fn = function() NULL,
        plot_ma_fn = function(...) NULL,
        plot_disp_fn = function(...) NULL,
        per_gene_q_fn = function(dxr) stats::setNames(c(0.2, 0.3), c("gene1", "gene2")),
        placeholder_pdf_fn = function(path, label) file.create(path)
    )

    contrast_dir <- file.path(argv$out_dir, "condition_treated_vs_control")
    dtu_tx <- utils::read.delim(
        file.path(contrast_dir, "results_dtu_transcript.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    dexseq <- utils::read.delim(
        file.path(contrast_dir, "results_dexseq.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_equal(captured$feature_ids, pipe_ids)
    testthat::expect_equal(dtu_tx$featureID, pipe_ids)
    testthat::expect_equal(dexseq$featureID, pipe_ids)
})

###
# CLI integration tests
#
# End-to-end test with real DESeq2/DRIMSeq (not mocked):
# - Run `supeRglue de_analysis` on synthetic count data
# - Verify output directories and files exist for downstream report generation

testthat::test_that("CLI integration produces expected outputs", {
    fixture_dir <- tempfile("de-cli-")
    dir.create(fixture_dir)
    bundle <- write_de_fixture_bundle(fixture_dir, levels = c("control", "treated"))
    out_dir <- file.path(fixture_dir, "out")

    result <- run_rscript(
        "supeRglue",
        c(
            "de_analysis",
            "--transcript_rds", bundle$transcript_rds,
            "--gene_rds", bundle$gene_rds,
            "--sample_sheet", bundle$sample_sheet,
            "--covariates", "batch",
            "--reference_level", "control",
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        result$status,
        0L,
        info = paste(result$output, collapse = "\n")
    )

    contrast_dir <- file.path(out_dir, "condition_treated_vs_control")
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dge.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_transcript.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_gene.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dexseq.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "samples_used.tsv")))

    dge <- utils::read.delim(file.path(contrast_dir, "results_dge.tsv"), check.names = FALSE)
    dtu_tx <- utils::read.delim(file.path(contrast_dir, "results_dtu_transcript.tsv"), check.names = FALSE)
    dexseq <- utils::read.delim(file.path(contrast_dir, "results_dexseq.tsv"), check.names = FALSE)

    testthat::expect_gt(nrow(dge), 0)
    testthat::expect_true(all(c("GENEID", "log2FoldChange", "padj") %in% names(dge)))
    testthat::expect_true(all(c("featureID", "groupID", "padj") %in% names(dtu_tx)))
    testthat::expect_true(all(c("featureID", "groupID", "padj") %in% names(dexseq)))
})
