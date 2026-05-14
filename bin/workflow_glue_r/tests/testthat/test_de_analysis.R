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

testthat::test_that("DTU transcript output renames contrast-specific log2fold column", {
    contrast_name <- "treated_control"
    dex_df <- data.frame(
        featureID = "tx1",
        groupID = "gene1",
        log2fold_treated_control = 1.5,
        pvalue = 0.01,
        padj = 0.05,
        exonBaseMean = 100,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    dtu_tx <- de_extract_dtu_transcript_table(dex_df, contrast_name)

    testthat::expect_equal(
        names(dtu_tx),
        c(
            "featureID",
            "groupID",
            "log2FoldChange",
            "pvalue",
            "padj",
            "exonBaseMean"
        )
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
        main_run_de_analysis(argv),
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
        suppressWarnings(main_run_de_analysis(argv)),
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
        main_run_de_analysis(argv),
        "requires at least two condition levels"
    )

    malformed_rds <- file.path(fixture_dir, "not-an-rds.rds")
    writeLines("not an RDS file", malformed_rds)
    argv$transcript_rds <- malformed_rds
    argv$sample_sheet <- sample_sheet
    argv$out_dir <- file.path(fixture_dir, "malformed-out")
    testthat::expect_error(
        main_run_de_analysis(argv)
    )
})

###
# Fallback helpers and metadata wiring
#
# Fixture-driven tests for DESeq2/DEXSeq helper behaviour and metadata.
# These avoid dependency injection and exercise the real package code paths.

testthat::test_that("de_choose_size_factor_type uses poscounts for zero-heavy matrices", {
    count_mat <- matrix(
        c(
            0, 4, 5,
            3, 0, 2,
            7, 1, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    testthat::expect_warning(
        sf_type <- de_choose_size_factor_type(count_mat, "test matrix"),
        "using DESeq2 size-factor estimation with sfType='poscounts'"
    )
    testthat::expect_equal(sf_type, "poscounts")
})

testthat::test_that("de_run_deseq_with_fallback returns structured metadata", {
    testthat::skip_if_not_installed("DESeq2")

    gene_se <- make_test_gene_se()
    sample_df <- data.frame(
        alias = colnames(gene_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    rownames(sample_df) <- sample_df$alias

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(gene_se, "counts"),
        colData = sample_df,
        design = ~ batch + condition
    )
    out_dir <- tempfile("deseq-fallback-")
    dir.create(out_dir)

    result <- suppressWarnings(de_run_deseq_with_fallback(
        dds = dds,
        contrast_name = "condition_treated_vs_control",
        out_dir = out_dir
    ))

    testthat::expect_true(!is.null(result$dds))
    testthat::expect_true(result$deseq2_dispersion_fallback$method_used %in% c(
        "parametric",
        "gene-wise"
    ))
    testthat::expect_true(is.logical(result$deseq2_dispersion_fallback$applied))
    if (isTRUE(result$deseq2_dispersion_fallback$applied)) {
        testthat::expect_true(file.exists(file.path(
            out_dir,
            result$deseq2_dispersion_fallback$diagnostic_file
        )))
    }
})

testthat::test_that("de_run_deseq_with_fallback rethrows non-recoverable errors", {
    testthat::skip_if_not_installed("DESeq2")

    out_dir <- tempfile("deseq-fallback-error-")
    dir.create(out_dir)

    testthat::expect_error(
        de_run_deseq_with_fallback(
            dds = list(not = "a DESeqDataSet"),
            contrast_name = "condition_treated_vs_control",
            out_dir = out_dir
        )
    )
})

testthat::test_that("de_estimate_dispersions_with_fallback reports method metadata", {
    testthat::skip_if_not_installed("DESeq2")

    gene_se <- make_test_gene_se()
    sample_df <- data.frame(
        alias = colnames(gene_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    rownames(sample_df) <- sample_df$alias
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(gene_se, "counts"),
        colData = sample_df,
        design = ~ batch + condition
    )
    dds <- DESeq2::estimateSizeFactors(dds)

    result <- suppressWarnings(suppressMessages(
        de_estimate_dispersions_with_fallback(dds, "DESeq2")
    ))

    testthat::expect_true(result$method_used %in% c(
        "parametric",
        "local",
        "mean",
        "gene-wise"
    ))
    testthat::expect_true(is.logical(result$fallback_applied))
    testthat::expect_true(!is.null(result$object))
})

testthat::test_that("de_is_recoverable_dexseq_error recognises expected messages", {
    testthat::expect_true(de_is_recoverable_dexseq_error(
        "all gene-wise dispersion estimates are within 2 orders of magnitude"
    ))
    testthat::expect_true(de_is_recoverable_dexseq_error("model matrix is not full rank"))
    testthat::expect_true(de_is_recoverable_dexseq_error("replacement has 1 row, data has 0"))
    testthat::expect_false(de_is_recoverable_dexseq_error("random unrelated failure"))
})

testthat::test_that("de_run_dexseq_result records rank-deficiency covariate drops", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    tx_se <- make_test_tx_se(
        sample_names = c(
            "control_rep1", "control_rep2", "control_rep3",
            "treated_rep1", "treated_rep2", "treated_rep3"
        )
    )
    tx_counts <- SummarizedExperiment::assay(tx_se, "counts")
    tx_meta <- as.data.frame(SummarizedExperiment::rowData(tx_se))
    coldata <- data.frame(
        alias = colnames(tx_counts),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("control", "treated"), each = 3),
        stringsAsFactors = FALSE
    )

    seen_messages <- character(0)
    result <- withCallingHandlers(
        tryCatch(
            de_run_dexseq_result(
                tx_counts,
                tx_meta,
                coldata,
                condition_column = "condition",
                covariates = c("batch")
            ),
            error = function(err) err
        ),
        message = function(m) {
            seen_messages <<- c(seen_messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )

    testthat::expect_true(any(grepl(
        "retrying without it",
        seen_messages,
        fixed = TRUE
    )))

    if (inherits(result, "error")) {
        testthat::expect_match(
            conditionMessage(result),
            "model matrix is not full rank",
            fixed = TRUE
        )
    } else {
        testthat::expect_equal(result$dexseq_covariates_dropped, "batch")
        testthat::expect_true(result$dexseq_dispersion_method %in% c(
            "parametric",
            "local",
            "mean",
            "gene-wise"
        ))
        testthat::expect_true(nrow(as.data.frame(result$dxr)) > 0)
    }
})

###
# Contrast planning and output writing
#
# Use small synthetic fixtures with the real DESeq2/DEXSeq path to verify
# workflow-level planning and output writing:
# - One contrast created per non-reference condition level (treated vs control, treated2 vs control)
# - Each contrast subsets to only reference + target samples (not all samples)
# - Output directories created with correct naming

# Multi-level design expands to multiple pairwise contrasts (all vs reference).
# Each contrast subsets samples to just reference + target level.
testthat::test_that("contrasts expanded and samples subsetted", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    fixture_dir <- tempfile("de-multi-")
    dir.create(fixture_dir)

    levels <- c("control", "treated", "treated2")
    bundle <- write_de_fixture_bundle(fixture_dir, levels = levels)

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )

    suppressWarnings(suppressMessages(main_run_de_analysis(argv)))

    treated_dir <- file.path(argv$out_dir, "condition_treated_vs_control")
    treated2_dir <- file.path(argv$out_dir, "condition_treated2_vs_control")

    treated_samples <- utils::read.delim(
        file.path(treated_dir, "samples_used.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    treated2_samples <- utils::read.delim(
        file.path(treated2_dir, "samples_used.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_true(dir.exists(treated_dir))
    testthat::expect_true(dir.exists(treated2_dir))
    testthat::expect_equal(
        sort(treated_samples$alias),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated_rep1", "treated_rep2", "treated_rep3"))
    )
    testthat::expect_equal(
        sort(treated2_samples$alias),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated2_rep1", "treated2_rep2", "treated2_rep3"))
    )
})

# Transcript IDs may contain '|' (e.g., Ensembl IDs like ENST00000123.4|ENSG00000456.7).
# TSV reading defaults to using '|' as separator - verify workflow preserves these IDs.
testthat::test_that("pipe characters in transcript IDs preserved", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

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

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )

    suppressWarnings(suppressMessages(main_run_de_analysis(argv)))

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
    de_qc <- jsonlite::read_json(
        file.path(argv$out_dir, "de_qc_stats.json"),
        simplifyVector = TRUE
    )
    contrast_qc <- de_qc$contrasts[["condition_treated_vs_control"]]

    if (identical(contrast_qc$dtu_status, "SUCCESS")) {
        testthat::expect_equal(dtu_tx$featureID, pipe_ids)
        testthat::expect_equal(dexseq$featureID, pipe_ids)
    } else {
        testthat::expect_true(file.exists(file.path(
            contrast_dir,
            "DTU_ANALYSIS_FAILED.txt"
        )))
        testthat::expect_true(nrow(dtu_tx) == 0)
        testthat::expect_true(nrow(dexseq) == 0)
    }
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
    de_qc <- jsonlite::read_json(
        file.path(out_dir, "de_qc_stats.json"),
        simplifyVector = TRUE
    )

    testthat::expect_gt(nrow(dge), 0)
    testthat::expect_true(all(c("GENEID", "log2FoldChange", "padj") %in% names(dge)))
    testthat::expect_true(all(c("featureID", "groupID", "padj") %in% names(dtu_tx)))
    testthat::expect_true(all(c("featureID", "groupID", "padj") %in% names(dexseq)))
    testthat::expect_true("analysis_fallbacks" %in% names(de_qc))
    testthat::expect_true("contrasts" %in% names(de_qc))
    contrast_qc <- de_qc$contrasts[["condition_treated_vs_control"]]
    testthat::expect_true("deseq2_size_factor_method" %in% names(contrast_qc))
    testthat::expect_true("deseq2_dispersion_fallback" %in% names(contrast_qc))
    testthat::expect_true("dexseq_size_factor_method" %in% names(contrast_qc))
    testthat::expect_true("dexseq_dispersion_method" %in% names(contrast_qc))
    testthat::expect_true("dexseq_covariates_dropped" %in% names(contrast_qc))
})
