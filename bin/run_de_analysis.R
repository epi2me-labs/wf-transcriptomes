#!/usr/bin/env Rscript

# Set seed for reproducibility
set.seed(42)

suppressPackageStartupMessages({
    library(argparser)
    library(DESeq2)
    library(DEXSeq)
    library(SummarizedExperiment)
    library(jsonlite)
})

parser <- arg_parser("Run DESeq2 and DEXSeq on bambu output.")
parser <- add_argument(parser, "--transcript_rds", help = "bambu transcript RDS.")
parser <- add_argument(parser, "--gene_rds", help = "bambu gene RDS.")
parser <- add_argument(parser, "--sample_sheet", help = "Sample sheet CSV.")
parser <- add_argument(parser, "--condition_column", help = "Primary condition column.", default = "condition")
parser <- add_argument(parser, "--covariates", help = "Comma-separated nuisance covariates.")
parser <- add_argument(parser, "--reference_level", help = "Reference level for the condition column.")
parser <- add_argument(parser, "--out_dir", help = "Output directory.", default = "de_analysis")
argv <- parse_args(parser)

arg_missing <- function(value) {
    if (is.null(value) || length(value) == 0 || all(is.na(value))) {
        return(TRUE)
    }
    if (is.character(value)) {
        return(all(!nzchar(value)))
    }
    FALSE
}

required_args <- c("transcript_rds", "gene_rds", "sample_sheet")
missing_args <- required_args[vapply(required_args, function(arg_name) {
    value <- argv[[arg_name]]
    arg_missing(value)
}, logical(1))]
if (length(missing_args) > 0) {
    stop(sprintf(
        "Missing required arguments: %s",
        paste(sprintf("--%s", missing_args), collapse = ", ")
    ))
}

dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

tx_se <- readRDS(argv$transcript_rds)
gene_se <- readRDS(argv$gene_rds)
sample_df <- read.csv(argv$sample_sheet, check.names = FALSE, stringsAsFactors = FALSE)

if (!"alias" %in% names(sample_df)) {
    stop("Sample sheet must contain an 'alias' column.")
}
if (!(argv$condition_column %in% names(sample_df))) {
    stop(sprintf("Sample sheet must contain the '%s' column.", argv$condition_column))
}

covariates <- character(0)
if (!arg_missing(argv$covariates)) {
    covariates <- trimws(strsplit(argv$covariates, ",", fixed = TRUE)[[1]])
    covariates <- covariates[nzchar(covariates)]
}
missing_covariates <- setdiff(covariates, names(sample_df))
if (length(missing_covariates) > 0) {
    stop(sprintf("Missing covariate columns: %s", paste(missing_covariates, collapse = ", ")))
}

sample_df <- sample_df[match(colnames(tx_se), sample_df$alias), , drop = FALSE]
if (any(is.na(sample_df$alias))) {
    stop("Sample sheet aliases do not match the bambu output sample names.")
}

condition_values <- unique(sample_df[[argv$condition_column]])
if (length(condition_values) < 2) {
    stop("Differential analysis requires at least two condition levels.")
}

reference_level <- argv$reference_level
if (arg_missing(reference_level)) {
    if ("control" %in% condition_values) {
        reference_level <- "control"
    } else {
        stop("Provide --reference_level when the condition column does not contain 'control'.")
    }
}
if (!(reference_level %in% condition_values)) {
    stop("The requested reference level is not present in the condition column.")
}

sample_df[[argv$condition_column]] <- factor(sample_df[[argv$condition_column]])
for (covariate in covariates) {
    sample_df[[covariate]] <- factor(sample_df[[covariate]])
}

run_deseq_with_fallback <- function(dds, contrast_name = "unknown") {
    tryCatch(
        DESeq(dds, quiet = TRUE),
        error = function(err) {
            if (!grepl(
                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                conditionMessage(err),
                fixed = TRUE
            )) {
                stop(err)
            }

            warning(
                "STATISTICAL POWER REDUCED: DESeq2 dispersion estimation failed for ", contrast_name, ".\n",
                "This usually indicates:\n",
                "  1. Too few replicates (recommend n>=3 per group)\n",
                "  2. High biological variability\n",
                "  3. Poor data quality\n",
                "Falling back to gene-wise dispersion (no information sharing).\n",
                "Results will have reduced power and wider confidence intervals."
            )

            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
            dispersions(dds) <- mcols(dds)$dispGeneEst

            # Write diagnostic file
            diag_content <- c(
                "DESeq2 Dispersion Estimation Fallback Applied",
                "==============================================",
                "",
                sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Contrast: %s", contrast_name),
                sprintf("Samples: %d", ncol(dds)),
                sprintf("Genes tested: %d", nrow(dds)),
                sprintf("Dispersion range: %.3f to %.3f", min(dispersions(dds)), max(dispersions(dds))),
                "",
                "WHAT HAPPENED:",
                "  Curve fitting failed. Using gene-wise dispersion estimates.",
                "",
                "IMPLICATIONS:",
                "  - No information sharing across genes",
                "  - Reduced statistical power",
                "  - Wider confidence intervals",
                "  - More conservative results (fewer discoveries)",
                "",
                "LIKELY CAUSES:",
                "  1. Too few replicates (recommend n>=3 per group)",
                "  2. High biological variability",
                "  3. Poor data quality or outlier samples",
                "",
                "RECOMMENDATIONS:",
                "  - Add more biological replicates if possible",
                "  - Check sample quality metrics",
                "  - Consider filtering low-count genes more stringently"
            )

            diag_file <- file.path(argv$out_dir, sprintf("DESeq2_dispersion_fallback_%s.txt", gsub("[^A-Za-z0-9_-]", "_", contrast_name)))
            writeLines(diag_content, diag_file)

            nbinomWaldTest(dds)
        }
    )
}

estimate_dispersions_with_fallback <- function(object, context_label, allow_gene_est = TRUE) {
    tryCatch(
        estimateDispersions(object),
        error = function(err) {
            if (!grepl(
                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                conditionMessage(err),
                fixed = TRUE
            )) {
                stop(err)
            }

            message(
                context_label,
                " dispersion fitting failed; ",
                "retrying with fitType='local'."
            )
            tryCatch(
                estimateDispersions(object, fitType = "local"),
                error = function(local_err) {
                    if (!grepl(
                        "all gene-wise dispersion estimates are within 2 orders of magnitude",
                        conditionMessage(local_err),
                        fixed = TRUE
                    )) {
                        stop(local_err)
                    }

                    message(
                        context_label,
                        " local-fit dispersion retry failed; ",
                        "retrying with fitType='mean'."
                    )
                    tryCatch(
                        estimateDispersions(object, fitType = "mean"),
                        error = function(mean_err) {
                            if (!grepl(
                                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                                conditionMessage(mean_err),
                                fixed = TRUE
                            )) {
                                stop(mean_err)
                            }
                            if (!allow_gene_est) {
                                stop(mean_err)
                            }

                            message(
                                context_label,
                                " mean-fit dispersion retry failed; ",
                                "falling back to gene-wise dispersion estimates."
                            )
                            object <- estimateDispersionsGeneEst(object)
                            dispersions(object) <- mcols(object)$dispGeneEst
                            object
                        }
                    )
                }
            )
        }
    )
}

normalise_tsv_value <- function(value) {
    if (length(value) == 0 || all(is.na(value))) {
        return(NA_character_)
    }
    if (is.list(value)) {
        value <- unlist(value, recursive = TRUE, use.names = FALSE)
    }
    if (length(value) == 0 || all(is.na(value))) {
        return(NA_character_)
    }
    paste(as.character(value), collapse = ";")
}

normalise_tsv_df <- function(df) {
    as.data.frame(
        lapply(df, function(column) {
            if (is.list(column)) {
                vapply(column, normalise_tsv_value, character(1))
            } else {
                column
            }
        }),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

is_recoverable_dexseq_error <- function(message_text) {
    grepl(
        "all gene-wise dispersion estimates are within 2 orders of magnitude",
        message_text,
        fixed = TRUE
    ) || grepl(
        "model matrix is not full rank",
        message_text,
        fixed = TRUE
    ) || grepl(
        "replacement has 1 row, data has 0",
        message_text,
        fixed = TRUE
    )
}

empty_tsv <- function(columns) {
    out <- as.data.frame(matrix(nrow = 0, ncol = length(columns)))
    names(out) <- columns
    out
}

write_placeholder_pdf <- function(path, label) {
    pdf(path)
    plot.new()
    text(0.5, 0.5, label, cex = 0.9)
    dev.off()
}

run_deseq2 <- function(count_mat, coldata, target_level, contrast_name) {
    design_terms <- c(covariates, argv$condition_column)
    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    dds <- DESeqDataSetFromMatrix(
        countData = round(count_mat),
        colData = coldata,
        design = design_formula
    )
    dds <- run_deseq_with_fallback(dds, contrast_name)
    results(dds, contrast = c(argv$condition_column, target_level, reference_level), independentFiltering = TRUE)
}

run_dexseq <- function(tx_counts, tx_meta, coldata, active_covariates = covariates) {
    coldata$sample <- factor(coldata$alias)
    coldata[[argv$condition_column]] <- factor(coldata[[argv$condition_column]])
    for (covariate in active_covariates) {
        coldata[[covariate]] <- factor(coldata[[covariate]])
    }

    covariate_exon_terms <- if (length(active_covariates) > 0) {
        paste0(active_covariates, ":exon")
    } else {
        character(0)
    }
    design_terms <- c("sample", "exon", covariate_exon_terms, paste0(argv$condition_column, ":exon"))
    reduced_terms <- c("sample", "exon", covariate_exon_terms)
    full_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))

    tryCatch({
        dxd <- DEXSeqDataSet(
            countData = round(tx_counts),
            sampleData = as.data.frame(coldata),
            design = full_formula,
            featureID = tx_meta$TXNAME,
            groupID = tx_meta$GENEID
        )
        dxd <- estimateSizeFactors(dxd)
        dxd <- estimate_dispersions_with_fallback(dxd, "DEXSeq", allow_gene_est = TRUE)
        dxd <- testForDEU(dxd, reducedModel = reduced_formula)
        dxd <- estimateExonFoldChanges(dxd, fitExpToVar = argv$condition_column)
        dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
        list(dxd = dxd, dxr = dxr)
    }, error = function(err) {
        if (length(active_covariates) == 0 || !grepl(
            "model matrix is not full rank",
            conditionMessage(err),
            fixed = TRUE
        )) {
            stop(err)
        }

        dropped_covariate <- tail(active_covariates, 1)
        kept_covariates <- head(active_covariates, -1)
        message(
            "DEXSeq design was not full rank with covariate '",
            dropped_covariate,
            "'; retrying without it."
        )
        run_dexseq(tx_counts, tx_meta, coldata, kept_covariates)
    })
}

tx_meta <- as.data.frame(rowData(tx_se))
if (!"TXNAME" %in% names(tx_meta)) {
    tx_meta$TXNAME <- rownames(tx_se)
}
if (!"GENEID" %in% names(tx_meta)) {
    stop("Transcript rowData must contain GENEID for DEXSeq.")
}

gene_meta <- as.data.frame(rowData(gene_se))
if (!"GENEID" %in% names(gene_meta)) {
    gene_meta$GENEID <- rownames(gene_se)
}

targets <- setdiff(as.character(condition_values), reference_level)

# Initialize QC statistics collector
de_qc_stats <- list()
de_qc_stats$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
de_qc_stats$total_samples <- nrow(sample_df)
de_qc_stats$condition_column <- argv$condition_column
de_qc_stats$reference_level <- reference_level
de_qc_stats$covariates <- if (length(covariates) > 0) covariates else "none"
de_qc_stats$num_contrasts <- length(targets)
de_qc_stats$contrasts <- list()

# Check sample sizes and warn if underpowered
n_per_group <- table(sample_df[[argv$condition_column]])
de_qc_stats$samples_per_group <- as.list(n_per_group)

sample_size_warnings <- c()
if (any(n_per_group < 3)) {
    warning(
        "WARNING: Some condition groups have fewer than 3 replicates.\n",
        "Recommended minimum for DGE: n=3 per group\n",
        "Current sample sizes: ", paste(names(n_per_group), "=", n_per_group, collapse=", "), "\n",
        "Results may have reduced statistical power."
    )
    sample_size_warnings <- c(sample_size_warnings, "Some groups have n<3 (recommended minimum)")
}

if (any(n_per_group < 2)) {
    stop("ERROR: Some condition groups have fewer than 2 replicates. Cannot perform statistical testing.")
}

de_qc_stats$sample_size_warnings <- if (length(sample_size_warnings) > 0) sample_size_warnings else "none"

# Multiple testing warning
if (length(targets) > 1) {
    fwer <- (1 - (1-0.05)^length(targets)) * 100
    mt_warning <- sprintf(
        "Multiple contrasts tested (%d). Per-contrast FDR < 0.05 yields family-wise error rate of ~%.1f%%",
        length(targets), fwer
    )
    message("WARNING: ", mt_warning)
    de_qc_stats$multiple_testing_note <- mt_warning

    mt_content <- c(
        "Multiple Testing Across Contrasts",
        "==================================",
        "",
        sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Number of contrasts tested: %d", length(targets)),
        sprintf("Contrasts: %s", paste(sprintf("%s vs %s", targets, reference_level), collapse=", ")),
        "",
        "PER-CONTRAST FDR THRESHOLD: 0.05",
        sprintf("FAMILY-WISE ERROR RATE: ~%.1f%%", fwer),
        "",
        "WHAT THIS MEANS:",
        "  Each contrast uses FDR < 0.05 independently.",
        "  When testing multiple contrasts, the overall false positive rate increases.",
        sprintf("  Expected: %.1f%% chance of at least one false positive across all contrasts", fwer),
        "",
        "RECOMMENDATIONS:",
        "  1. Use stricter per-contrast threshold:",
        sprintf("     Bonferroni correction: 0.05 / %d = %.4f", length(targets), 0.05/length(targets)),
        "  2. Focus on pre-specified contrasts of interest",
        "  3. Treat results as exploratory and validate key findings",
        "  4. Consider using hierarchical testing procedures",
        "",
        "INTERPRETATION:",
        "  - Results passing FDR < 0.05 in each contrast are discoveries for that contrast",
        "  - But the overall false discovery burden is higher than 5%",
        "  - Prioritize genes significant across multiple contrasts",
        "  - Validate top findings experimentally"
    )
    writeLines(mt_content, file.path(argv$out_dir, "MULTIPLE_TESTING_WARNING.txt"))
}

for (target_level in targets) {
    contrast_name <- sprintf("%s_%s_vs_%s", argv$condition_column, target_level, reference_level)
    contrast_dir <- file.path(argv$out_dir, contrast_name)
    dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

    keep_samples <- sample_df[[argv$condition_column]] %in% c(reference_level, target_level)
    contrast_samples <- droplevels(sample_df[keep_samples, , drop = FALSE])
    contrast_samples[[argv$condition_column]] <- relevel(
        factor(contrast_samples[[argv$condition_column]]),
        ref = reference_level
    )

    # Collect per-contrast QC stats
    contrast_qc <- list()
    contrast_qc$name <- contrast_name
    contrast_qc$target_level <- target_level
    contrast_qc$reference_level <- reference_level
    contrast_qc$n_samples <- nrow(contrast_samples)
    contrast_qc$n_target <- sum(contrast_samples[[argv$condition_column]] == target_level)
    contrast_qc$n_reference <- sum(contrast_samples[[argv$condition_column]] == reference_level)

    # DTU power warning
    if (nrow(contrast_samples) < 6) {
        dtu_warning <- sprintf(
            "DTU analysis may be underpowered (n=%d, recommend n>=6 with >=3 per group)",
            nrow(contrast_samples)
        )
        warning(dtu_warning)
        contrast_qc$dtu_power_warning <- dtu_warning
    }

    gene_counts <- assays(gene_se)$counts[, contrast_samples$alias, drop = FALSE]
    tx_counts <- assays(tx_se)$counts[, contrast_samples$alias, drop = FALSE]

    contrast_qc$genes_tested <- nrow(gene_counts)
    contrast_qc$transcripts_tested <- nrow(tx_counts)

    dge_res <- as.data.frame(run_deseq2(gene_counts, contrast_samples, target_level, contrast_name))
    dge_res$GENEID <- rownames(dge_res)
    dge_res <- merge(gene_meta, dge_res, by = "GENEID", all.y = TRUE, sort = FALSE)
    dge_res <- normalise_tsv_df(dge_res)

    # Collect DGE statistics
    contrast_qc$dge_total_genes <- nrow(dge_res)
    contrast_qc$dge_significant_fdr05 <- sum(dge_res$padj < 0.05, na.rm = TRUE)
    contrast_qc$dge_significant_fdr01 <- sum(dge_res$padj < 0.01, na.rm = TRUE)
    contrast_qc$dge_upregulated <- sum(dge_res$padj < 0.05 & dge_res$log2FoldChange > 0, na.rm = TRUE)
    contrast_qc$dge_downregulated <- sum(dge_res$padj < 0.05 & dge_res$log2FoldChange < 0, na.rm = TRUE)

    write.table(
        dge_res[order(dge_res$padj), ],
        file = file.path(contrast_dir, "results_dge.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    pdf(file.path(contrast_dir, "results_dge.pdf"))
    dds_plot <- DESeqDataSetFromMatrix(
        countData = round(gene_counts),
        colData = contrast_samples,
        design = as.formula(paste("~", paste(c(covariates, argv$condition_column), collapse = " + ")))
    )
    dds_plot <- run_deseq_with_fallback(dds_plot, contrast_name)
    plotMA(results(dds_plot, contrast = c(argv$condition_column, target_level, reference_level), independentFiltering = TRUE))
    dev.off()

    dex_res <- tryCatch(
        run_dexseq(tx_counts, tx_meta, contrast_samples),
        error = function(err) {
            message_text <- conditionMessage(err)
            if (!is_recoverable_dexseq_error(message_text)) {
                stop(err)
            }

            warning(
                "DEXSeq failed for contrast ", target_level, " vs ", reference_level, "\n",
                "Error: ", message_text
            )

            # Write explicit failure report
            failure_content <- c(
                "DTU Analysis Failed",
                "===================",
                "",
                sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Contrast: %s vs %s", target_level, reference_level),
                sprintf("Samples: %d (%d %s, %d %s)",
                        nrow(contrast_samples),
                        sum(contrast_samples[[argv$condition_column]] == target_level), target_level,
                        sum(contrast_samples[[argv$condition_column]] == reference_level), reference_level),
                sprintf("Transcripts: %d", nrow(tx_counts)),
                "",
                "ERROR MESSAGE:",
                sprintf("  %s", message_text),
                "",
                "DTU RESULTS CANNOT BE INTERPRETED",
                "",
                "This failure is likely due to:",
                "  1. Insufficient samples (need >=3 per group, recommend >=6 total for DTU)",
                "  2. Too few transcripts with sufficient counts",
                "  3. Design matrix not full rank (covariate confounding)",
                "  4. Extreme count distributions",
                "",
                "RECOMMENDATIONS:",
                "  - Use gene-level DGE results (less power required)",
                "  - Add more biological replicates",
                "  - Filter transcripts more stringently",
                "  - Simplify experimental design (remove problematic covariates)",
                "",
                "NOTE: Empty DTU result files indicate analysis failure, not 'no DTU detected'"
            )
            writeLines(failure_content, file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt"))

            NULL
        }
    )

    if (is.null(dex_res)) {
        dex_df <- empty_tsv(c(
            "featureID",
            "groupID",
            "log2fold",
            "pvalue",
            "padj",
            "exonBaseMean"
        ))
        tx_dtu <- dex_df
        gene_dtu <- empty_tsv(c("GENEID", "qval"))
        write_placeholder_pdf(
            file.path(contrast_dir, "results_dtu.pdf"),
            "DEXSeq did not converge for this contrast.\nSee DTU_ANALYSIS_FAILED.txt for details."
        )
        contrast_qc$dtu_status <- "FAILED"
        contrast_qc$dtu_significant_transcripts <- 0
        contrast_qc$dtu_significant_genes <- 0
    } else {
        dxr <- dex_res$dxr
        dxd <- dex_res$dxd
        dex_df <- as.data.frame(dxr)
        dex_df <- normalise_tsv_df(dex_df)
        tx_dtu <- dex_df[, intersect(
            c("featureID", "groupID", "log2fold", "pvalue", "padj", "exonBaseMean"),
            names(dex_df)
        ), drop = FALSE]
        tx_dtu <- normalise_tsv_df(tx_dtu)

        gene_q <- perGeneQValue(dxr)
        gene_dtu <- data.frame(
            GENEID = names(gene_q),
            qval = unname(gene_q),
            row.names = NULL
        )

        # Collect DTU statistics
        contrast_qc$dtu_status <- "SUCCESS"
        contrast_qc$dtu_significant_transcripts <- sum(tx_dtu$padj < 0.05, na.rm = TRUE)
        contrast_qc$dtu_significant_genes <- sum(gene_dtu$qval < 0.05, na.rm = TRUE)

        pdf(file.path(contrast_dir, "results_dtu.pdf"))
        plotMA(dxr, cex = 0.8, alpha = 0.05)
        plotDispEsts(dxd)
        dev.off()
    }

    write.table(
        dex_df,
        file = file.path(contrast_dir, "results_dexseq.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        tx_dtu[order(tx_dtu$padj), ],
        file = file.path(contrast_dir, "results_dtu_transcript.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        gene_dtu[order(gene_dtu$qval), ],
        file = file.path(contrast_dir, "results_dtu_gene.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        contrast_samples,
        file = file.path(contrast_dir, "samples_used.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Write per-contrast QC summary
    contrast_qc_summary <- c(
        sprintf("Contrast QC Summary: %s", contrast_name),
        paste(rep("=", 50), collapse = ""),
        "",
        "Sample Information:",
        sprintf("  Target level (%s): %d samples", target_level, contrast_qc$n_target),
        sprintf("  Reference level (%s): %d samples", reference_level, contrast_qc$n_reference),
        sprintf("  Total samples: %d", contrast_qc$n_samples),
        "",
        "DGE Results:",
        sprintf("  Genes tested: %d", contrast_qc$genes_tested),
        sprintf("  Significant (FDR < 0.05): %d", contrast_qc$dge_significant_fdr05),
        sprintf("  Significant (FDR < 0.01): %d", contrast_qc$dge_significant_fdr01),
        sprintf("  Upregulated: %d", contrast_qc$dge_upregulated),
        sprintf("  Downregulated: %d", contrast_qc$dge_downregulated),
        "",
        "DTU Results:",
        sprintf("  Status: %s", contrast_qc$dtu_status),
        sprintf("  Transcripts tested: %d", contrast_qc$transcripts_tested),
        if (contrast_qc$dtu_status == "SUCCESS") {
            c(
                sprintf("  Significant transcripts (FDR < 0.05): %d", contrast_qc$dtu_significant_transcripts),
                sprintf("  Genes with DTU (q < 0.05): %d", contrast_qc$dtu_significant_genes)
            )
        } else {
            "  See DTU_ANALYSIS_FAILED.txt for details"
        },
        if (!is.null(contrast_qc$dtu_power_warning)) paste0("  WARNING: ", contrast_qc$dtu_power_warning) else NULL,
        ""
    )
    writeLines(contrast_qc_summary, file.path(contrast_dir, "contrast_qc_summary.txt"))

    # Add to overall QC stats
    de_qc_stats$contrasts[[contrast_name]] <- contrast_qc
}

# Write overall DE/DTU QC statistics as JSON for HTML report
write_json(
    de_qc_stats,
    file.path(argv$out_dir, "de_qc_stats.json"),
    pretty = TRUE,
    auto_unbox = TRUE
)

# Write human-readable overall summary
overall_summary <- c(
    "Differential Expression/Usage Analysis Summary",
    paste(rep("=", 50), collapse = ""),
    "",
    sprintf("Timestamp: %s", de_qc_stats$timestamp),
    sprintf("Total samples: %d", de_qc_stats$total_samples),
    sprintf("Condition column: %s", de_qc_stats$condition_column),
    sprintf("Reference level: %s", de_qc_stats$reference_level),
    sprintf("Covariates: %s", paste(de_qc_stats$covariates, collapse = ", ")),
    "",
    "Sample Sizes:",
    sapply(names(de_qc_stats$samples_per_group), function(grp) {
        sprintf("  %s: %d samples", grp, de_qc_stats$samples_per_group[[grp]])
    }),
    if (de_qc_stats$sample_size_warnings != "none") paste0("  WARNING: ", de_qc_stats$sample_size_warnings) else NULL,
    "",
    sprintf("Number of contrasts tested: %d", de_qc_stats$num_contrasts),
    if (!is.null(de_qc_stats$multiple_testing_note)) paste0("  NOTE: ", de_qc_stats$multiple_testing_note) else NULL,
    "",
    "Per-Contrast Results:",
    sapply(names(de_qc_stats$contrasts), function(cname) {
        cqc <- de_qc_stats$contrasts[[cname]]
        c(
            "",
            sprintf("  %s:", cname),
            sprintf("    Samples: %d (%d vs %d)", cqc$n_samples, cqc$n_target, cqc$n_reference),
            sprintf("    DGE significant: %d genes (FDR<0.05)", cqc$dge_significant_fdr05),
            sprintf("    DTU status: %s", cqc$dtu_status),
            if (cqc$dtu_status == "SUCCESS") sprintf("    DTU significant: %d genes", cqc$dtu_significant_genes) else NULL
        )
    }),
    "",
    "For detailed per-contrast statistics, see:",
    "  - <contrast>/contrast_qc_summary.txt",
    "  - <contrast>/results_dge.tsv",
    "  - <contrast>/results_dtu_gene.tsv",
    ""
)
writeLines(unlist(overall_summary), file.path(argv$out_dir, "de_overall_summary.txt"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(argv$out_dir, "session_info.txt"))
message("QC statistics written to de_qc_stats.json and de_overall_summary.txt")
message("Session info saved for reproducibility")
