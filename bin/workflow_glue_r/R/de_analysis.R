de_analysis_arg_parser <- function() {
    parser <- argparser::arg_parser("Run DESeq2 and DEXSeq on bambu output.")
    parser <- argparser::add_argument(parser, "--transcript_rds", help = "bambu transcript RDS.")
    parser <- argparser::add_argument(parser, "--gene_rds", help = "bambu gene RDS.")
    parser <- argparser::add_argument(parser, "--sample_sheet", help = "Sample sheet CSV.")
    parser <- argparser::add_argument(
        parser,
        "--condition_column",
        help = "Primary condition column.",
        default = "condition"
    )
    parser <- argparser::add_argument(parser, "--covariates", help = "Comma-separated nuisance covariates.")
    parser <- argparser::add_argument(parser, "--reference_level", help = "Reference level for the condition column.")
    argparser::add_argument(parser, "--out_dir", help = "Output directory.", default = "de_analysis")
}

de_parse_covariates <- function(value) {
    workflow_glue_r_parse_csv_list(value)
}

de_validate_inputs <- function(tx_se, gene_se, sample_df, argv) {
    workflow_glue_r_require_args(argv, c("transcript_rds", "gene_rds", "sample_sheet"))

    covariates <- de_parse_covariates(argv$covariates)
    workflow_glue_r_validate_r_formula_names(
        c(argv$condition_column, covariates),
        label = "Design column"
    )

    if (!"alias" %in% names(sample_df)) {
        stop("Sample sheet must contain an 'alias' column.", call. = FALSE)
    }
    duplicate_sample_aliases <- unique(sample_df$alias[duplicated(sample_df$alias)])
    if (length(duplicate_sample_aliases) > 0) {
        stop(
            sprintf(
                "Sample sheet aliases must be unique; duplicated aliases: %s",
                paste(duplicate_sample_aliases, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    if (!(argv$condition_column %in% names(sample_df))) {
        stop(
            sprintf("Sample sheet must contain the '%s' column.", argv$condition_column),
            call. = FALSE
        )
    }

    missing_covariates <- setdiff(covariates, names(sample_df))
    if (length(missing_covariates) > 0) {
        stop(
            sprintf(
                "Missing covariate columns: %s",
                paste(missing_covariates, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (any(duplicated(colnames(tx_se)))) {
        stop("Transcript RDS sample names must be unique.", call. = FALSE)
    }
    if (any(duplicated(colnames(gene_se)))) {
        stop("Gene RDS sample names must be unique.", call. = FALSE)
    }
    if (!setequal(colnames(tx_se), colnames(gene_se))) {
        stop("Transcript and gene RDS sample names must match.", call. = FALSE)
    }

    tx_counts <- SummarizedExperiment::assays(tx_se)$counts
    gene_counts <- SummarizedExperiment::assays(gene_se)$counts
    if (is.null(tx_counts)) {
        stop("Transcript RDS must contain a 'counts' assay.", call. = FALSE)
    }
    if (is.null(gene_counts)) {
        stop("Gene RDS must contain a 'counts' assay.", call. = FALSE)
    }
    if (!is.numeric(tx_counts) || !is.numeric(gene_counts)) {
        stop("Count matrices must be numeric.", call. = FALSE)
    }
    if (anyNA(tx_counts) || anyNA(gene_counts)) {
        stop("Count matrices must not contain NA values.", call. = FALSE)
    }
    zero_count_samples <- unique(c(
        colnames(tx_counts)[colSums(tx_counts) == 0],
        colnames(gene_counts)[colSums(gene_counts) == 0]
    ))
    if (length(zero_count_samples) > 0) {
        stop(
            sprintf(
                "Count matrices contain samples with zero total counts: %s",
                paste(zero_count_samples, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    sample_df <- sample_df[match(colnames(tx_se), sample_df$alias), , drop = FALSE]
    if (any(is.na(sample_df$alias))) {
        stop("Sample sheet aliases do not match the bambu output sample names.", call. = FALSE)
    }

    condition_values <- unique(sample_df[[argv$condition_column]])
    if (length(condition_values) < 2) {
        stop("Differential analysis requires at least two condition levels.", call. = FALSE)
    }

    reference_level <- argv$reference_level
    if (workflow_glue_r_arg_missing(reference_level)) {
        if ("control" %in% condition_values) {
            reference_level <- "control"
        } else {
            stop(
                "Provide --reference_level when the condition column does not contain 'control'.",
                call. = FALSE
            )
        }
    }
    if (!(reference_level %in% condition_values)) {
        stop("The requested reference level is not present in the condition column.", call. = FALSE)
    }

    sample_df[[argv$condition_column]] <- factor(sample_df[[argv$condition_column]])
    for (covariate in covariates) {
        sample_df[[covariate]] <- factor(sample_df[[covariate]])
    }

    list(
        tx_se = tx_se,
        gene_se = gene_se,
        sample_df = sample_df,
        covariates = covariates,
        condition_values = condition_values,
        reference_level = reference_level
    )
}

de_build_contrast_name <- function(condition_column, target_level, reference_level) {
    sprintf("%s_%s_vs_%s", condition_column, target_level, reference_level)
}

de_set_dispersions <- function(object, value) {
    setter <- get("dispersions<-", envir = asNamespace("DESeq2"))
    setter(object, value = value)
}

de_extract_disp_gene_est <- function(object) {
    S4Vectors::mcols(object)$dispGeneEst
}

de_run_deseq_with_fallback <- function(
    dds,
    contrast_name,
    out_dir
) {
    fallback_info <- list(
        applied = FALSE,
        method_used = "parametric",
        reason = NULL,
        diagnostic_file = NULL
    )

    de_out <- tryCatch(
        DESeq2::DESeq(dds, quiet = TRUE),
        error = function(err) {
            if (!grepl(
                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                conditionMessage(err),
                fixed = TRUE
            )) {
                stop(err)
            }

            warning(
                "STATISTICAL POWER REDUCED: DESeq2 dispersion estimation failed for ",
                contrast_name,
                ".\n",
                "This usually indicates:\n",
                "  1. Too few replicates (recommend n>=3 per group)\n",
                "  2. High biological variability\n",
                "  3. Poor data quality\n",
                "Falling back to gene-wise dispersion (no information sharing).\n",
                "Results will have reduced power and wider confidence intervals."
            )

            dds <- DESeq2::estimateSizeFactors(dds)
            dds <- DESeq2::estimateDispersionsGeneEst(dds)
            dds <- de_set_dispersions(dds, de_extract_disp_gene_est(dds))

            dispersion_values <- suppressWarnings(as.numeric(DESeq2::dispersions(dds)))
            dispersion_values <- dispersion_values[is.finite(dispersion_values)]
            dispersion_range <- if (length(dispersion_values) > 0) {
                sprintf(
                    "Dispersion range: %.3f to %.3f",
                    min(dispersion_values),
                    max(dispersion_values)
                )
            } else {
                "Dispersion range: unavailable"
            }

            diag_content <- c(
                "DESeq2 Dispersion Estimation Fallback Applied",
                "==============================================",
                "",
                sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Contrast: %s", contrast_name),
                sprintf("Samples: %d", ncol(dds)),
                sprintf("Genes tested: %d", nrow(dds)),
                dispersion_range,
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

            diag_file <- file.path(
                out_dir,
                sprintf(
                    "DESeq2_dispersion_fallback_%s.txt",
                    gsub("[^A-Za-z0-9_-]", "_", contrast_name)
                )
            )
            writeLines(diag_content, diag_file)
            fallback_info <<- list(
                applied = TRUE,
                method_used = "gene-wise",
                reason = conditionMessage(err),
                diagnostic_file = basename(diag_file)
            )

            DESeq2::nbinomWaldTest(dds)
        }
    )
    list(dds = de_out, deseq2_dispersion_fallback = fallback_info)
}

de_estimate_dispersions_with_fallback <- function(
    object,
    context_label,
    allow_gene_est = TRUE
) {
    tryCatch(
        list(
            object = DESeq2::estimateDispersions(object),
            method_used = "parametric",
            fallback_applied = FALSE,
            reason = NULL
        ),
        error = function(err) {
            if (!grepl(
                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                conditionMessage(err),
                fixed = TRUE
            )) {
                stop(err)
            }

            primary_reason <- conditionMessage(err)
            message(
                context_label,
                " dispersion fitting failed; retrying with fitType='local'."
            )
            tryCatch(
                list(
                    object = DESeq2::estimateDispersions(object, fitType = "local"),
                    method_used = "local",
                    fallback_applied = TRUE,
                    reason = primary_reason
                ),
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
                        " local-fit dispersion retry failed; retrying with fitType='mean'."
                    )
                    tryCatch(
                        list(
                            object = DESeq2::estimateDispersions(object, fitType = "mean"),
                            method_used = "mean",
                            fallback_applied = TRUE,
                            reason = primary_reason
                        ),
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
                                " mean-fit dispersion retry failed; falling back to gene-wise dispersion estimates."
                            )
                            object <- DESeq2::estimateDispersionsGeneEst(object)
                            object <- de_set_dispersions(object, de_extract_disp_gene_est(object))
                            list(
                                object = object,
                                method_used = "gene-wise",
                                fallback_applied = TRUE,
                                reason = primary_reason
                            )
                        }
                    )
                }
            )
        }
    )
}

de_is_recoverable_dexseq_error <- function(message_text) {
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

de_write_placeholder_pdf <- function(path, label) {
    grDevices::pdf(path)
    graphics::plot.new()
    graphics::text(0.5, 0.5, label, cex = 0.9)
    grDevices::dev.off()
}

de_run_deseq2_result <- function(
    count_mat,
    coldata,
    target_level,
    reference_level,
    condition_column,
    covariates,
    out_dir,
    contrast_name
) {
    design_terms <- c(covariates, condition_column)
    design_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(count_mat),
        colData = coldata,
        design = design_formula
    )
    deseq_run <- de_run_deseq_with_fallback(dds, contrast_name, out_dir)
    dds <- deseq_run$dds
    deseq2_dispersion_fallback <- deseq_run$deseq2_dispersion_fallback
    result <- DESeq2::results(
        dds,
        contrast = c(condition_column, target_level, reference_level),
        independentFiltering = TRUE
    )
    list(
        dds = dds,
        result = result,
        deseq2_dispersion_fallback = deseq2_dispersion_fallback
    )
}

de_run_dexseq_result <- function(
    tx_counts,
    tx_meta,
    coldata,
    condition_column,
    covariates
) {
    coldata$sample <- factor(coldata$alias)
    coldata[[condition_column]] <- factor(coldata[[condition_column]])
    for (covariate in covariates) {
        coldata[[covariate]] <- factor(coldata[[covariate]])
    }
    dropped_covariates <- character(0)

    run_inner <- function(active_covariates) {
        covariate_exon_terms <- if (length(active_covariates) > 0) {
            paste0(active_covariates, ":exon")
        } else {
            character(0)
        }
        design_terms <- c("sample", "exon", covariate_exon_terms, paste0(condition_column, ":exon"))
        reduced_terms <- c("sample", "exon", covariate_exon_terms)
        full_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))
        reduced_formula <- stats::as.formula(paste("~", paste(reduced_terms, collapse = " + ")))

        tryCatch({
            dxd <- DEXSeq::DEXSeqDataSet(
                countData = round(tx_counts),
                sampleData = as.data.frame(coldata),
                design = full_formula,
                featureID = tx_meta$TXNAME,
                groupID = tx_meta$GENEID
            )
            dxd <- DESeq2::estimateSizeFactors(dxd)
            dispersion_result <- de_estimate_dispersions_with_fallback(
                dxd,
                "DEXSeq",
                allow_gene_est = TRUE
            )
            if (is.list(dispersion_result) && !is.null(dispersion_result$object)) {
                dxd <- dispersion_result$object
                dispersion_method <- dispersion_result$method_used
                dispersion_reason <- dispersion_result$reason
            } else {
                dxd <- dispersion_result
                dispersion_method <- "parametric"
                dispersion_reason <- NULL
            }
            dxd <- DEXSeq::testForDEU(dxd, reducedModel = reduced_formula)
            dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = condition_column)
            dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
            list(
                dxd = dxd,
                dxr = dxr,
                dexseq_dispersion_method = dispersion_method,
                dexseq_dispersion_reason = dispersion_reason
            )
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
            dropped_covariates <<- c(dropped_covariates, dropped_covariate)
            message(
                "DEXSeq design was not full rank with covariate '",
                dropped_covariate,
                "'; retrying without it."
            )
            run_inner(kept_covariates)
        })
    }

    result <- run_inner(covariates)
    result$dexseq_covariates_dropped <- dropped_covariates
    result
}

main_run_de_analysis <- function(argv) {
    set.seed(42)
    dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

    tx_se <- readRDS(argv$transcript_rds)
    gene_se <- readRDS(argv$gene_rds)
    sample_df <- workflow_glue_r_read_csv(argv$sample_sheet)
    validated <- de_validate_inputs(tx_se, gene_se, sample_df, argv)
    sample_df <- validated$sample_df
    covariates <- validated$covariates
    condition_values <- validated$condition_values
    reference_level <- validated$reference_level

    tx_meta <- as.data.frame(SummarizedExperiment::rowData(tx_se))
    if (!"TXNAME" %in% names(tx_meta)) {
        tx_meta$TXNAME <- rownames(tx_se)
    }
    if (!"GENEID" %in% names(tx_meta)) {
        stop("Transcript rowData must contain GENEID for DEXSeq.", call. = FALSE)
    }

    gene_meta <- as.data.frame(SummarizedExperiment::rowData(gene_se))
    if (!"GENEID" %in% names(gene_meta)) {
        gene_meta$GENEID <- rownames(gene_se)
    }

    targets <- setdiff(as.character(condition_values), reference_level)

    de_qc_stats <- list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        total_samples = nrow(sample_df),
        condition_column = argv$condition_column,
        reference_level = reference_level,
        covariates = if (length(covariates) > 0) covariates else "none",
        num_contrasts = length(targets),
        contrasts = list()
    )

    n_per_group <- table(sample_df[[argv$condition_column]])
    de_qc_stats$samples_per_group <- as.list(n_per_group)

    sample_size_warnings <- character(0)
    if (any(n_per_group < 3)) {
        warning(
            "WARNING: Some condition groups have fewer than 3 replicates.\n",
            "Recommended minimum for DGE: n=3 per group\n",
            "Current sample sizes: ",
            paste(names(n_per_group), "=", n_per_group, collapse = ", "),
            "\nResults may have reduced statistical power."
        )
        sample_size_warnings <- c(sample_size_warnings, "Some groups have n<3 (recommended minimum)")
    }
    if (any(n_per_group < 2)) {
        stop(
            "ERROR: Some condition groups have fewer than 2 replicates. Cannot perform statistical testing.",
            call. = FALSE
        )
    }
    de_qc_stats$sample_size_warnings <- if (length(sample_size_warnings) > 0) {
        sample_size_warnings
    } else {
        "none"
    }

    if (length(targets) > 1) {
        fwer <- (1 - (1 - 0.05)^length(targets)) * 100
        mt_warning <- sprintf(
            "Multiple contrasts tested (%d). Per-contrast FDR < 0.05 yields family-wise error rate of ~%.1f%%",
            length(targets),
            fwer
        )
        message("WARNING: ", mt_warning)
        de_qc_stats$multiple_testing_note <- mt_warning
        mt_content <- c(
            "Multiple Testing Across Contrasts",
            "==================================",
            "",
            sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            sprintf("Number of contrasts tested: %d", length(targets)),
            sprintf("Contrasts: %s", paste(sprintf("%s vs %s", targets, reference_level), collapse = ", ")),
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
            sprintf("  1. Bonferroni correction: 0.05 / %d = %.4f", length(targets), 0.05 / length(targets)),
            "  2. Focus on pre-specified contrasts of interest",
            "  3. Treat results as exploratory and validate key findings",
            "  4. Consider using hierarchical testing procedures",
            ""
        )
        writeLines(mt_content, file.path(argv$out_dir, "MULTIPLE_TESTING_WARNING.txt"))
    }

    for (target_level in targets) {
        contrast_name <- de_build_contrast_name(argv$condition_column, target_level, reference_level)
        contrast_dir <- file.path(argv$out_dir, contrast_name)
        dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

        keep_samples <- sample_df[[argv$condition_column]] %in% c(reference_level, target_level)
        contrast_samples <- droplevels(sample_df[keep_samples, , drop = FALSE])
        contrast_samples[[argv$condition_column]] <- stats::relevel(
            factor(contrast_samples[[argv$condition_column]]),
            ref = reference_level
        )

        contrast_qc <- list(
            name = contrast_name,
            target_level = target_level,
            reference_level = reference_level,
            n_samples = nrow(contrast_samples),
            n_target = sum(contrast_samples[[argv$condition_column]] == target_level),
            n_reference = sum(contrast_samples[[argv$condition_column]] == reference_level),
            deseq2_dispersion_fallback = list(
                applied = FALSE,
                method_used = "parametric",
                reason = NULL,
                diagnostic_file = NULL
            ),
            dexseq_dispersion_method = "parametric",
            dexseq_covariates_dropped = list()
        )

        if (nrow(contrast_samples) < 6) {
            contrast_qc$dtu_power_warning <- sprintf(
                "DTU analysis may be underpowered (n=%d, recommend n>=6 with >=3 per group)",
                nrow(contrast_samples)
            )
            warning(contrast_qc$dtu_power_warning)
        }

        gene_counts <- SummarizedExperiment::assays(gene_se)$counts[, contrast_samples$alias, drop = FALSE]
        tx_counts <- SummarizedExperiment::assays(tx_se)$counts[, contrast_samples$alias, drop = FALSE]
        contrast_qc$genes_tested <- nrow(gene_counts)
        contrast_qc$transcripts_tested <- nrow(tx_counts)

        dge_run <- de_run_deseq2_result(
            gene_counts,
            contrast_samples,
            target_level,
            reference_level,
            argv$condition_column,
            covariates,
            argv$out_dir,
            contrast_name
        )
        if (!is.null(dge_run$deseq2_dispersion_fallback)) {
            fallback <- dge_run$deseq2_dispersion_fallback
            fallback_applied <- isTRUE(fallback$applied)
            fallback_method <- fallback$method_used
            if (is.null(fallback_method) || identical(fallback_method, "")) {
                fallback_method <- if (fallback_applied) "gene-wise" else "parametric"
            }
            contrast_qc$deseq2_dispersion_fallback <- list(
                applied = fallback_applied,
                method_used = fallback_method,
                reason = fallback$reason,
                diagnostic_file = fallback$diagnostic_file
            )
        }
        dge_res <- as.data.frame(dge_run$result)
        dge_res$GENEID <- rownames(dge_res)
        dge_res <- merge(gene_meta, dge_res, by = "GENEID", all.y = TRUE, sort = FALSE)
        dge_res <- workflow_glue_r_normalise_tsv_df(dge_res)

        contrast_qc$dge_total_genes <- nrow(dge_res)
        contrast_qc$dge_significant_fdr05 <- sum(dge_res$padj < 0.05, na.rm = TRUE)
        contrast_qc$dge_significant_fdr01 <- sum(dge_res$padj < 0.01, na.rm = TRUE)
        contrast_qc$dge_upregulated <- sum(
            dge_res$padj < 0.05 & dge_res$log2FoldChange > 0,
            na.rm = TRUE
        )
        contrast_qc$dge_downregulated <- sum(
            dge_res$padj < 0.05 & dge_res$log2FoldChange < 0,
            na.rm = TRUE
        )

        utils::write.table(
            dge_res[order(dge_res$padj), ],
            file = file.path(contrast_dir, "results_dge.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        grDevices::pdf(file.path(contrast_dir, "results_dge.pdf"))
        DESeq2::plotMA(dge_run$result)
        grDevices::dev.off()

        dex_res <- tryCatch(
            de_run_dexseq_result(
                tx_counts,
                tx_meta,
                contrast_samples,
                argv$condition_column,
                covariates
            ),
            error = function(err) {
                message_text <- conditionMessage(err)
                if (!de_is_recoverable_dexseq_error(message_text)) {
                    stop(err)
                }

                warning(
                    "DEXSeq failed for contrast ",
                    target_level,
                    " vs ",
                    reference_level,
                    "\nError: ",
                    message_text
                )

                failure_content <- c(
                    "DTU Analysis Failed",
                    "===================",
                    "",
                    sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                    sprintf("Contrast: %s vs %s", target_level, reference_level),
                    sprintf(
                        "Samples: %d (%d %s, %d %s)",
                        nrow(contrast_samples),
                        sum(contrast_samples[[argv$condition_column]] == target_level),
                        target_level,
                        sum(contrast_samples[[argv$condition_column]] == reference_level),
                        reference_level
                    ),
                    sprintf("Transcripts: %d", nrow(tx_counts)),
                    "",
                    "ERROR MESSAGE:",
                    sprintf("  %s", message_text),
                    "",
                    "DTU RESULTS CANNOT BE INTERPRETED",
                    ""
                )
                writeLines(failure_content, file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt"))
                NULL
            }
        )

        if (is.null(dex_res)) {
            dex_df <- workflow_glue_r_empty_tsv(c(
                "featureID",
                "groupID",
                "log2fold",
                "pvalue",
                "padj",
                "exonBaseMean"
            ))
            tx_dtu <- dex_df
            gene_dtu <- workflow_glue_r_empty_tsv(c("GENEID", "qval"))
            de_write_placeholder_pdf(
                file.path(contrast_dir, "results_dtu.pdf"),
                "DEXSeq did not converge for this contrast.\nSee DTU_ANALYSIS_FAILED.txt for details."
            )
            contrast_qc$dtu_status <- "FAILED"
            contrast_qc$dtu_significant_transcripts <- 0
            contrast_qc$dtu_significant_genes <- 0
        } else {
            if (!is.null(dex_res$dexseq_dispersion_method)) {
                contrast_qc$dexseq_dispersion_method <- dex_res$dexseq_dispersion_method
            }
            if (!is.null(dex_res$dexseq_covariates_dropped)) {
                contrast_qc$dexseq_covariates_dropped <- as.list(dex_res$dexseq_covariates_dropped)
            }
            dex_df <- as.data.frame(dex_res$dxr)
            dex_df <- workflow_glue_r_normalise_tsv_df(dex_df)
            tx_dtu <- dex_df[, intersect(
                c("featureID", "groupID", "log2fold", "pvalue", "padj", "exonBaseMean"),
                names(dex_df)
            ), drop = FALSE]
            tx_dtu <- workflow_glue_r_normalise_tsv_df(tx_dtu)

            gene_q <- DEXSeq::perGeneQValue(dex_res$dxr)
            gene_dtu <- data.frame(
                GENEID = names(gene_q),
                qval = unname(gene_q),
                row.names = NULL
            )

            contrast_qc$dtu_status <- "SUCCESS"
            contrast_qc$dtu_significant_transcripts <- sum(tx_dtu$padj < 0.05, na.rm = TRUE)
            contrast_qc$dtu_significant_genes <- sum(gene_dtu$qval < 0.05, na.rm = TRUE)

            grDevices::pdf(file.path(contrast_dir, "results_dtu.pdf"))
            DESeq2::plotMA(dex_res$dxr, cex = 0.8, alpha = 0.05)
            DESeq2::plotDispEsts(dex_res$dxd)
            grDevices::dev.off()
        }

        utils::write.table(
            dex_df,
            file = file.path(contrast_dir, "results_dexseq.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        utils::write.table(
            tx_dtu[order(tx_dtu$padj), ],
            file = file.path(contrast_dir, "results_dtu_transcript.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        utils::write.table(
            gene_dtu[order(gene_dtu$qval), ],
            file = file.path(contrast_dir, "results_dtu_gene.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        utils::write.table(
            contrast_samples,
            file = file.path(contrast_dir, "samples_used.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

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
                    sprintf(
                        "  Significant transcripts (FDR < 0.05): %d",
                        contrast_qc$dtu_significant_transcripts
                    ),
                    sprintf("  Genes with DTU (q < 0.05): %d", contrast_qc$dtu_significant_genes)
                )
            } else {
                "  See DTU_ANALYSIS_FAILED.txt for details"
            },
            if (!is.null(contrast_qc$dtu_power_warning)) paste0("  WARNING: ", contrast_qc$dtu_power_warning) else NULL,
            ""
        )
        writeLines(contrast_qc_summary, file.path(contrast_dir, "contrast_qc_summary.txt"))
        de_qc_stats$contrasts[[contrast_name]] <- contrast_qc
    }

    deseq2_dispersion_fallbacks <- names(Filter(
        function(cqc) isTRUE(cqc$deseq2_dispersion_fallback$applied),
        de_qc_stats$contrasts
    ))
    deseq2_gene_wise <- names(Filter(
        function(cqc) identical(cqc$deseq2_dispersion_fallback$method_used, "gene-wise"),
        de_qc_stats$contrasts
    ))
    dexseq_non_parametric <- names(Filter(
        function(cqc) {
            method <- cqc$dexseq_dispersion_method
            !is.null(method) && !identical(method, "parametric")
        },
        de_qc_stats$contrasts
    ))
    dexseq_gene_wise <- names(Filter(
        function(cqc) identical(cqc$dexseq_dispersion_method, "gene-wise"),
        de_qc_stats$contrasts
    ))
    dexseq_covariate_drop <- names(Filter(
        function(cqc) length(cqc$dexseq_covariates_dropped) > 0,
        de_qc_stats$contrasts
    ))
    total_covariates_dropped <- sum(vapply(
        de_qc_stats$contrasts,
        function(cqc) length(cqc$dexseq_covariates_dropped),
        integer(1)
    ))
    de_qc_stats$analysis_fallbacks <- list(
        deseq2_dispersion_fallback_contrasts = length(deseq2_dispersion_fallbacks),
        deseq2_dispersion_fallback_contrast_names = as.list(deseq2_dispersion_fallbacks),
        deseq2_gene_wise_contrasts = length(deseq2_gene_wise),
        deseq2_gene_wise_contrast_names = as.list(deseq2_gene_wise),
        dexseq_non_parametric_dispersion_contrasts = length(dexseq_non_parametric),
        dexseq_non_parametric_dispersion_contrast_names = as.list(dexseq_non_parametric),
        dexseq_gene_wise_dispersion_contrasts = length(dexseq_gene_wise),
        dexseq_gene_wise_dispersion_contrast_names = as.list(dexseq_gene_wise),
        dexseq_covariate_drop_contrasts = length(dexseq_covariate_drop),
        dexseq_covariate_drop_contrast_names = as.list(dexseq_covariate_drop),
        total_covariates_dropped = total_covariates_dropped
    )

    jsonlite::write_json(
        de_qc_stats,
        file.path(argv$out_dir, "de_qc_stats.json"),
        pretty = TRUE,
        auto_unbox = TRUE
    )

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
    writeLines(capture.output(sessionInfo()), file.path(argv$out_dir, "session_info.txt"))

    invisible(list(qc = de_qc_stats))
}

run_de_analysis_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
    parsed <- argparser::parse_args(de_analysis_arg_parser(), argv = argv)
    main_run_de_analysis(parsed)
}
