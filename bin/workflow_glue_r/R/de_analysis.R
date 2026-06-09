de_analysis_arg_spec <- function() {
    list(
        list(
            name = "transcript_rds",
            flag = "--transcript_rds",
            help = "bambu transcript RDS.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "gene_rds",
            flag = "--gene_rds",
            help = "bambu gene RDS.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "sample_sheet",
            flag = "--sample_sheet",
            help = "Sample sheet CSV.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "condition_column",
            flag = "--condition_column",
            help = "Primary condition column.",
            type = "character",
            default = "condition"
        ),
        list(
            name = "covariates",
            flag = "--covariates",
            help = "Comma-separated nuisance covariates.",
            type = "character"
        ),
        list(
            name = "reference_level",
            flag = "--reference_level",
            help = "Reference level for the condition column.",
            type = "character"
        ),
        list(
            name = "out_dir",
            flag = "--out_dir",
            help = "Output directory.",
            type = "character",
            default = "de_analysis"
        )
    )
}

de_analysis_arg_parser <- function() {
    workflow_glue_r_arg_parser_from_spec(
        "Run DESeq2 and DEXSeq on bambu output.",
        de_analysis_arg_spec()
    )
}

de_validate_inputs <- function(tx_se, gene_se, sample_df, argv) {
    covariates <- workflow_glue_r_parse_csv_list(argv$covariates)

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
    if (is.null(reference_level)) {
        reference_level <- "control"
    }

    sample_df[[argv$condition_column]] <- factor(sample_df[[argv$condition_column]])
    for (covariate in covariates) {
        sample_df[[covariate]] <- factor(sample_df[[covariate]])
    }

    list(
        sample_df = sample_df,
        covariates = covariates,
        condition_values = condition_values,
        reference_level = reference_level
    )
}

# DESeq2's default geometric-mean size-factor estimator is undefined when
# every gene has at least one zero across samples, so use poscounts then.
de_choose_size_factor_type <- function(count_mat, context_label = "Count matrix") {
    if (all(rowSums(count_mat == 0) > 0)) {
        warning(
            context_label,
            " has every gene containing at least one zero; using DESeq2 size-factor estimation with sfType='poscounts'."
        )
        return("poscounts")
    }
    "ratio"
}
de_run_deseq_with_fallback <- function(
    dds,
    contrast_name,
    out_dir
) {
    sf_type <- de_choose_size_factor_type(
        DESeq2::counts(dds),
        context_label = paste0("DGE count matrix for ", contrast_name)
    )
    fallback_info <- list(
        applied = FALSE,
        method_used = "parametric",
        reason = NULL,
        diagnostic_file = NULL,
        size_factor_type = sf_type
    )

    de_out <- tryCatch(
        list(
            dds = DESeq2::DESeq(dds, quiet = TRUE, sfType = sf_type),
            deseq2_dispersion_fallback = fallback_info
        ),
        error = function(err) {
            msg <- conditionMessage(err)
            zero_geometric_means <- grepl(
                "every gene contains at least one zero",
                msg,
                fixed = TRUE
            )

            # Size factor estimation failure: all genes have at least one zero.
            # de_choose_size_factor_type should prevent this proactively, but as
            # a safety net switch to poscounts so the gene-wise fallback below
            # can call estimateSizeFactors without failing again.
            fallback_sf_type <- if (zero_geometric_means) "poscounts" else sf_type

            recoverable <- grepl(
                "all gene-wise dispersion estimates are within 2 orders of magnitude",
                msg,
                fixed = TRUE
            ) || grepl(
                "newsplit: out of vertex space",
                msg,
                fixed = TRUE
            ) || grepl(
                "every gene contains at least one zero",
                msg,
                fixed = TRUE
            )
            if (!recoverable) {
                stop(err)
            }

            if (zero_geometric_means) {
                warning(
                    "STATISTICAL POWER REDUCED: DESeq2 default size-factor estimation failed for ",
                    contrast_name,
                    " because every gene contains at least one zero.\n",
                    "Switching to sfType='poscounts' and using gene-wise dispersion estimates.\n",
                    "Results will have reduced power and wider confidence intervals."
                )
            } else {
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
            }

            dds <- DESeq2::estimateSizeFactors(dds, type = fallback_sf_type)
            dds <- DESeq2::estimateDispersionsGeneEst(dds)
            dispersions_setter <- get("dispersions<-", envir = asNamespace("DESeq2"))
            dds <- dispersions_setter(dds, value = S4Vectors::mcols(dds)$dispGeneEst)

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
                if (zero_geometric_means) {
                    "DESeq2 Poscounts and Gene-wise Fallback Applied"
                } else {
                    "DESeq2 Dispersion Estimation Fallback Applied"
                },
                if (zero_geometric_means) {
                    "============================================="
                } else {
                    "=============================================="
                },
                "",
                sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                sprintf("Contrast: %s", contrast_name),
                sprintf("Samples: %d", ncol(dds)),
                sprintf("Genes tested: %d", nrow(dds)),
                sprintf("Size factor method: %s", fallback_sf_type),
                dispersion_range,
                "",
                "WHAT HAPPENED:",
                if (zero_geometric_means) {
                    "  Every gene contained at least one zero, so DESeq2's default size-factor estimator was undefined."
                } else {
                    "  Curve fitting failed. Using gene-wise dispersion estimates."
                },
                if (zero_geometric_means) {
                    "  The analysis switched to sfType='poscounts' and used gene-wise dispersion estimates."
                } else {
                    NULL
                },
                "",
                "IMPLICATIONS:",
                "  - No information sharing across genes",
                "  - Reduced statistical power",
                "  - Wider confidence intervals",
                "  - More conservative results (fewer discoveries)",
                "",
                "LIKELY CAUSES:",
                if (zero_geometric_means) {
                    c(
                        "  1. Sparse count matrices with many zeros",
                        "  2. Low-expression features dominating the design",
                        "  3. Small or weakly separated groups"
                    )
                } else {
                    c(
                        "  1. Too few replicates (recommend n>=3 per group)",
                        "  2. High biological variability",
                        "  3. Poor data quality or outlier samples"
                    )
                },
                "",
                "RECOMMENDATIONS:",
                if (zero_geometric_means) {
                    c(
                        "  - Check whether the contrast is dominated by zeros or nearly identical samples",
                        "  - Add more biological replicates if possible",
                        "  - Consider filtering very low-count genes more stringently"
                    )
                } else {
                    c(
                        "  - Add more biological replicates if possible",
                        "  - Check sample quality metrics",
                        "  - Consider filtering low-count genes more stringently"
                    )
                }
            )

            diag_file <- file.path(
                out_dir,
                sprintf(
                    "DESeq2_dispersion_fallback_%s.txt",
                    gsub("[^A-Za-z0-9_-]", "_", contrast_name)
                )
            )
            writeLines(diag_content, diag_file)
            list(
                dds = DESeq2::nbinomWaldTest(dds),
                deseq2_dispersion_fallback = list(
                    applied = TRUE,
                    method_used = "gene-wise",
                    reason = conditionMessage(err),
                    diagnostic_file = basename(diag_file),
                    size_factor_type = fallback_sf_type
                )
            )
        }
    )
    de_out
}

de_estimate_dispersions_with_fallback <- function(
    object,
    context_label,
    allow_gene_est = TRUE
) {
    de_dispersion_recoverable <- function(msg) {
        grepl(
            "all gene-wise dispersion estimates are within 2 orders of magnitude",
            msg, fixed = TRUE
        ) || grepl(
            "newsplit: out of vertex space",
            msg, fixed = TRUE
        )
    }
    tryCatch(
        list(
            object = DESeq2::estimateDispersions(object),
            method_used = "parametric",
            fallback_applied = FALSE
        ),
        error = function(err) {
            if (!de_dispersion_recoverable(conditionMessage(err))) {
                stop(err)
            }

            message(
                context_label,
                " dispersion fitting failed; retrying with fitType='local'."
            )
            tryCatch(
                list(
                    object = DESeq2::estimateDispersions(object, fitType = "local"),
                    method_used = "local",
                    fallback_applied = TRUE
                ),
                error = function(local_err) {
                    if (!de_dispersion_recoverable(conditionMessage(local_err))) {
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
                            fallback_applied = TRUE
                        ),
                        error = function(mean_err) {
                            if (!de_dispersion_recoverable(conditionMessage(mean_err))) {
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
                            dispersions_setter <- get("dispersions<-", envir = asNamespace("DESeq2"))
                            object <- dispersions_setter(
                                object,
                                value = S4Vectors::mcols(object)$dispGeneEst
                            )
                            list(
                                object = object,
                                method_used = "gene-wise",
                                fallback_applied = TRUE
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
    ) || grepl(
        "newsplit: out of vertex space",
        message_text,
        fixed = TRUE
    ) || grepl(
        "every gene contains at least one zero",
        message_text,
        fixed = TRUE
    )
}

de_dexseq_failure_hint <- function(message_text) {
    if (grepl("nb of cols in 'assay'", message_text, fixed = TRUE)) {
        return(
            paste(
                "Probable cause: degenerate or near-identical sample groups can",
                "drive DEXSeq into an internal object-shape mismatch during model",
                "fitting. This is often seen when the contrast has little real",
                "separation or too few informative non-zero features."
            )
        )
    }
    if (grepl("'x' must be an array of at least two dimensions", message_text, fixed = TRUE)) {
        return(
            paste(
                "Probable cause: degenerate or near-identical sample groups can",
                "leave DEXSeq with too little informative structure for downstream",
                "gene-level DTU summarisation. This is often seen when the contrast",
                "has little real separation or too few informative non-zero features."
            )
        )
    }
    NULL
}

de_write_placeholder_pdf <- function(path, label) {
    grDevices::pdf(path)
    graphics::plot.new()
    graphics::text(0.5, 0.5, label, cex = 0.9)
    grDevices::dev.off()
}

de_write_dtu_failure_outputs <- function(
    contrast_dir,
    target_level,
    reference_level,
    contrast_samples,
    condition_column,
    tx_counts,
    message_text,
    failure_hint
) {
    failure_content <- c(
        "DTU Analysis Failed",
        "===================",
        "",
        sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Contrast: %s vs %s", target_level, reference_level),
        sprintf(
            "Samples: %d (%d %s, %d %s)",
            nrow(contrast_samples),
            sum(contrast_samples[[condition_column]] == target_level),
            target_level,
            sum(contrast_samples[[condition_column]] == reference_level),
            reference_level
        ),
        sprintf("Transcripts: %d", nrow(tx_counts)),
        "",
        "ERROR MESSAGE:",
        sprintf("  %s", message_text),
        if (!is.null(failure_hint)) {
            c(
                "",
                "PROBABLE CAUSE:",
                sprintf("  %s", failure_hint)
            )
        } else {
            NULL
        },
        "",
        "DTU RESULTS CANNOT BE INTERPRETED",
        ""
    )
    writeLines(failure_content, file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt"))
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
        result = result,
        deseq2_dispersion_fallback = deseq2_dispersion_fallback
    )
}

de_estimate_size_factors_for_dexseq <- function(dxd, count_mat, coldata, sf_type) {
    sf_dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(count_mat),
        colData = coldata,
        design = ~ 1
    )
    sf_dds <- DESeq2::estimateSizeFactors(sf_dds, type = sf_type)
    DESeq2::sizeFactors(dxd) <- DESeq2::sizeFactors(sf_dds)
    dxd
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
    run_inner <- function(active_covariates, dropped_covariates = character(0)) {
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
            dexseq_sf_type <- de_choose_size_factor_type(
                round(tx_counts),
                context_label = "DEXSeq transcript count matrix"
            )
            dxd <- de_estimate_size_factors_for_dexseq(
                dxd,
                tx_counts,
                coldata,
                dexseq_sf_type
            )
            dispersion_result <- de_estimate_dispersions_with_fallback(
                dxd,
                "DEXSeq",
                allow_gene_est = TRUE
            )
            dxd <- dispersion_result$object
            dispersion_method <- dispersion_result$method_used
            dxd <- DEXSeq::testForDEU(dxd, reducedModel = reduced_formula)
            dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = condition_column)
            dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
            list(
                dxd = dxd,
                dxr = dxr,
                dexseq_dispersion_method = dispersion_method,
                dexseq_size_factor_type = dexseq_sf_type,
                dexseq_covariates_dropped = dropped_covariates
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
            message(
                "DEXSeq design was not full rank with covariate '",
                dropped_covariate,
                "'; retrying without it."
            )
            run_inner(
                kept_covariates,
                c(dropped_covariates, dropped_covariate)
            )
        })
    }

    run_inner(covariates)
}

de_dge_columns <- c("GENEID", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

de_dtu_transcript_columns <- c(
    "featureID",
    "groupID",
    "gene_name",
    "transcript_name",
    "log2FoldChange",
    "pvalue",
    "padj",
    "exonBaseMean"
)

#' Extract transcript-level DTU columns for TSV output.
#'
#' Renames the contrast-specific DEXSeq fold-change column to
#' `log2FoldChange`, normalizes data for TSV output
#' and returns only the transcript output columns.
#'
#' @param dex_df DEXSeq results as a data frame.
#' @param contrast_name Contrast suffix used in the DEXSeq fold-change column.
#'
#' @return A normalized data frame ready for `results_dtu_transcript.tsv`.
de_extract_dtu_transcript_table <- function(dex_df, contrast_name) {
    log2fold_column <- paste0("log2fold_", contrast_name)
    if (log2fold_column %in% names(dex_df)) {
        names(dex_df)[names(dex_df) == log2fold_column] <- "log2FoldChange"
    }
    tx_dtu <- dex_df[, intersect(
        de_dtu_transcript_columns,
        names(dex_df)
    ), drop = FALSE]
    workflow_glue_r_normalise_tsv_df(tx_dtu)
}

de_postprocess_dexseq_result <- function(
    dex_res,
    contrast_dir,
    target_level,
    reference_level,
    contrast_samples,
    condition_column,
    tx_counts,
    contrast_name,
    per_gene_qvalue_fn = DEXSeq::perGeneQValue,
    plot_writer = function(dex_res, contrast_dir) {
        grDevices::pdf(file.path(contrast_dir, "results_dtu.pdf"))
        DESeq2::plotMA(dex_res$dxr, cex = 0.8, alpha = 0.05)
        DESeq2::plotDispEsts(dex_res$dxd)
        grDevices::dev.off()
    }
) {
    tryCatch(
        {
            dex_df <- as.data.frame(dex_res$dxr)
            dex_df <- workflow_glue_r_normalise_tsv_df(dex_df)
            tx_dtu <- de_extract_dtu_transcript_table(dex_df, contrast_name)
            gene_q <- per_gene_qvalue_fn(dex_res$dxr)
            gene_dtu <- data.frame(
                GENEID = names(gene_q),
                qval = unname(gene_q),
                row.names = NULL
            )
            plot_writer(dex_res, contrast_dir)
            list(
                ok = TRUE,
                dex_df = dex_df,
                tx_dtu = tx_dtu,
                gene_dtu = gene_dtu,
                failure_hint = NULL
            )
        },
        error = function(err) {
            message_text <- conditionMessage(err)
            failure_hint <- de_dexseq_failure_hint(message_text)

            warning(
                "DEXSeq post-processing failed for contrast ",
                target_level,
                " vs ",
                reference_level,
                if (!is.null(failure_hint)) {
                    paste0("\nProbable cause: ", failure_hint)
                } else {
                    ""
                },
                "\nError: ",
                message_text
            )
            de_write_dtu_failure_outputs(
                contrast_dir = contrast_dir,
                target_level = target_level,
                reference_level = reference_level,
                contrast_samples = contrast_samples,
                condition_column = condition_column,
                tx_counts = tx_counts,
                message_text = message_text,
                failure_hint = failure_hint
            )
            list(
                ok = FALSE,
                dex_df = NULL,
                tx_dtu = NULL,
                gene_dtu = NULL,
                failure_hint = failure_hint
            )
        }
    )
}

main_run_de_analysis <- function(args) {
    set.seed(42)
    dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

    tx_se <- readRDS(args$transcript_rds)
    gene_se <- readRDS(args$gene_rds)
    sample_df <- workflow_glue_r_read_sample_sheet(args$sample_sheet)
    validated <- de_validate_inputs(tx_se, gene_se, sample_df, args)
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
        condition_column = args$condition_column,
        reference_level = reference_level,
        covariates = if (length(covariates) > 0) covariates else "none",
        num_contrasts = length(targets),
        contrasts = list()
    )

    n_per_group <- table(sample_df[[args$condition_column]])
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
        writeLines(mt_content, file.path(args$out_dir, "MULTIPLE_TESTING_WARNING.txt"))
    }

    for (target_level in targets) {
        contrast_name <- sprintf(
            "%s_%s_vs_%s",
            args$condition_column,
            target_level,
            reference_level
        )
        contrast_dir <- file.path(args$out_dir, contrast_name)
        dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

        keep_samples <- sample_df[[args$condition_column]] %in% c(reference_level, target_level)
        contrast_samples <- droplevels(sample_df[keep_samples, , drop = FALSE])
        contrast_samples[[args$condition_column]] <- stats::relevel(
            factor(contrast_samples[[args$condition_column]]),
            ref = reference_level
        )

        contrast_qc <- list(
            name = contrast_name,
            target_level = target_level,
            reference_level = reference_level,
            n_samples = nrow(contrast_samples),
            n_target = sum(contrast_samples[[args$condition_column]] == target_level),
            n_reference = sum(contrast_samples[[args$condition_column]] == reference_level),
            deseq2_size_factor_method = "ratio",
            deseq2_dispersion_fallback = list(
                applied = FALSE,
                method_used = "parametric",
                reason = NULL,
                diagnostic_file = NULL
            ),
            dexseq_size_factor_method = "ratio",
            dexseq_dispersion_method = "parametric",
            dexseq_covariates_dropped = list(),
            dtu_failure_hint = NULL,
            dge_status = "SUCCESS"
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

        dge_run <- tryCatch(
            de_run_deseq2_result(
                gene_counts,
                contrast_samples,
                target_level,
                reference_level,
                args$condition_column,
                covariates,
                args$out_dir,
                contrast_name
            ),
            error = function(err) {
                warning(
                    "DESeq2 DGE failed for contrast ",
                    target_level, " vs ", reference_level,
                    "\nError: ", conditionMessage(err)
                )
                failure_content <- c(
                    "DGE Analysis Failed",
                    "===================",
                    "",
                    sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                    sprintf("Contrast: %s vs %s", target_level, reference_level),
                    sprintf(
                        "Samples: %d (%d %s, %d %s)",
                        nrow(contrast_samples),
                        sum(contrast_samples[[args$condition_column]] == target_level),
                        target_level,
                        sum(contrast_samples[[args$condition_column]] == reference_level),
                        reference_level
                    ),
                    sprintf("Genes: %d", nrow(gene_counts)),
                    "",
                    "ERROR MESSAGE:",
                    sprintf("  %s", conditionMessage(err)),
                    "",
                    "DGE RESULTS CANNOT BE INTERPRETED",
                    ""
                )
                writeLines(failure_content, file.path(contrast_dir, "DGE_ANALYSIS_FAILED.txt"))
                NULL
            }
        )

        if (is.null(dge_run)) {
            dge_res <- workflow_glue_r_empty_tsv(de_dge_columns)
            contrast_qc$dge_status <- "FAILED"
            contrast_qc$dge_total_genes <- 0
            contrast_qc$dge_significant_fdr05 <- 0
            contrast_qc$dge_significant_fdr01 <- 0
            contrast_qc$dge_upregulated <- 0
            contrast_qc$dge_downregulated <- 0
            de_write_placeholder_pdf(
                file.path(contrast_dir, "results_dge.pdf"),
                "DESeq2 did not converge for this contrast.\nSee DGE_ANALYSIS_FAILED.txt for details."
            )
        } else {
            if (!is.null(dge_run$deseq2_dispersion_fallback)) {
                fallback <- dge_run$deseq2_dispersion_fallback
                fallback_applied <- isTRUE(fallback$applied)
                fallback_method <- fallback$method_used
                if (is.null(fallback_method) || identical(fallback_method, "")) {
                    fallback_method <- if (fallback_applied) "gene-wise" else "parametric"
                }
                if (!is.null(fallback$size_factor_type) && !identical(fallback$size_factor_type, "")) {
                    contrast_qc$deseq2_size_factor_method <- fallback$size_factor_type
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
            contrast_qc$dge_status <- "SUCCESS"
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
            grDevices::pdf(file.path(contrast_dir, "results_dge.pdf"))
            DESeq2::plotMA(dge_run$result)
            grDevices::dev.off()
        }

        utils::write.table(
            dge_res[order(dge_res$padj), ],
            file = file.path(contrast_dir, "results_dge.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        dex_run <- tryCatch(
            list(
                ok = TRUE,
                result = de_run_dexseq_result(
                    tx_counts,
                    tx_meta,
                    contrast_samples,
                    args$condition_column,
                    covariates
                ),
                failure_hint = NULL
            ),
            error = function(err) {
                message_text <- conditionMessage(err)
                is_known_failure <- de_is_recoverable_dexseq_error(message_text)
                failure_hint <- de_dexseq_failure_hint(message_text)

                warning(
                    "DEXSeq failed for contrast ",
                    target_level,
                    " vs ",
                    reference_level,
                    if (!is_known_failure) {
                        "\nThis error did not match a known statistical fallback pattern."
                    } else {
                        ""
                    },
                    if (!is.null(failure_hint)) {
                        paste0("\nProbable cause: ", failure_hint)
                    } else {
                        ""
                    },
                    "\nError: ",
                    message_text
                )
                de_write_dtu_failure_outputs(
                    contrast_dir = contrast_dir,
                    target_level = target_level,
                    reference_level = reference_level,
                    contrast_samples = contrast_samples,
                    condition_column = args$condition_column,
                    tx_counts = tx_counts,
                    message_text = message_text,
                    failure_hint = failure_hint
                )
                list(ok = FALSE, result = NULL, failure_hint = failure_hint)
            }
        )

        if (!isTRUE(dex_run$ok)) {
            dex_res <- NULL
            dtu_failure_hint <- dex_run$failure_hint
        } else {
            dex_res <- dex_run$result
            dtu_failure_hint <- NULL
        }

        if (is.null(dex_res)) {
            dex_df <- workflow_glue_r_empty_tsv(de_dtu_transcript_columns)
            tx_dtu <- dex_df
            gene_dtu <- workflow_glue_r_empty_tsv(c("GENEID", "qval"))
            de_write_placeholder_pdf(
                file.path(contrast_dir, "results_dtu.pdf"),
                "DEXSeq did not converge for this contrast.\nSee DTU_ANALYSIS_FAILED.txt for details."
            )
            contrast_qc$dtu_status <- "FAILED"
            contrast_qc$dtu_failure_hint <- dtu_failure_hint
            contrast_qc$dtu_significant_transcripts <- 0
            contrast_qc$dtu_significant_genes <- 0
        } else {
            if (!is.null(dex_res$dexseq_size_factor_type) && !identical(dex_res$dexseq_size_factor_type, "")) {
                contrast_qc$dexseq_size_factor_method <- dex_res$dexseq_size_factor_type
            }
            if (!is.null(dex_res$dexseq_dispersion_method)) {
                contrast_qc$dexseq_dispersion_method <- dex_res$dexseq_dispersion_method
            }
            if (!is.null(dex_res$dexseq_covariates_dropped)) {
                contrast_qc$dexseq_covariates_dropped <- as.list(dex_res$dexseq_covariates_dropped)
            }

            postprocess_run <- de_postprocess_dexseq_result(
                dex_res = dex_res,
                contrast_dir = contrast_dir,
                target_level = target_level,
                reference_level = reference_level,
                contrast_samples = contrast_samples,
                condition_column = args$condition_column,
                tx_counts = tx_counts,
                contrast_name = paste(target_level, reference_level, sep = "_")
            )

            if (!isTRUE(postprocess_run$ok)) {
                dex_df <- workflow_glue_r_empty_tsv(de_dtu_transcript_columns)
                tx_dtu <- dex_df
                gene_dtu <- workflow_glue_r_empty_tsv(c("GENEID", "qval"))
                de_write_placeholder_pdf(
                    file.path(contrast_dir, "results_dtu.pdf"),
                    "DEXSeq did not converge for this contrast.\nSee DTU_ANALYSIS_FAILED.txt for details."
                )
                contrast_qc$dtu_status <- "FAILED"
                contrast_qc$dtu_failure_hint <- postprocess_run$failure_hint
                contrast_qc$dtu_significant_transcripts <- 0
                contrast_qc$dtu_significant_genes <- 0
            } else {
                dex_df <- postprocess_run$dex_df
                tx_dtu <- postprocess_run$tx_dtu
                if ("featureID" %in% names(tx_dtu)) {
                    tx_match <- match(as.character(tx_dtu$featureID), as.character(tx_meta$TXNAME))
                    if ("gene_name" %in% names(tx_meta)) {
                        tx_dtu$gene_name <- tx_meta$gene_name[tx_match]
                    }
                    if ("transcript_name" %in% names(tx_meta)) {
                        tx_dtu$transcript_name <- tx_meta$transcript_name[tx_match]
                    }
                }
                gene_dtu <- postprocess_run$gene_dtu

                contrast_qc$dtu_status <- "SUCCESS"
                contrast_qc$dtu_significant_transcripts <- sum(tx_dtu$padj < 0.05, na.rm = TRUE)
                contrast_qc$dtu_significant_genes <- sum(gene_dtu$qval < 0.05, na.rm = TRUE)
            }
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
            sprintf("  Status: %s", contrast_qc$dge_status),
            sprintf("  Genes tested: %d", contrast_qc$genes_tested),
            sprintf("  Size factor method: %s", contrast_qc$deseq2_size_factor_method),
            if (contrast_qc$dge_status == "SUCCESS") {
                c(
                    sprintf("  Significant (FDR < 0.05): %d", contrast_qc$dge_significant_fdr05),
                    sprintf("  Significant (FDR < 0.01): %d", contrast_qc$dge_significant_fdr01),
                    sprintf("  Upregulated: %d", contrast_qc$dge_upregulated),
                    sprintf("  Downregulated: %d", contrast_qc$dge_downregulated)
                )
            } else {
                c(
                    "  Significant (FDR < 0.05): N/A",
                    "  Significant (FDR < 0.01): N/A",
                    "  Upregulated: N/A",
                    "  Downregulated: N/A",
                    "  See DGE_ANALYSIS_FAILED.txt for details"
                )
            },
            "",
            "DTU Results:",
            sprintf("  Status: %s", contrast_qc$dtu_status),
            sprintf("  Transcripts tested: %d", contrast_qc$transcripts_tested),
            sprintf("  Size factor method: %s", contrast_qc$dexseq_size_factor_method),
            if (contrast_qc$dtu_status == "SUCCESS") {
                c(
                    sprintf(
                        "  Significant transcripts (FDR < 0.05): %d",
                        contrast_qc$dtu_significant_transcripts
                    ),
                    sprintf("  Genes with DTU (q < 0.05): %d", contrast_qc$dtu_significant_genes)
                )
            } else {
                c(
                    if (!is.null(contrast_qc$dtu_failure_hint)) {
                        paste0("  Probable cause: ", contrast_qc$dtu_failure_hint)
                    } else {
                        NULL
                    },
                    "  See DTU_ANALYSIS_FAILED.txt for details"
                )
            },
            if (!is.null(contrast_qc$dtu_power_warning)) paste0("  WARNING: ", contrast_qc$dtu_power_warning) else NULL,
            ""
        )
        writeLines(contrast_qc_summary, file.path(contrast_dir, "contrast_qc_summary.txt"))
        de_qc_stats$contrasts[[contrast_name]] <- contrast_qc
    }

    failed_dge_contrasts <- names(Filter(
        function(cqc) identical(cqc$dge_status, "FAILED"),
        de_qc_stats$contrasts
    ))
    failed_dtu_contrasts <- names(Filter(
        function(cqc) identical(cqc$dtu_status, "FAILED"),
        de_qc_stats$contrasts
    ))

    deseq2_dispersion_fallbacks <- names(Filter(
        function(cqc) isTRUE(cqc$deseq2_dispersion_fallback$applied),
        de_qc_stats$contrasts
    ))
    deseq2_gene_wise <- names(Filter(
        function(cqc) identical(cqc$deseq2_dispersion_fallback$method_used, "gene-wise"),
        de_qc_stats$contrasts
    ))
    deseq2_poscounts <- names(Filter(
        function(cqc) identical(cqc$deseq2_size_factor_method, "poscounts"),
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
    dexseq_poscounts <- names(Filter(
        function(cqc) identical(cqc$dexseq_size_factor_method, "poscounts"),
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
        failed_dge_contrasts = length(failed_dge_contrasts),
        failed_dge_contrast_names = as.list(failed_dge_contrasts),
        failed_dtu_contrasts = length(failed_dtu_contrasts),
        failed_dtu_contrast_names = as.list(failed_dtu_contrasts),
        deseq2_dispersion_fallback_contrasts = length(deseq2_dispersion_fallbacks),
        deseq2_dispersion_fallback_contrast_names = as.list(deseq2_dispersion_fallbacks),
        deseq2_gene_wise_contrasts = length(deseq2_gene_wise),
        deseq2_gene_wise_contrast_names = as.list(deseq2_gene_wise),
        deseq2_poscounts_contrasts = length(deseq2_poscounts),
        deseq2_poscounts_contrast_names = as.list(deseq2_poscounts),
        dexseq_non_parametric_dispersion_contrasts = length(dexseq_non_parametric),
        dexseq_non_parametric_dispersion_contrast_names = as.list(dexseq_non_parametric),
        dexseq_gene_wise_dispersion_contrasts = length(dexseq_gene_wise),
        dexseq_gene_wise_dispersion_contrast_names = as.list(dexseq_gene_wise),
        dexseq_poscounts_contrasts = length(dexseq_poscounts),
        dexseq_poscounts_contrast_names = as.list(dexseq_poscounts),
        dexseq_covariate_drop_contrasts = length(dexseq_covariate_drop),
        dexseq_covariate_drop_contrast_names = as.list(dexseq_covariate_drop),
        total_covariates_dropped = total_covariates_dropped
    )

    jsonlite::write_json(
        de_qc_stats,
        file.path(args$out_dir, "de_qc_stats.json"),
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
                sprintf("    DGE status: %s", cqc$dge_status),
                if (cqc$dge_status == "SUCCESS") {
                    sprintf("    DGE significant: %d genes (FDR<0.05)", cqc$dge_significant_fdr05)
                } else {
                    "    DGE significant: N/A (analysis failed)"
                },
                sprintf("    DESeq2 size factors: %s", cqc$deseq2_size_factor_method),
                if (cqc$dge_status != "SUCCESS") "    See DGE_ANALYSIS_FAILED.txt for details" else NULL,
                sprintf("    DTU status: %s", cqc$dtu_status),
                sprintf("    DEXSeq size factors: %s", cqc$dexseq_size_factor_method),
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
    writeLines(unlist(overall_summary), file.path(args$out_dir, "de_overall_summary.txt"))
    writeLines(capture.output(sessionInfo()), file.path(args$out_dir, "session_info.txt"))

    invisible(list(qc = de_qc_stats))
}

run_de_analysis_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
    parsed <- argparser::parse_args(de_analysis_arg_parser(), argv = argv)
    args <- workflow_glue_r_normalise_args(parsed, de_analysis_arg_spec(), raw_argv = argv)
    main_run_de_analysis(args)
}
