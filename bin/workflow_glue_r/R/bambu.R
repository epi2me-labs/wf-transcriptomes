bambu_arg_spec <- function() {
    list(
        list(
            name = "bams",
            flag = "--bams",
            help = "Comma-separated BAM paths.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "aliases",
            flag = "--aliases",
            help = "Comma-separated aliases for --bams.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "sample_sheet",
            flag = "--sample_sheet",
            help = "Optional sample sheet CSV.",
            type = "character"
        ),
        list(
            name = "annotation",
            flag = "--annotation",
            help = "Reference annotation GTF/GFF.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "genome",
            flag = "--genome",
            help = "Reference genome FASTA.",
            type = "character",
            required = TRUE
        ),
        list(
            name = "transcriptome_mode",
            flag = "--transcriptome_mode",
            help = "discover or fixed_annotation.",
            type = "character",
            default = "discover",
            choices = c("discover", "fixed_annotation")
        ),
        list(
            name = "threads",
            flag = "--threads",
            help = "Number of worker threads.",
            type = "integer",
            default = 1L,
            min = 1L
        ),
        list(
            name = "ndr",
            flag = "--ndr",
            help = "Optional novel discovery rate.",
            type = "numeric",
            min = 0,
            max = 1,
            value_error = "NDR (Novel Discovery Rate) must be between 0 and 1"
        ),
        list(
            name = "out_dir",
            flag = "--out_dir",
            help = "Output directory.",
            type = "character",
            required = TRUE
        )
    )
}

bambu_arg_parser <- function() {
    workflow_glue_r_arg_parser_from_spec(
        "Run bambu transcript discovery and quantification.",
        bambu_arg_spec()
    )
}

bambu_resolve_inputs <- function(args) {
    sample_df <- NULL
    if (!is.null(args$sample_sheet)) {
        sample_df <- utils::read.csv(
            args$sample_sheet,
            check.names = FALSE,
            stringsAsFactors = FALSE
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
    }

    bam_paths <- workflow_glue_r_parse_csv_list(args$bams)
    aliases <- workflow_glue_r_parse_csv_list(args$aliases)

    if (length(bam_paths) < 1) {
        stop("No BAM files were provided in --bams.", call. = FALSE)
    }
    if (length(aliases) != length(bam_paths)) {
        stop("Provide one alias per BAM in --bams.", call. = FALSE)
    }

    duplicate_bam_aliases <- unique(aliases[duplicated(aliases)])
    if (length(duplicate_bam_aliases) > 0) {
        stop(
            sprintf(
                "BAM aliases must be unique; duplicated aliases: %s",
                paste(duplicate_bam_aliases, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (!is.null(sample_df)) {
        missing_aliases <- setdiff(aliases, sample_df$alias)
        if (length(missing_aliases) > 0) {
            stop(
                sprintf(
                    "Sample sheet is missing alias rows for BAM files: %s",
                    paste(missing_aliases, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        sample_df <- sample_df[match(aliases, sample_df$alias), , drop = FALSE]
    } else {
        sample_df <- data.frame(alias = aliases, stringsAsFactors = FALSE)
    }

    reads <- if (length(bam_paths) == 1) {
        bam_paths
    } else {
        Rsamtools::BamFileList(bam_paths, yieldSize = 250000L)
    }

    list(
        bam_paths = bam_paths,
        aliases = aliases,
        sample_df = sample_df,
        reads = reads
    )
}

bambu_discovery_enabled <- function(args) {
    identical(args$transcriptome_mode, "discover")
}

bambu_build_args <- function(args, reads, annotation_obj) {
    bambu_args <- list(
        reads = reads,
        annotations = annotation_obj,
        genome = args$genome,
        ncore = args$threads,
        discovery = bambu_discovery_enabled(args),
        lowMemory = TRUE,
        yieldSize = 250000L,
        verbose = TRUE
    )

    if (bambu_discovery_enabled(args) && !is.null(args$ndr)) {
        bambu_args$NDR <- args$ndr
    }

    bambu_args
}

bambu_effective_threads <- function(args, bam_count) {
    # bambu's low-memory mode can have issues with multiple BAMs and
    # parallel threads due to BiocFileCache writes, 
    # so we enforce single-threading in that case.
    threads <- as.integer(args$threads)
    if (bam_count > 1 && threads > 1L) {
        warning(
            paste(
                "Low-memory mode with multiple BAMs can fail in bambu due to",
                "parallel BiocFileCache writes; forcing threads=1."
            ),
            call. = FALSE
        )
        return(1L)
    }
    threads
}

bambu_filter_transcripts <- function(se) {
    assays <- SummarizedExperiment::assays(se)
    counts_mat <- assays$counts

    row_data <- SummarizedExperiment::rowData(se)
    if (!"GENEID" %in% names(row_data)) {
        stop("rowData(se) does not contain required column 'GENEID'.")
    }
    if (is.null(counts_mat)) {
        stop("counts assay not found in bambu output.")
    }

    gene_ids <- row_data$GENEID
    qc_stats <- list(
        total_transcripts_before_filter = nrow(se),
        total_genes_before_filter = length(unique(gene_ids)),
        samples = ncol(se)
    )

    keep_idx <- rowSums(counts_mat) > 0
    if (!any(keep_idx)) {
        stop(
            "All transcripts have zero counts after filtering. ",
            "This suggests a problem with the bambu analysis or input data. ",
            "Check bambu logs and verify input quality."
        )
    }

    qc_stats$transcripts_filtered <- sum(!keep_idx)
    se <- se[keep_idx, ]

    # Keep incompatible gene-level counts in sync with the filtered transcript set.
    # transcriptToGeneExpression() expects incompatibleCounts GENEIDs to be a
    # subset of rowData(se)$GENEID after filtering.
    sample_names <- colnames(se)
    empty_incompatible_counts <- function() {
        cols <- c(
            list(GENEID = character(0)),
            stats::setNames(rep(list(numeric(0)), length(sample_names)), sample_names)
        )
        data.table::as.data.table(cols)
    }

    incompatible_counts <- S4Vectors::metadata(se)$incompatibleCounts
    if (is.null(incompatible_counts)) {
        incompatible_counts <- empty_incompatible_counts()
    } else {
        if (!"GENEID" %in% names(incompatible_counts)) {
            incompatible_counts <- empty_incompatible_counts()
        } else {
            kept_genes <- unique(SummarizedExperiment::rowData(se)$GENEID)
            incompatible_gene_ids <- incompatible_counts$GENEID

            rows_to_keep <- incompatible_gene_ids %in% kept_genes
            incompatible_counts <- incompatible_counts[rows_to_keep, , drop = FALSE]
            data.table::set(
                incompatible_counts,
                j = "GENEID",
                value = incompatible_gene_ids[rows_to_keep]
            )

            for (sample_name in sample_names) {
                if (!sample_name %in% names(incompatible_counts)) {
                    data.table::set(
                        incompatible_counts,
                        j = sample_name,
                        value = numeric(nrow(incompatible_counts))
                    )
                }
            }

            keep_cols <- c("GENEID", sample_names)
            incompatible_counts <- incompatible_counts[, keep_cols, with = FALSE]
        }
    }
    S4Vectors::metadata(se)$incompatibleCounts <- incompatible_counts

    qc_stats$total_transcripts_after_filter <- nrow(se)
    qc_stats$total_genes_after_filter <- length(unique(SummarizedExperiment::rowData(se)$GENEID))

    list(se = se, qc_stats = qc_stats)
}

bambu_matrix_to_df <- function(se_obj, assay_name, id_col, meta_df) {
    assay_df <- as.data.frame(SummarizedExperiment::assays(se_obj)[[assay_name]])
    assay_df[[id_col]] <- rownames(se_obj)
    assay_df <- assay_df[, c(id_col, setdiff(names(assay_df), id_col)), drop = FALSE]
    merge(meta_df, assay_df, by.x = id_col, by.y = id_col, all.y = TRUE, sort = FALSE)
}

bambu_format_count <- function(value) {
    if (length(value) == 0 || all(is.na(value))) {
        return("NA")
    }
    format(
        round(as.numeric(value), 0),
        scientific = FALSE,
        trim = TRUE,
        big.mark = ","
    )
}

bambu_write_matrix_tsv <- function(se_obj, assay_name, id_col, meta_df, output_path) {
    table_df <- bambu_matrix_to_df(se_obj, assay_name, id_col, meta_df)
    utils::write.table(
        table_df,
        file = output_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    rm(table_df)
    invisible(gc(verbose = FALSE))
}

bambu_write_outputs <- function(se, gene_se, sample_df, args, qc_stats) {
    bambu::writeToGTF(
        SummarizedExperiment::rowRanges(se),
        file = file.path(args$out_dir, "transcripts.gtf")
    )

    saveRDS(se, file.path(args$out_dir, "bambu_transcripts.rds"))
    saveRDS(gene_se, file.path(args$out_dir, "bambu_genes.rds"))
    utils::write.csv(
        sample_df,
        file.path(args$out_dir, "samples.csv"),
        row.names = FALSE,
        quote = FALSE
    )

    tx_meta <- as.data.frame(SummarizedExperiment::rowData(se))
    if (!"TXNAME" %in% names(tx_meta)) {
        tx_meta$TXNAME <- rownames(se)
    }
    if (!"GENEID" %in% names(tx_meta)) {
        tx_meta$GENEID <- NA_character_
    }

    gene_meta <- as.data.frame(SummarizedExperiment::rowData(gene_se))
    if (!"GENEID" %in% names(gene_meta)) {
        gene_meta$GENEID <- rownames(gene_se)
    }

    tx_meta <- workflow_glue_r_normalise_tsv_df(tx_meta)
    gene_meta <- workflow_glue_r_normalise_tsv_df(gene_meta)

    utils::write.table(
        tx_meta,
        file = file.path(args$out_dir, "transcript_metadata.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    utils::write.table(
        gene_meta,
        file = file.path(args$out_dir, "gene_metadata.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    bambu_write_matrix_tsv(
        se,
        "counts",
        "TXNAME",
        tx_meta,
        file.path(args$out_dir, "transcript_counts.tsv")
    )
    bambu_write_matrix_tsv(
        se,
        "CPM",
        "TXNAME",
        tx_meta,
        file.path(args$out_dir, "transcript_cpm.tsv")
    )
    bambu_write_matrix_tsv(
        gene_se,
        "counts",
        "GENEID",
        gene_meta,
        file.path(args$out_dir, "gene_counts.tsv")
    )
    bambu_write_matrix_tsv(
        gene_se,
        "CPM",
        "GENEID",
        gene_meta,
        file.path(args$out_dir, "gene_cpm.tsv")
    )

    qc_stats$transcriptome_mode <- args$transcriptome_mode
    qc_stats$ndr_used <- if (bambu_discovery_enabled(args)) {
        if (is.null(args$ndr)) "automatic" else args$ndr
    } else {
        "N/A"
    }
    qc_stats$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    jsonlite::write_json(
        qc_stats,
        file.path(args$out_dir, "bambu_qc_stats.json"),
        pretty = TRUE,
        auto_unbox = TRUE
    )

    qc_summary <- c(
        "Bambu Quantification QC Summary",
        "================================",
        "",
        sprintf("Timestamp: %s", qc_stats$timestamp),
        sprintf("Mode: %s", args$transcriptome_mode),
        if (bambu_discovery_enabled(args)) {
            if (is.null(args$ndr)) {
                "NDR: automatic (bambu-selected)"
            } else {
                sprintf("NDR: %.3f", args$ndr)
            }
        } else NULL,
        "",
        "Sample Statistics:",
        sprintf("  Samples analyzed: %s", bambu_format_count(qc_stats$samples)),
        sprintf("  Median library size: %s reads", bambu_format_count(qc_stats$median_library_size)),
        sprintf(
            "  Library size range: %s - %s reads",
            bambu_format_count(qc_stats$min_library_size),
            bambu_format_count(qc_stats$max_library_size)
        ),
        if (!is.null(qc_stats$library_size_warning)) sprintf("  WARNING: %s", qc_stats$library_size_warning) else NULL,
        "",
        "Transcript Discovery:",
        sprintf("  Transcripts before filtering: %s", bambu_format_count(qc_stats$total_transcripts_before_filter)),
        sprintf("  Transcripts after filtering: %s", bambu_format_count(qc_stats$total_transcripts_after_filter)),
        sprintf("  Transcripts removed: %s", bambu_format_count(qc_stats$transcripts_filtered)),
        sprintf(
            "  Median transcripts detected per sample: %s",
            bambu_format_count(qc_stats$median_transcripts_detected)
        ),
        "",
        "Gene-Level Summary:",
        sprintf("  Unique genes (before filter): %s", bambu_format_count(qc_stats$total_genes_before_filter)),
        sprintf("  Unique genes (after filter): %s", bambu_format_count(qc_stats$total_genes_after_filter)),
        ""
    )

    writeLines(qc_summary, file.path(args$out_dir, "bambu_qc_summary.txt"))
    writeLines(capture.output(sessionInfo()), file.path(args$out_dir, "session_info.txt"))
}

main_run_bambu <- function(args) {
    set.seed(42)
    # bambu's parallel worker code may rely on these being attached for generics
    # such as seqlengths().
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(Rsamtools))

    dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

    inputs <- bambu_resolve_inputs(args)
    args$threads <- bambu_effective_threads(args, length(inputs$bam_paths))
    annotation_obj <- bambu::prepareAnnotations(args$annotation)

    if (!is.null(args$ndr)) {
        message(sprintf("Using user-specified NDR = %.3f", args$ndr))
    } else if (bambu_discovery_enabled(args)) {
        message("Using bambu automatic NDR selection.")
    }

    if (bambu_discovery_enabled(args)) {
        message("Novel Discovery Rate (NDR) controls transcript discovery stringency:")
        message("  Lower NDR (e.g., 0.05) = fewer false positive transcripts, may miss real ones")
        message("  Higher NDR (e.g., 0.2) = more sensitive discovery, more false positives")
        if (is.null(args$ndr)) {
            message("  Current NDR = automatic (selected by bambu from the data)")
        } else {
            message(sprintf("  Current NDR = %.3f balances precision and recall", args$ndr))
        }
    }
    if (length(inputs$bam_paths) > 1) {
        message("Using BamFileList yieldSize = 250000")
    }
    message(sprintf("Running bambu with threads = %d", args$threads))

    message("Running bambu...")
    se <- do.call(bambu::bambu, bambu_build_args(args, inputs$reads, annotation_obj))
    message("Bambu completed successfully")
    colnames(se) <- inputs$aliases

    filtered <- bambu_filter_transcripts(se)
    se <- filtered$se
    qc_stats <- filtered$qc_stats
    gene_se <- bambu::transcriptToGeneExpression(se)
    colnames(gene_se) <- inputs$aliases
    message(
        sprintf(
            "Filtering: keeping %d / %d transcripts",
            qc_stats$total_transcripts_after_filter,
            qc_stats$total_transcripts_before_filter
        )
    )

    lib_sizes <- colSums(SummarizedExperiment::assays(se)$counts)
    qc_stats$library_sizes <- as.list(lib_sizes)
    qc_stats$min_library_size <- min(lib_sizes)
    qc_stats$max_library_size <- max(lib_sizes)
    qc_stats$median_library_size <- stats::median(lib_sizes)

    if (length(lib_sizes) > 1) {
        lib_size_ratio <- max(lib_sizes) / min(lib_sizes)
        qc_stats$library_size_ratio <- lib_size_ratio
        if (lib_size_ratio > 3) {
            warning(
                sprintf(
                    paste0(
                        "Large library size variation detected (%.1fx difference).\n",
                        "  Min: %d, Max: %d reads.\n",
                        "  CPM normalization may not be appropriate for such variation."
                    ),
                    lib_size_ratio,
                    min(lib_sizes),
                    max(lib_sizes)
                )
            )
            qc_stats$library_size_warning <- sprintf("%.1fx variation (>3x threshold)", lib_size_ratio)
        }
    }

    detected_per_sample <- colSums(SummarizedExperiment::assays(se)$counts > 0)
    qc_stats$transcripts_detected_per_sample <- as.list(detected_per_sample)
    qc_stats$median_transcripts_detected <- stats::median(detected_per_sample)

    bambu_write_outputs(
        se,
        gene_se,
        inputs$sample_df,
        args,
        qc_stats
    )

    invisible(
        list(
            se = se,
            gene_se = gene_se,
            sample_df = inputs$sample_df,
            qc_stats = qc_stats
        )
    )
}

run_bambu_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
    parsed <- argparser::parse_args(bambu_arg_parser(), argv = argv)
    args <- workflow_glue_r_normalise_args(parsed, bambu_arg_spec(), raw_argv = argv)
    main_run_bambu(args)
}
