bambu_arg_parser <- function() {
    parser <- argparser::arg_parser("Run bambu transcript discovery and quantification.")
    parser <- argparser::add_argument(parser, "--bams", help = "Comma-separated BAM paths.")
    parser <- argparser::add_argument(parser, "--aliases", help = "Comma-separated aliases for --bams.")
    parser <- argparser::add_argument(parser, "--sample_sheet", help = "Optional sample sheet CSV.")
    parser <- argparser::add_argument(parser, "--annotation", help = "Reference annotation GTF/GFF.")
    parser <- argparser::add_argument(parser, "--genome", help = "Reference genome FASTA.")
    parser <- argparser::add_argument(
        parser,
        "--transcriptome_mode",
        help = "discover or fixed_annotation.",
        default = "discover"
    )
    parser <- argparser::add_argument(
        parser,
        "--threads",
        help = "Number of worker threads.",
        type = "numeric",
        default = 1
    )
    parser <- argparser::add_argument(
        parser,
        "--ndr",
        help = "Optional novel discovery rate.",
        type = "numeric"
    )
    argparser::add_argument(parser, "--out_dir", help = "Output directory.")
}

bambu_validate_args <- function(argv) {
    workflow_glue_r_require_args(argv, c("annotation", "genome", "out_dir"))

    if (workflow_glue_r_arg_missing(argv$bams)) {
        stop("Missing required arguments: --bams", call. = FALSE)
    }
    if (workflow_glue_r_arg_missing(argv$aliases)) {
        stop("Missing required arguments: --aliases", call. = FALSE)
    }

    if (!argv$transcriptome_mode %in% c("discover", "fixed_annotation")) {
        stop(
            sprintf(
                "transcriptome_mode must be one of: %s",
                paste(c("discover", "fixed_annotation"), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (!workflow_glue_r_arg_missing(argv$ndr) && (argv$ndr < 0 || argv$ndr > 1)) {
        stop("NDR (Novel Discovery Rate) must be between 0 and 1", call. = FALSE)
    }

    invisible(argv)
}

bambu_resolve_inputs <- function(
    argv,
    bamfile_list_ctor = Rsamtools::BamFileList
) {
    sample_df <- NULL
    if (!workflow_glue_r_arg_missing(argv$sample_sheet)) {
        sample_df <- workflow_glue_r_read_csv(argv$sample_sheet)
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

    bam_paths <- workflow_glue_r_parse_csv_list(argv$bams)
    aliases <- workflow_glue_r_parse_csv_list(argv$aliases)

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
        bamfile_list_ctor(bam_paths, yieldSize = 1000000)
    }

    list(
        bam_paths = bam_paths,
        aliases = aliases,
        sample_df = sample_df,
        reads = reads
    )
}

bambu_discovery_enabled <- function(argv) {
    identical(argv$transcriptome_mode, "discover")
}

bambu_resolve_ndr <- function(argv, default_ndr = 0.1) {
    if (workflow_glue_r_arg_missing(argv$ndr)) {
        default_ndr
    } else {
        as.numeric(argv$ndr)
    }
}

bambu_build_args <- function(argv, reads, annotation_obj) {
    bambu_args <- list(
        reads = reads,
        annotations = annotation_obj,
        genome = argv$genome,
        ncore = as.integer(argv$threads),
        discovery = bambu_discovery_enabled(argv)
    )

    if (bambu_discovery_enabled(argv)) {
        bambu_args$NDR <- bambu_resolve_ndr(argv)
    }

    bambu_args
}

bambu_filter_transcripts <- function(se) {
    counts_mat <- SummarizedExperiment::assays(se)$counts
    full_length_mat <- SummarizedExperiment::assays(se)$fullLengthCounts

    gene_ids <- SummarizedExperiment::rowData(se)$GENEID
    qc_stats <- list(
        total_transcripts_before_filter = nrow(se),
        total_genes_before_filter = length(unique(gene_ids)),
        samples = ncol(se)
    )

    if (is.null(full_length_mat)) {
        keep_idx <- rowSums(counts_mat) > 0
    } else {
        keep_idx <- rowSums(full_length_mat) > 0
    }
    if (!any(keep_idx)) {
        stop(
            "All transcripts have zero counts after filtering. ",
            "This suggests a problem with the bambu analysis or input data. ",
            "Check bambu logs and verify input quality."
        )
    }

    qc_stats$transcripts_filtered <- sum(!keep_idx)
    se <- se[keep_idx, ]
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

bambu_extract_gtf_attribute <- function(attr_field, key) {
    match <- regexec(sprintf('%s "([^"]*)";', key), attr_field, perl = TRUE)
    captures <- regmatches(attr_field, match)[[1]]
    if (length(captures) < 2) {
        return(NULL)
    }
    captures[2]
}

bambu_normalise_gtf_attribute_value <- function(value) {
    if (is.null(value)) {
        return(NULL)
    }
    value <- gsub('[";]', "", value)
    value <- trimws(gsub("\\s+", " ", value))
    if (!nzchar(value)) {
        return(NULL)
    }
    value
}

bambu_sanitise_gtf_file <- function(path) {
    lines <- readLines(path, warn = FALSE)
    cleaned_lines <- vapply(lines, function(line) {
        if (!nzchar(line) || startsWith(line, "#")) {
            return(line)
        }

        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) < 9) {
            return(line)
        }

        attr_field <- fields[9]
        transcript_id <- bambu_normalise_gtf_attribute_value(
            bambu_extract_gtf_attribute(attr_field, "transcript_id")
        )
        gene_id <- bambu_normalise_gtf_attribute_value(
            bambu_extract_gtf_attribute(attr_field, "gene_id")
        )

        if (!is.null(gene_id) && identical(gene_id, "transcript_id")) {
            warning(
                sprintf(
                    "Replaced malformed gene_id 'transcript_id' with transcript_id '%s'.",
                    transcript_id
                ),
                call. = FALSE
            )
            gene_id <- transcript_id
        }
        if (is.null(gene_id)) {
            gene_id <- transcript_id
        }

        if (!is.null(gene_id)) {
            attr_field <- sub(
                'gene_id "([^"]*)";',
                sprintf('gene_id "%s";', gene_id),
                attr_field,
                perl = TRUE
            )
        }
        if (!is.null(transcript_id)) {
            attr_field <- sub(
                'transcript_id "([^"]*)";',
                sprintf('transcript_id "%s";', transcript_id),
                attr_field,
                perl = TRUE
            )
        }

        fields[9] <- attr_field
        paste(fields, collapse = "\t")
    }, character(1))

    writeLines(cleaned_lines, path)
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

bambu_write_outputs <- function(se, gene_se, sample_df, argv, qc_stats, write_gtf_fn = bambu::writeToGTF) {
    write_gtf_fn(
        SummarizedExperiment::rowRanges(se),
        file = file.path(argv$out_dir, "transcripts.gtf")
    )

    saveRDS(se, file.path(argv$out_dir, "bambu_transcripts.rds"))
    saveRDS(gene_se, file.path(argv$out_dir, "bambu_genes.rds"))
    utils::write.csv(
        sample_df,
        file.path(argv$out_dir, "samples.csv"),
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
        file = file.path(argv$out_dir, "transcript_metadata.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    utils::write.table(
        gene_meta,
        file = file.path(argv$out_dir, "gene_metadata.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    tx_counts <- bambu_matrix_to_df(se, "counts", "TXNAME", tx_meta)
    tx_cpm <- bambu_matrix_to_df(se, "CPM", "TXNAME", tx_meta)
    gene_counts <- bambu_matrix_to_df(gene_se, "counts", "GENEID", gene_meta)
    gene_cpm <- bambu_matrix_to_df(gene_se, "CPM", "GENEID", gene_meta)

    utils::write.table(
        tx_counts,
        file = file.path(argv$out_dir, "transcript_counts.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    utils::write.table(
        tx_cpm,
        file = file.path(argv$out_dir, "transcript_cpm.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    utils::write.table(
        gene_counts,
        file = file.path(argv$out_dir, "gene_counts.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    utils::write.table(
        gene_cpm,
        file = file.path(argv$out_dir, "gene_cpm.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    bambu_sanitise_gtf_file(file.path(argv$out_dir, "transcripts.gtf"))

    qc_stats$transcriptome_mode <- argv$transcriptome_mode
    qc_stats$ndr_used <- if (bambu_discovery_enabled(argv)) {
        bambu_resolve_ndr(argv)
    } else {
        "N/A"
    }
    qc_stats$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    jsonlite::write_json(
        qc_stats,
        file.path(argv$out_dir, "bambu_qc_stats.json"),
        pretty = TRUE,
        auto_unbox = TRUE
    )

    qc_summary <- c(
        "Bambu Quantification QC Summary",
        "================================",
        "",
        sprintf("Timestamp: %s", qc_stats$timestamp),
        sprintf("Mode: %s", argv$transcriptome_mode),
        if (bambu_discovery_enabled(argv)) sprintf("NDR: %.3f", bambu_resolve_ndr(argv)) else NULL,
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

    writeLines(qc_summary, file.path(argv$out_dir, "bambu_qc_summary.txt"))
    writeLines(capture.output(sessionInfo()), file.path(argv$out_dir, "session_info.txt"))
}

main_run_bambu <- function(
    argv,
    analysis_fn = bambu::bambu,
    prepare_annotations_fn = bambu::prepareAnnotations,
    gene_expression_fn = bambu::transcriptToGeneExpression,
    write_gtf_fn = bambu::writeToGTF,
    bamfile_list_ctor = Rsamtools::BamFileList
) {
    set.seed(42)
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(Rsamtools)
    })

    bambu_validate_args(argv)
    dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

    inputs <- bambu_resolve_inputs(
        argv,
        bamfile_list_ctor = bamfile_list_ctor
    )
    annotation_obj <- prepare_annotations_fn(argv$annotation)
    ndr_value <- bambu_resolve_ndr(argv)

    if (!workflow_glue_r_arg_missing(argv$ndr)) {
        message(sprintf("Using user-specified NDR = %.3f", ndr_value))
    } else {
        message(sprintf("Using default NDR = %.3f", ndr_value))
    }

    if (bambu_discovery_enabled(argv)) {
        message("Novel Discovery Rate (NDR) controls transcript discovery stringency:")
        message("  Lower NDR (e.g., 0.05) = fewer false positive transcripts, may miss real ones")
        message("  Higher NDR (e.g., 0.2) = more sensitive discovery, more false positives")
        message(sprintf("  Current NDR = %.3f balances precision and recall", ndr_value))
    }

    message("Running bambu...")
    se <- do.call(analysis_fn, bambu_build_args(argv, inputs$reads, annotation_obj))
    message("Bambu completed successfully")
    colnames(se) <- inputs$aliases

    filtered <- bambu_filter_transcripts(se)
    se <- filtered$se
    qc_stats <- filtered$qc_stats
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

    gene_se <- gene_expression_fn(se)
    colnames(gene_se) <- inputs$aliases

    bambu_write_outputs(
        se,
        gene_se,
        inputs$sample_df,
        argv,
        qc_stats,
        write_gtf_fn = write_gtf_fn
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
    main_run_bambu(parsed)
}
