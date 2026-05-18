bambu_arg_spec <- function() {
    list(
        list(
            name = "mode",
            flag = "--mode",
            help = "Run mode: discover, quant, collate, or empty.",
            type = "character",
            required = TRUE,
            choices = c("discover", "quant", "collate", "empty")
        ),
        list(
            name = "bams",
            flag = "--bams",
            help = "Comma-separated BAM paths.",
            type = "character"
        ),
        list(
            name = "aliases",
            flag = "--aliases",
            help = "Comma-separated aliases for --bams.",
            type = "character"
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
            type = "character"
        ),
        list(
            name = "chunk_rds",
            flag = "--chunk_rds",
            help = "Path to a chunked rcFile bundle RDS for quant mode.",
            type = "character"
        ),
        list(
            name = "discovered_annotation_rds",
            flag = "--discovered_annotation_rds",
            help = "Path to a discovered annotation RDS for quant mode.",
            type = "character"
        ),
        list(
            name = "chunk_dirs",
            flag = "--chunk_dirs",
            help = "Comma-separated chunk quantification directories for collate mode.",
            type = "character"
        ),
        list(
            name = "genome",
            flag = "--genome",
            help = "Reference genome FASTA.",
            type = "character"
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

bambu_default_yield_size <- 250000L
bambu_default_seed <- 42L

bambu_missing <- function(value) {
    is.null(value) ||
        length(value) == 0 ||
        (length(value) == 1 && is.na(value)) ||
        (is.character(value) && length(value) == 1 && !nzchar(value))
}

bambu_validate_args <- function(args) {
    mode <- args$mode
    if (bambu_missing(mode)) {
        stop("Missing required arguments: --mode", call. = FALSE)
    }
    if (!mode %in% c("discover", "quant", "collate", "empty")) {
        stop("mode must be one of: discover, quant, collate, empty", call. = FALSE)
    }
    if (bambu_missing(args$out_dir)) {
        stop("Missing required arguments: --out_dir", call. = FALSE)
    }

    if (identical(mode, "discover")) {
        missing_args <- character(0)
        if (bambu_missing(args$bams)) missing_args <- c(missing_args, "--bams")
        if (bambu_missing(args$aliases)) missing_args <- c(missing_args, "--aliases")
        if (bambu_missing(args$annotation)) missing_args <- c(missing_args, "--annotation")
        if (bambu_missing(args$genome)) missing_args <- c(missing_args, "--genome")
        if (length(missing_args) > 0) {
            stop(
                sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")),
                call. = FALSE
            )
        }
    }

    if (identical(mode, "quant")) {
        missing_args <- character(0)
        if (bambu_missing(args$chunk_rds)) missing_args <- c(missing_args, "--chunk_rds")
        if (bambu_missing(args$discovered_annotation_rds)) {
            missing_args <- c(missing_args, "--discovered_annotation_rds")
        }
        if (bambu_missing(args$genome)) missing_args <- c(missing_args, "--genome")
        if (length(missing_args) > 0) {
            stop(
                sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")),
                call. = FALSE
            )
        }
    }

    if (identical(mode, "collate") && bambu_missing(args$chunk_dirs)) {
        stop("Missing required arguments: --chunk_dirs", call. = FALSE)
    }

    if (identical(mode, "empty") && bambu_missing(args$aliases)) {
        stop("Missing required arguments: --aliases", call. = FALSE)
    }

    invisible(args)
}

bambu_resolve_chunk_dirs <- function(args) {
    chunk_dirs <- workflow_glue_r_parse_csv_list(args$chunk_dirs)
    if (length(chunk_dirs) < 1) {
        stop("No chunk quantification directories were provided in --chunk_dirs.", call. = FALSE)
    }
    chunk_dirs
}

bambu_resolve_inputs <- function(args, bamfile_list_ctor = Rsamtools::BamFileList) {
    sample_df <- NULL
    if (!bambu_missing(args$sample_sheet)) {
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
        bamfile_list_ctor(bam_paths, yieldSize = bambu_default_yield_size)
    }
    if (is.list(reads) && length(bam_paths) > 1) {
        names(reads) <- aliases
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

bambu_build_args <- function(args, reads, annotation_obj, discovery, quant) {
    bambu_args <- list(
        reads = reads,
        annotations = annotation_obj,
        genome = args$genome,
        ncore = as.integer(args$threads),
        discovery = discovery,
        quant = quant,
        lowMemory = TRUE,
        yieldSize = bambu_default_yield_size,
        verbose = TRUE
    )

    if (discovery && !is.null(args$ndr)) {
        bambu_args$NDR <- args$ndr
    }

    if (quant) {
        # Chunked quant re-estimates degradation bias per chunk, which changes
        # the EM inputs and breaks equivalence with unchunked bambu quant.
        bambu_args$opt.em <- list(degradationBias = FALSE)
    }

    bambu_args
}

bambu_message_ndr <- function(args) {
    if (!is.null(args$ndr)) {
        message(sprintf("Using user-specified NDR = %.3f", args$ndr))
    } else {
        message("Using bambu automatic NDR selection.")
    }
    message("Novel Discovery Rate (NDR) controls transcript discovery stringency:")
    message("  Lower NDR (e.g., 0.05) = fewer false positive transcripts, may miss real ones")
    message("  Higher NDR (e.g., 0.2) = more sensitive discovery, more false positives")
    if (is.null(args$ndr)) {
        message("  Current NDR = automatic (selected by bambu from the data)")
    } else {
        message(sprintf("  Current NDR = %.3f balances precision and recall", args$ndr))
    }
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

bambu_normalise_rc_file_list <- function(rc_files, aliases = NULL) {
    if (!is.list(rc_files)) {
        rc_files <- list(rc_files)
    }
    if (!is.null(aliases) && is.null(names(rc_files)) && length(rc_files) == length(aliases)) {
        names(rc_files) <- aliases
    }
    rc_files
}

bambu_chunk_id_for_seqname <- function(seqname) {
    seqname <- as.character(seqname)
    seqname <- gsub("[^A-Za-z0-9._-]+", "_", seqname)
    seqname <- sub("^_+", "", seqname)
    seqname <- sub("_+$", "", seqname)
    if (!nzchar(seqname)) {
        seqname <- "chunk"
    }
    seqname
}

bambu_chunk_rc_files <- function(rc_files, aliases, sample_df) {
    rc_files <- bambu_normalise_rc_file_list(rc_files, aliases = aliases)
    seqnames <- unique(unlist(lapply(rc_files, function(rcf) {
        as.character(SummarizedExperiment::rowData(rcf)$chr.rc)
    })))
    seqnames <- seqnames[!is.na(seqnames) & nzchar(seqnames)]
    chunk_ids <- make.unique(vapply(seqnames, bambu_chunk_id_for_seqname, character(1)), sep = "_")

    Map(function(seqname, chunk_id) {
        rc_chunk <- lapply(rc_files, function(rcf) {
            idx <- as.character(SummarizedExperiment::rowData(rcf)$chr.rc) == seqname
            rcf[idx, , drop = FALSE]
        })
        if (!is.null(names(rc_files))) {
            names(rc_chunk) <- names(rc_files)
        }
        list(
            chunk_id = chunk_id,
            seqname = seqname,
            aliases = aliases,
            sample_df = sample_df,
            rc_files = rc_chunk
        )
    }, seqnames, chunk_ids)
}

bambu_write_discovery_outputs <- function(out_dir, rc_files, discovered_annotations, chunk_bundles, sample_df) {
    saveRDS(rc_files, file.path(out_dir, "bambu_rcfiles.rds"))
    saveRDS(discovered_annotations, file.path(out_dir, "bambu_discovered_annotations.rds"))
    utils::write.csv(
        sample_df,
        file.path(out_dir, "samples.csv"),
        row.names = FALSE,
        quote = FALSE
    )

    annotation_tx_counts <- bambu_annotation_tx_counts_by_seqname(discovered_annotations)
    chunks_dir <- file.path(out_dir, "chunks")
    dir.create(chunks_dir, showWarnings = FALSE, recursive = TRUE)

    manifest <- do.call(rbind, lapply(chunk_bundles, function(bundle) {
        annotation_tx_count <- unname(annotation_tx_counts[bundle$seqname])
        if (length(annotation_tx_count) != 1L || is.na(annotation_tx_count)) {
            annotation_tx_count <- 0L
        }
        bundle$annotation_tx_count <- annotation_tx_count
        chunk_path <- file.path(chunks_dir, sprintf("%s.rds", bundle$chunk_id))
        saveRDS(bundle, chunk_path)
        data.frame(
            chunk_id = bundle$chunk_id,
            seqname = bundle$seqname,
            annotation_tx_count = annotation_tx_count,
            rds_path = chunk_path,
            stringsAsFactors = FALSE
        )
    }))

    utils::write.table(
        manifest,
        file = file.path(out_dir, "chunk_manifest.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

bambu_run_discover_mode <- function(args, analysis_fn, prepare_annotations_fn, bamfile_list_ctor) {
    inputs <- bambu_resolve_inputs(args, bamfile_list_ctor = bamfile_list_ctor)
    args$threads <- bambu_effective_threads(args, length(inputs$bam_paths))
    annotation_obj <- prepare_annotations_fn(args$annotation)

    if (length(inputs$bam_paths) > 1) {
        message(sprintf("Using BamFileList yieldSize = %d", bambu_default_yield_size))
    }
    message(sprintf("Running bambu discover setup with threads = %d", args$threads))
    message("Generating bambu rcFiles...")
    rc_files <- bambu_call_analysis(
        analysis_fn,
        bambu_build_args(
            args,
            inputs$reads,
            annotation_obj,
            discovery = FALSE,
            quant = FALSE
        )
    )
    rc_files <- bambu_normalise_rc_file_list(rc_files, aliases = inputs$aliases)

    discovered_annotations <- if (bambu_discovery_enabled(args)) {
        bambu_message_ndr(args)
        message("Running global bambu discovery from rcFiles...")
        bambu_call_analysis(
            analysis_fn,
            bambu_build_args(
                args,
                rc_files,
                annotation_obj,
                discovery = TRUE,
                quant = FALSE
            )
        )
    } else {
        message("Fixed annotation mode: using prepared annotation without bambu discovery.")
        annotation_obj
    }

    chunk_bundles <- bambu_chunk_rc_files(rc_files, inputs$aliases, inputs$sample_df)
    bambu_write_discovery_outputs(
        args$out_dir,
        rc_files,
        discovered_annotations,
        chunk_bundles,
        inputs$sample_df
    )

    invisible(
        list(
            rc_files = rc_files,
            discovered_annotations = discovered_annotations,
            sample_df = inputs$sample_df,
            chunk_bundles = chunk_bundles
        )
    )
}

bambu_run_quant_mode <- function(args, analysis_fn) {
    chunk_bundle <- readRDS(args$chunk_rds)
    discovered_annotations <- readRDS(args$discovered_annotation_rds)
    sample_names <- bambu_chunk_sample_names(chunk_bundle)

    annotation_tx_count <- chunk_bundle$annotation_tx_count
    if (length(annotation_tx_count) != 1L || is.na(annotation_tx_count)) {
        annotation_tx_count <- bambu_annotation_tx_count_for_seqname(
            discovered_annotations,
            chunk_bundle$seqname
        )
    }

    if (!is.na(annotation_tx_count) && annotation_tx_count < 1L) {
        warning(
            sprintf(
                paste(
                    "Skipping bambu quant for chunk '%s' (%s):",
                    "the discovered annotation contains no transcripts on this seqname."
                ),
                chunk_bundle$chunk_id,
                chunk_bundle$seqname
            ),
            call. = FALSE
        )
        se <- bambu_empty_quant_se(
            discovered_annotations = discovered_annotations,
            sample_names = sample_names
        )
    } else {
        se <- bambu_call_analysis(
            analysis_fn,
            bambu_build_args(
                args,
                chunk_bundle$rc_files,
                discovered_annotations,
                discovery = FALSE,
                quant = TRUE
            )
        )
    }

    if (length(sample_names) > 0 && ncol(se) == length(sample_names)) {
        colnames(se) <- sample_names
    }
    S4Vectors::metadata(se)$incompatibleCounts <- bambu_normalise_incompatible_counts(
        S4Vectors::metadata(se)$incompatibleCounts,
        sample_names = colnames(se)
    )

    saveRDS(se, file.path(args$out_dir, "bambu_transcripts.rds"))
    utils::write.csv(
        chunk_bundle$sample_df,
        file.path(args$out_dir, "samples.csv"),
        row.names = FALSE,
        quote = FALSE
    )
    qc_stats <- list(
        chunk_id = chunk_bundle$chunk_id,
        seqname = chunk_bundle$seqname,
        samples = ncol(se),
        total_transcripts_before_filter = nrow(se),
        total_genes_before_filter = length(unique(SummarizedExperiment::rowData(se)$GENEID))
    )
    jsonlite::write_json(
        qc_stats,
        file.path(args$out_dir, "bambu_qc_stats.json"),
        pretty = TRUE,
        auto_unbox = TRUE
    )

    invisible(
        list(
            se = se,
            sample_df = chunk_bundle$sample_df,
            qc_stats = qc_stats
        )
    )
}

bambu_run_collate_mode <- function(args, gene_expression_fn, write_gtf_fn) {
    chunk_dirs <- bambu_resolve_chunk_dirs(args)
    invisible(
        bambu_collate_chunk_outputs(
            chunk_dirs,
            out_dir = args$out_dir,
            transcriptome_mode = args$transcriptome_mode,
            ndr = args$ndr,
            gene_expression_fn = gene_expression_fn,
            write_gtf_fn = write_gtf_fn
        )
    )
}

bambu_run_empty_mode <- function(args, write_gtf_fn) {
    aliases <- workflow_glue_r_parse_csv_list(args$aliases)
    if (length(aliases) < 1) {
        stop("No sample aliases were provided in --aliases.", call. = FALSE)
    }

    invisible(
        bambu_write_empty_outputs(
            sample_aliases = aliases,
            out_dir = args$out_dir,
            transcriptome_mode = args$transcriptome_mode,
            ndr = args$ndr,
            write_gtf_fn = write_gtf_fn
        )
    )
}

bambu_call_analysis <- function(analysis_fn, analysis_args, seed = bambu_default_seed) {
    set.seed(seed)
    do.call(analysis_fn, analysis_args)
}

bambu_chunk_sample_names <- function(chunk_bundle) {
    if (!is.null(chunk_bundle$aliases)) {
        return(as.character(chunk_bundle$aliases))
    }
    if (!is.null(chunk_bundle$sample_df$alias)) {
        return(as.character(chunk_bundle$sample_df$alias))
    }

    rc_files <- chunk_bundle$rc_files
    if (!is.null(names(rc_files))) {
        return(names(rc_files))
    }

    paste0("sample", seq_along(rc_files))
}

bambu_annotation_tx_counts_by_seqname <- function(discovered_annotations) {
    if (methods::is(discovered_annotations, "GenomicRangesList")) {
        unlisted_annotations <- unlist(discovered_annotations, use.names = FALSE)
        if (length(unlisted_annotations) < 1L) {
            return(integer())
        }
        partitioning <- IRanges::PartitioningByEnd(discovered_annotations)
        tx_seqnames <- as.character(
            GenomeInfoDb::seqnames(unlisted_annotations)[IRanges::start(partitioning)]
        )
    } else if (methods::is(discovered_annotations, "GenomicRanges")) {
        tx_seqnames <- as.character(GenomeInfoDb::seqnames(discovered_annotations))
    } else {
        return(integer())
    }

    tx_seqnames <- tx_seqnames[!is.na(tx_seqnames) & nzchar(tx_seqnames)]
    if (length(tx_seqnames) < 1L) {
        return(integer())
    }

    table(tx_seqnames)
}

bambu_annotation_tx_count_for_seqname <- function(discovered_annotations, seqname) {
    annotation_tx_counts <- bambu_annotation_tx_counts_by_seqname(discovered_annotations)
    annotation_tx_count <- unname(annotation_tx_counts[seqname])
    if (length(annotation_tx_count) != 1L || is.na(annotation_tx_count)) {
        return(0L)
    }
    as.integer(annotation_tx_count)
}

bambu_empty_quant_se <- function(discovered_annotations, sample_names) {
    empty_row_ranges <- if (methods::is(discovered_annotations, "GenomicRangesList")) {
        discovered_annotations[0]
    } else if (methods::is(discovered_annotations, "GenomicRanges")) {
        discovered_annotations[0]
    } else {
        GenomicRanges::GRanges()
    }

    empty_matrix <- function() {
        matrix(
            numeric(0),
            nrow = 0,
            ncol = length(sample_names),
            dimnames = list(character(0), sample_names)
        )
    }

    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = empty_matrix(),
            CPM = empty_matrix(),
            fullLengthCounts = empty_matrix(),
            uniqueCounts = empty_matrix()
        ),
        rowRanges = empty_row_ranges,
        colData = S4Vectors::DataFrame(row.names = sample_names),
        metadata = list(
            incompatibleCounts = bambu_empty_incompatible_counts(sample_names),
            warnings = character(0)
        )
    )
}

bambu_empty_gene_se <- function(sample_names) {
    empty_matrix <- matrix(
        numeric(0),
        nrow = 0,
        ncol = length(sample_names),
        dimnames = list(character(0), sample_names)
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = empty_matrix,
            CPM = empty_matrix
        ),
        rowData = S4Vectors::DataFrame(GENEID = character(0)),
        colData = S4Vectors::DataFrame(row.names = sample_names)
    )
}

bambu_write_empty_outputs <- function(
    sample_aliases,
    out_dir,
    transcriptome_mode = "discover",
    ndr = NULL,
    write_gtf_fn = bambu::writeToGTF
) {
    sample_df <- data.frame(alias = sample_aliases, stringsAsFactors = FALSE)
    tx_se <- bambu_empty_quant_se(
        discovered_annotations = GenomicRanges::GRangesList(),
        sample_names = sample_aliases
    )
    gene_se <- bambu_empty_gene_se(sample_aliases)
    qc_stats <- list(
        samples = length(sample_aliases),
        total_transcripts_before_filter = 0L,
        total_genes_before_filter = 0L,
        total_transcripts_after_filter = 0L,
        total_genes_after_filter = 0L,
        transcripts_filtered = 0L,
        chunk_count = 0L,
        empty_output = TRUE
    )
    qc_stats <- bambu_add_library_qc_stats(tx_se, qc_stats)

    args <- list(
        out_dir = out_dir,
        transcriptome_mode = transcriptome_mode,
        ndr = ndr
    )
    bambu_write_outputs(
        tx_se,
        gene_se,
        sample_df,
        args,
        qc_stats,
        write_gtf_fn = write_gtf_fn,
        write_rds = TRUE
    )

    invisible(
        list(
            se = tx_se,
            gene_se = gene_se,
            sample_df = sample_df,
            qc_stats = qc_stats
        )
    )
}

bambu_empty_incompatible_counts <- function(sample_names) {
    cols <- c(
        list(GENEID = character(0)),
        stats::setNames(rep(list(numeric(0)), length(sample_names)), sample_names)
    )
    data.table::as.data.table(cols)
}

bambu_normalise_incompatible_counts <- function(incompatible_counts, sample_names) {
    if (is.null(incompatible_counts)) {
        return(bambu_empty_incompatible_counts(sample_names))
    }

    incompatible_counts <- data.table::as.data.table(incompatible_counts)
    if (!"GENEID" %in% names(incompatible_counts)) {
        return(bambu_empty_incompatible_counts(sample_names))
    }

    value_cols <- setdiff(names(incompatible_counts), c("GENEID", "TXNAME"))
    if (length(sample_names) > 0 &&
        !all(sample_names %in% value_cols) &&
        length(value_cols) == length(sample_names)) {
        data.table::setnames(incompatible_counts, value_cols, sample_names)
        value_cols <- sample_names
    }

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
    incompatible_counts[, keep_cols, with = FALSE]
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
    incompatible_counts <- bambu_normalise_incompatible_counts(
        S4Vectors::metadata(se)$incompatibleCounts,
        sample_names = sample_names
    )
    kept_genes <- unique(SummarizedExperiment::rowData(se)$GENEID)
    incompatible_gene_ids <- incompatible_counts$GENEID

    rows_to_keep <- incompatible_gene_ids %in% kept_genes
    incompatible_counts <- incompatible_counts[rows_to_keep, , drop = FALSE]
    data.table::set(
        incompatible_counts,
        j = "GENEID",
        value = incompatible_gene_ids[rows_to_keep]
    )
    S4Vectors::metadata(se)$incompatibleCounts <- incompatible_counts

    qc_stats$total_transcripts_after_filter <- nrow(se)
    qc_stats$total_genes_after_filter <- length(unique(SummarizedExperiment::rowData(se)$GENEID))

    list(se = se, qc_stats = qc_stats)
}

bambu_matrix_to_df <- function(se_obj, assay_name, id_col, meta_df) {
    assay_df <- as.data.frame(SummarizedExperiment::assays(se_obj)[[assay_name]])
    if (nrow(assay_df) < 1) {
        assay_df <- data.frame(stringsAsFactors = FALSE)
        assay_df[[id_col]] <- character(0)
        for (sample_name in colnames(se_obj)) {
            assay_df[[sample_name]] <- numeric(0)
        }
    } else {
        assay_df[[id_col]] <- rownames(se_obj)
        assay_df <- assay_df[, c(id_col, setdiff(names(assay_df), id_col)), drop = FALSE]
    }

    if (nrow(meta_df) < 1 || nrow(assay_df) < 1) {
        output_cols <- unique(c(names(meta_df), names(assay_df)))
        out_df <- data.frame(stringsAsFactors = FALSE)
        for (col_name in output_cols) {
            if (col_name %in% colnames(se_obj)) {
                out_df[[col_name]] <- numeric(0)
            } else {
                out_df[[col_name]] <- character(0)
            }
        }
        return(out_df[, output_cols, drop = FALSE])
    }

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

bambu_add_library_qc_stats <- function(se, qc_stats = list()) {
    lib_sizes <- colSums(SummarizedExperiment::assays(se)$counts)
    qc_stats$library_sizes <- as.list(lib_sizes)
    qc_stats$min_library_size <- min(lib_sizes)
    qc_stats$max_library_size <- max(lib_sizes)
    qc_stats$median_library_size <- stats::median(lib_sizes)

    if (length(lib_sizes) > 1 && min(lib_sizes) > 0) {
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
    } else if (length(lib_sizes) > 1) {
        qc_stats$library_size_ratio <- NA_real_
    }

    detected_per_sample <- colSums(SummarizedExperiment::assays(se)$counts > 0)
    qc_stats$transcripts_detected_per_sample <- as.list(detected_per_sample)
    qc_stats$median_transcripts_detected <- stats::median(detected_per_sample)
    qc_stats
}

bambu_sum_incompatible_counts <- function(tx_ses) {
    incompatible_counts <- lapply(tx_ses, function(se) {
        bambu_normalise_incompatible_counts(
            S4Vectors::metadata(se)$incompatibleCounts,
            sample_names = colnames(se)
        )
    })
    incompatible_counts <- incompatible_counts[!vapply(incompatible_counts, is.null, logical(1))]
    if (length(incompatible_counts) < 1) {
        return(NULL)
    }

    sample_cols <- setdiff(colnames(incompatible_counts[[1]]), c("GENEID", "TXNAME"))
    gene_ids <- unique(unlist(lapply(incompatible_counts, function(ic) {
        as.character(ic$GENEID)
    })))

    combined <- data.frame(
        GENEID = gene_ids,
        stringsAsFactors = FALSE
    )
    for (sample_col in sample_cols) {
        combined[[sample_col]] <- numeric(length(gene_ids))
    }

    for (ic in incompatible_counts) {
        current_sample_cols <- setdiff(colnames(ic), c("GENEID", "TXNAME"))
        if (!identical(current_sample_cols, sample_cols)) {
            stop(
                "Chunk quantification outputs have mismatched incompatible count columns.",
                call. = FALSE
            )
        }

        idx <- match(as.character(ic$GENEID), combined$GENEID)
        for (sample_col in sample_cols) {
            values <- ic[[sample_col]]
            values[is.na(values)] <- 0
            combined[[sample_col]][idx] <- combined[[sample_col]][idx] + values
        }
    }

    data.table::as.data.table(combined)
}

bambu_unique_row_ranges <- function(tx_ses, tx_names) {
    combined_row_ranges <- do.call(c, lapply(tx_ses, SummarizedExperiment::rowRanges))
    unique_row_ranges <- combined_row_ranges[!duplicated(names(combined_row_ranges))]
    unique_row_ranges[match(tx_names, names(unique_row_ranges))]
}

bambu_sum_chunk_assay <- function(tx_ses, assay_name, tx_names) {
    sample_names <- colnames(tx_ses[[1]])
    combined <- matrix(
        0,
        nrow = length(tx_names),
        ncol = length(sample_names),
        dimnames = list(tx_names, sample_names)
    )

    for (se in tx_ses) {
        assay_mat <- SummarizedExperiment::assay(se, assay_name)
        collapsed <- rowsum(assay_mat, group = rownames(se), reorder = FALSE)
        idx <- match(rownames(collapsed), tx_names)
        combined[idx, ] <- combined[idx, , drop = FALSE] + collapsed
    }

    combined
}

bambu_combine_transcript_chunks <- function(tx_ses) {
    if (length(tx_ses) < 1) {
        stop("No chunk quantification results were provided for collation.", call. = FALSE)
    }
    if (length(tx_ses) == 1) {
        return(tx_ses[[1]])
    }

    assay_names <- SummarizedExperiment::assayNames(tx_ses[[1]])
    sample_names <- colnames(tx_ses[[1]])
    col_data <- SummarizedExperiment::colData(tx_ses[[1]])
    object_metadata <- S4Vectors::metadata(tx_ses[[1]])

    for (se in tx_ses[-1]) {
        if (!identical(SummarizedExperiment::assayNames(se), assay_names)) {
            stop("Chunk quantification outputs have mismatched assay sets.", call. = FALSE)
        }
        if (!identical(colnames(se), sample_names)) {
            stop("Chunk quantification outputs have mismatched sample columns.", call. = FALSE)
        }
    }

    tx_names <- unique(unlist(lapply(tx_ses, rownames)))
    combined_assays <- setNames(
        lapply(assay_names, function(assay_name) {
            bambu_sum_chunk_assay(tx_ses, assay_name, tx_names)
        }),
        assay_names
    )
    if ("counts" %in% assay_names && "CPM" %in% assay_names) {
        counts_mat <- combined_assays[["counts"]]
        combined_assays[["CPM"]] <- t(t(counts_mat) / pmax(colSums(counts_mat), 1)) * 1e6
    }

    combined_row_ranges <- bambu_unique_row_ranges(tx_ses, tx_names)

    object_metadata$incompatibleCounts <- bambu_sum_incompatible_counts(tx_ses)

    SummarizedExperiment::SummarizedExperiment(
        assays = combined_assays,
        rowRanges = combined_row_ranges,
        colData = col_data,
        metadata = object_metadata
    )
}

bambu_collate_chunk_outputs <- function(
    chunk_dirs,
    out_dir,
    transcriptome_mode = "discover",
    ndr = NULL,
    gene_expression_fn = bambu::transcriptToGeneExpression,
    write_gtf_fn = bambu::writeToGTF
) {
    if (length(chunk_dirs) < 1) {
        stop("No chunk quantification directories were provided for collation.", call. = FALSE)
    }

    tx_ses <- lapply(chunk_dirs, function(chunk_dir) {
        readRDS(file.path(chunk_dir, "bambu_transcripts.rds"))
    })
    raw_se <- bambu_combine_transcript_chunks(tx_ses)

    sample_df <- utils::read.csv(
        file.path(chunk_dirs[[1]], "samples.csv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    gene_se <- gene_expression_fn(raw_se)
    filtered <- bambu_filter_transcripts(raw_se)
    se <- filtered$se
    qc_stats <- filtered$qc_stats
    message(
        sprintf(
            "Filtering: keeping %d / %d transcripts",
            qc_stats$total_transcripts_after_filter,
            qc_stats$total_transcripts_before_filter
        )
    )

    keep_gene_ids <- unique(as.character(SummarizedExperiment::rowData(se)$GENEID))
    gene_keep_idx <- rownames(gene_se) %in% keep_gene_ids
    gene_se <- gene_se[gene_keep_idx, ]

    qc_stats$chunk_count <- length(chunk_dirs)
    qc_stats <- bambu_add_library_qc_stats(se, qc_stats)
    colnames(se) <- sample_df$alias
    colnames(gene_se) <- sample_df$alias

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    args <- list(
        out_dir = out_dir,
        transcriptome_mode = transcriptome_mode,
        ndr = ndr
    )
    bambu_write_outputs(
        se,
        gene_se,
        sample_df,
        args,
        qc_stats,
        write_gtf_fn = write_gtf_fn
    )

    invisible(
        list(
            se = se,
            gene_se = gene_se,
            sample_df = sample_df,
            qc_stats = qc_stats
        )
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

bambu_write_outputs <- function(
    se,
    gene_se,
    sample_df,
    args,
    qc_stats,
    write_gtf_fn = bambu::writeToGTF,
    write_rds = TRUE
) {
    row_ranges <- SummarizedExperiment::rowRanges(se)
    if (length(row_ranges) > 0) {
        write_gtf_fn(
            row_ranges,
            file = file.path(args$out_dir, "transcripts.gtf")
        )
    } else {
        dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
        writeLines(
            c(
                "##gff-version 2",
                "# Empty GTF generated by supeRglue bambu"
            ),
            con = file.path(args$out_dir, "transcripts.gtf")
        )
    }

    if (isTRUE(write_rds)) {
        saveRDS(se, file.path(args$out_dir, "bambu_transcripts.rds"))
        saveRDS(gene_se, file.path(args$out_dir, "bambu_genes.rds"))
    }
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
        tx_meta$GENEID <- rep(NA_character_, nrow(tx_meta))
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
    qc_stats$ndr_used <- if (identical(args$transcriptome_mode, "discover")) {
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
        if (identical(args$transcriptome_mode, "discover")) {
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

main_run_bambu <- function(
    args,
    analysis_fn = bambu::bambu,
    prepare_annotations_fn = bambu::prepareAnnotations,
    gene_expression_fn = bambu::transcriptToGeneExpression,
    write_gtf_fn = bambu::writeToGTF,
    bamfile_list_ctor = Rsamtools::BamFileList
) {
    set.seed(bambu_default_seed)
    # bambu's parallel worker code may rely on these being attached for generics
    # such as seqlengths().
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(Rsamtools)
    })

    bambu_validate_args(args)
    dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

    if (identical(args$mode, "discover")) {
        return(invisible(bambu_run_discover_mode(
            args,
            analysis_fn = analysis_fn,
            prepare_annotations_fn = prepare_annotations_fn,
            bamfile_list_ctor = bamfile_list_ctor
        )))
    }

    if (identical(args$mode, "collate")) {
        return(invisible(bambu_run_collate_mode(
            args,
            gene_expression_fn = gene_expression_fn,
            write_gtf_fn = write_gtf_fn
        )))
    }

    if (identical(args$mode, "empty")) {
        return(invisible(bambu_run_empty_mode(
            args,
            write_gtf_fn = write_gtf_fn
        )))
    }

    invisible(bambu_run_quant_mode(
        args,
        analysis_fn = analysis_fn
    ))
}

run_bambu_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
    if (length(argv) >= 1 && !startsWith(argv[[1]], "-")) {
        if (argv[[1]] %in% c("discover", "quant", "collate", "empty")) {
            argv <- c("--mode", argv[[1]], argv[-1])
        }
    }
    parsed <- argparser::parse_args(bambu_arg_parser(), argv = argv)
    args <- workflow_glue_r_normalise_args(parsed, bambu_arg_spec(), raw_argv = argv)
    main_run_bambu(args)
}
