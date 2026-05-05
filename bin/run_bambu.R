#!/usr/bin/env Rscript

# Set seed for reproducibility
set.seed(42)

suppressPackageStartupMessages({
    library(argparser)
    library(bambu)
    library(Rsamtools)
    library(SummarizedExperiment)
    library(jsonlite)
})

parser <- arg_parser("Run bambu transcript discovery and quantification.")
parser <- add_argument(parser, "--bam_dir", help = "Directory containing BAM files.")
parser <- add_argument(parser, "--bam_path", help = "Path to a single BAM file.")
parser <- add_argument(parser, "--sample_alias", help = "Alias to use for a single BAM file.")
parser <- add_argument(parser, "--sample_sheet", help = "Optional sample sheet CSV.")
parser <- add_argument(parser, "--annotation", help = "Reference annotation GTF/GFF.")
parser <- add_argument(parser, "--genome", help = "Reference genome FASTA.")
parser <- add_argument(parser, "--transcriptome_mode", help = "discover or fixed_annotation.", default = "discover")
parser <- add_argument(parser, "--threads", help = "Number of worker threads.", type = "numeric", default = 1)
parser <- add_argument(parser, "--ndr", help = "Optional novel discovery rate.", type = "numeric")
parser <- add_argument(parser, "--out_dir", help = "Output directory.")
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

required_args <- c("annotation", "genome", "out_dir")
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

if (arg_missing(argv$bam_dir) == arg_missing(argv$bam_path)) {
    stop("Provide exactly one of --bam_dir or --bam_path.")
}

dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

sample_df <- NULL
if (!arg_missing(argv$sample_sheet)) {
    sample_df <- read.csv(argv$sample_sheet, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"alias" %in% names(sample_df)) {
        stop("Sample sheet must contain an 'alias' column.")
    }
}

strip_alias <- function(path) {
    name <- basename(path)
    name <- sub("\\.aligned\\.sorted\\.bam$", "", name)
    name <- tools::file_path_sans_ext(name)
    name
}

if (!arg_missing(argv$bam_dir)) {
    bam_paths <- sort(list.files(argv$bam_dir, pattern = "\\.bam$", full.names = TRUE))
    if (length(bam_paths) < 1) {
        stop("No BAM files were found in bam_dir.")
    }
    aliases <- vapply(bam_paths, strip_alias, character(1))
} else {
    bam_paths <- argv$bam_path
    aliases <- if (!arg_missing(argv$sample_alias)) argv$sample_alias else strip_alias(argv$bam_path)
}

if (!is.null(sample_df)) {
    missing_aliases <- setdiff(aliases, sample_df$alias)
    if (length(missing_aliases) > 0) {
        stop(sprintf(
            "Sample sheet is missing alias rows for BAM files: %s",
            paste(missing_aliases, collapse = ", ")
        ))
    }
    sample_df <- sample_df[match(aliases, sample_df$alias), , drop = FALSE]
} else {
    sample_df <- data.frame(alias = aliases, stringsAsFactors = FALSE)
}

annotation_obj <- prepareAnnotations(argv$annotation)
reads <- if (length(bam_paths) == 1) bam_paths else BamFileList(bam_paths, yieldSize = 1000000)

# Handle NDR parameter with validation and documentation
default_ndr <- 0.1
ndr_value <- default_ndr

if (!arg_missing(argv$ndr)) {
    if (argv$ndr < 0 || argv$ndr > 1) {
        stop("NDR (Novel Discovery Rate) must be between 0 and 1")
    }
    ndr_value <- argv$ndr
    message(sprintf("Using user-specified NDR = %.3f", ndr_value))
} else {
    message(sprintf("Using default NDR = %.3f", default_ndr))
}

if (identical(argv$transcriptome_mode, "discover")) {
    message("Novel Discovery Rate (NDR) controls transcript discovery stringency:")
    message("  Lower NDR (e.g., 0.05) = fewer false positive transcripts, may miss real ones")
    message("  Higher NDR (e.g., 0.2) = more sensitive discovery, more false positives")
    message(sprintf("  Current NDR = %.3f balances precision and recall", ndr_value))
}

bambu_args <- list(
    reads = reads,
    annotations = annotation_obj,
    genome = argv$genome,
    ncore = as.integer(argv$threads),
    discovery = identical(argv$transcriptome_mode, "discover")
)

if (identical(argv$transcriptome_mode, "discover")) {
    bambu_args$NDR <- ndr_value
}

message("Running bambu...")
se <- do.call(bambu, bambu_args)
message("Bambu completed successfully")
colnames(se) <- aliases

counts_mat <- assays(se)$counts
full_length_mat <- assays(se)$fullLengthCounts

# Collect QC statistics before filtering
qc_stats <- list()
qc_stats$total_transcripts_before_filter <- nrow(se)
qc_stats$total_genes_before_filter <- length(unique(rowData(se)$GENEID))
qc_stats$samples <- ncol(se)

# Filter low-count transcripts
if (is.null(full_length_mat)) {
    keep_idx <- rowSums(counts_mat) > 0
} else {
    keep_idx <- rowSums(full_length_mat) > 0
}
if (!any(keep_idx)) {
    keep_idx <- rowSums(counts_mat) >= 0
}

qc_stats$transcripts_filtered <- sum(!keep_idx)
message(sprintf("Filtering: keeping %d / %d transcripts", sum(keep_idx), length(keep_idx)))

se <- se[keep_idx, ]

# Library size statistics and warnings
lib_sizes <- colSums(assays(se)$counts)
qc_stats$library_sizes <- as.list(lib_sizes)
qc_stats$min_library_size <- min(lib_sizes)
qc_stats$max_library_size <- max(lib_sizes)
qc_stats$median_library_size <- median(lib_sizes)

if (length(lib_sizes) > 1) {
    lib_size_ratio <- max(lib_sizes) / min(lib_sizes)
    qc_stats$library_size_ratio <- lib_size_ratio

    if (lib_size_ratio > 3) {
        warning(sprintf(
            "Large library size variation detected (%.1fx difference).\n  Min: %d, Max: %d reads.\n  CPM normalization may not be appropriate for such variation.",
            lib_size_ratio, min(lib_sizes), max(lib_sizes)
        ))
        qc_stats$library_size_warning <- sprintf("%.1fx variation (>3x threshold)", lib_size_ratio)
    }
}

# Per-sample detection statistics
qc_stats$transcripts_detected_per_sample <- as.list(colSums(assays(se)$counts > 0))
qc_stats$median_transcripts_detected <- median(colSums(assays(se)$counts > 0))

qc_stats$total_transcripts_after_filter <- nrow(se)
qc_stats$total_genes_after_filter <- length(unique(rowData(se)$GENEID))

row_ranges <- rowRanges(se)
writeToGTF(row_ranges, file = file.path(argv$out_dir, "transcripts.gtf"))

gene_se <- transcriptToGeneExpression(se)
colnames(gene_se) <- aliases

saveRDS(se, file.path(argv$out_dir, "bambu_transcripts.rds"))
saveRDS(gene_se, file.path(argv$out_dir, "bambu_genes.rds"))
write.csv(sample_df, file.path(argv$out_dir, "samples.csv"), row.names = FALSE, quote = FALSE)

tx_meta <- as.data.frame(rowData(se))
if (!"TXNAME" %in% names(tx_meta)) {
    tx_meta$TXNAME <- rownames(se)
}
if (!"GENEID" %in% names(tx_meta)) {
    tx_meta$GENEID <- NA_character_
}
gene_meta <- as.data.frame(rowData(gene_se))
if (!"GENEID" %in% names(gene_meta)) {
    gene_meta$GENEID <- rownames(gene_se)
}

matrix_to_df <- function(se_obj, assay_name, id_col, meta_df) {
    assay_df <- as.data.frame(assays(se_obj)[[assay_name]])
    assay_df[[id_col]] <- rownames(se_obj)
    assay_df <- assay_df[, c(id_col, setdiff(names(assay_df), id_col)), drop = FALSE]
    merge(meta_df, assay_df, by.x = id_col, by.y = id_col, all.y = TRUE, sort = FALSE)
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

extract_gtf_attribute <- function(attr_field, key) {
    match <- regexec(sprintf('%s "([^"]*)";', key), attr_field, perl = TRUE)
    captures <- regmatches(attr_field, match)[[1]]
    if (length(captures) < 2) {
        return(NULL)
    }
    captures[2]
}

normalise_gtf_attribute_value <- function(value) {
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

sanitise_gtf_file <- function(path) {
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
        transcript_id <- normalise_gtf_attribute_value(
            extract_gtf_attribute(attr_field, "transcript_id")
        )
        gene_id <- normalise_gtf_attribute_value(
            extract_gtf_attribute(attr_field, "gene_id")
        )

        if (!is.null(gene_id) && grepl("\\btranscript_id\\b", gene_id)) {
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

tx_meta <- normalise_tsv_df(tx_meta)
gene_meta <- normalise_tsv_df(gene_meta)

write.table(
    tx_meta,
    file = file.path(argv$out_dir, "transcript_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
write.table(
    gene_meta,
    file = file.path(argv$out_dir, "gene_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

tx_counts <- matrix_to_df(se, "counts", "TXNAME", tx_meta)
tx_cpm <- matrix_to_df(se, "CPM", "TXNAME", tx_meta)
gene_counts <- matrix_to_df(gene_se, "counts", "GENEID", gene_meta)
gene_cpm <- matrix_to_df(gene_se, "CPM", "GENEID", gene_meta)

write.table(
    tx_counts,
    file = file.path(argv$out_dir, "transcript_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
write.table(
    tx_cpm,
    file = file.path(argv$out_dir, "transcript_cpm.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
write.table(
    gene_counts,
    file = file.path(argv$out_dir, "gene_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
write.table(
    gene_cpm,
    file = file.path(argv$out_dir, "gene_cpm.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

sanitise_gtf_file(file.path(argv$out_dir, "transcripts.gtf"))

# Write QC statistics as JSON for HTML report
qc_stats$transcriptome_mode <- argv$transcriptome_mode
qc_stats$ndr_used <- if (identical(argv$transcriptome_mode, "discover")) ndr_value else "N/A"
qc_stats$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write_json(
    qc_stats,
    file.path(argv$out_dir, "bambu_qc_stats.json"),
    pretty = TRUE,
    auto_unbox = TRUE
)

format_count <- function(value) {
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

# Write human-readable QC summary
qc_summary <- c(
    "Bambu Quantification QC Summary",
    "================================",
    "",
    sprintf("Timestamp: %s", qc_stats$timestamp),
    sprintf("Mode: %s", argv$transcriptome_mode),
    if (identical(argv$transcriptome_mode, "discover")) sprintf("NDR: %.3f", ndr_value) else NULL,
    "",
    "Sample Statistics:",
    sprintf("  Samples analyzed: %s", format_count(qc_stats$samples)),
    sprintf(
        "  Median library size: %s reads",
        format_count(qc_stats$median_library_size)
    ),
    sprintf(
        "  Library size range: %s - %s reads",
        format_count(qc_stats$min_library_size),
        format_count(qc_stats$max_library_size)
    ),
    if (!is.null(qc_stats$library_size_warning)) sprintf("  WARNING: %s", qc_stats$library_size_warning) else NULL,
    "",
    "Transcript Discovery:",
    sprintf(
        "  Transcripts before filtering: %s",
        format_count(qc_stats$total_transcripts_before_filter)
    ),
    sprintf(
        "  Transcripts after filtering: %s",
        format_count(qc_stats$total_transcripts_after_filter)
    ),
    sprintf(
        "  Transcripts removed: %s",
        format_count(qc_stats$transcripts_filtered)
    ),
    sprintf(
        "  Median transcripts detected per sample: %s",
        format_count(qc_stats$median_transcripts_detected)
    ),
    "",
    "Gene-Level Summary:",
    sprintf(
        "  Unique genes (before filter): %s",
        format_count(qc_stats$total_genes_before_filter)
    ),
    sprintf(
        "  Unique genes (after filter): %s",
        format_count(qc_stats$total_genes_after_filter)
    ),
    ""
)

writeLines(qc_summary, file.path(argv$out_dir, "bambu_qc_summary.txt"))
message("QC statistics written to bambu_qc_stats.json and bambu_qc_summary.txt")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(argv$out_dir, "session_info.txt"))
message("Session info saved for reproducibility")
