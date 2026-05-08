#' These tests cover the validation logic owned by supeRglue bambu before bambu
#' itself is invoked: mutually exclusive BAM inputs, alias derivation, sample
#' sheet alignment, transcriptome mode selection, and NDR handling.
#'
#' NOTE: Annotation/reference preparation is handled by Python
#' (bin/workflow_glue/prepare_annotation_reference.py) with pytest coverage.
#' These tests focus on bambu-specific validation and integration.

# Workflow must accept exactly one of --bam_dir or --bam_path, not both, not neither.
# Fail fast with clear error rather than passing invalid inputs to bambu.
testthat::test_that("exactly one BAM input required", {
    args <- list(
        annotation = "annotation.gtf",
        genome = "genome.fa",
        out_dir = tempfile("bambu-out-"),
        bam_dir = NULL,
        bam_path = NULL,
        transcriptome_mode = "discover",
        ndr = NULL
    )

    testthat::expect_error(
        bambu_validate_args(args),
        "Provide exactly one of --bam_dir or --bam_path"
    )

    args$bam_dir <- tempfile("bam-dir-")
    args$bam_path <- "sample.bam"
    testthat::expect_error(
        bambu_validate_args(args),
        "Provide exactly one of --bam_dir or --bam_path"
    )

    args$bam_path <- NULL
    testthat::expect_silent(bambu_validate_args(args))
})

# transcriptome_mode must be "discover" or "fixed_annotation".
# NDR (Novel Discovery Rate) must be between 0 and 1 when provided.
testthat::test_that("invalid discovery settings rejected", {
    args <- list(
        annotation = "annotation.gtf",
        genome = "genome.fa",
        out_dir = tempfile("bambu-out-"),
        bam_dir = tempfile("bam-dir-"),
        bam_path = NULL,
        transcriptome_mode = "novel",
        ndr = NULL
    )

    testthat::expect_error(
        bambu_validate_args(args),
        "transcriptome_mode must be one of"
    )

    args$transcriptome_mode <- "discover"
    args$ndr <- -0.01
    testthat::expect_error(
        bambu_validate_args(args),
        "NDR .* must be between 0 and 1"
    )

    args$ndr <- 1.01
    testthat::expect_error(
        bambu_validate_args(args),
        "NDR .* must be between 0 and 1"
    )

    args$ndr <- 0
    testthat::expect_silent(bambu_validate_args(args))
    args$ndr <- 1
    testthat::expect_silent(bambu_validate_args(args))
})

# Derive clean sample aliases from BAM filenames by removing workflow-specific suffixes.
testthat::test_that("BAM suffixes stripped from aliases", {
    testthat::expect_equal(
        bambu_strip_alias("/tmp/sample.aligned.sorted.bam"),
        "sample"
    )
    testthat::expect_equal(
        bambu_strip_alias("/tmp/sample.bam"),
        "sample"
    )
})

# Fail fast if --bam_dir contains no BAM files rather than passing empty input to bambu.
testthat::test_that("empty BAM directory rejected", {
    bam_dir <- tempfile("empty-bam-dir-")
    dir.create(bam_dir)
    args <- list(
        bam_dir = bam_dir,
        bam_path = NULL,
        sample_alias = NULL,
        sample_sheet = NULL
    )

    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "No BAM files were found in bam_dir"
    )
})

# Sample aliases must be unique and explicitly specified.
# Derived aliases (from filenames) can collide: sample.bam vs sample.aligned.sorted.bam.
# Sample sheets must have an 'alias' column with unique values that match BAM inputs.
testthat::test_that("unique sample aliases required", {
    bam_dir <- tempfile("ambiguous-bam-dir-")
    dir.create(bam_dir)
    file.create(file.path(bam_dir, "sample.bam"))
    file.create(file.path(bam_dir, "sample.aligned.sorted.bam"))

    args <- list(
        bam_dir = bam_dir,
        bam_path = NULL,
        sample_alias = NULL,
        sample_sheet = NULL
    )
    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "BAM aliases must be unique"
    )

    sample_bam <- tempfile(fileext = ".bam")
    file.create(sample_bam)
    missing_alias_sheet <- tempfile(fileext = ".csv")
    writeLines(
        paste(
            "condition",
            "control",
            sep = "\n"
        ),
        missing_alias_sheet
    )
    args <- list(
        bam_dir = NULL,
        bam_path = sample_bam,
        sample_alias = "sampleA",
        sample_sheet = missing_alias_sheet
    )
    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "Sample sheet must contain an 'alias' column"
    )

    duplicate_alias_sheet <- tempfile(fileext = ".csv")
    writeLines(
        paste(
            "alias,condition",
            "sampleA,control",
            "sampleA,treated",
            sep = "\n"
        ),
        duplicate_alias_sheet
    )
    args$sample_sheet <- duplicate_alias_sheet
    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "Sample sheet aliases must be unique"
    )
})

# Sample sheet rows must align with BAM file order.
# If sheet lists samples in different order than filesystem, reorder sheet to match.
# If sheet is missing aliases found in BAMs, fail.
testthat::test_that("sample sheet reordered to match BAMs", {
    bam_dir <- tempfile("bam-dir-")
    dir.create(bam_dir)
    file.create(file.path(bam_dir, "sampleA.aligned.sorted.bam"))
    file.create(file.path(bam_dir, "sampleB.bam"))

    sample_sheet <- tempfile(fileext = ".csv")
    writeLines(
        paste(
            "alias,condition",
            "sampleB,treated",
            "sampleA,control",
            sep = "\n"
        ),
        sample_sheet
    )

    args <- list(
        bam_dir = bam_dir,
        bam_path = NULL,
        sample_alias = NULL,
        sample_sheet = sample_sheet
    )
    resolved <- bambu_resolve_inputs(
        args,
        bamfile_list_ctor = function(paths, yieldSize) paths
    )

    testthat::expect_equal(resolved$aliases, c("sampleA", "sampleB"))
    testthat::expect_equal(resolved$sample_df$alias, c("sampleA", "sampleB"))

    bad_sheet <- tempfile(fileext = ".csv")
    writeLines(
        paste(
            "alias,condition",
            "sampleA,control",
            sep = "\n"
        ),
        bad_sheet
    )
    args$sample_sheet <- bad_sheet

    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "Sample sheet is missing alias rows"
    )
})

# transcriptome_mode="discover" → bambu(discovery=TRUE, NDR=value)
# transcriptome_mode="fixed_annotation" → bambu(discovery=FALSE, no NDR param)
testthat::test_that("transcriptome mode mapped to bambu args", {
    annotation_obj <- structure(list(annotation = TRUE), class = "mockAnnotation")

    discover_args <- list(
        genome = "genome.fa",
        threads = 3,
        transcriptome_mode = "discover",
        ndr = 0.2
    )
    discover <- bambu_build_args(discover_args, reads = "sample.bam", annotation_obj = annotation_obj)
    testthat::expect_true(discover$discovery)
    testthat::expect_equal(discover$NDR, 0.2)
    testthat::expect_equal(discover$ncore, 3L)

    fixed_args <- list(
        genome = "genome.fa",
        threads = 1,
        transcriptome_mode = "fixed_annotation",
        ndr = NULL
    )
    fixed <- bambu_build_args(fixed_args, reads = "sample.bam", annotation_obj = annotation_obj)
    testthat::expect_false(fixed$discovery)
    testthat::expect_false("NDR" %in% names(fixed))
})

# End-to-end unit test with mocked bambu analysis function.
# Verifies --bam_dir input resolution, sample sheet reordering, and discovery settings.
testthat::test_that("bam_dir input with discovery mode", {
    fixture_dir <- tempfile("bambu-discover-")
    dir.create(fixture_dir)
    bam_dir <- file.path(fixture_dir, "bams")
    dir.create(bam_dir)
    sample_a <- file.path(bam_dir, "sampleA.aligned.sorted.bam")
    sample_b <- file.path(bam_dir, "sampleB.bam")
    file.create(sample_b)
    file.create(sample_a)

    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    writeLines(
        paste(
            "alias,condition",
            "sampleB,treated",
            "sampleA,control",
            sep = "\n"
        ),
        sample_sheet
    )

    captured <- new.env(parent = emptyenv())
    fake_bamfile_list <- function(paths, yieldSize) {
        captured$bamfile_paths <- paths
        captured$yield_size <- yieldSize
        structure(paths, class = "mockBamFileList")
    }
    fake_analysis <- function(reads, annotations, genome, ncore, discovery, NDR) {
        captured$reads <- reads
        captured$annotations <- annotations
        captured$genome <- genome
        captured$ncore <- ncore
        captured$discovery <- discovery
        captured$NDR <- NDR
        make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    }

    argv <- list(
        annotation = "annotation.gtf",
        genome = "genome.fa",
        out_dir = file.path(fixture_dir, "out"),
        bam_dir = bam_dir,
        bam_path = NULL,
        sample_alias = NULL,
        sample_sheet = sample_sheet,
        transcriptome_mode = "discover",
        ndr = 0.25,
        threads = 2
    )

    result <- suppressMessages(main_run_bambu(
        argv,
        analysis_fn = fake_analysis,
        prepare_annotations_fn = function(annotation) {
            captured$annotation_path <- annotation
            structure(list(path = annotation), class = "mockAnnotation")
        },
        gene_expression_fn = function(se) make_test_gene_se(sample_names = colnames(se)),
        write_gtf_fn = function(row_ranges, file) {
            writeLines(
                'chr1\tsim\texon\t1\t50\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";',
                file
            )
        },
        bamfile_list_ctor = fake_bamfile_list
    ))

    testthat::expect_equal(captured$annotation_path, "annotation.gtf")
    testthat::expect_equal(captured$genome, "genome.fa")
    testthat::expect_equal(captured$ncore, 2L)
    testthat::expect_true(captured$discovery)
    testthat::expect_equal(captured$NDR, 0.25)
    testthat::expect_equal(captured$yield_size, 1000000)
    testthat::expect_equal(captured$bamfile_paths, c(sample_a, sample_b))
    testthat::expect_equal(result$sample_df$alias, c("sampleA", "sampleB"))
    testthat::expect_equal(result$sample_df$condition, c("control", "treated"))
    testthat::expect_true(file.exists(file.path(argv$out_dir, "bambu_qc_stats.json")))
})

###
# Output serialization
#
# bambu returns Bioconductor objects with metadata columns that may be list-like
# or range-backed. Check that these are flattened or simplified appropriately
# when written to TSV.

# List-valued metadata columns (e.g., eqClassById) must be collapsed to strings for TSV.
testthat::test_that("list columns flattened for TSV output", {
    df <- data.frame(name = c("a", "b"), stringsAsFactors = FALSE)
    df$list_col <- I(list(c("x", "y"), "z"))

    normalised <- workflow_glue_r_normalise_tsv_df(df)

    testthat::expect_equal(normalised$list_col, c("x;y", "z"))
})

# NCBI annotations may have gene_id="transcript_id" (literal string, not value reference).
# Replace with gene_id=<actual transcript_id value>, but leave gene_id="MYTRANSCRIPT_ID" alone.
testthat::test_that("malformed gene_id values sanitized", {
    gtf_path <- tempfile(fileext = ".gtf")
    writeLines(
        c(
            'chr1\tsim\texon\t1\t50\t.\t+\t.\tgene_id "MYTRANSCRIPT_ID"; transcript_id "tx_keep";',
            'chr1\tsim\texon\t101\t150\t.\t+\t.\tgene_id "transcript_id"; transcript_id "tx_replace";'
        ),
        gtf_path
    )

    testthat::expect_warning(
        bambu_sanitise_gtf_file(gtf_path),
        "Replaced malformed gene_id"
    )

    lines <- readLines(gtf_path, warn = FALSE)
    testthat::expect_match(lines[[1]], 'gene_id "MYTRANSCRIPT_ID";', fixed = TRUE)
    testthat::expect_match(lines[[2]], 'gene_id "tx_replace";', fixed = TRUE)
})

# Verify all expected output files are created with correct structure.
testthat::test_that("bambu outputs written correctly", {
    out_dir <- tempfile("bambu-write-")
    dir.create(out_dir)

    sample_names <- c("sampleA", "sampleB")
    se <- make_test_tx_se(sample_names = sample_names)
    gene_se <- make_test_gene_se(sample_names = sample_names)
    sample_df <- data.frame(alias = sample_names, stringsAsFactors = FALSE)
    argv <- list(out_dir = out_dir, transcriptome_mode = "discover", ndr = 0.15)
    qc_stats <- list(
        samples = 2,
        total_transcripts_before_filter = 4,
        total_genes_before_filter = 2,
        total_transcripts_after_filter = 4,
        total_genes_after_filter = 2,
        transcripts_filtered = 0,
        median_library_size = 400,
        min_library_size = 390,
        max_library_size = 410,
        median_transcripts_detected = 4
    )

    bambu_write_outputs(
        se,
        gene_se,
        sample_df,
        argv,
        qc_stats,
        write_gtf_fn = function(row_ranges, file) {
            writeLines(
                'chr1\tsim\texon\t1\t50\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";',
                file
            )
        }
    )

    testthat::expect_true(file.exists(file.path(out_dir, "transcripts.gtf")))
    testthat::expect_true(file.exists(file.path(out_dir, "transcript_counts.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "gene_counts.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "bambu_qc_stats.json")))

    tx_meta <- utils::read.delim(
        file.path(out_dir, "transcript_metadata.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    tx_counts <- utils::read.delim(
        file.path(out_dir, "transcript_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_true("eqClassById" %in% names(tx_meta))
    testthat::expect_equal(tx_meta$eqClassById[[1]], "1;2")
    testthat::expect_true(all(c("TXNAME", "sampleA", "sampleB") %in% names(tx_counts)))
})

# Transcript filtering edge cases
testthat::test_that("zero full-length count transcripts filtered", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # Set tx1 to have full-length counts, tx2 to have none (will be filtered)
    counts <- SummarizedExperiment::assay(se, "counts")
    counts[1, ] <- c(10, 8)
    counts[2, ] <- c(5, 3)
    counts[3, ] <- c(0, 0)
    counts[4, ] <- c(0, 0)
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    full_length <- matrix(0, nrow = 4, ncol = 2)
    full_length[1, ] <- c(5, 4)
    SummarizedExperiment::assays(se, withDimnames = FALSE)[["fullLengthCounts"]] <- full_length

    result <- bambu_filter_transcripts(se)

    testthat::expect_equal(nrow(result$se), 1)
    testthat::expect_equal(result$qc_stats$transcripts_filtered, 3)
})

testthat::test_that("filters on counts when no fullLengthCounts assay", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # Set some transcripts with counts
    counts <- SummarizedExperiment::assay(se, "counts")
    counts[1, ] <- c(10, 8)
    counts[2, ] <- c(5, 3)
    counts[3, ] <- c(0, 0)
    counts[4, ] <- c(0, 0)
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    # Remove fullLengthCounts assay entirely so it falls back to counts
    assay_list <- SummarizedExperiment::assays(se)
    assay_list[["fullLengthCounts"]] <- NULL
    SummarizedExperiment::assays(se, withDimnames = FALSE) <- assay_list

    result <- bambu_filter_transcripts(se)

    testthat::expect_equal(nrow(result$se), 2)
    testthat::expect_equal(result$qc_stats$transcripts_filtered, 2)
})

testthat::test_that("error when all transcripts filtered", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # All transcripts have zero counts
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- matrix(0, nrow = 4, ncol = 2)
    SummarizedExperiment::assays(se, withDimnames = FALSE)[["fullLengthCounts"]] <- matrix(0, nrow = 4, ncol = 2)

    testthat::expect_error(
        bambu_filter_transcripts(se),
        "All transcripts have zero counts after filtering"
    )
})

testthat::test_that("QC stats match filtered results", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # Set one transcript with counts, rest without
    counts <- matrix(0, nrow = 4, ncol = 2)
    counts[1, ] <- c(10, 8)
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    full_length <- matrix(0, nrow = 4, ncol = 2)
    full_length[1, ] <- c(5, 4)
    SummarizedExperiment::assays(se, withDimnames = FALSE)[["fullLengthCounts"]] <- full_length

    result <- bambu_filter_transcripts(se)

    testthat::expect_equal(result$qc_stats$total_transcripts_after_filter, nrow(result$se))
    testthat::expect_equal(result$qc_stats$transcripts_filtered, 3)
})

###
# CLI integration tests
#
# End-to-end tests with real bambu library (not mocked):
# - Build BAMs from committed fixtures (reference.fa, annotation.gtf, reads.fastq)
# - Run `supeRglue bambu` with both --bam_path and --bam_dir inputs
# - Verify output files exist and contain data for downstream workflow steps

testthat::test_that("CLI single BAM input with fixed annotation", {
    fixture_dir <- tempfile("bambu-cli-")
    dir.create(fixture_dir)

    reference <- workflow_glue_r_fixture("bambu", "reference.fa")
    annotation <- workflow_glue_r_fixture("bambu", "annotation.gtf")
    reads <- workflow_glue_r_fixture("bambu", "reads.fastq")
    sample_sheet <- workflow_glue_r_fixture("bambu", "sample_sheet.csv")
    bam_path <- expect_bam_fixture_built(reference, reads, fixture_dir, alias = "sampleA")
    out_dir <- file.path(fixture_dir, "out")

    result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "--bam_path", bam_path,
            "--sample_alias", "sampleA",
            "--sample_sheet", sample_sheet,
            "--annotation", annotation,
            "--genome", reference,
            "--transcriptome_mode", "fixed_annotation",
            "--threads", "1",
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        result$status,
        0L,
        info = paste(result$output, collapse = "\n")
    )
    testthat::expect_true(file.exists(file.path(out_dir, "transcripts.gtf")))
    testthat::expect_true(file.exists(file.path(out_dir, "transcript_counts.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "gene_counts.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "samples.csv")))

    tx_counts <- utils::read.delim(
        file.path(out_dir, "transcript_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    gene_counts <- utils::read.delim(
        file.path(out_dir, "gene_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_gt(nrow(tx_counts), 0)
    testthat::expect_gt(nrow(gene_counts), 0)
})

testthat::test_that("CLI bam_dir input preserves sample order", {
    fixture_dir <- tempfile("bambu-cli-dir-")
    dir.create(fixture_dir)
    bam_dir <- file.path(fixture_dir, "bams")
    dir.create(bam_dir)

    reference <- workflow_glue_r_fixture("bambu", "reference.fa")
    annotation <- workflow_glue_r_fixture("bambu", "annotation.gtf")
    reads <- workflow_glue_r_fixture("bambu", "reads.fastq")
    expect_bam_fixture_built(reference, reads, bam_dir, alias = "sampleA")
    expect_bam_fixture_built(reference, reads, bam_dir, alias = "sampleB")
    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    writeLines(
        paste(
            "alias",
            "sampleB",
            "sampleA",
            sep = "\n"
        ),
        sample_sheet
    )
    out_dir <- file.path(fixture_dir, "out")

    result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "--bam_dir", bam_dir,
            "--sample_sheet", sample_sheet,
            "--annotation", annotation,
            "--genome", reference,
            "--transcriptome_mode", "fixed_annotation",
            "--threads", "1",
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        result$status,
        0L,
        info = paste(result$output, collapse = "\n")
    )

    samples <- utils::read.csv(
        file.path(out_dir, "samples.csv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    tx_counts <- utils::read.delim(
        file.path(out_dir, "transcript_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_equal(samples$alias, c("sampleA", "sampleB"))
    testthat::expect_true(all(c("sampleA", "sampleB") %in% names(tx_counts)))
    testthat::expect_gt(nrow(tx_counts), 0)
})

# Ensembl GTF annotations include version numbers in IDs (e.g., ENST000001.7).
# Verify bambu doesn't strip versions during prepareAnnotations().
testthat::test_that("Ensembl versioned identifiers preserved", {
    fixture_dir <- tempfile("annotation-ensembl-")
    dir.create(fixture_dir)

    gtf <- file.path(fixture_dir, "ensembl_versions.gtf")
    writeLines(
        c(
            paste(
                "chr1", "Ensembl", "transcript", "1", "100", ".", "+", ".",
                'gene_id "ENSG000001.16"; transcript_id "ENST000001.7"; gene_name "GENEA";',
                sep = "\t"
            ),
            paste(
                "chr1", "Ensembl", "exon", "1", "100", ".", "+", ".",
                'gene_id "ENSG000001.16"; transcript_id "ENST000001.7"; exon_number "1";',
                sep = "\t"
            )
        ),
        gtf
    )

    prepared <- bambu::prepareAnnotations(gtf)
    testthat::expect_s4_class(prepared, "GRangesList")
    testthat::expect_true("ENST000001.7" %in% names(prepared))
})
