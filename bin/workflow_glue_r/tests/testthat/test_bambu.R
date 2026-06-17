#' These tests cover the validation logic owned by supeRglue bambu before bambu
#' itself is invoked: mode-specific inputs, explicit BAM aliases, sample-sheet
#' alignment, transcriptome mode selection, NDR handling, and chunk collation.
#'
#' NOTE: Annotation/reference preparation is handled by Python
#' (bin/workflow_glue/prepare_annotation_reference.py) with pytest coverage.
#' These tests focus on bambu-specific validation and integration.

# Workflow requires an explicit mode plus the appropriate inputs for that mode.
# Fail fast with clear errors rather than passing invalid inputs to bambu.
testthat::test_that("mode-specific bambu inputs required", {
    args <- workflow_glue_r_normalise_args(
        list(mode = "discover", out_dir = tempfile("bambu-out-")),
        bambu_arg_spec()
    )
    testthat::expect_error(
        bambu_validate_args(args),
        "Missing required arguments: --bams, --aliases, --annotation, --genome"
    )

    args$bams <- "sampleA.bam"
    testthat::expect_error(
        bambu_validate_args(args),
        "Missing required arguments: --aliases, --annotation, --genome"
    )

    args$aliases <- "sampleA"
    args$annotation <- "annotation.gtf"
    args$genome <- "genome.fa"
    testthat::expect_silent(bambu_validate_args(args))

    quant_args <- workflow_glue_r_normalise_args(
        list(mode = "quant", out_dir = tempfile("bambu-out-")),
        bambu_arg_spec()
    )
    testthat::expect_error(
        bambu_validate_args(quant_args),
        "Missing required arguments: --chunk_rds, --discovered_annotation_rds, --genome"
    )
    quant_args$chunk_rds <- "chunk.rds"
    quant_args$discovered_annotation_rds <- "annotations.rds"
    quant_args$genome <- "genome.fa"
    testthat::expect_silent(bambu_validate_args(quant_args))

    collate_args <- workflow_glue_r_normalise_args(
        list(mode = "collate", out_dir = tempfile("bambu-out-")),
        bambu_arg_spec()
    )
    testthat::expect_error(
        bambu_validate_args(collate_args),
        "Missing required arguments: --chunk_dirs"
    )
    collate_args$chunk_dirs <- "chunkA,chunkB"
    testthat::expect_silent(bambu_validate_args(collate_args))

    empty_args <- workflow_glue_r_normalise_args(
        list(mode = "empty", out_dir = tempfile("bambu-out-")),
        bambu_arg_spec()
    )
    testthat::expect_error(
        bambu_validate_args(empty_args),
        "Missing required arguments: --aliases"
    )
    empty_args$aliases <- "sampleA,sampleB"
    testthat::expect_silent(bambu_validate_args(empty_args))
})

# transcriptome_mode must be "discover" or "fixed_annotation".
# NDR (Novel Discovery Rate) must be between 0 and 1 when provided.
testthat::test_that("invalid discovery settings rejected", {
    args <- list(
        mode = "discover",
        annotation = "annotation.gtf",
        genome = "genome.fa",
        out_dir = tempfile("bambu-out-"),
        bams = "sampleA.bam",
        aliases = "sampleA",
        transcriptome_mode = "novel",
        ndr = NULL
    )

    testthat::expect_error(
        workflow_glue_r_normalise_args(args, bambu_arg_spec()),
        "transcriptome_mode must be one of"
    )

    args$transcriptome_mode <- "discover"
    args$ndr <- -0.01
    testthat::expect_error(
        workflow_glue_r_normalise_args(args, bambu_arg_spec()),
        "NDR .* must be between 0 and 1"
    )

    args$ndr <- 1.01
    testthat::expect_error(
        workflow_glue_r_normalise_args(args, bambu_arg_spec()),
        "NDR .* must be between 0 and 1"
    )

    args$ndr <- 0
    testthat::expect_silent(workflow_glue_r_normalise_args(args, bambu_arg_spec()))
    args$ndr <- 1
    testthat::expect_silent(workflow_glue_r_normalise_args(args, bambu_arg_spec()))
})

# Fail fast if --bams is empty rather than passing empty input to bambu.
testthat::test_that("empty BAM list rejected", {
    args <- list(
        bams = "",
        aliases = "",
        sample_sheet = NULL
    )

    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "No BAM files were provided in --bams"
    )
})

# Sample aliases must be unique and explicitly specified.
# Sample sheets must have an 'alias' column with unique values that match BAM inputs.
testthat::test_that("unique sample aliases required", {
    args <- list(
        bams = "sampleA.bam,sampleB.bam",
        aliases = "sampleA,sampleA",
        sample_sheet = NULL
    )
    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "BAM aliases must be unique"
    )

    args$aliases <- "sampleA"
    testthat::expect_error(
        bambu_resolve_inputs(args, bamfile_list_ctor = function(paths, yieldSize) paths),
        "Provide one alias per BAM in --bams"
    )
})

# Sample sheet rows must align with BAM file order.
# If sheet lists samples in different order than CLI aliases, reorder sheet to match.
# If sheet is missing aliases found in BAMs, fail.
testthat::test_that("sample sheet reordered to match BAMs", {
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
        bams = "sampleA.aligned.sorted.bam,sampleB.bam",
        aliases = "sampleA,sampleB",
        sample_sheet = sample_sheet
    )
    resolved <- bambu_resolve_inputs(
        args,
        bamfile_list_ctor = function(paths, yieldSize) paths
    )

    testthat::expect_equal(resolved$aliases, c("sampleA", "sampleB"))
    testthat::expect_equal(resolved$sample_df$alias, c("sampleA", "sampleB"))
    testthat::expect_equal(resolved$sample_df$condition, c("control", "treated"))

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

testthat::test_that("numeric alias and sample_id values are preserved as strings", {
    sample_sheet <- tempfile(fileext = ".csv")
    writeLines(
        paste(
            "barcode,sample_id,alias,condition",
            "barcode01,01,01,control",
            "barcode02,02,02,treated",
            sep = "\n"
        ),
        sample_sheet
    )

    args <- list(
        bams = "sample1.bam,sample2.bam",
        aliases = "01,02",
        sample_sheet = sample_sheet
    )
    resolved <- bambu_resolve_inputs(
        args,
        bamfile_list_ctor = function(paths, yieldSize) paths
    )

    testthat::expect_equal(resolved$aliases, c("01", "02"))
    testthat::expect_equal(resolved$sample_df$alias, c("01", "02"))
    testthat::expect_equal(resolved$sample_df$sample_id, c("01", "02"))
    testthat::expect_type(resolved$sample_df$alias, "character")
    testthat::expect_type(resolved$sample_df$sample_id, "character")
})

# Explicit discovery/quant flags are passed through to bambu consistently.
# NDR is only passed during discovery and omitted when automatic selection is wanted.
testthat::test_that("bambu args include requested discovery and quant flags", {
    annotation_obj <- structure(list(annotation = TRUE), class = "mockAnnotation")

    args <- list(
        genome = "genome.fa",
        transcriptome_mode = "discover",
        ndr = 0.2
    )
    discover <- bambu_build_args(
        args,
        reads = "sample.bam",
        annotation_obj = annotation_obj,
        discovery = TRUE,
        quant = FALSE
    )
    testthat::expect_true(discover$discovery)
    testthat::expect_false(discover$quant)
    testthat::expect_equal(discover$NDR, 0.2)
    testthat::expect_equal(discover$ncore, 1L)
    testthat::expect_true(discover$lowMemory)
    testthat::expect_equal(discover$yieldSize, 250000L)

    args$ndr <- NULL
    auto_ndr <- bambu_build_args(
        args,
        reads = "sample.bam",
        annotation_obj = annotation_obj,
        discovery = TRUE,
        quant = FALSE
    )
    testthat::expect_false("NDR" %in% names(auto_ndr))

    quant <- bambu_build_args(
        args,
        reads = "sample.bam",
        annotation_obj = annotation_obj,
        discovery = FALSE,
        quant = TRUE
    )
    testthat::expect_false(quant$discovery)
    testthat::expect_true(quant$quant)
    testthat::expect_false("NDR" %in% names(quant))
})

# Discover mode should write reusable rcFiles, discovered annotations, and chunk bundles.
# This is the scatter source for later per-chromosome quantification.
testthat::test_that("discover mode writes chunked rc outputs", {
    fixture_dir <- tempfile("bambu-discover-mode-")
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
        structure(paths, names = c("sampleA", "sampleB"), class = "mockBamFileList")
    }
    make_rc_sample <- function(alias) {
        rcf <- make_test_tx_se(sample_names = alias)
        S4Vectors::mcols(SummarizedExperiment::rowRanges(rcf))$chr.rc <- c("chr1", "chr1", "chr2", "chr2")
        rcf
    }
    fake_analysis <- function(
        reads,
        annotations,
        genome,
        ncore,
        discovery,
        quant,
        lowMemory,
        yieldSize,
        verbose,
        NDR = NULL
    ) {
        captured$calls <- c(captured$calls, list(list(
            reads = reads,
            annotations = annotations,
            genome = genome,
            ncore = ncore,
            discovery = discovery,
            quant = quant,
            lowMemory = lowMemory,
            yieldSize = yieldSize,
            verbose = verbose,
            NDR = NDR
        )))
        if (!discovery && !quant) {
            return(list(
                sampleA = make_rc_sample("sampleA"),
                sampleB = make_rc_sample("sampleB")
            ))
        }
        structure(list(discovered = TRUE), class = "mockDiscoveredAnnotation")
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "discover",
            annotation = "annotation.gtf",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            bams = paste(c(sample_a, sample_b), collapse = ","),
            aliases = "sampleA,sampleB",
            sample_sheet = sample_sheet,
            transcriptome_mode = "discover",
            ndr = 0.25
        ),
        bambu_arg_spec()
    )

    result <- suppressMessages(main_run_bambu(
        args,
        analysis_fn = fake_analysis,
        prepare_annotations_fn = function(annotation) {
            captured$annotation_path <- annotation
            structure(list(path = annotation), class = "mockAnnotation")
        },
        bamfile_list_ctor = fake_bamfile_list
    ))

    testthat::expect_equal(captured$annotation_path, "annotation.gtf")
    testthat::expect_equal(captured$yield_size, 250000L)
    testthat::expect_equal(captured$bamfile_paths, c(sample_a, sample_b))
    testthat::expect_equal(result$sample_df$alias, c("sampleA", "sampleB"))
    testthat::expect_equal(result$sample_df$condition, c("control", "treated"))
    testthat::expect_length(captured$calls, 2)
    testthat::expect_false(captured$calls[[1]]$discovery)
    testthat::expect_false(captured$calls[[1]]$quant)
    testthat::expect_true(captured$calls[[1]]$lowMemory)
    testthat::expect_equal(captured$calls[[1]]$yieldSize, 250000L)
    testthat::expect_equal(captured$calls[[1]]$ncore, 1L)
    testthat::expect_true(captured$calls[[2]]$discovery)
    testthat::expect_false(captured$calls[[2]]$quant)
    testthat::expect_equal(captured$calls[[2]]$NDR, 0.25)
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_rcfiles.rds")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_discovered_annotations.rds")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "chunk_manifest.tsv")))

    manifest <- utils::read.delim(
        file.path(args$out_dir, "chunk_manifest.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    testthat::expect_equal(sort(manifest$seqname), c("chr1", "chr2"))
    chunk_bundle <- readRDS(manifest$rds_path[[1]])
    testthat::expect_true("annotation_tx_count" %in% names(manifest))
    testthat::expect_true(all(c(
        "chunk_id", "seqname", "aliases", "sample_df", "rc_files", "annotation_tx_count"
    ) %in% names(chunk_bundle)))
    testthat::expect_equal(manifest$annotation_tx_count, c(0L, 0L))
    testthat::expect_equal(chunk_bundle$aliases, c("sampleA", "sampleB"))
})

testthat::test_that("cw-7338 chunking accepts path-backed rc outputs", {
    fixture_dir <- tempfile("bambu-path-backed-rc-")
    dir.create(fixture_dir)

    make_rc_sample <- function(alias) {
        rcf <- make_test_tx_se(sample_names = alias)
        S4Vectors::mcols(SummarizedExperiment::rowRanges(rcf))$chr.rc <- c("chr1", "chr1", "chr2", "chr2")
        rcf
    }

    rc_paths <- file.path(fixture_dir, paste0(c("sampleA", "sampleB"), "_readClassSe.rds"))
    saveRDS(make_rc_sample("sampleA"), rc_paths[[1]])
    saveRDS(make_rc_sample("sampleB"), rc_paths[[2]])
    names(rc_paths) <- c("sampleA", "sampleB")

    chunk_bundles <- bambu_chunk_rc_files(
        as.list(rc_paths),
        aliases = c("sampleA", "sampleB"),
        sample_df = data.frame(alias = c("sampleA", "sampleB"), stringsAsFactors = FALSE)
    )

    testthat::expect_equal(unname(sort(vapply(chunk_bundles, `[[`, character(1), "seqname"))), c("chr1", "chr2"))
    testthat::expect_equal(chunk_bundles[[1]]$aliases, c("sampleA", "sampleB"))
    testthat::expect_s4_class(chunk_bundles[[1]]$rc_files$sampleA, "RangedSummarizedExperiment")
    testthat::expect_s4_class(chunk_bundles[[1]]$rc_files$sampleB, "RangedSummarizedExperiment")
})

testthat::test_that("annotation tx counts are computed once per seqname", {
    discovered_annotations <- GenomicRanges::GRangesList(
        tx1 = GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(c(1, 10), width = 5)),
        tx2 = GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(20, width = 5)),
        tx3 = GenomicRanges::GRanges(seqnames = "chr2", ranges = IRanges::IRanges(c(30, 40), width = 5))
    )

    counts <- bambu_annotation_tx_counts_by_seqname(discovered_annotations)

    testthat::expect_equal(as.integer(counts[c("chr1", "chr2")]), c(2L, 1L))
    testthat::expect_equal(bambu_annotation_tx_count_for_seqname(discovered_annotations, "chr1"), 2L)
    testthat::expect_equal(bambu_annotation_tx_count_for_seqname(discovered_annotations, "chr3"), 0L)
})

# Quant mode should consume one chunk bundle and write only raw chunk quantification
# artifacts. Filtering and gene aggregation happen during collate.
testthat::test_that("quant mode writes chunk quantification outputs", {
    fixture_dir <- tempfile("bambu-quant-mode-")
    dir.create(fixture_dir)

    make_rc_sample <- function(alias) {
        rcf <- make_test_tx_se(sample_names = alias)
        S4Vectors::mcols(SummarizedExperiment::rowRanges(rcf))$chr.rc <- c("chr1", "chr1", "chr1", "chr1")
        rcf
    }

    chunk_bundle <- list(
        chunk_id = "chr1",
        seqname = "chr1",
        annotation_tx_count = 1L,
        aliases = c("sampleA", "sampleB"),
        sample_df = data.frame(alias = c("sampleA", "sampleB"), stringsAsFactors = FALSE),
        rc_files = list(sampleA = make_rc_sample("sampleA"), sampleB = make_rc_sample("sampleB"))
    )
    chunk_rds <- file.path(fixture_dir, "chr1.rds")
    saveRDS(chunk_bundle, chunk_rds)
    discovered_annotation_rds <- file.path(fixture_dir, "annotations.rds")
    saveRDS(structure(list(discovered = TRUE), class = "mockDiscoveredAnnotation"), discovered_annotation_rds)

    captured <- new.env(parent = emptyenv())
    fake_analysis <- function(
        reads,
        annotations,
        genome,
        ncore,
        discovery,
        quant,
        lowMemory,
        yieldSize,
        verbose,
        ...
    ) {
        captured$reads <- reads
        captured$annotations <- annotations
        captured$genome <- genome
        captured$ncore <- ncore
        captured$discovery <- discovery
        captured$quant <- quant
        captured$lowMemory <- lowMemory
        captured$yieldSize <- yieldSize
        captured$verbose <- verbose
        make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "quant",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            chunk_rds = chunk_rds,
            discovered_annotation_rds = discovered_annotation_rds,
            transcriptome_mode = "discover",
            ndr = NULL
        ),
        bambu_arg_spec()
    )

    result <- suppressMessages(main_run_bambu(
        args,
        analysis_fn = fake_analysis
    ))

    testthat::expect_equal(captured$genome, "genome.fa")
    testthat::expect_false(captured$discovery)
    testthat::expect_true(captured$quant)
    testthat::expect_true(captured$lowMemory)
    testthat::expect_equal(captured$yieldSize, 250000L)
    testthat::expect_equal(result$sample_df$alias, c("sampleA", "sampleB"))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_transcripts.rds")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "samples.csv")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_qc_stats.json")))
    qc <- jsonlite::read_json(file.path(args$out_dir, "bambu_qc_stats.json"), simplifyVector = TRUE)
    testthat::expect_equal(qc$chunk_id, "chr1")
    testthat::expect_equal(qc$seqname, "chr1")
})

testthat::test_that("quant mode skips chunks with no discovered annotations on the seqname", {
    fixture_dir <- tempfile("bambu-quant-empty-chunk-")
    dir.create(fixture_dir)

    make_rc_sample <- function(alias) {
        rcf <- make_test_tx_se(sample_names = alias)
        S4Vectors::mcols(SummarizedExperiment::rowRanges(rcf))$chr.rc <- rep(
            "chr3_GL000221v1_random",
            nrow(rcf)
        )
        rcf
    }

    chunk_bundle <- list(
        chunk_id = "chr3_GL000221v1_random",
        seqname = "chr3_GL000221v1_random",
        aliases = c("sampleA", "sampleB"),
        sample_df = data.frame(alias = c("sampleA", "sampleB"), stringsAsFactors = FALSE),
        rc_files = list(sampleA = make_rc_sample("sampleA"), sampleB = make_rc_sample("sampleB"))
    )
    chunk_rds <- file.path(fixture_dir, "chr3_GL000221v1_random.rds")
    saveRDS(chunk_bundle, chunk_rds)

    discovered_annotation_rds <- file.path(fixture_dir, "annotations.rds")
    saveRDS(make_test_bambu_row_ranges(fixture_dir), discovered_annotation_rds)

    analysis_called <- FALSE
    fake_analysis <- function(...) {
        analysis_called <<- TRUE
        stop("analysis_fn should not be called for chunks with no matching annotations")
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "quant",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            chunk_rds = chunk_rds,
            discovered_annotation_rds = discovered_annotation_rds,
            transcriptome_mode = "discover",
            ndr = NULL
        ),
        bambu_arg_spec()
    )

    result <- testthat::expect_warning(
        suppressMessages(main_run_bambu(
            args,
            analysis_fn = fake_analysis
        )) ,
        "contains no transcripts on this seqname"
    )

    testthat::expect_false(analysis_called)
    testthat::expect_equal(result$sample_df$alias, c("sampleA", "sampleB"))

    se <- readRDS(file.path(args$out_dir, "bambu_transcripts.rds"))
    testthat::expect_equal(nrow(se), 0)
    testthat::expect_equal(colnames(se), c("sampleA", "sampleB"))
    testthat::expect_equal(
        SummarizedExperiment::assayNames(se),
        c("counts", "CPM", "fullLengthCounts", "uniqueCounts")
    )
    testthat::expect_equal(
        colnames(S4Vectors::metadata(se)$incompatibleCounts),
        c("GENEID", "sampleA", "sampleB")
    )

    qc <- jsonlite::read_json(file.path(args$out_dir, "bambu_qc_stats.json"), simplifyVector = TRUE)
    testthat::expect_equal(qc$chunk_id, "chr3_GL000221v1_random")
    testthat::expect_equal(qc$seqname, "chr3_GL000221v1_random")
    testthat::expect_equal(qc$total_transcripts_before_filter, 0)
    testthat::expect_equal(qc$total_genes_before_filter, 0)
})

testthat::test_that("quant mode catches known uniqueStartLengthQuery/min-Inf edge case and writes empty outputs", {
    fixture_dir <- tempfile("bambu-quant-known-edge-")
    dir.create(fixture_dir)

    chunk_bundle <- list(
        chunk_id = "chr1",
        seqname = "chr1",
        aliases = c("sampleA"),
        sample_df = data.frame(alias = c("sampleA"), stringsAsFactors = FALSE),
        rc_files = list(sampleA = make_test_tx_se(sample_names = "sampleA")),
        annotation_tx_count = 1L
    )
    chunk_rds <- file.path(fixture_dir, "chr1.rds")
    saveRDS(chunk_bundle, chunk_rds)

    discovered_annotation_rds <- file.path(fixture_dir, "annotations.rds")
    saveRDS(make_test_bambu_row_ranges(fixture_dir), discovered_annotation_rds)

    analysis_called <- FALSE
    fake_analysis <- function(...) {
        analysis_called <<- TRUE
        stop(
            paste(
                "Error in filter(., (uniqueStartLengthQuery <= primarySecondaryDistStartEnd & : In argument: `==...`.",
                "Caused by warning in `min()`: no non-missing arguments to min; returning Inf"
            ),
            call. = FALSE
        )
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "quant",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            chunk_rds = chunk_rds,
            discovered_annotation_rds = discovered_annotation_rds,
            transcriptome_mode = "discover",
            ndr = NULL,
            threads = 2
        ),
        bambu_arg_spec()
    )

    testthat::expect_warning(
        suppressMessages(main_run_bambu(args, analysis_fn = fake_analysis)),
        "known bambu chunk edge case"
    )

    testthat::expect_true(analysis_called)

    se <- readRDS(file.path(args$out_dir, "bambu_transcripts.rds"))
    testthat::expect_equal(nrow(se), 0)
    testthat::expect_equal(colnames(se), c("sampleA"))
})

testthat::test_that("quant mode catches known eqClassById incompatible-type edge case and writes empty outputs", {
    fixture_dir <- tempfile("bambu-quant-known-eqclass-edge-")
    dir.create(fixture_dir)

    chunk_bundle <- list(
        chunk_id = "chr2",
        seqname = "chr2",
        aliases = c("sampleA"),
        sample_df = data.frame(alias = c("sampleA"), stringsAsFactors = FALSE),
        rc_files = list(sampleA = make_test_tx_se(sample_names = "sampleA")),
        annotation_tx_count = 1L
    )
    chunk_rds <- file.path(fixture_dir, "chr2.rds")
    saveRDS(chunk_bundle, chunk_rds)

    discovered_annotation_rds <- file.path(fixture_dir, "annotations.rds")
    saveRDS(make_test_bambu_row_ranges(fixture_dir), discovered_annotation_rds)

    analysis_called <- FALSE
    fake_analysis <- function(...) {
        analysis_called <<- TRUE
        stop(
            "Can't join `x$eqClassById` with `y$eqClassById` due to incompatible types.",
            call. = FALSE
        )
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "quant",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            chunk_rds = chunk_rds,
            discovered_annotation_rds = discovered_annotation_rds,
            transcriptome_mode = "discover",
            ndr = NULL,
            threads = 2
        ),
        bambu_arg_spec()
    )

    testthat::expect_warning(
        suppressMessages(main_run_bambu(args, analysis_fn = fake_analysis)),
        "known bambu chunk edge case"
    )

    testthat::expect_true(analysis_called)

    se <- readRDS(file.path(args$out_dir, "bambu_transcripts.rds"))
    testthat::expect_equal(nrow(se), 0)
    testthat::expect_equal(colnames(se), c("sampleA"))
})

testthat::test_that("quant mode catches known txid unspecified edge case and writes empty outputs", {
    fixture_dir <- tempfile("bambu-quant-known-txid-edge-")
    dir.create(fixture_dir)

    chunk_bundle <- list(
        chunk_id = "chr2",
        seqname = "chr2",
        aliases = c("sampleA"),
        sample_df = data.frame(alias = c("sampleA"), stringsAsFactors = FALSE),
        rc_files = list(sampleA = make_test_tx_se(sample_names = "sampleA")),
        annotation_tx_count = 1L
    )
    chunk_rds <- file.path(fixture_dir, "chr2.rds")
    saveRDS(chunk_bundle, chunk_rds)

    discovered_annotation_rds <- file.path(fixture_dir, "annotations.rds")
    saveRDS(make_test_bambu_row_ranges(fixture_dir), discovered_annotation_rds)

    analysis_called <- FALSE
    fake_analysis <- function(...) {
        analysis_called <<- TRUE
        stop(
            paste0(
                "Can't convert `x$txid` <vctrs_unspecified> ",
                "to match type of `txid` <integer>."
            ),
            call. = FALSE
        )
    }

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "quant",
            genome = "genome.fa",
            out_dir = file.path(fixture_dir, "out"),
            chunk_rds = chunk_rds,
            discovered_annotation_rds = discovered_annotation_rds,
            transcriptome_mode = "discover",
            ndr = NULL,
            threads = 2
        ),
        bambu_arg_spec()
    )

    testthat::expect_warning(
        suppressMessages(main_run_bambu(args, analysis_fn = fake_analysis)),
        "known bambu chunk edge case"
    )

    testthat::expect_true(analysis_called)

    se <- readRDS(file.path(args$out_dir, "bambu_transcripts.rds"))
    testthat::expect_equal(nrow(se), 0)
    testthat::expect_equal(colnames(se), c("sampleA"))
})

testthat::test_that("empty mode writes valid empty outputs including bambu rds files", {
    fixture_dir <- tempfile("bambu-empty-mode-")
    dir.create(fixture_dir)

    args <- workflow_glue_r_normalise_args(
        list(
            mode = "empty",
            aliases = "sampleA,sampleB",
            out_dir = file.path(fixture_dir, "out"),
            transcriptome_mode = "fixed_annotation"
        ),
        bambu_arg_spec()
    )

    result <- suppressMessages(main_run_bambu(args))

    testthat::expect_equal(result$sample_df$alias, c("sampleA", "sampleB"))
    testthat::expect_true(file.exists(file.path(args$out_dir, "transcripts.gtf")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "transcript_counts.tsv")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "gene_counts.tsv")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "samples.csv")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_qc_stats.json")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_transcripts.rds")))
    testthat::expect_true(file.exists(file.path(args$out_dir, "bambu_genes.rds")))

    tx_counts <- utils::read.delim(
        file.path(args$out_dir, "transcript_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    gene_counts <- utils::read.delim(
        file.path(args$out_dir, "gene_counts.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    qc <- jsonlite::read_json(file.path(args$out_dir, "bambu_qc_stats.json"), simplifyVector = TRUE)
    gtf_lines <- readLines(file.path(args$out_dir, "transcripts.gtf"), warn = FALSE)
    tx_rds <- readRDS(file.path(args$out_dir, "bambu_transcripts.rds"))
    gene_rds <- readRDS(file.path(args$out_dir, "bambu_genes.rds"))

    testthat::expect_equal(nrow(tx_counts), 0)
    testthat::expect_equal(nrow(gene_counts), 0)
    testthat::expect_true(all(c("sampleA", "sampleB") %in% names(tx_counts)))
    testthat::expect_true(all(c("sampleA", "sampleB") %in% names(gene_counts)))
    testthat::expect_equal(gtf_lines[[1]], "##gff-version 2")
    testthat::expect_true(any(grepl("^#", gtf_lines)))
    testthat::expect_s4_class(tx_rds, "RangedSummarizedExperiment")
    testthat::expect_s4_class(gene_rds, "SummarizedExperiment")
    testthat::expect_equal(nrow(tx_rds), 0)
    testthat::expect_equal(nrow(gene_rds), 0)
    testthat::expect_true(isTRUE(qc$empty_output))
    testthat::expect_equal(qc$chunk_count, 0)
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

testthat::test_that("annotation name maps are extracted from GTF", {
    fixture_dir <- tempfile("gene-name-map-")
    dir.create(fixture_dir)

    gtf <- file.path(fixture_dir, "annotation.gtf")
    writeLines(
        c(
            paste(
                "chr1", "sim", "transcript", "1", "100", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1"; gene_name "GENEA"; transcript_name "TXA";',
                sep = "\t"
            ),
            paste(
                "chr1", "sim", "exon", "1", "100", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1"; gene_name "GENEA"; transcript_name "TXA";',
                sep = "\t"
            ),
            paste(
                "chr1", "sim", "transcript", "201", "300", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx2"; gene_name "GENEB"; transcript_name "TXB";',
                sep = "\t"
            )
        ),
        gtf
    )

    maps <- workflow_glue_r_annotation_name_maps(gtf)
    gene_name_map <- maps$gene
    transcript_name_map <- maps$transcript

    testthat::expect_equal(names(gene_name_map), c("GENEID", "gene_name"))
    testthat::expect_equal(nrow(gene_name_map), 2)
    testthat::expect_equal(gene_name_map$gene_name[match("gene1", gene_name_map$GENEID)], "GENEA")
    testthat::expect_equal(gene_name_map$gene_name[match("gene2", gene_name_map$GENEID)], "GENEB")

    testthat::expect_equal(names(transcript_name_map), c("TXNAME", "transcript_name"))
    testthat::expect_equal(nrow(transcript_name_map), 2)
    testthat::expect_equal(
        transcript_name_map$transcript_name[match("tx1", transcript_name_map$TXNAME)],
        "TXA"
    )
    testthat::expect_equal(
        transcript_name_map$transcript_name[match("tx2", transcript_name_map$TXNAME)],
        "TXB"
    )
})

# Verify all expected output files are created with correct structure.
testthat::test_that("bambu outputs written correctly", {
    out_dir <- tempfile("bambu-write-")
    dir.create(out_dir)
    fixture_dir <- tempfile("bambu-write-fixture-")
    dir.create(fixture_dir)

    sample_names <- c("sampleA", "sampleB")
    base_se <- make_test_tx_se(sample_names = sample_names)
    row_ranges <- make_test_bambu_row_ranges(out_dir)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(base_se),
        rowRanges = row_ranges
    )
    gene_se <- make_test_gene_se(sample_names = sample_names)
    sample_df <- data.frame(alias = sample_names, stringsAsFactors = FALSE)
    annotation <- file.path(fixture_dir, "gene_names.gtf")
    writeLines(
        c(
            paste(
                "chr1", "test", "transcript", "1", "50", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1"; gene_name "GENEA"; transcript_name "TXA";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "201", "250", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx3"; gene_name "GENEB"; transcript_name "TXC";',
                sep = "\t"
            )
        ),
        annotation
    )

    args <- list(
        out_dir = out_dir,
        transcriptome_mode = "discover",
        ndr = 0.15,
        annotation = annotation
    )
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
        args,
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
    testthat::expect_true("gene_name" %in% names(tx_meta))
    testthat::expect_true("transcript_name" %in% names(tx_meta))
    testthat::expect_equal(
        unique(stats::na.omit(tx_meta$gene_name[tx_meta$GENEID == "gene1"])),
        "GENEA"
    )
    testthat::expect_equal(
        tx_meta$transcript_name[match("tx1", tx_meta$TXNAME)],
        "TXA"
    )
    testthat::expect_true(all(c("TXNAME", "sampleA", "sampleB") %in% names(tx_counts)))

    tx_rds <- readRDS(file.path(out_dir, "bambu_transcripts.rds"))
    gene_rds <- readRDS(file.path(out_dir, "bambu_genes.rds"))
    tx_rds_meta <- as.data.frame(SummarizedExperiment::rowData(tx_rds))
    gene_rds_meta <- as.data.frame(SummarizedExperiment::rowData(gene_rds))
    testthat::expect_true("gene_name" %in% names(tx_rds_meta))
    testthat::expect_true("transcript_name" %in% names(tx_rds_meta))
    testthat::expect_true("gene_name" %in% names(gene_rds_meta))
    testthat::expect_equal(
        unique(stats::na.omit(tx_rds_meta$gene_name[tx_rds_meta$GENEID == "gene2"])),
        "GENEB"
    )
    testthat::expect_equal(
        tx_rds_meta$transcript_name[match("tx3", tx_rds_meta$TXNAME)],
        "TXC"
    )
    testthat::expect_equal(
        gene_rds_meta$gene_name[match("gene1", gene_rds_meta$GENEID)],
        "GENEA"
    )
})

testthat::test_that("collate combines multiple chunk quantification outputs", {
    fixture_dir <- tempfile("bambu-collate-multi-")
    dir.create(fixture_dir)

    sample_df <- data.frame(alias = c("sampleA", "sampleB"), stringsAsFactors = FALSE)
    base_se <- make_test_tx_se(sample_names = sample_df$alias)
    tx_se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(base_se),
        rowRanges = make_test_bambu_row_ranges(fixture_dir)
    )

    chunk_dirs <- c(file.path(fixture_dir, "chunk1"), file.path(fixture_dir, "chunk2"))
    dir.create(chunk_dirs[[1]])
    dir.create(chunk_dirs[[2]])
    saveRDS(tx_se[1:2, ], file.path(chunk_dirs[[1]], "bambu_transcripts.rds"))
    saveRDS(tx_se[3:4, ], file.path(chunk_dirs[[2]], "bambu_transcripts.rds"))
    utils::write.csv(sample_df, file.path(chunk_dirs[[1]], "samples.csv"), row.names = FALSE, quote = FALSE)
    utils::write.csv(sample_df, file.path(chunk_dirs[[2]], "samples.csv"), row.names = FALSE, quote = FALSE)

    mock_gene_expression <- function(se) {
        counts <- rowsum(
            SummarizedExperiment::assay(se, "counts"),
            group = as.character(SummarizedExperiment::rowData(se)$GENEID),
            reorder = FALSE
        )
        cpm <- t(t(counts) / colSums(counts)) * 1e6
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = counts, CPM = cpm),
            rowData = S4Vectors::DataFrame(GENEID = rownames(counts))
        )
    }

    out_dir <- file.path(fixture_dir, "collated")
    suppressWarnings(suppressMessages(
        bambu_collate_chunk_outputs(
            chunk_dirs,
            out_dir = out_dir,
            transcriptome_mode = "fixed_annotation",
            gene_expression_fn = mock_gene_expression,
            write_gtf_fn = function(row_ranges, file) {
                writeLines(
                    'chr1\tsim\texon\t1\t50\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";',
                    file
                )
            }
        )
    ))

    collated_tx <- readRDS(file.path(out_dir, "bambu_transcripts.rds"))
    collated_gene <- readRDS(file.path(out_dir, "bambu_genes.rds"))
    testthat::expect_equal(nrow(collated_tx), 4)
    testthat::expect_equal(nrow(collated_gene), 2)
    testthat::expect_equal(colnames(collated_tx), sample_df$alias)
    testthat::expect_true(file.exists(file.path(out_dir, "transcript_counts.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "gene_counts.tsv")))
})

testthat::test_that("chunk combiner sums duplicate transcript rows across chunks", {
    sample_df <- data.frame(alias = "sampleA", stringsAsFactors = FALSE)
    base_se <- make_test_tx_se(sample_names = sample_df$alias)
    fixture_dir <- tempfile("bambu-combine-")
    dir.create(fixture_dir)
    tx_se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(base_se),
        rowRanges = make_test_bambu_row_ranges(fixture_dir)
    )

    chunk1 <- tx_se
    chunk2 <- tx_se
    counts1 <- SummarizedExperiment::assay(chunk1, "counts")
    counts2 <- SummarizedExperiment::assay(chunk2, "counts")
    counts1[3:4, ] <- 0
    counts2[1:2, ] <- 0
    SummarizedExperiment::assay(chunk1, "counts", withDimnames = FALSE) <- counts1
    SummarizedExperiment::assay(chunk2, "counts", withDimnames = FALSE) <- counts2

    cpm1 <- t(t(counts1) / pmax(colSums(counts1), 1)) * 1e6
    cpm2 <- t(t(counts2) / pmax(colSums(counts2), 1)) * 1e6
    SummarizedExperiment::assay(chunk1, "CPM", withDimnames = FALSE) <- cpm1
    SummarizedExperiment::assay(chunk2, "CPM", withDimnames = FALSE) <- cpm2

    S4Vectors::metadata(chunk1)$incompatibleCounts <- data.frame(
        GENEID = c("gene1", "gene2"),
        `01` = c(10, 0),
        stringsAsFactors = FALSE
    )
    S4Vectors::metadata(chunk2)$incompatibleCounts <- data.frame(
        GENEID = c("gene1", "gene2"),
        `01` = c(0, 20),
        stringsAsFactors = FALSE
    )

    combined <- bambu_combine_transcript_chunks(list(chunk1, chunk2))

    testthat::expect_equal(nrow(combined), 4)
    testthat::expect_equal(rownames(combined), rownames(tx_se))
    testthat::expect_equal(
        SummarizedExperiment::assay(combined, "counts"),
        SummarizedExperiment::assay(tx_se, "counts")
    )
    expected_cpm <- t(t(SummarizedExperiment::assay(tx_se, "counts")) /
        colSums(SummarizedExperiment::assay(tx_se, "counts"))) * 1e6
    testthat::expect_equal(SummarizedExperiment::assay(combined, "CPM"), expected_cpm)

    incompatible <- S4Vectors::metadata(combined)$incompatibleCounts
    testthat::expect_equal(names(incompatible), c("GENEID", "sampleA"))
    testthat::expect_equal(incompatible$sampleA, c(10, 20))
})

testthat::test_that("streaming chunk directory combiner matches in-memory combiner", {
    sample_df <- data.frame(alias = "sampleA", stringsAsFactors = FALSE)
    base_se <- make_test_tx_se(sample_names = sample_df$alias)
    fixture_dir <- tempfile("bambu-combine-dirs-")
    dir.create(fixture_dir)
    tx_se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(base_se),
        rowRanges = make_test_bambu_row_ranges(fixture_dir)
    )

    chunk1 <- tx_se
    chunk2 <- tx_se
    counts1 <- SummarizedExperiment::assay(chunk1, "counts")
    counts2 <- SummarizedExperiment::assay(chunk2, "counts")
    counts1[3:4, ] <- 0
    counts2[1:2, ] <- 0
    SummarizedExperiment::assay(chunk1, "counts", withDimnames = FALSE) <- counts1
    SummarizedExperiment::assay(chunk2, "counts", withDimnames = FALSE) <- counts2
    SummarizedExperiment::assay(chunk1, "CPM", withDimnames = FALSE) <-
        t(t(counts1) / pmax(colSums(counts1), 1)) * 1e6
    SummarizedExperiment::assay(chunk2, "CPM", withDimnames = FALSE) <-
        t(t(counts2) / pmax(colSums(counts2), 1)) * 1e6

    S4Vectors::metadata(chunk1)$incompatibleCounts <- data.frame(
        GENEID = c("gene1", "gene2"),
        `01` = c(10, 0),
        stringsAsFactors = FALSE
    )
    S4Vectors::metadata(chunk2)$incompatibleCounts <- data.frame(
        GENEID = c("gene1", "gene2"),
        `01` = c(0, 20),
        stringsAsFactors = FALSE
    )

    chunk_dirs <- file.path(fixture_dir, c("chunk1", "chunk2"))
    dir.create(chunk_dirs[[1]])
    dir.create(chunk_dirs[[2]])
    saveRDS(chunk1, file.path(chunk_dirs[[1]], "bambu_transcripts.rds"))
    saveRDS(chunk2, file.path(chunk_dirs[[2]], "bambu_transcripts.rds"))

    in_memory <- bambu_combine_transcript_chunks(list(chunk1, chunk2))
    streaming <- bambu_combine_transcript_chunk_dirs(chunk_dirs)

    testthat::expect_equal(rownames(streaming), rownames(in_memory))
    testthat::expect_equal(SummarizedExperiment::rowData(streaming), SummarizedExperiment::rowData(in_memory))
    testthat::expect_equal(
        SummarizedExperiment::assay(streaming, "counts"),
        SummarizedExperiment::assay(in_memory, "counts")
    )
    testthat::expect_equal(
        SummarizedExperiment::assay(streaming, "CPM"),
        SummarizedExperiment::assay(in_memory, "CPM")
    )
    testthat::expect_equal(
        S4Vectors::metadata(streaming)$incompatibleCounts,
        S4Vectors::metadata(in_memory)$incompatibleCounts
    )
})

# Transcript filtering edge cases
testthat::test_that("filters on counts when fullLengthCounts disagrees", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # tx2 has counts but zero full-length counts: must be kept.
    # tx3 has zero counts but non-zero full-length counts: must be filtered.
    counts <- SummarizedExperiment::assay(se, "counts")
    counts[1, ] <- c(10, 8)
    counts[2, ] <- c(5, 3)
    counts[3, ] <- c(0, 0)
    counts[4, ] <- c(0, 0)
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    full_length <- matrix(0, nrow = 4, ncol = 2)
    full_length[1, ] <- c(5, 4)
    full_length[3, ] <- c(6, 6)
    SummarizedExperiment::assays(se, withDimnames = FALSE)[["fullLengthCounts"]] <- full_length

    result <- bambu_filter_transcripts(se)

    testthat::expect_equal(rownames(result$se), c("tx1", "tx2"))
    testthat::expect_equal(nrow(result$se), 2)
    testthat::expect_equal(result$qc_stats$transcripts_filtered, 2)
})

testthat::test_that("filters on counts when no fullLengthCounts assay", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # Set some transcripts with counts.
    counts <- SummarizedExperiment::assay(se, "counts")
    counts[1, ] <- c(10, 8)
    counts[2, ] <- c(5, 3)
    counts[3, ] <- c(0, 0)
    counts[4, ] <- c(0, 0)
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    # Remove fullLengthCounts assay entirely so it falls back to counts.
    assay_list <- SummarizedExperiment::assays(se)
    assay_list[["fullLengthCounts"]] <- NULL
    SummarizedExperiment::assays(se, withDimnames = FALSE) <- assay_list

    result <- bambu_filter_transcripts(se)

    testthat::expect_equal(nrow(result$se), 2)
    testthat::expect_equal(result$qc_stats$transcripts_filtered, 2)
})

testthat::test_that("error when all transcripts filtered", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # All transcripts have zero counts.
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- matrix(0, nrow = 4, ncol = 2)

    testthat::expect_error(
        bambu_filter_transcripts(se),
        "All transcripts have zero counts after filtering"
    )
})

testthat::test_that("QC stats match filtered results", {
    se <- make_test_tx_se(sample_names = c("sampleA", "sampleB"))
    # Set one transcript with counts, rest without.
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

# Contract test: after transcript filtering, the filtered object must remain
# acceptable input for gene-level aggregation.
testthat::test_that("bambu_filter_transcripts contract with transcriptToGeneExpression", {
    fixture_dir <- tempfile("bambu-filter-contract-")
    dir.create(fixture_dir)

    base_se <- make_test_tx_se(sample_names = "sampleA")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(base_se),
        rowRanges = make_test_bambu_row_ranges(fixture_dir)
    )

    counts <- SummarizedExperiment::assays(se)$counts
    counts["tx1", ] <- 10
    counts["tx2", ] <- 8
    counts["tx3", ] <- 0
    counts["tx4", ] <- 0
    SummarizedExperiment::assay(se, "counts", withDimnames = FALSE) <- counts

    full_length <- matrix(
        c(5, 4, 6, 0),
        nrow = nrow(se),
        ncol = ncol(se),
        dimnames = dimnames(SummarizedExperiment::assays(se)$counts)
    )
    SummarizedExperiment::assays(se, withDimnames = FALSE)[["fullLengthCounts"]] <- full_length

    S4Vectors::metadata(se)$incompatibleCounts <- data.table::data.table(
        GENEID = "gene2",
        sampleA = 7L
    )

    filtered <- bambu_filter_transcripts(se)
    testthat::expect_true(all(SummarizedExperiment::rowData(filtered$se)$GENEID == "gene1"))
    testthat::expect_equal(
        unique(S4Vectors::metadata(filtered$se)$incompatibleCounts$GENEID),
        character(0)
    )

    testthat::expect_error(bambu::transcriptToGeneExpression(filtered$se), NA)
})

testthat::test_that("bambu_filter_transcripts renames generic incompatibleCounts columns", {
    se <- make_test_tx_se(sample_names = "sampleA")
    S4Vectors::metadata(se)$incompatibleCounts <- data.table::data.table(
        GENEID = c("gene1", "gene2"),
        `01` = c(4, 9)
    )

    filtered <- bambu_filter_transcripts(se)
    incompatible <- S4Vectors::metadata(filtered$se)$incompatibleCounts

    testthat::expect_equal(names(incompatible), c("GENEID", "sampleA"))
    testthat::expect_equal(incompatible$sampleA, c(4, 9))
})

###
# CLI integration tests
#
# End-to-end tests with real bambu library (not mocked):
# - Build BAMs from committed fixtures (reference.fa, annotation.gtf, reads.fastq)
# - Run `supeRglue bambu discover` to produce rcFiles/chunks
# - Run `supeRglue bambu quant` on one emitted chunk
# - Run `supeRglue bambu collate` to produce downstream workflow outputs

testthat::test_that("CLI discover writes reusable chunk artifacts", {
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
            "discover",
            "--bams", bam_path,
            "--aliases", "sampleA",
            "--sample_sheet", sample_sheet,
            "--annotation", annotation,
            "--genome", reference,
            "--transcriptome_mode", "fixed_annotation",
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        result$status,
        0L,
        info = paste(result$output, collapse = "\n")
    )
    testthat::expect_true(file.exists(file.path(out_dir, "bambu_rcfiles.rds")))
    testthat::expect_true(file.exists(file.path(out_dir, "bambu_discovered_annotations.rds")))
    testthat::expect_true(file.exists(file.path(out_dir, "chunk_manifest.tsv")))
    testthat::expect_true(file.exists(file.path(out_dir, "samples.csv")))

    samples <- utils::read.csv(
        file.path(out_dir, "samples.csv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    manifest <- utils::read.delim(
        file.path(out_dir, "chunk_manifest.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    chunk_bundle <- readRDS(manifest$rds_path[[1]])

    testthat::expect_equal(samples$alias, "sampleA")
    testthat::expect_gt(nrow(manifest), 0)
    testthat::expect_true(all(c("chunk_id", "seqname", "annotation_tx_count", "rds_path") %in% names(manifest)))
    testthat::expect_true("annotation_tx_count" %in% names(chunk_bundle))
    testthat::expect_true(all(manifest$annotation_tx_count >= 0))
    testthat::expect_equal(chunk_bundle$aliases, "sampleA")
})

testthat::test_that("CLI quant consumes a discover chunk", {
    fixture_dir <- tempfile("bambu-cli-dir-")
    dir.create(fixture_dir)

    reference <- workflow_glue_r_fixture("bambu", "reference.fa")
    annotation <- workflow_glue_r_fixture("bambu", "annotation.gtf")
    reads <- workflow_glue_r_fixture("bambu", "reads.fastq")
    bam_path <- expect_bam_fixture_built(reference, reads, fixture_dir, alias = "sampleA")
    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    writeLines(
        paste(
            "alias",
            "sampleA",
            sep = "\n"
        ),
        sample_sheet
    )
    discover_out_dir <- file.path(fixture_dir, "discover")
    discover_result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "discover",
            "--bams", bam_path,
            "--aliases", "sampleA",
            "--sample_sheet", sample_sheet,
            "--annotation", annotation,
            "--genome", reference,
            "--transcriptome_mode", "fixed_annotation",
            "--out_dir", discover_out_dir
        )
    )

    testthat::expect_equal(
        discover_result$status,
        0L,
        info = paste(discover_result$output, collapse = "\n")
    )

    manifest <- utils::read.delim(
        file.path(discover_out_dir, "chunk_manifest.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    out_dir <- file.path(fixture_dir, "quant")
    quant_result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "quant",
            "--chunk_rds", manifest$rds_path[[1]],
            "--discovered_annotation_rds", file.path(discover_out_dir, "bambu_discovered_annotations.rds"),
            "--genome", reference,
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        quant_result$status,
        0L,
        info = paste(quant_result$output, collapse = "\n")
    )

    samples <- utils::read.csv(
        file.path(out_dir, "samples.csv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    tx_se <- readRDS(file.path(out_dir, "bambu_transcripts.rds"))
    qc <- jsonlite::read_json(
        file.path(out_dir, "bambu_qc_stats.json"),
        simplifyVector = TRUE
    )

    testthat::expect_equal(samples$alias, "sampleA")
    testthat::expect_gt(nrow(tx_se), 0)
    testthat::expect_equal(qc$chunk_id, manifest$chunk_id[[1]])
    testthat::expect_equal(qc$seqname, manifest$seqname[[1]])
})

testthat::test_that("CLI collate consumes quant chunk directories", {
    fixture_dir <- tempfile("bambu-cli-collate-")
    dir.create(fixture_dir)

    reference <- workflow_glue_r_fixture("bambu", "reference.fa")
    annotation <- workflow_glue_r_fixture("bambu", "annotation.gtf")
    reads <- workflow_glue_r_fixture("bambu", "reads.fastq")
    bam_path <- expect_bam_fixture_built(reference, reads, fixture_dir, alias = "sampleA")
    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    writeLines(
        paste(
            "alias",
            "sampleA",
            sep = "\n"
        ),
        sample_sheet
    )

    discover_out_dir <- file.path(fixture_dir, "discover")
    discover_result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "discover",
            "--bams", bam_path,
            "--aliases", "sampleA",
            "--sample_sheet", sample_sheet,
            "--annotation", annotation,
            "--genome", reference,
            "--transcriptome_mode", "fixed_annotation",
            "--out_dir", discover_out_dir
        )
    )
    testthat::expect_equal(
        discover_result$status,
        0L,
        info = paste(discover_result$output, collapse = "\n")
    )

    manifest <- utils::read.delim(
        file.path(discover_out_dir, "chunk_manifest.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    chunk_out_dir <- file.path(fixture_dir, "chunk1")
    quant_result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "quant",
            "--chunk_rds", manifest$rds_path[[1]],
            "--discovered_annotation_rds", file.path(discover_out_dir, "bambu_discovered_annotations.rds"),
            "--genome", reference,
            "--out_dir", chunk_out_dir
        )
    )
    testthat::expect_equal(
        quant_result$status,
        0L,
        info = paste(quant_result$output, collapse = "\n")
    )

    collate_out_dir <- file.path(fixture_dir, "sampleA")
    collate_result <- run_rscript(
        "supeRglue",
        c(
            "bambu",
            "collate",
            "--chunk_dirs", chunk_out_dir,
            "--transcriptome_mode", "fixed_annotation",
            "--out_dir", collate_out_dir
        )
    )
    testthat::expect_equal(
        collate_result$status,
        0L,
        info = paste(collate_result$output, collapse = "\n")
    )

    testthat::expect_true(file.exists(file.path(collate_out_dir, "bambu_transcripts.rds")))
    testthat::expect_true(file.exists(file.path(collate_out_dir, "bambu_genes.rds")))
    testthat::expect_true(file.exists(file.path(collate_out_dir, "transcript_counts.tsv")))
    testthat::expect_true(file.exists(file.path(collate_out_dir, "gene_counts.tsv")))
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
