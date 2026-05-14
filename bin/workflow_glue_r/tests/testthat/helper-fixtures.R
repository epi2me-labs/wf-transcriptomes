# Helpers

#' Get workflow_glue_r repository root directory from environment variable.
#' @return Normalized path to workflow_glue_r repository root directory.
#' @export
workflow_glue_r_repo_root <- function() {
    normalizePath(Sys.getenv("WORKFLOW_GLUE_R_REPO_ROOT"))
}


#' Get path to a test fixture file within the workflow_glue_r test data directory.
#' @param ... Path components within the test data directory.
#' @return Normalized path to the specified test fixture file.
#' @export
workflow_glue_r_fixture <- function(...) {
    file.path(Sys.getenv("TEST_DATA"), "workflow_glue_r", ...)
}

#' Run a system command and capture its output and status.
#' @param command Command to run.
#' @param args Arguments to the command.
#' @param stdout Whether to capture standard output (TRUE or file path).
#' @param stderr Whether to capture standard error (TRUE or file path).
#' @return A list with 'status' (exit status) and 'output' (captured output).
#' @export
run_system2 <- function(command, args, stdout = TRUE, stderr = TRUE) {
    output <- system2(command, args, stdout = stdout, stderr = stderr)
    if (is.character(output)) {
        status <- attr(output, "status")
        if (is.null(status)) {
            status <- 0L
        }
        return(list(status = status, output = output))
    }

    captured <- character(0)
    if (is.character(stderr) && file.exists(stderr)) {
        captured <- readLines(stderr, warn = FALSE)
    }

    list(status = as.integer(output), output = captured)
}
#' Run an R script using Rscript.
#' @param script_name Name of the R script to run.
#' @param args Arguments to pass to the R script.
#' @return A list with 'status' (exit status) and 'output' (captured output).
#' @export
run_rscript <- function(script_name, args) {
    command <- file.path(R.home("bin"), "Rscript")
    script <- file.path(workflow_glue_r_repo_root(), "bin", script_name)
    run_system2(command, c(script, args), stdout = TRUE, stderr = TRUE)
}

#' Create a test transcript-level SummarizedExperiment with synthetic count data.
#'
#' SummarizedExperiment is a Bioconductor container that holds:
#' - assays: matrices of counts/CPM (rows=transcripts, cols=samples)
#' - rowRanges: genomic coordinates and metadata for each transcript
#' - colData: sample metadata (not used here)
#'
#' This fixture creates realistic differential expression patterns:
#' - tx1: high in control, low in treated (downregulated)
#' - tx2: low in control, high in treated (upregulated)
#' - tx3: similar across conditions (not DE)
#' - tx4: low counts in both (filtered in real DE analysis)
#'
#' @param include_geneid Whether to include GENEID metadata column in rowRanges.
#' @param sample_names Optional vector of sample names to use as column names.
#' @return A SummarizedExperiment object with synthetic transcript-level data.
#' @export
make_test_tx_se <- function(include_geneid = TRUE, sample_names = NULL) {
    if (is.null(sample_names)) {
        sample_names <- c(
            "control_rep1",
            "control_rep2",
            "control_rep3",
            "treated_rep1",
            "treated_rep2",
            "treated_rep3"
        )
    }

    counts <- vapply(sample_names, function(sample_name) {
        condition <- sub("_rep[0-9]+$", "", sample_name)
        replicate_id <- suppressWarnings(as.integer(sub("^.*_rep", "", sample_name)))
        if (is.na(replicate_id)) {
            replicate_id <- match(sample_name, sample_names)
        }
        offset <- ((replicate_id - 1L) %% 3L) - 1L

        base_counts <- switch(
            condition,
            control = c(120, 18, 80, 15),
            treated = c(25, 95, 76, 14),
            treated2 = c(35, 85, 60, 40),
            baseline = c(120, 18, 80, 15),
            c(80, 40, 70, 20)
        )
        pmax(base_counts + c(offset, -offset, offset, 0), 1)
    }, numeric(4))
    dimnames(counts) <- list(c("tx1", "tx2", "tx3", "tx4"), sample_names)
    cpm <- t(t(counts) / colSums(counts)) * 1e6

    row_ranges <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(1, 101, 201, 301), width = 50),
        strand = rep("+", 4),
        TXNAME = rownames(counts),
        eqClassById = IRanges::CharacterList(list(c("1", "2"), "3", "4", "5"))
    )
    if (include_geneid) {
        S4Vectors::mcols(row_ranges)$GENEID <- c("gene1", "gene1", "gene2", "gene2")
    }

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts, CPM = cpm),
        rowRanges = row_ranges
    )
}

#' Create bambu-style transcript row ranges for output-writing tests.
#'
#' `bambu::writeToGTF()` expects a `GRangesList` shaped like the object
#' returned by `bambu::prepareAnnotations()`. This helper writes a minimal GTF
#' and returns that object with metadata columns used by test assertions.
#'
#' @param out_dir Directory where temporary annotation fixture is written.
#' @return A GRangesList with tx1-tx4 transcript entries and metadata.
#' @export
make_test_bambu_row_ranges <- function(out_dir) {
    gtf <- file.path(out_dir, "annotation.gtf")
    writeLines(
        c(
            paste(
                "chr1", "test", "transcript", "1", "50", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "1", "50", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "101", "150", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx2";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "101", "150", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx2"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "201", "250", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx3";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "201", "250", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx3"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "301", "350", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx4";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "301", "350", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx4"; exon_number "1";',
                sep = "\t"
            )
        ),
        gtf
    )

    row_ranges <- bambu::prepareAnnotations(gtf)[c("tx1", "tx2", "tx3", "tx4")]
    S4Vectors::mcols(row_ranges)$TXNAME <- names(row_ranges)
    S4Vectors::mcols(row_ranges)$GENEID <- c("gene1", "gene1", "gene2", "gene2")
    S4Vectors::mcols(row_ranges)$eqClassById <- IRanges::CharacterList(list(c("1", "2"), "3", "4", "5"))

    row_ranges
}

#' Create a test gene-level SummarizedExperiment by aggregating a transcript-level SummarizedExperiment.
#' @param sample_names Optional vector of sample names to use as column names.
#' @return A SummarizedExperiment object with synthetic gene-level data.
#' @export
make_test_gene_se <- function(sample_names = NULL) {
    if (is.null(sample_names)) {
        sample_names <- c(
            "control_rep1",
            "control_rep2",
            "control_rep3",
            "treated_rep1",
            "treated_rep2",
            "treated_rep3"
        )
    }

    tx_se <- make_test_tx_se(sample_names = sample_names)
    tx_counts <- SummarizedExperiment::assays(tx_se)$counts
    counts <- rbind(
        gene1 = tx_counts["tx1", ] + tx_counts["tx2", ],
        gene2 = tx_counts["tx3", ] + tx_counts["tx4", ]
    )
    cpm <- t(t(counts) / colSums(counts)) * 1e6

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts, CPM = cpm),
        rowData = S4Vectors::DataFrame(GENEID = rownames(counts))
    )
}

#' Create a test DE analysis input bundle with synthetic transcript and gene
#' SummarizedExperiments and a sample sheet.
#' @param out_dir Directory to write the output files.
#' @param levels Optional vector of condition levels to use in the sample sheet.
#' @return A list with paths to the generated transcript RDS, gene RDS, and sample sheet CSV.
#' @export
write_de_fixture_bundle <- function(out_dir, levels = c("control", "treated")) {
    sample_names <- unlist(lapply(levels, function(level) {
        sprintf("%s_rep%d", level, seq_len(3))
    }))

    tx_se <- make_test_tx_se(sample_names = sample_names)
    gene_se <- make_test_gene_se(sample_names = sample_names)

    tx_path <- file.path(out_dir, "transcripts.rds")
    gene_path <- file.path(out_dir, "genes.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(out_dir, "sample_sheet.csv")
    sample_rows <- do.call(
        rbind,
        lapply(levels, function(level) {
            data.frame(
                alias = sprintf("%s_rep%d", level, seq_len(3)),
                condition = level,
                batch = c("b1", "b2", "b1"),
                stringsAsFactors = FALSE
            )
        })
    )
    utils::write.csv(sample_rows, sample_sheet, row.names = FALSE, quote = FALSE)

    list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet
    )
}

# Expecters

#' Expect a system command to succeed.
#' @param command Command to run.
#' @param args Arguments to the command.
#' @param stdout Whether to capture standard output (TRUE or file path).
#' @param stderr Whether to capture standard error (TRUE or file path).
#' @return Invisibly returns the result of the command.
#' @export
expect_command_success <- function(command, args, stdout = TRUE, stderr = TRUE) {
    result <- run_system2(command, args, stdout = stdout, stderr = stderr)
    testthat::expect_equal(
        result$status,
        0L,
        info = paste(c(command, args, result$output), collapse = "\n")
    )
    invisible(result)
}

#' Expect a BAM fixture can be built with minimap2 and samtools.
#' @param reference Path to the reference FASTA file.
#' @param reads Path to the FASTQ file containing reads to align.
#' @param out_dir Directory to write the output BAM file and intermediate files.
#' @param alias Sample alias to use in output file names.
#' @return Path to the generated BAM file.
#' @export
expect_bam_fixture_built <- function(reference, reads, out_dir, alias = "sampleA") {
    sam_path <- file.path(out_dir, sprintf("%s.sam", alias))
    bam_path <- file.path(out_dir, sprintf("%s.aligned.sorted.bam", alias))
    minimap2_stderr <- file.path(out_dir, sprintf("%s.minimap2.stderr.txt", alias))
    sort_stderr <- file.path(out_dir, sprintf("%s.samtools-sort.stderr.txt", alias))
    index_stderr <- file.path(out_dir, sprintf("%s.samtools-index.stderr.txt", alias))

    expect_command_success(
        "minimap2",
        c("-ax", "splice", "-uf", reference, reads),
        stdout = sam_path,
        stderr = minimap2_stderr
    )
    testthat::expect_true(file.exists(sam_path))
    testthat::expect_gt(file.info(sam_path)$size, 0)

    expect_command_success(
        "samtools",
        c("sort", "-o", bam_path, sam_path),
        stderr = sort_stderr
    )
    expect_command_success(
        "samtools",
        c("index", bam_path),
        stderr = index_stderr
    )

    bam_path
}
