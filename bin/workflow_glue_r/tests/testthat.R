args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) < 1) {
    stop("Unable to determine testthat.R path", call. = FALSE)
}

test_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
pkg_dir <- normalizePath(file.path(test_dir, ".."))
repo_root <- normalizePath(file.path(pkg_dir, "..", ".."))

source(file.path(pkg_dir, "load.R"))
workflow_glue_r_load(pkg_dir)

Sys.setenv(
    WORKFLOW_GLUE_R_PACKAGE_DIR = pkg_dir,
    WORKFLOW_GLUE_R_REPO_ROOT = repo_root,
    TEST_DATA = Sys.getenv("TEST_DATA", unset = file.path(repo_root, "test_data"))
)

testthat_dir <- file.path(pkg_dir, "tests", "testthat")
test_files <- if (dir.exists(testthat_dir)) {
    list.files(testthat_dir, pattern = "\\.[Rr]$", full.names = TRUE)
} else {
    character(0)
}

if (length(test_files) < 1) {
    message("No workflow-local R test files found; skipping.")
    quit(save = "no", status = 0)
}

testthat::test_dir(
    testthat_dir,
    reporter = "summary",
    stop_on_failure = TRUE
)
