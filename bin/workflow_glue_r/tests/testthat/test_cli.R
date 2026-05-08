#' Tests covering the CLI introspection logic of workflow_glue_r_cli.

testthat::test_that("workflow_glue_r_components discovers command entrypoints", {
    components <- workflow_glue_r_components()

    testthat::expect_true("bambu" %in% names(components))
    testthat::expect_true("de_analysis" %in% names(components))
    testthat::expect_equal(components$bambu$parser_name, "bambu_arg_parser")
    testthat::expect_equal(components$bambu$runner_name, "run_bambu_cli")
    testthat::expect_equal(components$de_analysis$parser_name, "de_analysis_arg_parser")
    testthat::expect_equal(components$de_analysis$runner_name, "run_de_analysis_cli")
})

testthat::test_that("workflow_glue_r_cli routes to an introspected command", {
    env <- new.env(parent = emptyenv())
    env$toy_arg_parser <- function() TRUE
    env$run_toy_cli <- function(argv) paste(argv, collapse = ",")

    testthat::expect_equal(
        workflow_glue_r_cli(c("toy", "alpha", "beta"), env = env),
        "alpha,beta"
    )
})

testthat::test_that("workflow_glue_r_cli rejects unknown commands", {
    testthat::expect_error(
        workflow_glue_r_cli(c("missing")),
        "Unknown supeRglue command"
    )
})
