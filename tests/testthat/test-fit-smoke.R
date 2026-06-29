test_that("eDNA_dmm runs end-to-end on a tiny simulated dataset", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  # This compiles Stan on first run; can be slow. Skip unless explicitly enabled.
  testthat::skip_if(
    !identical(Sys.getenv("EDNASTRUCTURE_RUN_STAN_TESTS"), "true"),
    "Set EDNASTRUCTURE_RUN_STAN_TESTS=true to run Stan sampling tests."
  )

  sim <- simulate_eDNA_survey(n_communities = 2, n_species = 10,
                              samples_per_community = 5, seed = 7)
  fit <- eDNA_dmm(
    counts     = sim$counts,
    covariates = sim$covariates[, c("Depth", "Distance_shore")],
    K          = 2,
    iter       = 300, warmup = 150,
    verbose    = FALSE
  )
  expect_s3_class(fit, "edna_dmm_fit")
  expect_equal(fit$K, 2)
  expect_equal(nrow(fit$pi_mean), 2)
})
