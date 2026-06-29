test_that("validate_counts rejects non-numeric and bad shapes", {
  expect_error(eDNAstructure:::validate_counts("not a matrix"))
  expect_error(eDNAstructure:::validate_counts(matrix(1:2, nrow = 1)))  # too few rows
})

test_that("validate_K enforces 2 <= K < N", {
  expect_error(eDNAstructure:::validate_K(1, N = 10))
  expect_error(eDNAstructure:::validate_K(10, N = 10))
  expect_equal(eDNAstructure:::validate_K(3, N = 10), 3L)
})

test_that("simulate_eDNA_survey returns a count matrix of expected shape", {
  sim <- simulate_eDNA_survey(
    n_communities = 2, n_species = 12,
    samples_per_community = 4, seed = 1
  )
  expect_true(is.matrix(sim$counts))
  expect_equal(nrow(sim$counts), 8)          # 2 communities x 4 samples
  expect_true(all(sim$counts >= 0))
  expect_true(is.integer(sim$counts))
})
