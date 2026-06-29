# =============================================================================
# data-raw/generate_example_data.R
# =============================================================================
# Regenerates the bundled example dataset shipped at data/example_edna.rda and
# returned by get_example_data(). Run this once (and whenever the simulation
# code changes) from the package root:
#
#   source("data-raw/generate_example_data.R")
#
# Requires the package to be loaded/installed so simulate_eDNA_survey() is
# available. Easiest: devtools::load_all(".") first, then source this file.
# =============================================================================

# If running interactively without the package loaded, load it:
if (!exists("simulate_eDNA_survey", mode = "function")) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".")
  } else {
    stop("Load the package first: devtools::load_all('.')")
  }
}

set.seed(2026)

# Matches the @format documented in R/data.R:
#   20 samples (STN_001..STN_020), 40 taxa (Sp_1..Sp_40), K = 4 communities,
#   separated by Depth and Distance_shore, spillover = 0.15.
example_edna <- simulate_eDNA_survey(
  n_communities         = 4,
  n_species             = 40,
  samples_per_community = 5,    # 4 communities x 5 = 20 samples
  spillover             = 0.15,
  mean_read_depth       = 10000,
  seed                  = 2026
)

usethis::use_data(example_edna, overwrite = TRUE)

message("Wrote data/example_edna.rda (", nrow(example_edna$counts),
        " samples x ", ncol(example_edna$counts), " taxa).")
