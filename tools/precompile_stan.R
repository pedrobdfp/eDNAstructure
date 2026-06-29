#!/usr/bin/env Rscript
# =============================================================================
# tools/precompile_stan.R
# =============================================================================
# Run this ONCE on your machine (and re-run whenever you change dmm.stan or
# upgrade rstan/StanHeaders) to generate the precompiled Stan model that ships
# inside the package at inst/stan/dmm_model.rds.
#
# Users on a matching rstan/StanHeaders/OS/arch will load this object instantly
# with no compilation. Users on a mismatched setup fall back automatically to
# compiling dmm.stan on first use (see R/stanmodels.R) — so shipping this file
# is a pure speed-up, never a correctness requirement.
#
# USAGE (from the package root directory):
#   Rscript tools/precompile_stan.R
#
# or inside R:
#   source("tools/precompile_stan.R")
# =============================================================================

message("== eDNAstructure: precompiling Stan model ==")

if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("rstan is required to precompile the model. install.packages('rstan')")
}

# Locate dmm.stan relative to this script / the package root.
stan_src <- NULL
candidates <- c(
  file.path("inst", "stan", "dmm.stan"),  # run from package root
  file.path("stan", "dmm.stan"),          # run from inst/
  "dmm.stan"                               # run from inst/stan/
)
for (cand in candidates) {
  if (file.exists(cand)) { stan_src <- cand; break }
}
if (is.null(stan_src)) {
  stop(
    "Could not find inst/stan/dmm.stan. ",
    "Run this script from the package root: Rscript tools/precompile_stan.R"
  )
}
message("Found Stan source: ", normalizePath(stan_src))

out_path <- file.path(dirname(stan_src), "dmm_model.rds")

# Report the toolchain this object will be tied to, for reproducibility.
message("rstan version:       ", as.character(utils::packageVersion("rstan")))
if (requireNamespace("StanHeaders", quietly = TRUE)) {
  message("StanHeaders version: ", as.character(utils::packageVersion("StanHeaders")))
}
message("R version:           ", paste0(R.version$major, ".", R.version$minor))
message("Platform:            ", R.version$platform)

message("\nCompiling (this takes ~30-60 seconds)...")
t0 <- Sys.time()
model <- rstan::stan_model(file = stan_src, model_name = "dmm")
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message("Compiled in ", elapsed, "s.")

saveRDS(model, out_path)
message("Wrote precompiled model to: ", normalizePath(out_path))
message(
  "\nNext steps:\n",
  "  1. git add inst/stan/dmm_model.rds\n",
  "  2. git commit -m 'Add precompiled Stan model'\n",
  "  3. git push\n",
  "Users will then load it instantly (with automatic fallback to compile-",
  "from-source where the precompiled object is not usable)."
)
