#' eDNAstructure: Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding
#'
#' @description
#' `eDNAstructure` fits Bayesian Dirichlet-Multinomial Mixture (DMM) models
#' to eDNA read count data from metabarcoding surveys. The model simultaneously
#' estimates:
#'
#' - How many latent ecological communities are present
#' - The taxonomic composition (relative frequencies) of each community
#' - How environmental covariates (depth, latitude, etc.) drive community membership
#'
#' Input can be raw ASV tables, taxonomically-annotated taxon tables, or any
#' sample by feature count matrix.
#'
#' @section Core workflow:
#' 1. Fit the model: [eDNA_dmm()]
#' 2. Visualize assignments: [eDNA_dmm_structure()]
#' 3. Ordinate communities: [eDNA_dmm_nmds()]
#' 4. Inspect covariate effects: [eDNA_dmm_beta()]
#' 5. Get example data: [get_example_data()]
#'
#' @section Getting started:
#' ```r
#' data <- get_example_data()
#' fit  <- eDNA_dmm(counts = data$counts, covariates = data$covariates, K = 2)
#' eDNA_dmm_structure(fit)
#' ```
#'
#' @section Stan model compilation:
#' The Stan model is provided precompiled where possible and otherwise compiled
#' once on your machine and cached, so repeated calls to [eDNA_dmm()] go
#' straight to sampling. See [eDNA_clear_stan_cache()] if you ever need to
#' force a recompile (for example after upgrading `rstan`).
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rstan sampling stan_model extract rstan_options
#'   check_hmc_diagnostics
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_line geom_hline
#'   geom_vline geom_density geom_ribbon geom_text stat_ellipse
#'   scale_fill_manual scale_color_manual scale_alpha_manual
#'   scale_x_continuous scale_y_continuous scale_x_discrete
#'   scale_size_continuous facet_wrap facet_grid labs theme theme_bw
#'   element_text element_blank guides guide_legend unit ggsave
#' @importFrom dplyr select group_by summarise across all_of arrange mutate
#'   distinct case_when
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom vegan vegdist metaMDS scores
#' @importFrom posterior ess_bulk as_draws_matrix
#' @importFrom loo loo loo_compare extract_log_lik
#' @importFrom scales percent
#' @importFrom stats median quantile rnorm var sd setNames reformulate runif
#'   rbeta rgamma rnbinom rlnorm rmultinom
#' @importFrom grDevices hcl
#' @importFrom methods is slot
#' @importFrom tools R_user_dir
#' @importFrom utils data packageVersion globalVariables
#' @importFrom rlang .data abort warn inform caller_env
## usethis namespace: end
NULL

# Internal null-coalescing operator used across the package.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

.onLoad <- function(libname, pkgname) {
  # Cache compiled Stan programs that rstan itself writes (harmless, speeds up
  # any auxiliary rstan use). The DMM model has its own caching in
  # R/stanmodels.R and does not rely on a precompiled Rcpp module.
  rstan::rstan_options(auto_write = TRUE)
}
