# =============================================================================
# eDNA_loo() — Leave-one-out cross-validation for K selection
# =============================================================================

#' Compare DMM models across values of K using LOO cross-validation
#'
#' @description
#' Fits DMM models for a range of K values and compares them using
#' Leave-One-Out cross-validation (LOO-CV) via the `loo` package. Returns
#' both a comparison table and an elbow plot of LOO-ELPD vs K.
#'
#' LOO-ELPD (expected log predictive density) measures how well each model
#' predicts held-out observations. Higher values are better. The "elbow" —
#' the point where additional communities yield diminishing LOO-ELPD gains —
#' is a useful heuristic for choosing K.
#'
#' @param counts A count matrix (same as passed to [eDNA_dmm()]).
#' @param covariates A covariate matrix or data frame, or `NULL`.
#' @param K_range An integer vector of K values to evaluate. Default is `2:5`.
#'   Example: `K_range = 2:7`.
#' @param scale_covariates Logical. Default `TRUE`. See [eDNA_dmm()].
#' @param chains Integer. Number of chains per fit. Default `1`.
#' @param iter Integer. Iterations per chain. Default `4000`.
#' @param warmup Integer. Warmup iterations. Default `2000`.
#' @param adapt_delta Number. Default `0.95`.
#' @param seed Integer. Default `13`.
#' @param conc Number. Default `0.5`. See [eDNA_dmm()].
#' @param alpha_shape Number. Default `5`. See [eDNA_dmm()].
#' @param alpha_rate Number. Default `2`. See [eDNA_dmm()].
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`loo_table`}{A data frame with LOO-ELPD estimates and SE per K.}
#'   \item{`loo_compare`}{The output of [loo::loo_compare()].}
#'   \item{`plot`}{A [ggplot2::ggplot()] elbow plot of LOO-ELPD vs K.}
#'   \item{`fits`}{A named list of `edna_dmm_fit` objects, one per K.}
#' }
#'
#' @seealso [eDNA_dmm()]
#' @export
eDNA_loo <- function(
    counts,
    covariates      = NULL,
    K_range         = 2:5,
    scale_covariates = TRUE,
    chains          = 1,
    iter            = 4000,
    warmup          = 2000,
    adapt_delta     = 0.95,
    seed            = 13,
    conc            = 0.5,
    alpha_shape     = 5,
    alpha_rate      = 2,
    verbose         = TRUE
) {
  if (!all(K_range == round(K_range)) || any(K_range < 2)) {
    rlang::abort("`K_range` must be a vector of integers all >= 2.")
  }
  K_range <- as.integer(K_range)

  fits      <- vector("list", length(K_range))
  names(fits) <- paste0("K", K_range)
  loo_list  <- vector("list", length(K_range))

  for (i in seq_along(K_range)) {
    k <- K_range[i]
    if (verbose) message(sprintf("\n===== Fitting K = %d =====", k))

    fits[[i]] <- eDNA_dmm(
      counts           = counts,
      covariates       = covariates,
      K                = k,
      scale_covariates = scale_covariates,
      chains           = chains,
      iter             = iter,
      warmup           = warmup,
      adapt_delta      = adapt_delta,
      seed             = seed,
      conc             = conc,
      alpha_shape      = alpha_shape,
      alpha_rate       = alpha_rate,
      verbose          = verbose
    )
    ll_mat       <- loo::extract_log_lik(fits[[i]]$stan_fit, parameter_name = "log_lik")
    loo_list[[i]] <- loo::loo(ll_mat)
    if (verbose) message(sprintf("K = %d: LOO-ELPD = %.1f (SE = %.1f)",
                                 k,
                                 loo_list[[i]]$estimates["elpd_loo", "Estimate"],
                                 loo_list[[i]]$estimates["elpd_loo", "SE"]))
  }

  names(loo_list) <- paste0("K", K_range)

  # ── LOO comparison ────────────────────────────────────────────────────────────
  loo_compare_result <- do.call(loo::loo_compare, unname(loo_list))

  loo_df <- data.frame(
    K    = K_range,
    elpd = vapply(loo_list, function(l) l$estimates["elpd_loo", "Estimate"], numeric(1)),
    se   = vapply(loo_list, function(l) l$estimates["elpd_loo", "SE"],       numeric(1))
  )

  # ── Elbow plot ────────────────────────────────────────────────────────────────
  p_loo <- ggplot2::ggplot(loo_df, ggplot2::aes(x = .data$K, y = .data$elpd)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$elpd - .data$se, ymax = .data$elpd + .data$se),
      alpha = 0.2, fill = "#457B9D"
    ) +
    ggplot2::geom_line(linewidth = 0.9, color = "#457B9D") +
    ggplot2::geom_point(size = 3.5, color = "#457B9D") +
    ggplot2::scale_x_continuous(breaks = K_range) +
    ggplot2::labs(
      x        = "K (number of communities)",
      y        = "LOO-ELPD (higher is better)",
      title    = "K selection — LOO cross-validation",
      subtitle = "Shaded band = \u00b11 SE  |  Elbow = point of diminishing returns"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40")
    )

  list(
    loo_table   = loo_df,
    loo_compare = loo_compare_result,
    plot        = p_loo,
    fits        = fits
  )
}


# =============================================================================
# get_example_data() — Built-in example dataset
# =============================================================================

#' Load the bundled example dataset
#'
#' @description
#' Returns the package's built-in simulated eDNA metabarcoding dataset: a small,
#' self-contained survey with known ground-truth community structure, generated
#' by [simulate_eDNA_survey()] (seed `2026`). Use it to learn the expected input
#' format for [eDNA_dmm()], follow the tutorial, and verify your installation.
#'
#' For the full field-by-field description of the returned object, see
#' [example_edna].
#'
#' @section Data format:
#' A named list. The two elements you need for [eDNA_dmm()] are:
#' \describe{
#'   \item{`counts`}{Integer matrix (20 samples x 40 taxa). Row names are
#'     stations (`STN_001`-`STN_020`); column names are taxa (`Sp_1`-`Sp_40`).}
#'   \item{`covariates`}{Data frame (20 rows) with columns `sample_id`,
#'     `TrueCommunity`, `Depth` (m), and `Distance_shore`.}
#' }
#' The list also carries the true community compositions and simulation
#' metadata; see [example_edna].
#'
#' @section Formatting your own data:
#' The `counts` matrix format is the required input to [eDNA_dmm()]:
#' - **Rows** = samples (one row per station, individual, or replicate)
#' - **Columns** = taxa (or ASVs - taxonomic annotation is not required)
#' - **Values** = non-negative integer read counts
#'
#' If your count table is in long format (sample, taxon, count in three
#' columns), convert it to wide format with [tidyr::pivot_wider()]:
#' ```r
#' library(tidyr)
#' count_matrix <- pivot_wider(
#'   long_data,
#'   names_from  = taxon,
#'   values_from = reads,
#'   values_fill = 0
#' ) |>
#'   tibble::column_to_rownames("sample_id") |>
#'   as.matrix()
#' ```
#'
#' @return A named list; see the **Data format** section and [example_edna].
#'
#' @seealso [example_edna], [simulate_eDNA_survey()], [eDNA_dmm()]
#'
#' @examples
#' data <- get_example_data()
#'
#' # Inspect the count matrix
#' dim(data$counts)           # 20 samples x 40 taxa
#' head(data$counts[, 1:5])   # first 5 taxa
#'
#' # Inspect the covariate data frame
#' head(data$covariates)
#'
#' \dontrun{
#' fit <- eDNA_dmm(
#'   counts     = data$counts,
#'   covariates = data$covariates[, c("Depth", "Distance_shore")],
#'   K          = 2
#' )
#' print(fit)
#' eDNA_dmm_structure(fit)
#' }
#'
#' @export
get_example_data <- function() {
  # example_edna is a lazy-loaded dataset (data/example_edna.rda). If the
  # package was installed before the dataset was generated, give a clear,
  # actionable message instead of a cryptic "object not found".
  env <- new.env(parent = emptyenv())
  loaded <- tryCatch({
    utils::data("example_edna", package = "eDNAstructure", envir = env)
    exists("example_edna", envir = env, inherits = FALSE)
  }, warning = function(w) FALSE, error = function(e) FALSE)

  if (!isTRUE(loaded)) {
    rlang::abort(
      c(
        "The bundled example dataset is not available in this installation.",
        i = "It is generated by data-raw/generate_example_data.R.",
        i = "Alternatively, simulate your own: simulate_eDNA_survey(seed = 2026)."
      )
    )
  }
  get("example_edna", envir = env, inherits = FALSE)
}
