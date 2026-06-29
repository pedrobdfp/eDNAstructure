# =============================================================================
# eDNA_dmm() — Fit a Dirichlet-Multinomial Mixture model
# =============================================================================

#' Fit a Dirichlet-Multinomial Mixture model to eDNA count data
#'
#' @description
#' `eDNA_dmm()` fits a Bayesian Dirichlet-Multinomial Mixture (DMM) model
#' to a sample × taxon read count matrix. It identifies `K` latent ecological
#' communities, estimates their taxonomic compositions, and models how
#' environmental covariates (e.g., depth, latitude) drive community membership
#' through a softmax regression.
#'
#' Inference is performed with Stan via [rstan::stan()]. The result is an
#' `edna_dmm_fit` object that can be passed to the visualization functions
#' [eDNA_dmm_structure()], [eDNA_dmm_nmds()], and [eDNA_dmm_beta()].
#'
#' @section Input format:
#' `counts` should be a **sample × taxon** matrix of non-negative integer read
#' counts. Rows are samples (stations, replicates, etc.) and columns are taxa
#' (or ASVs — the model does not require taxonomic annotation). Row names
#' are used as sample identifiers in plots; column names are used as taxon
#' labels.
#'
#' `covariates` should be a **sample × covariate** matrix or data frame. It
#' must have the same number of rows as `counts`, in the same order. Covariates
#' are Z-score standardized internally by default (`scale_covariates = TRUE`),
#' which is strongly recommended for interpretable beta coefficients and good
#' MCMC mixing. Pass `NULL` to fit an intercept-only model.
#'
#' @section The model:
#' For each sample `i`, the model marginalizes over a latent community
#' assignment `z_i`:
#'
#' - Community compositions: `pi_k ~ Dirichlet(conc * 1_S)` for each
#'   community `k`
#' - Community membership: `P(z_i = k) = softmax(beta_0k + beta_covariates)`
#' - Observed counts: `x_i | z_i = k ~ DirichletMultinomial(N_i, alpha * pi_k)`
#'
#' `alpha` is a global overdispersion scalar estimated from the data. Larger
#' `alpha` means lower overdispersion (counts are more multinomial-like).
#' Typical eDNA data have `alpha` in the range 1–5.
#'
#' @section Single-chain recommendation:
#' By default, `eDNA_dmm()` fits a **single MCMC chain** (`chains = 1`).
#' This is intentional. Mixture models suffer from *label switching*: across
#' multiple chains, "Community 1" may refer to different groups, making
#' multi-chain Rhat diagnostics uninformative (they will always look bad, even
#' when each chain converges perfectly). Within a single long chain, label
#' switching is extremely rare given good initialization. Use the within-chain
#' ESS reported in the output to assess convergence instead.
#'
#' @section Output object:
#' Returns an `edna_dmm_fit` object, which is also a standard R `list` with
#' elements:
#' \describe{
#'   \item{`stan_fit`}{The raw [rstan::stanfit-class] object.}
#'   \item{`sample_info`}{Data frame: one row per sample with posterior
#'     community membership probabilities (`prob_comm1`, `prob_comm2`, ...),
#'     MAP assignment (`z_hat`), and assignment certainty.}
#'   \item{`pi_mean`}{Matrix [K × S]: posterior mean community compositions.
#'     `pi_mean[k, j]` = posterior mean relative frequency of taxon `j` in
#'     community `k`.}
#'   \item{`beta_summary`}{Data frame: posterior summaries for each softmax
#'     regression coefficient (mean, 90% CI, P(direction), ESS).}
#'   \item{`alpha_mean`}{Scalar: posterior mean of the overdispersion parameter.}
#'   \item{`K`}{Number of communities fitted.}
#'   \item{`N`}{Number of samples.}
#'   \item{`S`}{Number of taxa.}
#'   \item{`taxa_names`}{Character vector of taxon names (column names of `counts`).}
#'   \item{`covariate_names`}{Character vector of covariate names.}
#'   \item{`scale_info`}{List with `center` and `scale` vectors used for
#'     Z-scoring (NULL if `scale_covariates = FALSE`).}
#'   \item{`counts`}{The filtered count matrix as passed to Stan.}
#'   \item{`stan_data`}{The list passed to Stan (for advanced diagnostics).}
#'   \item{`call`}{The matched function call.}
#' }
#'
#' @param counts A numeric matrix or data frame of **non-negative integer** read
#'   counts. Rows are samples; columns are taxa or ASVs. Row names (if present)
#'   are used as sample IDs in plots. Column names are used as taxon labels.
#' @param covariates A numeric matrix or data frame of environmental covariates
#'   (rows = samples, columns = covariates). Must have the same number of rows
#'   as `counts`, in the same order. Set to `NULL` (the default) to fit an
#'   intercept-only model with no covariate effects on community membership.
#' @param K A single positive integer ≥ 2: the number of latent communities.
#'   If you are unsure, start with `K = 2` and increase. The function
#'   [eDNA_loo()] can help compare models across values of K.
#' @param scale_covariates Logical. If `TRUE` (the default), covariates are
#'   Z-score standardized (mean = 0, SD = 1) before fitting. This is
#'   **strongly recommended** because it (a) improves MCMC mixing,
#'   (b) makes beta coefficients directly comparable across covariates, and
#'   (c) ensures the Normal(0, 1) prior on beta is weakly informative. Set
#'   to `FALSE` only if you have already standardized your covariates manually.
#' @param chains A single positive integer: the number of MCMC chains.
#'   Default is `1`. See the **Single-chain recommendation** section above.
#'   If you want multi-chain runs (e.g., for sensitivity checks), increase
#'   this, but interpret Rhat values with caution in the context of mixture
#'   models.
#' @param iter A single positive integer: total number of MCMC iterations per
#'   chain, **including** warmup. Default is `4000`. Increase for complex models
#'   (large K, many covariates) or when ESS is low. A chain of length 4000 with
#'   2000 warmup gives 2000 post-warmup draws.
#' @param warmup A single positive integer: the number of warmup (burn-in)
#'   iterations per chain. Default is `2000` (half of `iter`). Must be less
#'   than `iter`. During warmup, Stan adapts the step size; these draws are
#'   discarded.
#' @param adapt_delta A number in (0, 1): the target average acceptance
#'   probability for the NUTS sampler. Default is `0.95`. Increase toward
#'   `0.99` if you see divergent transitions in the diagnostics. Higher values
#'   slow sampling but reduce divergences.
#' @param max_treedepth A positive integer: the maximum tree depth for the
#'   NUTS sampler. Default is `12`. Increase to `14` or `15` if you see
#'   "maximum treedepth exceeded" warnings.
#' @param seed A single integer: the random seed for Stan. Default is `13`.
#'   Set a fixed seed for reproducibility across runs.
#' @param conc A positive number: the Dirichlet concentration parameter for
#'   the prior on community compositions (`pi`). Default is `0.5`.
#'   - `conc < 1` (e.g., 0.5): sparse compositions — each community is
#'     dominated by a few taxa. Usually appropriate for eDNA data.
#'   - `conc = 1`: flat (symmetric Dirichlet) prior — all compositions
#'     equally likely. Uninformative.
#'   - `conc > 1`: concentrates compositions toward uniform. Use only
#'     if you expect all taxa to be equally abundant in each community.
#' @param alpha_shape A positive number: the shape parameter of the Gamma
#'   prior on the overdispersion parameter `alpha`. Default is `5`. The
#'   prior mean of `alpha` is `alpha_shape / alpha_rate`. Larger `alpha`
#'   means less overdispersion (counts closer to multinomial).
#' @param alpha_rate A positive number: the rate parameter of the Gamma
#'   prior on `alpha`. Default is `2`. Prior mean = `alpha_shape / alpha_rate`
#'   = 2.5. Adjust if you have strong prior information about overdispersion.
#' @param verbose Logical. If `TRUE` (the default), print Stan compilation and
#'   sampling progress. Set to `FALSE` for silent fitting (useful in loops
#'   over multiple K values).
#'
#' @return An `edna_dmm_fit` object (a list). See the **Output object** section
#'   for full documentation of all elements.
#'
#' @seealso
#' - [eDNA_dmm_structure()] for structure (STRUCTURE-like bar) plots
#' - [eDNA_dmm_nmds()] for NMDS ordination colored by community
#' - [eDNA_dmm_beta()] for covariate coefficient plots with prior/posterior comparison
#' - [eDNA_loo()] for LOO-CV model comparison across K values
#' - [get_example_data()] for a built-in example dataset and tutorial
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data <- get_example_data()
#'
#' # Fit K=2, with depth and latitude as covariates
#' fit <- eDNA_dmm(
#'   counts     = data$counts,
#'   covariates = data$covariates[, c("depth", "latitude")],
#'   K          = 2
#' )
#'
#' # Summary of the fit
#' print(fit)
#' summary(fit)
#'
#' # Intercept-only model (no covariates)
#' fit_null <- eDNA_dmm(
#'   counts     = data$counts,
#'   covariates = NULL,
#'   K          = 2
#' )
#' }
#'
#' @export
eDNA_dmm <- function(
    counts,
    covariates      = NULL,
    K               = 2,
    scale_covariates = TRUE,
    chains          = 1,
    iter            = 4000,
    warmup          = 2000,
    adapt_delta     = 0.95,
    max_treedepth   = 12,
    seed            = 13,
    conc            = 0.5,
    alpha_shape     = 5,
    alpha_rate      = 2,
    verbose         = TRUE
) {
  cl <- match.call()
  
  # ── Input validation ────────────────────────────────────────────────────────
  counts     <- validate_counts(counts)
  K          <- validate_K(K, nrow(counts))
  
  # Drop all-zero taxa silently after warning in validate_counts
  counts <- counts[, colSums(counts) > 0, drop = FALSE]
  
  covariates <- validate_covariates(covariates, counts, scale_covariates)
  
  # ── Argument checks ──────────────────────────────────────────────────────────
  if (!is.numeric(chains) || chains < 1 || chains != round(chains)) {
    rlang::abort("`chains` must be a positive integer.")
  }
  if (!is.numeric(iter) || iter < 100) {
    rlang::abort("`iter` must be a positive integer >= 100.")
  }
  if (!is.numeric(warmup) || warmup < 1 || warmup >= iter) {
    rlang::abort("`warmup` must be a positive integer less than `iter`.")
  }
  if (!is.numeric(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    rlang::abort("`adapt_delta` must be a number strictly between 0 and 1 (e.g., 0.95).")
  }
  if (!is.numeric(conc) || conc <= 0) {
    rlang::abort("`conc` must be a positive number.")
  }
  if (!is.numeric(alpha_shape) || alpha_shape <= 0) {
    rlang::abort("`alpha_shape` must be a positive number.")
  }
  if (!is.numeric(alpha_rate) || alpha_rate <= 0) {
    rlang::abort("`alpha_rate` must be a positive number.")
  }
  
  N <- nrow(counts)
  S <- ncol(counts)
  P <- ncol(covariates)
  
  # ── Recover scaling info before we lose attributes ───────────────────────────
  scale_info <- NULL
  if (scale_covariates && P > 0) {
    scale_info <- list(
      center = attr(covariates, "scale_center"),
      scale  = attr(covariates, "scale_scale")
    )
    # Strip custom attributes so Stan doesn't complain
    attr(covariates, "scale_center") <- NULL
    attr(covariates, "scale_scale")  <- NULL
  }
  
  taxa_names      <- colnames(counts)
  covariate_names <- if (P > 0) colnames(covariates) else character(0)
  sample_ids      <- rownames(counts)
  
  # ── Build Stan data list ─────────────────────────────────────────────────────
  cov_matrix <- if (P > 0) covariates else matrix(numeric(0), nrow = N, ncol = 0)
  
  stan_data <- list(
    N           = N,
    S           = S,
    K           = K,
    P           = P,
    X           = counts,
    covariates  = cov_matrix,
    conc        = conc,
    alpha_shape = alpha_shape,
    alpha_rate  = alpha_rate
  )
  
  # ── Fit using lazily-compiled, disk-cached Stan model ────────────────────────
  # .get_dmm_stanmodel() (see R/stanmodels.R) compiles the model from
  # inst/stan/dmm.stan on first use and caches it to disk. No precompiled
  # Rcpp Module is shipped in the package's compiled code.
  if (verbose) {
    message(sprintf(
      "Fitting DMM: N=%d samples, S=%d taxa, K=%d communities, P=%d covariates",
      N, S, K, P
    ))
    message(sprintf(
      "MCMC: %d chain(s), %d iterations (%d warmup, %d sampling)",
      chains, iter, warmup, iter - warmup
    ))
    if (chains == 1) {
      message("Note: Using single chain (recommended for mixture models; see ?eDNA_dmm).")
    }
  }
  
  stan_fit <- rstan::sampling(
    object  = .get_dmm_stanmodel(),
    data    = stan_data,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    cores   = 1L,
    seed    = seed,
    verbose = FALSE,
    refresh = if (verbose) max(1, (iter - warmup) %/% 10) else 0,
    control = list(
      adapt_delta   = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
  
  # ── Extract posterior summaries ───────────────────────────────────────────────
  n_post <- iter - warmup
  
  cp_draws  <- rstan::extract(stan_fit, pars = "community_probs")$community_probs
  # cp_draws is [draws, N, K]; take chain 1 draws only (first n_post rows)
  cp_mean   <- apply(cp_draws[seq_len(n_post), , , drop = FALSE], c(2, 3), mean)
  
  pi_draws  <- rstan::extract(stan_fit, pars = "pi")$pi
  # pi_draws is [draws, K, S]
  pi_mean   <- apply(pi_draws[seq_len(n_post), , , drop = FALSE], c(2, 3), mean)
  colnames(pi_mean) <- taxa_names
  
  alpha_draws <- rstan::extract(stan_fit, pars = "alpha")$alpha
  alpha_mean  <- mean(alpha_draws[seq_len(n_post)])
  
  # ── Sample assignment data frame ──────────────────────────────────────────────
  prob_df    <- as.data.frame(cp_mean)
  colnames(prob_df) <- paste0("prob_comm", seq_len(K))
  
  sample_info <- data.frame(
    sample_id           = if (!is.null(sample_ids)) sample_ids else paste0("Sample_", seq_len(N)),
    stringsAsFactors    = FALSE
  )
  sample_info <- cbind(sample_info, prob_df)
  sample_info$z_hat <- factor(
    apply(cp_mean, 1, which.max),
    levels = seq_len(K),
    labels = paste0("Community ", seq_len(K))
  )
  sample_info$assignment_certainty <- apply(cp_mean, 1, max)
  
  # ── Beta summaries ─────────────────────────────────────────────────────────────
  beta_summary <- .extract_beta_summary(stan_fit, K, P, covariate_names, n_post)
  
  # ── Diagnostics summary ───────────────────────────────────────────────────────
  if (verbose) {
    rstan::check_hmc_diagnostics(stan_fit)
    pi_mat  <- as.matrix(stan_fit, pars = "pi")[seq_len(n_post), ]
    pi_ess  <- posterior::ess_bulk(posterior::as_draws_matrix(pi_mat))
    message(sprintf(
      "Within-chain ESS (pi): min=%.0f, median=%.0f",
      min(pi_ess), stats::median(pi_ess)
    ))
    message(sprintf("Posterior mean alpha: %.2f", alpha_mean))
  }
  
  # ── Return ─────────────────────────────────────────────────────────────────────
  structure(
    list(
      stan_fit        = stan_fit,
      sample_info     = sample_info,
      pi_mean         = pi_mean,
      beta_summary    = beta_summary,
      alpha_mean      = alpha_mean,
      K               = K,
      N               = N,
      S               = S,
      taxa_names      = taxa_names,
      covariate_names = covariate_names,
      scale_info      = scale_info,
      counts          = counts,
      stan_data       = stan_data,
      call            = cl
    ),
    class = c("edna_dmm_fit", "list")
  )
}


# Internal helper: extract beta posterior summaries
#' @keywords internal
.extract_beta_summary <- function(stan_fit, K, P, covariate_names, n_post) {
  if (K < 2) return(data.frame())
  
  cov_labels <- c("intercept", covariate_names)
  n_beta_cols <- K - 1   # communities 1..(K-1) vs reference K
  
  beta_mat <- as.matrix(stan_fit, pars = "beta")
  if (nrow(beta_mat) >= n_post) beta_mat <- beta_mat[seq_len(n_post), , drop = FALSE]
  
  results <- vector("list", (K - 1) * (P + 1))
  idx <- 1L
  for (comm_i in seq_len(K - 1)) {
    for (cov_j in seq_len(P + 1)) {
      pname <- sprintf("beta[%d,%d]", comm_i, cov_j)
      if (!pname %in% colnames(beta_mat)) next
      draws   <- beta_mat[, pname]
      ess_val <- as.numeric(
        posterior::ess_bulk(
          posterior::as_draws_matrix(matrix(draws, ncol = 1,
                                            dimnames = list(NULL, pname)))
        )
      )
      results[[idx]] <- data.frame(
        community     = paste0("Community ", comm_i),
        reference     = paste0("Community ", K, " (reference)"),
        covariate     = cov_labels[cov_j],
        mean          = mean(draws),
        median        = stats::median(draws),
        ci_5          = stats::quantile(draws, 0.05),
        ci_95         = stats::quantile(draws, 0.95),
        ci_10         = stats::quantile(draws, 0.10),
        ci_90         = stats::quantile(draws, 0.90),
        prob_positive = mean(draws > 0),
        prob_negative = mean(draws < 0),
        ess           = ess_val,
        reliability   = dplyr::case_when(
          ess_val > 400 ~ "trustworthy",
          ess_val > 100 ~ "cautious",
          TRUE          ~ "unreliable"
        ),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, Filter(Negate(is.null), results))
}


# =============================================================================
# S3 methods for edna_dmm_fit
# =============================================================================

#' @export
print.edna_dmm_fit <- function(x, ...) {
  cat("eDNA Dirichlet-Multinomial Mixture Model\n")
  cat("=========================================\n")
  cat(sprintf("  K (communities)  : %d\n", x$K))
  cat(sprintf("  N (samples)      : %d\n", x$N))
  cat(sprintf("  S (taxa)         : %d\n", x$S))
  cat(sprintf("  Covariates       : %s\n",
              if (length(x$covariate_names) > 0)
                paste(x$covariate_names, collapse = ", ")
              else "(none — intercept-only)"))
  cat(sprintf("  Mean alpha       : %.2f\n", x$alpha_mean))
  cat("\nCommunity sizes (MAP assignment):\n")
  tbl <- table(x$sample_info$z_hat)
  for (nm in names(tbl)) {
    cat(sprintf("  %-15s: %d samples (%.1f%%)\n",
                nm, tbl[[nm]], 100 * tbl[[nm]] / x$N))
  }
  cat("\nUse summary() for convergence diagnostics, or pass this object to:\n")
  cat("  eDNA_dmm_structure()  — structure bar plots\n")
  cat("  eDNA_dmm_nmds()       — NMDS ordination\n")
  cat("  eDNA_dmm_beta()       — covariate coefficient plots\n")
  invisible(x)
}

#' @export
summary.edna_dmm_fit <- function(object, ...) {
  cat("eDNA DMM — Fit Summary\n")
  cat("======================\n\n")
  
  cat("Call:\n  ")
  print(object$call)
  cat("\n")
  
  # Model dimensions
  cat(sprintf("Dimensions: N=%d samples, S=%d taxa, K=%d communities\n\n",
              object$N, object$S, object$K))
  
  # Alpha (overdispersion)
  cat("Overdispersion (alpha):\n")
  alpha_draws <- rstan::extract(object$stan_fit, pars = "alpha")$alpha
  cat(sprintf("  Mean = %.2f, Median = %.2f, 90%% CI = [%.2f, %.2f]\n",
              mean(alpha_draws), stats::median(alpha_draws),
              stats::quantile(alpha_draws, 0.05),
              stats::quantile(alpha_draws, 0.95)))
  cat("  (alpha >> 1: low overdispersion / near-multinomial;\n")
  cat("   alpha ~  1: high overdispersion / typical eDNA)\n\n")
  
  # Assignment certainty
  cert <- object$sample_info$assignment_certainty
  cat("Assignment certainty across samples:\n")
  cat(sprintf("  Mean = %.2f, Min = %.2f, Max = %.2f\n",
              mean(cert), min(cert), max(cert)))
  cat(sprintf("  Decisive (>=80%% certainty): %d/%d samples (%.0f%%)\n\n",
              sum(cert >= 0.8), object$N, 100 * mean(cert >= 0.8)))
  
  # Beta summary (non-intercept terms only)
  if (nrow(object$beta_summary) > 0 && length(object$covariate_names) > 0) {
    beta_show <- object$beta_summary[object$beta_summary$covariate != "intercept", ]
    if (nrow(beta_show) > 0) {
      cat("Covariate effects (beta coefficients, 90% CI):\n")
      fmt <- "  %-18s vs ref  |  %s: %6.2f [%5.2f, %5.2f]  P(dir)=%.0f%%  [%s]\n"
      for (i in seq_len(nrow(beta_show))) {
        r <- beta_show[i, ]
        p_dir <- max(r$prob_positive, r$prob_negative)
        cat(sprintf(fmt,
                    r$community, r$covariate, r$mean, r$ci_5, r$ci_95,
                    100 * p_dir, r$reliability))
      }
      cat("\n  Reliability: trustworthy (ESS>400), cautious (100-400), unreliable (<100)\n\n")
    }
  }
  
  # HMC diagnostics
  cat("HMC Diagnostics:\n")
  rstan::check_hmc_diagnostics(object$stan_fit)
  
  # ESS (within-chain, single long chain recommended for mixtures)
  pi_mat   <- as.matrix(object$stan_fit, pars = "pi")
  pi_ess   <- posterior::ess_bulk(posterior::as_draws_matrix(pi_mat))
  cat(sprintf("\nWithin-chain ESS (pi): min=%.0f, median=%.0f, max=%.0f\n",
              min(pi_ess), stats::median(pi_ess), max(pi_ess)))
  cat("  ESS > 400 = trustworthy; ESS 100-400 = cautious; ESS < 100 = unreliable\n")
  
  invisible(object)
}