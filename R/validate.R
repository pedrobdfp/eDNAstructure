# =============================================================================
# Input validation helpers
# All user-facing checks live here. The philosophy: catch mistakes early,
# explain what went wrong AND how to fix it.
# =============================================================================

#' @keywords internal
validate_counts <- function(counts, call = rlang::caller_env()) {
  # ── Type check ──────────────────────────────────────────────────────────────
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    rlang::abort(
      c(
        "`counts` must be a numeric matrix or data frame.",
        i = paste0("You supplied an object of class: ", paste(class(counts), collapse = ", ")),
        i = "Expected format: rows = samples, columns = taxa (or ASVs).",
        i = "Example: a matrix where counts[i, j] = reads of taxon j in sample i."
      ),
      call = call
    )
  }

  if (is.data.frame(counts)) {
    # Check for non-numeric columns (e.g. a stray Sample ID column)
    non_num <- names(counts)[!vapply(counts, is.numeric, logical(1))]
    if (length(non_num) > 0) {
      rlang::abort(
        c(
          "`counts` data frame contains non-numeric columns.",
          i = paste0("Non-numeric columns found: ", paste(non_num, collapse = ", ")),
          i = "All columns must be numeric read counts.",
          i = "If one column is a sample ID, move it to `rownames(counts)` first.",
          i = 'Example: `rownames(counts) <- counts$SampleID; counts$SampleID <- NULL`'
        ),
        call = call
      )
    }
    counts <- as.matrix(counts)
  }

  # ── Dimension check ─────────────────────────────────────────────────────────
  if (nrow(counts) < 3) {
    rlang::abort(
      c(
        paste0("`counts` has only ", nrow(counts), " row(s); at least 3 samples are required."),
        i = "The DMM cannot meaningfully partition community structure with fewer than 3 samples.",
        i = "Each row of `counts` should be one sample (station, replicate, etc.)."
      ),
      call = call
    )
  }

  if (ncol(counts) < 2) {
    rlang::abort(
      c(
        paste0("`counts` has only ", ncol(counts), " column(s); at least 2 taxa are required."),
        i = "Each column of `counts` should be one taxon or ASV."
      ),
      call = call
    )
  }

  # ── Value checks ─────────────────────────────────────────────────────────────
  if (any(counts < 0, na.rm = TRUE)) {
    rlang::abort(
      c(
        "`counts` contains negative values.",
        i = "Read counts must be zero or positive integers.",
        i = paste0("Negative values found at: ",
                   paste(which(counts < 0, arr.ind = TRUE)[1:min(3, sum(counts < 0)), ,
                               drop = FALSE] |>
                           apply(1, function(r) paste0("[", r[1], ",", r[2], "]")),
                         collapse = ", "))
      ),
      call = call
    )
  }

  if (any(is.na(counts))) {
    n_na <- sum(is.na(counts))
    rlang::abort(
      c(
        paste0("`counts` contains ", n_na, " NA value(s)."),
        i = "All count values must be non-missing.",
        i = "Replace NAs with 0 if a taxon was simply not detected: `counts[is.na(counts)] <- 0`"
      ),
      call = call
    )
  }

  # ── Zero-row samples ─────────────────────────────────────────────────────────
  row_sums <- rowSums(counts)
  empty_samples <- which(row_sums == 0)
  if (length(empty_samples) > 0) {
    sample_names <- if (!is.null(rownames(counts))) rownames(counts)[empty_samples] else
      paste0("row ", empty_samples)
    rlang::abort(
      c(
        paste0(length(empty_samples), " sample(s) have zero total reads and cannot be modelled."),
        i = paste0("Empty sample(s): ", paste(sample_names[seq_len(min(5, length(sample_names)))],
                                              collapse = ", ")),
        i = "Remove these rows before fitting: `counts <- counts[rowSums(counts) > 0, ]`"
      ),
      call = call
    )
  }

  # ── Warn about suspiciously low read counts ───────────────────────────────────
  low_samples <- which(row_sums < 100)
  if (length(low_samples) > 0) {
    rlang::warn(
      c(
        paste0(length(low_samples), " sample(s) have fewer than 100 total reads."),
        i = "Very low read counts produce unreliable composition estimates.",
        i = "Consider filtering: `counts <- counts[rowSums(counts) >= 100, ]`"
      )
    )
  }

  # ── All-zero taxa ─────────────────────────────────────────────────────────────
  col_sums <- colSums(counts)
  zero_taxa <- which(col_sums == 0)
  if (length(zero_taxa) > 0) {
    taxa_names <- if (!is.null(colnames(counts))) colnames(counts)[zero_taxa] else
      paste0("column ", zero_taxa)
    rlang::warn(
      c(
        paste0(length(zero_taxa), " taxon/taxa column(s) have zero reads across all samples."),
        i = paste0("Zero-count column(s): ",
                   paste(taxa_names[seq_len(min(5, length(taxa_names)))], collapse = ", ")),
        i = "These are dropped automatically before fitting.",
        i = "To suppress this warning, remove them yourself: `counts <- counts[, colSums(counts) > 0]`"
      )
    )
  }

  # ── Storage mode ──────────────────────────────────────────────────────────────
  if (!is.integer(counts)) {
    if (!all(counts == floor(counts), na.rm = TRUE)) {
      rlang::warn(
        c(
          "`counts` contains non-integer values. These will be coerced to integers via `round()`.",
          i = "Stan requires integer count data.",
          i = "If your data are already integers stored as doubles, this is expected and safe.",
          i = "If your data are normalized read proportions, you need raw counts."
        )
      )
    }
    storage.mode(counts) <- "integer"
  }

  counts
}

#' @keywords internal
validate_covariates <- function(covariates, counts, scale_covariates, call = rlang::caller_env()) {
  N <- nrow(counts)

  # ── NULL means intercept-only model ──────────────────────────────────────────
  if (is.null(covariates)) {
    return(matrix(numeric(0), nrow = N, ncol = 0))
  }

  # ── Type ──────────────────────────────────────────────────────────────────────
  if (!is.matrix(covariates) && !is.data.frame(covariates)) {
    rlang::abort(
      c(
        "`covariates` must be a numeric matrix or data frame (or NULL for an intercept-only model).",
        i = paste0("You supplied an object of class: ", paste(class(covariates), collapse = ", ")),
        i = "Expected format: rows = samples (same order as `counts`), columns = covariates.",
        i = "Example: a data frame with columns `depth` and `latitude`, one row per sample."
      ),
      call = call
    )
  }

  if (is.data.frame(covariates)) {
    non_num <- names(covariates)[!vapply(covariates, is.numeric, logical(1))]
    if (length(non_num) > 0) {
      rlang::abort(
        c(
          "`covariates` data frame contains non-numeric columns.",
          i = paste0("Non-numeric columns: ", paste(non_num, collapse = ", ")),
          i = "All covariate columns must be numeric.",
          i = "For categorical covariates, create indicator (dummy) columns manually."
        ),
        call = call
      )
    }
    covariates <- as.matrix(covariates)
  }

  # ── Dimension alignment ───────────────────────────────────────────────────────
  if (nrow(covariates) != N) {
    rlang::abort(
      c(
        paste0("`covariates` has ", nrow(covariates), " rows but `counts` has ", N, " rows."),
        i = "Both must have one row per sample, in the same order.",
        i = "Check that you haven't filtered one object without filtering the other."
      ),
      call = call
    )
  }

  # ── Missing values ────────────────────────────────────────────────────────────
  if (any(is.na(covariates))) {
    n_na <- sum(is.na(covariates))
    na_cols <- colnames(covariates)[apply(covariates, 2, anyNA)]
    rlang::abort(
      c(
        paste0("`covariates` contains ", n_na, " NA value(s)."),
        i = paste0("Columns with NAs: ", paste(na_cols, collapse = ", ")),
        i = "Stan cannot handle missing covariate values.",
        i = "Options: (1) impute NAs, (2) drop rows with NAs from both `counts` and `covariates`,",
        i = "         (3) set `covariates = NULL` to fit an intercept-only model."
      ),
      call = call
    )
  }

  # ── Near-zero variance covariates ─────────────────────────────────────────────
  col_vars <- apply(covariates, 2, var)
  zero_var <- which(col_vars < .Machine$double.eps * 100)
  if (length(zero_var) > 0) {
    zero_names <- if (!is.null(colnames(covariates))) colnames(covariates)[zero_var] else
      paste0("column ", zero_var)
    rlang::abort(
      c(
        paste0("Covariate(s) have (near-)zero variance and cannot be used: ",
               paste(zero_names, collapse = ", ")),
        i = "A constant covariate provides no information about community membership.",
        i = "Remove it from `covariates` or check that your data are correct."
      ),
      call = call
    )
  }

  # ── Scaling ───────────────────────────────────────────────────────────────────
  if (scale_covariates) {
    scaled <- scale(covariates)
    attr(scaled, "scale_center") <- attr(scaled, "scaled:center")
    attr(scaled, "scale_scale")  <- attr(scaled, "scaled:scale")
    attr(scaled, "scaled:center") <- NULL
    attr(scaled, "scaled:scale")  <- NULL
    covariates <- scaled
  } else {
    # Warn if covariates look unscaled and the user said not to scale
    col_ranges <- apply(covariates, 2, function(x) diff(range(x)))
    big_range  <- col_ranges > 100
    if (any(big_range)) {
      rlang::warn(
        c(
          "Some covariates have a large numeric range and `scale_covariates = FALSE`.",
          i = paste0("Wide-range columns: ",
                     paste(colnames(covariates)[big_range], collapse = ", ")),
          i = "Unscaled covariates can cause slow mixing and poor convergence.",
          i = "Strongly recommended: set `scale_covariates = TRUE` (the default)."
        )
      )
    }
  }

  covariates
}

#' @keywords internal
validate_K <- function(K, N, call = rlang::caller_env()) {
  if (!is.numeric(K) || length(K) != 1 || K != round(K)) {
    rlang::abort(
      c(
        "`K` must be a single positive integer (the number of communities to fit).",
        i = paste0("You supplied: K = ", paste(K, collapse = ", ")),
        i = "Example: `K = 3` fits a model with 3 latent communities."
      ),
      call = call
    )
  }
  K <- as.integer(K)
  if (K < 2) {
    rlang::abort(
      c(
        "`K` must be at least 2.",
        i = "K = 1 is a single-community model (no mixture), which is not meaningful here.",
        i = "Start with K = 2 and increase if needed."
      ),
      call = call
    )
  }
  if (K >= N) {
    rlang::abort(
      c(
        paste0("`K` (", K, ") must be less than the number of samples N (", N, ")."),
        i = "You cannot have more communities than samples.",
        i = paste0("With ", N, " samples, maximum K is ", N - 1, ".",
                   " Typical values are K = 2 to 6.")
      ),
      call = call
    )
  }
  if (K > 10) {
    rlang::warn(
      c(
        paste0("K = ", K, " is unusually large."),
        i = "Models with many communities are prone to label switching and slow convergence.",
        i = "Consider fitting K = 2 through K = 6 and using LOO-CV to select."
      )
    )
  }
  K
}

#' @keywords internal
validate_covariate_names <- function(covariates, covariate_names, call = rlang::caller_env()) {
  if (is.null(covariate_names)) return(invisible(NULL))

  current_names <- colnames(covariates)
  if (is.null(current_names)) {
    rlang::warn(
      c(
        "`covariates` has no column names. The `covariate_names` argument will be used as labels.",
        i = "To avoid this, set colnames on your covariates matrix: `colnames(covariates) <- c('depth', 'lat')`"
      )
    )
    return(invisible(NULL))
  }

  bad <- setdiff(covariate_names, current_names)
  if (length(bad) > 0) {
    # Suggest near-matches (typos)
    suggestions <- vapply(bad, function(b) {
      dists   <- adist(b, current_names, ignore.case = TRUE)
      nearest <- current_names[which.min(dists)]
      if (min(dists) <= 2) paste0("Did you mean '", nearest, "'?") else ""
    }, character(1))
    sugg_msgs <- suggestions[nchar(suggestions) > 0]

    rlang::abort(
      c(
        paste0("Covariate name(s) not found in `covariates`: ",
               paste(bad, collapse = ", ")),
        i = paste0("Available names: ", paste(current_names, collapse = ", ")),
        if (length(sugg_msgs) > 0) i = paste(sugg_msgs, collapse = "; ") else NULL
      ),
      call = call
    )
  }
}

#' @keywords internal
check_stan_fit_object <- function(x, fn_name = "this function", call = rlang::caller_env()) {
  if (!inherits(x, "edna_dmm_fit")) {
    rlang::abort(
      c(
        paste0("`fit` must be an `edna_dmm_fit` object (the output of `eDNA_dmm()`)."),
        i = paste0("You supplied an object of class: ", paste(class(x), collapse = ", ")),
        i = paste0("Run `fit <- eDNA_dmm(counts, covariates, K = 2)` first, then pass `fit` to `", fn_name, "`.")
      ),
      call = call
    )
  }
}

#' @keywords internal
make_community_colors <- function(k) {
  hues <- seq(15, 375, length.out = k + 1)[seq_len(k)]
  setNames(grDevices::hcl(h = hues, c = 80, l = 55),
           paste0("Community ", seq_len(k)))
}

#' @keywords internal
make_taxa_colors <- function(taxa) {
  taxa_sorted <- sort(taxa[taxa != "Other"])
  n           <- length(taxa_sorted)
  if (n == 0) return(c("Other" = "grey70"))
  n_shades   <- ceiling(n / 7)
  base_hues  <- seq(15, 375, length.out = 8)[seq_len(7)]
  lum_vals   <- seq(75, 40, length.out = n_shades)
  color_grid <- outer(lum_vals, base_hues,
                      function(l, h) grDevices::hcl(h = h, c = 80, l = l))
  c(setNames(as.vector(color_grid)[seq_len(n)], taxa_sorted),
    "Other" = "grey70")
}
