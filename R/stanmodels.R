# =============================================================================
# Stan model loading: bundled precompiled object + automatic fallback
# =============================================================================
#
# This package does NOT use the rstantools precompiled-module mechanism
# (RCPP_MODULE / stanExports_*.cc baked into the package DLL). That mechanism
# is fragile on Windows: the module boot symbol can be dropped from the DLL at
# link time, breaking `library()` with a "no such symbol
# _rcpp_module_boot_*" error. There is therefore no Stan C++ in src/ and no
# configure script.
#
# Instead, the compiled Stan model is loaded by the following priority:
#
#   1. In-session cache (.eDNA_stanmodels_cache$dmm) — instant.
#   2. Bundled precompiled object (inst/stan/dmm_model.rds), if present AND
#      it loads cleanly on THIS machine. Generated once by the package author
#      via tools/precompile_stan.R. Instant for users on a matching
#      rstan/StanHeaders/OS.
#   3. User-level disk cache (R_user_dir cache), from a previous compile on
#      this machine — instant after the first compile.
#   4. Compile from inst/stan/dmm.stan via rstan::stan_model(), then cache to
#      the user-level cache. ~30-60s, happens at most once per machine.
#
# A precompiled rstan model embeds a pointer to compiled C++ built for a
# specific rstan/StanHeaders/OS/arch. Loading such an object on a mismatched
# machine yields an invalid pointer. .stanmodel_is_usable() detects that and
# triggers the compile-from-source fallback, so the package works everywhere.
# =============================================================================

.eDNA_stanmodels_cache <- new.env(parent = emptyenv())

# Internal: path to the per-machine user cache for the compiled model.
#' @keywords internal
.dmm_cache_path <- function() {
  file.path(tools::R_user_dir("eDNAstructure", which = "cache"), "dmm_model.rds")
}

# Internal: path to the precompiled object shipped inside the package (if any).
#' @keywords internal
.dmm_bundled_path <- function() {
  p <- system.file("stan", "dmm_model.rds", package = "eDNAstructure")
  if (nzchar(p)) p else ""
}

# Internal: path to the Stan source shipped inside the package.
#' @keywords internal
.dmm_stan_source <- function() {
  system.file("stan", "dmm.stan", package = "eDNAstructure")
}

# Internal: verify a loaded object is a usable stanmodel on THIS machine.
# A stanmodel whose compiled C++ pointer is stale (wrong rstan/OS/arch) will
# typically fail when its module is touched; we probe it defensively.
#' @keywords internal
.stanmodel_is_usable <- function(model) {
  if (is.null(model)) return(FALSE)
  if (!methods::is(model, "stanmodel")) return(FALSE)
  ok <- tryCatch({
    # Touching the model name / DSO forces rstan to validate the backing
    # module. On a mismatched machine this errors; we catch and fall back.
    invisible(model@model_name)
    dso <- methods::slot(model, "dso")
    # If a DSO is present, confirm it can report its loaded state without error.
    if (!is.null(dso)) invisible(dso@dso_saved)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

# Internal: read an .rds and return a usable stanmodel, or NULL.
#' @keywords internal
.try_load_stanmodel <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  model <- tryCatch(readRDS(path), error = function(e) NULL)
  if (.stanmodel_is_usable(model)) model else NULL
}

#' Internal accessor for the compiled DMM Stan model
#'
#' Returns a compiled `stanmodel`, loading from (in order) the in-session
#' cache, the bundled precompiled object, the per-machine disk cache, or by
#' compiling from source. See the top of `R/stanmodels.R` for the rationale.
#'
#' @return A [rstan::stanmodel-class] object.
#' @keywords internal
.get_dmm_stanmodel <- function() {
  # 1. In-session cache.
  if (!is.null(.eDNA_stanmodels_cache$dmm)) {
    return(.eDNA_stanmodels_cache$dmm)
  }

  # 2. Bundled precompiled object (author-generated), if usable here.
  model <- .try_load_stanmodel(.dmm_bundled_path())

  # 3. Per-machine disk cache from a previous compile.
  if (is.null(model)) {
    model <- .try_load_stanmodel(.dmm_cache_path())
  }

  # 4. Compile from source and cache.
  if (is.null(model)) {
    stan_file <- .dmm_stan_source()
    if (!nzchar(stan_file) || !file.exists(stan_file)) {
      rlang::abort(
        c(
          "Could not find the Stan model source (dmm.stan) in the installed package.",
          i = "Try reinstalling: remotes::install_github('pedrobdfp/eDNA_structure', force = TRUE)"
        )
      )
    }

    message("Compiling the DMM Stan model (first use on this machine).")
    message("This takes ~30-60 seconds and is cached, so it only happens once.")

    model <- rstan::stan_model(file = stan_file, model_name = "dmm")

    cache_path <- .dmm_cache_path()
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      saveRDS(model, cache_path),
      error = function(e) {
        message(
          "Note: could not cache the compiled model to disk (",
          conditionMessage(e), "). It will be recompiled next session."
        )
      }
    )
  }

  .eDNA_stanmodels_cache$dmm <- model
  model
}

#' Clear the cached compiled Stan model
#'
#' @description
#' Removes the in-session and per-machine disk caches of the compiled DMM Stan
#' model, forcing a fresh compile on next use. Useful after updating the
#' package, after upgrading `rstan`/`StanHeaders`, or if the cached model ever
#' becomes stale or corrupted.
#'
#' This does not remove the precompiled object bundled inside the package; it
#' only clears caches created on your machine.
#'
#' @return Invisibly `NULL`. Called for its side effect.
#'
#' @examples
#' \dontrun{
#' eDNA_clear_stan_cache()
#' }
#'
#' @export
eDNA_clear_stan_cache <- function() {
  cache_path <- .dmm_cache_path()
  if (file.exists(cache_path)) {
    file.remove(cache_path)
    message("Cleared cached Stan model: ", cache_path)
  } else {
    message("No per-machine cached Stan model found.")
  }
  if (exists("dmm", envir = .eDNA_stanmodels_cache, inherits = FALSE)) {
    rm("dmm", envir = .eDNA_stanmodels_cache)
  }
  invisible(NULL)
}
