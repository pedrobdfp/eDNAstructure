# =============================================================================
# eDNA_dmm_nmds() — NMDS ordination colored by community assignment
# =============================================================================

#' NMDS ordination of community composition colored by community assignment
#'
#' @description
#' Runs a non-metric multidimensional scaling (NMDS) ordination on the sample
#' community composition and plots the result, with points colored by their
#' MAP community assignment from the DMM. Point size reflects assignment
#' certainty: larger points are more confidently assigned.
#'
#' NMDS is performed using [vegan::metaMDS()] on Bray-Curtis dissimilarities
#' computed from the eDNA index transformation. The eDNA index (Kelly et al.
#' 2019) normalizes read counts to relative abundances within each sample and
#' then divides by the across-sample maximum for each taxon, which reduces
#' the influence of dominant taxa on distance calculations.
#'
#' The function returns a [ggplot2::ggplot()] object that can be further
#' customized.
#'
#' @section Axes displayed:
#' `nmds_axes` controls which two NMDS axes to display. For `k = 2`
#' dimensions (the default), only axes 1 and 2 exist. For `k = 3`, you can
#' pass `nmds_axes = c(1, 3)` to view NMDS1 vs NMDS3. Use multiple calls to
#' inspect all axis pairs when fitting in 3D.
#'
#' @section NMDS stress:
#' The stress value (a measure of how well 2D positions represent the true
#' dissimilarities) is shown in the subtitle. Rough guidelines:
#' - Stress < 0.05: excellent representation
#' - Stress 0.05–0.10: good
#' - Stress 0.10–0.20: adequate, but interpret with care
#' - Stress > 0.20: poor; consider `k = 3` and examining multiple axis pairs
#'
#' @section Ellipses:
#' If `show_ellipse = TRUE` (the default), a 95% confidence ellipse is drawn
#' for each community using [ggplot2::stat_ellipse()]. Ellipses require at
#' least 3 samples per community to be drawn; communities with fewer samples
#' are silently skipped.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param k A positive integer: the number of NMDS dimensions. Default is `2`.
#'   Increasing to `3` can reduce stress but requires inspecting multiple
#'   axis pairs.
#' @param nmds_axes An integer vector of length 2: which NMDS axes to plot.
#'   Default is `c(1, 2)` (NMDS1 vs NMDS2). Change to `c(1, 3)` or `c(2, 3)`
#'   to view other axis combinations when `k >= 3`.
#' @param distance A string: the dissimilarity index passed to
#'   [vegan::vegdist()]. Default is `"bray"` (Bray-Curtis). Other options
#'   include `"jaccard"`, `"euclidean"`, `"kulczynski"`. See
#'   [vegan::vegdist()] for all options.
#' @param use_edna_index Logical. If `TRUE` (the default), apply the eDNA
#'   index transformation before computing distances. If `FALSE`, raw
#'   relative abundances (row sums normalized to 1) are used.
#' @param trymax A positive integer: the maximum number of random starts for
#'   NMDS. Default is `100`. More starts reduce the risk of local optima but
#'   increase computation time.
#' @param seed A single integer: the random seed for NMDS. Default is `42`.
#'   Set for reproducibility.
#' @param show_ellipse Logical. If `TRUE` (the default), draw a 95% confidence
#'   ellipse per community.
#' @param ellipse_type A string: the type of ellipse passed to
#'   [ggplot2::stat_ellipse()]. Default is `"t"` (t-distribution based, more
#'   robust with few points). Use `"norm"` for a normal-distribution based
#'   ellipse.
#' @param community_colors A named character vector mapping community names to
#'   hex colors. Pass `NULL` (the default) for the automatic HCL palette.
#' @param size_range A numeric vector of length 2: the minimum and maximum
#'   point sizes, representing assignment certainties of 50% and 100%
#'   respectively. Default is `c(1.5, 5)`.
#' @param alpha A number in (0, 1]: the transparency of points. Default `0.85`.
#' @param base_size A positive number: the base font size for the plot theme.
#'   Default is `13`.
#' @param title A string: the plot title. Default is auto-generated.
#' @param subtitle A string: the plot subtitle (shown below the title, includes
#'   stress value by default). Default is auto-generated.
#' @param legend_position A string: legend placement. One of `"right"`
#'   (default), `"bottom"`, `"left"`, `"top"`, or `"none"`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`plot`}{A [ggplot2::ggplot()] object showing the NMDS ordination.}
#'   \item{`nmds`}{The [vegan::metaMDS()] result object, for advanced use
#'     (e.g., adding species scores with `vegan::scores()`).}
#' }
#'
#' @seealso [eDNA_dmm()], [eDNA_dmm_structure()], [eDNA_dmm_beta()],
#'   [vegan::metaMDS()]
#'
#' @examples
#' \dontrun{
#' data <- get_example_data()
#' fit  <- eDNA_dmm(data$counts, data$covariates, K = 3)
#'
#' # Default NMDS (2D, Bray-Curtis, eDNA index)
#' result <- eDNA_dmm_nmds(fit)
#' result$plot
#'
#' # Access stress value
#' message("NMDS stress: ", round(result$nmds$stress, 3))
#'
#' # 3D NMDS, view axes 1 and 3
#' result3d <- eDNA_dmm_nmds(fit, k = 3, nmds_axes = c(1, 3))
#' result3d$plot
#'
#' # Custom colors and no ellipses
#' eDNA_dmm_nmds(
#'   fit,
#'   show_ellipse     = FALSE,
#'   community_colors = c("Community 1" = "#E63946",
#'                        "Community 2" = "#1D3557",
#'                        "Community 3" = "#A8DADC")
#' )$plot
#'
#' # Save the plot
#' ggplot2::ggsave("nmds.png", result$plot, width = 7, height = 6, dpi = 300)
#' }
#'
#' @export
eDNA_dmm_nmds <- function(
    fit,
    k                = 2,
    nmds_axes        = c(1, 2),
    distance         = "bray",
    use_edna_index   = TRUE,
    trymax           = 100,
    seed             = 42,
    show_ellipse     = TRUE,
    ellipse_type     = "t",
    community_colors = NULL,
    size_range       = c(1.5, 5),
    alpha            = 0.85,
    base_size        = 13,
    title            = NULL,
    subtitle         = NULL,
    legend_position  = "right"
) {
  check_stan_fit_object(fit, "eDNA_dmm_nmds")

  K  <- fit$K
  si <- fit$sample_info

  # ── Argument checks ───────────────────────────────────────────────────────────
  if (!is.numeric(k) || k < 2 || k != round(k)) {
    rlang::abort("`k` must be a positive integer >= 2 (the number of NMDS dimensions).")
  }
  if (length(nmds_axes) != 2 || !all(nmds_axes %in% seq_len(k))) {
    rlang::abort(
      c(
        paste0("`nmds_axes` must be a length-2 vector of axis indices within 1:", k, "."),
        i = paste0("With k=", k, " dimensions, valid axes are 1 through ", k, "."),
        i = "Example: `nmds_axes = c(1, 2)` for NMDS1 vs NMDS2."
      )
    )
  }

  valid_distances <- c("bray", "jaccard", "euclidean", "kulczynski", "canberra",
                       "clark", "gower", "horn", "morisita", "raup")
  if (!distance %in% valid_distances) {
    rlang::warn(
      c(
        paste0("`distance = '", distance, "'` is not one of the commonly tested options."),
        i = paste0("Commonly used: ", paste(valid_distances, collapse = ", ")),
        i = "Passing to vegan::vegdist() anyway."
      )
    )
  }

  # ── Community colors ──────────────────────────────────────────────────────────
  if (is.null(community_colors)) {
    community_colors <- make_community_colors(K)
  }

  # ── eDNA index transform or relative abundance ────────────────────────────────
  counts_mat <- fit$counts
  if (use_edna_index) {
    rel      <- counts_mat / rowSums(counts_mat)
    col_max  <- apply(rel, 2, max)
    col_max[col_max == 0] <- 1
    ord_mat  <- sweep(rel, 2, col_max, "/")
  } else {
    ord_mat  <- counts_mat / rowSums(counts_mat)
  }

  # ── NMDS ──────────────────────────────────────────────────────────────────────
  dist_mat <- vegan::vegdist(ord_mat, method = distance)
  set.seed(seed)
  nmds_obj <- vegan::metaMDS(
    ord_mat,
    distance = distance,
    k        = k,
    trymax   = trymax,
    trace    = FALSE
  )

  stress <- nmds_obj$stress

  # ── Scores ────────────────────────────────────────────────────────────────────
  sc        <- vegan::scores(nmds_obj, display = "sites")
  ax_labels <- paste0("NMDS", nmds_axes)
  plot_df   <- data.frame(
    NMDS_x = sc[, nmds_axes[1]],
    NMDS_y = sc[, nmds_axes[2]]
  )
  plot_df$z_hat                <- si$z_hat
  plot_df$assignment_certainty <- si$assignment_certainty

  # ── Stress quality message ────────────────────────────────────────────────────
  stress_label <- dplyr::case_when(
    stress < 0.05  ~ "excellent",
    stress < 0.10  ~ "good",
    stress < 0.20  ~ "adequate",
    TRUE           ~ "poor — consider k=3"
  )

  # ── Build plot ────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("NMDS  (K = %d communities)", K)
  subtitle_str <- subtitle %||% sprintf(
    "Bray-Curtis %s | Stress = %.3f (%s) | N = %d samples",
    if (use_edna_index) "(eDNA index)" else "",
    stress, stress_label, nrow(plot_df)
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$NMDS_x, y = .data$NMDS_y, color = .data$z_hat)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$assignment_certainty),
      alpha = alpha
    ) +
    ggplot2::scale_color_manual(values = community_colors, name = "Community") +
    ggplot2::scale_size_continuous(
      name   = "Assignment\ncertainty",
      range  = size_range,
      breaks = c(0.5, 0.75, 1.0),
      labels = c("50%", "75%", "100%")
    ) +
    ggplot2::labs(
      x        = ax_labels[1],
      y        = ax_labels[2],
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = legend_position,
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    )

  # ── Ellipses ──────────────────────────────────────────────────────────────────
  if (show_ellipse) {
    # Only draw ellipses for communities with >= 3 samples
    comm_counts <- table(plot_df$z_hat)
    comms_with_ellipse <- names(comm_counts)[comm_counts >= 3]
    if (length(comms_with_ellipse) > 0) {
      ellipse_df <- plot_df[plot_df$z_hat %in% comms_with_ellipse, ]
      p <- p +
        ggplot2::stat_ellipse(
          data     = ellipse_df,
          ggplot2::aes(group = .data$z_hat),
          type     = ellipse_type,
          level    = 0.95,
          linewidth = 0.6,
          linetype = "dashed"
        )
    }
    skipped <- setdiff(levels(plot_df$z_hat), comms_with_ellipse)
    if (length(skipped) > 0) {
      rlang::inform(
        paste0("Ellipses skipped for communities with < 3 samples: ",
               paste(skipped, collapse = ", "))
      )
    }
  }

  list(plot = p, nmds = nmds_obj)
}
