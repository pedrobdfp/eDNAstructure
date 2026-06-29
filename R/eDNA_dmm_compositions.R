# =============================================================================
# eDNA_dmm_compositions() — posterior community composition bar plots
# =============================================================================

#' Plot posterior mean species composition for each community
#'
#' @description
#' Produces a stacked bar plot of the **posterior mean species composition**
#' for each fitted community (the pi matrix). One bar per community, colored
#' by species using the same structured palette as [plot_true_compositions()],
#' so species colors are directly comparable between the two plots.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param top_n Number of most abundant taxa to show individually; the rest
#'   are collapsed into `"Other"`. Default `20`.
#' @param base_size Base font size. Default `13`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default auto-generated.
#' @param legend_position Legend position. Default `"right"`.
#' @param bar_width Bar width. Default `0.7`.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [plot_true_compositions()], [eDNA_dmm_structure()]
#' @export
eDNA_dmm_compositions <- function(
    fit,
    top_n            = 20,
    base_size        = 13,
    title            = NULL,
    subtitle         = NULL,
    legend_position  = "right",
    bar_width        = 0.7
) {
  check_stan_fit_object(fit, "eDNA_dmm_compositions")
  
  K          <- fit$K
  pi_mean    <- fit$pi_mean          # K x S matrix, posterior mean compositions
  taxa_names <- colnames(pi_mean)
  
  if (is.null(taxa_names))
    taxa_names <- paste0("Sp_", seq_len(ncol(pi_mean)))
  
  # ── Top N taxa ─────────────────────────────────────────────────────────────
  # Use total weight across communities to pick top taxa — same logic as
  # plot_true_compositions() so colors stay consistent.
  taxon_totals <- colSums(pi_mean)
  top_taxa     <- names(sort(taxon_totals, decreasing = TRUE))[seq_len(min(top_n, ncol(pi_mean)))]
  other_taxa   <- setdiff(taxa_names, top_taxa)
  
  # ── Structured color palette — identical to plot_true_compositions() ───────
  named_taxa <- sort(top_taxa)
  n_named    <- length(named_taxa)
  n_shades   <- ceiling(n_named / 7)
  base_hues  <- seq(15, 375, length.out = 8)[1:7]
  lum_vals   <- seq(75, 40, length.out = n_shades)
  color_grid <- outer(lum_vals, base_hues,
                      function(l, h) grDevices::hcl(h = h, c = 80, l = l))
  tax_colors <- c(
    stats::setNames(as.vector(color_grid)[seq_len(n_named)], named_taxa),
    if (length(other_taxa) > 0) c(Other = "grey70") else NULL
  )
  
  taxon_levels <- c(sort(top_taxa), if (length(other_taxa) > 0) "Other")
  
  # ── Build long data frame ──────────────────────────────────────────────────
  comp_df <- as.data.frame(pi_mean)
  colnames(comp_df) <- taxa_names
  comp_df$Community <- as.character(seq_len(K))
  
  if (length(other_taxa) > 0) {
    comp_df$Other <- rowSums(comp_df[, other_taxa, drop = FALSE])
  }
  
  plot_long <- tidyr::pivot_longer(
    comp_df,
    cols      = dplyr::all_of(c(top_taxa, if (length(other_taxa) > 0) "Other")),
    names_to  = "taxon",
    values_to = "proportion"
  )
  plot_long$taxon     <- factor(plot_long$taxon, levels = taxon_levels)
  plot_long$Community <- factor(plot_long$Community, levels = as.character(seq_len(K)))
  
  # ── Titles ─────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("Posterior community compositions  (K = %d)", K)
  subtitle_str <- subtitle %||% sprintf(
    "Posterior mean π  |  top %d taxa  |  colors match observed composition plot",
    min(top_n, ncol(pi_mean))
  )
  
  # ── Plot ───────────────────────────────────────────────────────────────────
  ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = .data$Community, y = .data$proportion, fill = .data$taxon)
  ) +
    ggplot2::geom_bar(stat = "identity", width = bar_width, color = NA) +
    ggplot2::scale_fill_manual(values = tax_colors, name = "Taxon", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(
      x        = "Community",
      y        = "Posterior mean proportion",
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = legend_position,
      strip.background = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend())
}