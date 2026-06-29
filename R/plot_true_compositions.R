# =============================================================================
# plot_true_compositions() — visualize observed species composition per sample
# =============================================================================

#' Plot observed species composition per sample (pre-fit diagnostic)
#'
#' @description
#' Produces a stacked bar plot of **observed species frequencies** — one bar
#' per sample, colored by species. This is a pre-fitting diagnostic that shows
#' the raw community signal in your data before any model is applied. It uses
#' the same visual format as [eDNA_dmm_structure()] so the two plots can be
#' compared directly: the top panel shows what species are present, the bottom
#' shows what community they were assigned to.
#'
#' When `facet_row_var` is supplied, one sub-plot is built per row-level and
#' they are assembled vertically with cowplot. This mirrors the original
#' make_bar_figure() approach that gives each depth stratum its own
#' proportional-width panels — the only architecture that produces gapless
#' bars with correct proportional panel widths via space = "free_x".
#'
#' @param counts A numeric matrix or data frame of read counts (samples × taxa).
#' @param metadata An optional data frame with one row per sample.
#' @param sample_id_col Column in `metadata` matching `rownames(counts)`.
#'   Default `"sample_id"`.
#' @param facet_var Column for column facets (e.g. `"year"`). Default `NULL`.
#' @param facet_row_var Column for row facets (e.g. `"depth_bin"`). When
#'   supplied, one sub-plot is built per level and stacked with cowplot.
#'   Default `NULL`.
#' @param sort_var Column to sort samples within panels. Default `NULL`.
#' @param top_n Number of most abundant taxa shown individually; rest collapsed
#'   to `"Other"`. Default `20`.
#' @param base_size Base font size. Default `11`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default auto-generated.
#' @param legend_position Legend position. Default `"bottom"`.
#' @param vline_var Numeric column for vertical reference line. Default `NULL`.
#' @param vline_value Threshold value in `vline_var` space. Default `NULL`.
#' @param vline_color Line color. Default `"black"`.
#' @param vline_linetype Line type. Default `"dashed"`.
#' @param vline_linewidth Line width. Default `0.7`.
#'
#' @return A ggplot2 object (no `facet_row_var`) or a cowplot grid object.
#'
#' @seealso [eDNA_dmm_structure()], [simulate_eDNA_survey()]
#' @export
plot_true_compositions <- function(
    counts,
    metadata         = NULL,
    sample_id_col    = "sample_id",
    facet_var        = NULL,
    facet_row_var    = NULL,
    sort_var         = NULL,
    top_n            = 20,
    base_size        = 11,
    title            = NULL,
    subtitle         = NULL,
    show_legend      = FALSE,
    legend_position  = "bottom",
    vline_var        = NULL,
    vline_value      = NULL,
    vline_color      = "black",
    vline_linetype   = "dashed",
    vline_linewidth  = 0.7
) {
  # ── Coerce to matrix ───────────────────────────────────────────────────────
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  N          <- nrow(counts)
  S          <- ncol(counts)
  taxa_names <- if (!is.null(colnames(counts))) colnames(counts) else paste0("Sp_", seq_len(S))
  sample_ids <- if (!is.null(rownames(counts))) rownames(counts) else paste0("S", seq_len(N))
  
  # ── Relative frequencies ──────────────────────────────────────────────────
  row_tots <- rowSums(counts)
  row_tots[row_tots == 0] <- 1
  freq_mat <- counts / row_tots
  
  # ── Top N taxa ─────────────────────────────────────────────────────────────
  mean_freq  <- colMeans(freq_mat)
  top_taxa   <- names(sort(mean_freq, decreasing = TRUE))[seq_len(min(top_n, S))]
  other_taxa <- setdiff(taxa_names, top_taxa)
  
  # ── Long data frame ────────────────────────────────────────────────────────
  plot_df <- as.data.frame(freq_mat)
  colnames(plot_df) <- taxa_names
  plot_df$sample_id <- sample_ids
  
  if (length(other_taxa) > 0) {
    plot_df$Other <- rowSums(plot_df[, other_taxa, drop = FALSE])
    plot_df       <- plot_df[, c("sample_id", top_taxa, "Other")]
  } else {
    plot_df <- plot_df[, c("sample_id", top_taxa)]
  }
  
  # ── Merge metadata ─────────────────────────────────────────────────────────
  si <- plot_df
  if (!is.null(metadata)) {
    if (!sample_id_col %in% names(metadata))
      rlang::abort(paste0("Column '", sample_id_col, "' not found in `metadata`."))
    si <- merge(si, metadata, by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  }
  
  # ── Colors — structured palette matching make_structured_colors() ──────────
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
  
  value_cols   <- c(top_taxa, if (length(other_taxa) > 0) "Other")
  taxon_levels <- c(sort(top_taxa), if (length(other_taxa) > 0) "Other")
  
  # ── Titles ─────────────────────────────────────────────────────────────────
  title_str    <- title    %||% "Observed species composition"
  subtitle_str <- subtitle %||% sprintf("%d samples  |  top %d taxa shown individually",
                                        N, min(top_n, S))
  
  # ── Core single-panel builder ──────────────────────────────────────────────
  # Builds one ggplot for a given subset of si. This is called once per
  # facet_row_var level (or once for the whole dataset if no row facet).
  # x = reorder(sample_id, sort_var) makes x discrete and per-subset, which
  # is the only reliable way to get gapless bars with proportional panel widths.
  
  build_panel <- function(dat, row_label = NULL, show_legend = FALSE,
                          show_title = FALSE) {
    # Sort within this subset
    if (!is.null(sort_var) && sort_var %in% names(dat)) {
      dat <- dat[order(dat[[sort_var]]), ]
    }
    # Ordered factor — levels only from this subset so bars fill the panel
    dat$x_label <- factor(dat$sample_id, levels = unique(dat$sample_id))
    
    # Vline: find position within each facet_var group
    vline_df <- NULL
    if (!is.null(vline_var) && !is.null(vline_value) && vline_var %in% names(dat)) {
      if (!is.null(facet_var) && facet_var %in% names(dat)) {
        vline_df <- dat |>
          dplyr::select(dplyr::all_of(c(facet_var, vline_var, "x_label"))) |>
          dplyr::distinct() |>
          dplyr::group_by(dplyr::across(dplyr::all_of(facet_var))) |>
          dplyr::summarise(
            cape_x = as.integer(x_label[which.min(abs(.data[[vline_var]] - vline_value))]) + 0.5,
            .groups = "drop"
          )
      } else {
        closest <- which.min(abs(dat[[vline_var]] - vline_value))
        vline_df <- data.frame(cape_x = as.integer(dat$x_label[closest]) + 0.5)
      }
    }
    
    # Pivot long
    plot_long <- tidyr::pivot_longer(dat, cols = dplyr::all_of(value_cols),
                                     names_to = "taxon", values_to = "frequency")
    plot_long$taxon <- factor(plot_long$taxon, levels = taxon_levels)
    
    p <- ggplot2::ggplot(
      plot_long,
      ggplot2::aes(x = .data$x_label, y = .data$frequency, fill = .data$taxon)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = NA) +
      ggplot2::scale_fill_manual(values = tax_colors, name = NULL, drop = FALSE) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                                  breaks = c(0, 0.5, 1)) +
      ggplot2::labs(
        x     = NULL,
        y     = if (!is.null(row_label)) row_label else "Species frequency",
        title = NULL
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position  = if (show_legend) legend_position else "none",
        axis.text.x      = ggplot2::element_blank(),
        axis.ticks.x     = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(face = "bold"),
        panel.spacing.x  = ggplot2::unit(0.3, "lines"),
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(color = "grey40",
                                                 size = base_size - 2)
      )
    
    # Vline
    if (!is.null(vline_df)) {
      p <- p + ggplot2::geom_vline(
        data        = vline_df,
        ggplot2::aes(xintercept = .data$cape_x),
        color       = vline_color,
        linetype    = vline_linetype,
        linewidth   = vline_linewidth,
        inherit.aes = FALSE
      )
    }
    
    # Column facet
    if (!is.null(facet_var) && facet_var %in% names(dat)) {
      p <- p + ggplot2::facet_grid(
        reformulate(facet_var),
        scales = "free_x",
        space  = "free_x"
      )
    }
    
    p
  }
  
  # ── No row facet: single plot ──────────────────────────────────────────────
  if (is.null(facet_row_var)) {
    p <- build_panel(si, show_legend = show_legend, show_title = TRUE)
    p <- p + ggplot2::labs(title = title_str, subtitle = subtitle_str,
                           y = "Species frequency")
    return(p)
  }
  
  # ── Row facet: one panel per level, assembled with cowplot ─────────────────
  if (!facet_row_var %in% names(si))
    rlang::abort(paste0("`facet_row_var = '", facet_row_var, "'` not found."))
  
  row_levels <- if (is.factor(si[[facet_row_var]])) {
    levels(si[[facet_row_var]])
  } else {
    sort(unique(si[[facet_row_var]]))
  }
  
  panels <- lapply(seq_along(row_levels), function(i) {
    lv  <- row_levels[i]
    dat <- si[si[[facet_row_var]] == lv, ]
    if (nrow(dat) == 0) return(NULL)
    build_panel(
      dat,
      row_label   = as.character(lv),
      show_legend = FALSE,
      show_title  = FALSE
    )
  })
  panels <- Filter(Negate(is.null), panels)
  
  # Title row
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_str, fontface = "bold", size = base_size + 1,
                        x = 0.02, hjust = 0) +
    cowplot::draw_label(subtitle_str, size = base_size - 1, color = "grey40",
                        x = 0.02, y = 0.25, hjust = 0)
  
  # Assemble with or without legend
  stacked <- cowplot::plot_grid(plotlist = panels, ncol = 1,
                                labels = letters[seq_along(panels)],
                                label_size = 14)
  
  if (show_legend) {
    # Extract legend from a dummy plot using the requested position
    legend_plot <- build_panel(si, show_legend = TRUE) +
      ggplot2::theme(legend.position = legend_position) +
      ggplot2::guides(fill = ggplot2::guide_legend())
    
    if (legend_position == "right") {
      # Legend to the right of the stacked panels
      legend_grob <- cowplot::get_plot_component(legend_plot,
                                                 "guide-box-right",
                                                 return_all = TRUE)
      inner <- cowplot::plot_grid(stacked, legend_grob,
                                  nrow = 1, rel_widths = c(1, 0.2))
      cowplot::plot_grid(title_grob, inner,
                         ncol = 1, rel_heights = c(0.06, 1))
      
    } else if (legend_position == "left") {
      legend_grob <- cowplot::get_plot_component(legend_plot,
                                                 "guide-box-left",
                                                 return_all = TRUE)
      inner <- cowplot::plot_grid(legend_grob, stacked,
                                  nrow = 1, rel_widths = c(0.2, 1))
      cowplot::plot_grid(title_grob, inner,
                         ncol = 1, rel_heights = c(0.06, 1))
      
    } else {
      # bottom or top
      legend_grob <- cowplot::get_plot_component(legend_plot,
                                                 "guide-box-bottom",
                                                 return_all = TRUE)
      if (legend_position == "top") {
        cowplot::plot_grid(title_grob, legend_grob, stacked,
                           ncol = 1, rel_heights = c(0.06, 0.1, 1))
      } else {
        cowplot::plot_grid(title_grob, stacked, legend_grob,
                           ncol = 1, rel_heights = c(0.06, 1, 0.1))
      }
    }
  } else {
    cowplot::plot_grid(title_grob, stacked,
                       ncol = 1, rel_heights = c(0.06, 1))
  }
}

