# =============================================================================
# eDNA_dmm_structure() — STRUCTURE-like community assignment bar plots
# =============================================================================

#' Plot posterior community assignments as STRUCTURE-like bar charts
#'
#' @description
#' Produces a stacked bar plot where each vertical bar is one sample, and the
#' bar is divided into colored segments whose heights represent the posterior
#' probability that the sample belongs to each community.
#'
#' When `facet_row_var` is supplied, one sub-plot is built per row-level and
#' assembled vertically with cowplot — the only architecture that produces
#' gapless bars with correct proportional panel widths via space = "free_x".
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param metadata An optional data frame with one row per sample.
#' @param sample_id_col Column in `metadata` matching sample IDs. Default `"sample_id"`.
#' @param facet_var Column for column facets (e.g. `"year"`). Default `NULL`.
#' @param facet_row_var Column for row facets (e.g. `"depth_bin"`). When
#'   supplied, one sub-plot is built per level and stacked with cowplot.
#'   Default `NULL`.
#' @param sort_var Column to sort samples within panels. Default `NULL`.
#' @param community_colors Named character vector mapping community names to
#'   colors. Default `NULL` (automatic HCL palette).
#' @param x_text Logical. Show sample labels on x-axis. Default `FALSE`.
#' @param base_size Base font size. Default `11`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default auto-generated.
#' @param show_legend Logical. Show legend. Default `TRUE`.
#' @param legend_position Legend position. Default `"bottom"`.
#' @param vline_var Numeric column for vertical reference line. Default `NULL`.
#' @param vline_value Threshold value in `vline_var` space. Default `NULL`.
#' @param vline_color Line color. Default `"black"`.
#' @param vline_linetype Line type. Default `"dashed"`.
#' @param vline_linewidth Line width. Default `0.7`.
#'
#' @return A ggplot2 object (no `facet_row_var`) or a cowplot grid object.
#'
#' @seealso [eDNA_dmm()], [eDNA_dmm_nmds()], [eDNA_dmm_beta()]
#' @export
eDNA_dmm_structure <- function(
    fit,
    metadata          = NULL,
    sample_id_col     = "sample_id",
    facet_var         = NULL,
    facet_row_var     = NULL,
    sort_var          = NULL,
    community_colors  = NULL,
    x_text            = FALSE,
    base_size         = 11,
    title             = NULL,
    subtitle          = NULL,
    show_legend       = TRUE,
    legend_position   = "bottom",
    vline_var         = NULL,
    vline_value       = NULL,
    vline_color       = "black",
    vline_linetype    = "dashed",
    vline_linewidth   = 0.7
) {
  check_stan_fit_object(fit, "eDNA_dmm_structure")
  
  K         <- fit$K
  si        <- fit$sample_info
  prob_cols <- paste0("prob_comm", seq_len(K))
  
  # ── Community colors ──────────────────────────────────────────────────────
  if (is.null(community_colors)) {
    community_colors <- make_community_colors(K)
  } else {
    expected      <- paste0("Community ", seq_len(K))
    missing_comms <- setdiff(expected, names(community_colors))
    if (length(missing_comms) > 0)
      rlang::abort(c("`community_colors` missing names.",
                     i = paste0("Missing: ", paste(missing_comms, collapse = ", "))))
  }
  
  # ── Merge metadata ────────────────────────────────────────────────────────
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata))
      rlang::abort("`metadata` must be a data frame.")
    if (!sample_id_col %in% names(metadata)) {
      near <- names(metadata)[which.min(adist(sample_id_col, names(metadata),
                                              ignore.case = TRUE))]
      rlang::abort(c(paste0("Column '", sample_id_col, "' not found in `metadata`."),
                     i = paste0("Did you mean '", near, "'?")))
    }
    si <- merge(si, metadata, by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  }
  
  # ── Validate ──────────────────────────────────────────────────────────────
  for (v in c(facet_var, facet_row_var, sort_var, vline_var)) {
    if (!is.null(v) && !v %in% names(si))
      rlang::abort(paste0("'", v, "' not found in merged data. Available: ",
                          paste(names(si), collapse = ", ")))
  }
  if (!is.null(vline_var) && !is.numeric(si[[vline_var]]))
    rlang::abort("`vline_var` must be numeric.")
  
  # ── Titles ────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("Posterior Community Assignments  (K = %d)", K)
  subtitle_str <- subtitle %||% sprintf(
    "%d samples  |  bar height = posterior membership probability", nrow(si))
  
  # ── Core single-panel builder ─────────────────────────────────────────────
  # Uses sample_id as discrete x within each subset — gapless bars guaranteed.
  build_panel <- function(dat, row_label = NULL, show_leg = FALSE) {
    if (!is.null(sort_var)) dat <- dat[order(dat[[sort_var]]), ]
    
    # Ordered factor from this subset only — key to gapless bars
    dat$x_label <- factor(dat$sample_id, levels = unique(dat$sample_id))
    
    # Vline position within facet_var panels
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
        closest  <- which.min(abs(dat[[vline_var]] - vline_value))
        vline_df <- data.frame(cape_x = as.integer(dat$x_label[closest]) + 0.5)
      }
    }
    
    plot_long <- tidyr::pivot_longer(dat, cols = dplyr::all_of(prob_cols),
                                     names_to = "community", values_to = "probability")
    plot_long$community <- factor(plot_long$community, levels = prob_cols,
                                  labels = paste0("Community ", seq_len(K)))
    # ensure x_label factor levels are preserved after pivot
    plot_long$x_label <- factor(plot_long$x_label, levels = levels(dat$x_label))
    
    p <- ggplot2::ggplot(
      plot_long,
      ggplot2::aes(x = .data$x_label, y = .data$probability, fill = .data$community)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = NA) +
      ggplot2::scale_fill_manual(values = community_colors, name = NULL) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                                  breaks = c(0, 0.5, 1)) +
      ggplot2::labs(
        x = NULL,
        y = if (!is.null(row_label)) row_label else "Membership probability"
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position  = if (show_leg) legend_position else "none",
        panel.spacing.x  = ggplot2::unit(0.3, "lines"),
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(face = "bold"),
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
      )
    
    if (!x_text) {
      p <- p + ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                              axis.ticks.x = ggplot2::element_blank())
    }
    
    if (!is.null(vline_df)) {
      p <- p + ggplot2::geom_vline(
        data = vline_df,
        ggplot2::aes(xintercept = .data$cape_x),
        color     = vline_color, linetype  = vline_linetype,
        linewidth = vline_linewidth, inherit.aes = FALSE
      )
    }
    
    if (!is.null(facet_var) && facet_var %in% names(dat)) {
      p <- p + ggplot2::facet_grid(reformulate(facet_var),
                                   scales = "free_x", space = "free_x")
    }
    
    p
  }
  
  # ── No row facet: single ggplot ───────────────────────────────────────────
  if (is.null(facet_row_var)) {
    p <- build_panel(si, show_leg = show_legend)
    return(p + ggplot2::labs(title = title_str, subtitle = subtitle_str,
                             y = "Membership probability"))
  }
  
  # ── Row facet: cowplot assembly ───────────────────────────────────────────
  row_levels <- if (is.factor(si[[facet_row_var]])) levels(si[[facet_row_var]]) else
    sort(unique(si[[facet_row_var]]))
  
  panels <- Filter(Negate(is.null), lapply(row_levels, function(lv) {
    dat <- si[si[[facet_row_var]] == lv, ]
    if (nrow(dat) == 0) return(NULL)
    build_panel(dat, row_label = as.character(lv), show_leg = FALSE)
  }))
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_str, fontface = "bold", size = base_size + 1,
                        x = 0.02, hjust = 0) +
    cowplot::draw_label(subtitle_str, size = base_size - 1, color = "grey40",
                        x = 0.02, y = 0.25, hjust = 0)
  
  stacked <- cowplot::plot_grid(plotlist = panels, ncol = 1,
                                labels = letters[seq_along(panels)],
                                label_size = 14)
  
  if (!show_legend) {
    return(cowplot::plot_grid(title_grob, stacked,
                              ncol = 1, rel_heights = c(0.06, 1)))
  }
  
  legend_plot <- build_panel(si, show_leg = TRUE) +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::guides(fill = ggplot2::guide_legend())
  
  if (legend_position == "right") {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-right", return_all = TRUE)
    inner <- cowplot::plot_grid(stacked, legend_grob, nrow = 1, rel_widths = c(1, 0.15))
    cowplot::plot_grid(title_grob, inner, ncol = 1, rel_heights = c(0.06, 1))
  } else if (legend_position == "left") {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-left", return_all = TRUE)
    inner <- cowplot::plot_grid(legend_grob, stacked, nrow = 1, rel_widths = c(0.15, 1))
    cowplot::plot_grid(title_grob, inner, ncol = 1, rel_heights = c(0.06, 1))
  } else if (legend_position == "top") {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)
    cowplot::plot_grid(title_grob, legend_grob, stacked,
                       ncol = 1, rel_heights = c(0.06, 0.08, 1))
  } else {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)
    cowplot::plot_grid(title_grob, stacked, legend_grob,
                       ncol = 1, rel_heights = c(0.06, 1, 0.08))
  }
}

