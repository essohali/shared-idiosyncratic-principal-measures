# ============================================================================
# 08_publication_figure.R
#
# Figure 4 -- Temporal localisation of shared and idiosyncratic
#             principal measures
#
# Publication layout matching the manuscript:
#   A. Dominant shared principal measure
#   B. First idiosyncratic principal measures
#
# IMPORTANT:
# - No statistical estimation is performed or modified here.
# - Principal measures come from principal_measure_table().
# - Temporal profiles come exclusively from pm_tv_profile().
# - Peak 50-ms intervals are detected automatically.
# - No Allen peak interval is hard-coded.
# ============================================================================


export_temporal_decomposition_figure <- function(
    fit,
    areas,
    out_root,
    bin_width = 0.05
) {
  
  # ==========================================================================
  # 1. CHECK DEPENDENCIES
  # ==========================================================================
  
  if (!exists("principal_measure_table")) {
    stop(
      "Function 'principal_measure_table' is not available. ",
      "Source 05_principal_measures.R first."
    )
  }
  
  if (!exists("pm_tv_profile")) {
    stop(
      "Function 'pm_tv_profile' is not available. ",
      "Source 07_temporal_interpretation.R first."
    )
  }
  
  if (length(areas) != 5L) {
    stop(
      "This publication layout expects exactly five cortical areas."
    )
  }
  
  
  # ==========================================================================
  # 2. OUTPUT DIRECTORIES
  # ==========================================================================
  
  figdir <- file.path(out_root, "figures")
  tabdir <- file.path(out_root, "tables")
  
  dir.create(
    figdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    tabdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  outfile <- file.path(
    figdir,
    "Figure_temporal_decomposition.pdf"
  )
  
  peakfile <- file.path(
    tabdir,
    "Figure4_temporal_peak_intervals.csv"
  )
  
  
  # ==========================================================================
  # 3. PRINCIPAL MEASURES
  # ==========================================================================
  
  # Shared principal measures
  pmS <- principal_measure_table(
    fit$grid,
    fit$eigS,
    "shared",
    "shared"
  )
  
  # Idiosyncratic principal measures
  pmI <- do.call(
    rbind,
    lapply(
      seq_along(areas),
      function(r) {
        
        principal_measure_table(
          fit$grid,
          fit$eigI[[r]],
          "idiosyncratic",
          areas[r]
        )
      }
    )
  )
  
  
  # ==========================================================================
  # 4. TOTAL-VARIATION TEMPORAL PROFILES
  # ==========================================================================
  
  # First shared principal measure
  profS <- pm_tv_profile(
    pmS[pmS$k == 1L, , drop = FALSE],
    bin_width = bin_width
  )
  
  # First idiosyncratic principal measure for each cortical area
  profI <- pm_tv_profile(
    pmI[pmI$k == 1L, , drop = FALSE],
    bin_width = bin_width
  )
  
  
  # ==========================================================================
  # 5. AUTOMATIC PEAK DETECTION
  # ==========================================================================
  
  peak_info <- function(df) {
    
    df <- df[order(df$bin_left), , drop = FALSE]
    
    i <- which.max(df$tv_share)
    
    data.frame(
      peak_left = df$bin_left[i],
      peak_right = df$bin_right[i],
      peak_mid = (
        df$bin_left[i] +
          df$bin_right[i]
      ) / 2,
      peak_share = df$tv_share[i],
      stringsAsFactors = FALSE
    )
  }
  
  
  peakS <- peak_info(profS)
  
  peakI <- lapply(
    areas,
    function(a) {
      
      g <- profI[
        profI$area == a,
        ,
        drop = FALSE
      ]
      
      peak_info(g)
    }
  )
  
  names(peakI) <- areas
  
  
  # ==========================================================================
  # 6. TABLE OF VALUES USED IN FIGURE 4
  # ==========================================================================
  
  peak_table <- data.frame(
    component = c(
      "shared",
      rep("idiosyncratic", length(areas))
    ),
    
    area = c(
      "shared",
      areas
    ),
    
    peak_50ms_left = c(
      peakS$peak_left,
      vapply(
        peakI,
        function(x) x$peak_left,
        numeric(1)
      )
    ),
    
    peak_50ms_right = c(
      peakS$peak_right,
      vapply(
        peakI,
        function(x) x$peak_right,
        numeric(1)
      )
    ),
    
    peak_50ms_tv_share = c(
      peakS$peak_share,
      vapply(
        peakI,
        function(x) x$peak_share,
        numeric(1)
      )
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    peak_table,
    peakfile,
    row.names = FALSE
  )
  
  
  # ==========================================================================
  # 7. COLOUR PALETTE
  # ==========================================================================
  
  # Palette close to the manuscript Figure 4.
  # These choices are purely graphical.
  
  col_shared <- "#225EA8"
  
  col_area <- c(
    VISp  = "#E41A1C",
    VISl  = "#FF7F00",
    VISrl = "#1A9850",
    VISal = "#7A35C5",
    VISam = "#16A6B6"
  )
  
  
  # Fallback if area names differ
  missing_areas <- setdiff(
    areas,
    names(col_area)
  )
  
  if (length(missing_areas) > 0L) {
    
    extra_cols <- grDevices::hcl.colors(
      length(missing_areas),
      palette = "Dark 3"
    )
    
    names(extra_cols) <- missing_areas
    
    col_area <- c(
      col_area,
      extra_cols
    )
  }
  
  
  # ==========================================================================
  # 8. HELPER -- DRAW FILLED TEMPORAL PROFILE
  # ==========================================================================
  
  draw_profile <- function(
    df,
    col,
    ylim,
    show_x_axis = TRUE,
    show_y_axis = TRUE,
    ylab = "",
    xlab = "",
    annotation_cex = 0.82,
    annotation_side = "right"
  ) {
    
    df <- df[
      order(df$bin_left),
      ,
      drop = FALSE
    ]
    
    mid <- (
      df$bin_left +
        df$bin_right
    ) / 2
    
    y <- 100 * df$tv_share
    
    pk <- peak_info(df)
    
    peak_left <- pk$peak_left
    peak_right <- pk$peak_right
    peak_mid <- pk$peak_mid
    peak_y <- 100 * pk$peak_share
    
    
    # ------------------------------------------------------------------------
    # Empty plotting region
    # ------------------------------------------------------------------------
    
    plot(
      mid,
      y,
      type = "n",
      xlim = c(0, 2),
      ylim = ylim,
      xaxs = "i",
      yaxs = "i",
      axes = FALSE,
      xlab = "",
      ylab = "",
      bty = "n"
    )
    
    
    # ------------------------------------------------------------------------
    # Peak 50-ms shaded interval
    # ------------------------------------------------------------------------
    
    rect(
      xleft = peak_left,
      ybottom = ylim[1],
      xright = peak_right,
      ytop = ylim[2],
      col = adjustcolor(
        col,
        alpha.f = 0.12
      ),
      border = NA
    )
    
    
    # ------------------------------------------------------------------------
    # Area under the profile
    # ------------------------------------------------------------------------
    
    polygon(
      x = c(
        mid,
        rev(mid)
      ),
      y = c(
        y,
        rep(0, length(y))
      ),
      col = adjustcolor(
        col,
        alpha.f = 0.18
      ),
      border = NA
    )
    
    
    # ------------------------------------------------------------------------
    # Profile curve
    # ------------------------------------------------------------------------
    
    lines(
      mid,
      y,
      col = col,
      lwd = 2.1
    )
    
    
    # ------------------------------------------------------------------------
    # Dashed line at centre of maximal 50-ms interval
    # ------------------------------------------------------------------------
    
    abline(
      v = peak_mid,
      col = col,
      lty = 2,
      lwd = 1.25
    )
    
    
    # ------------------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------------------
    
    if (show_x_axis) {
      
      axis(
        side = 1,
        at = seq(
          0,
          2,
          by = 0.5
        ),
        labels = format(
          seq(
            0,
            2,
            by = 0.5
          ),
          nsmall = 1
        ),
        cex.axis = 0.82
      )
      
    } else {
      
      axis(
        side = 1,
        at = seq(
          0,
          2,
          by = 0.5
        ),
        labels = FALSE,
        tck = -0.015
      )
    }
    
    
    if (show_y_axis) {
      
      axis(
        side = 2,
        las = 1,
        cex.axis = 0.80
      )
    }
    
    
    # Left and bottom axes only
    abline(
      h = 0,
      col = "black",
      lwd = 0.8
    )
    
    
    # ------------------------------------------------------------------------
    # Axis labels
    # ------------------------------------------------------------------------
    
    if (nzchar(ylab)) {
      
      mtext(
        ylab,
        side = 2,
        line = 2.6,
        cex = 0.90
      )
    }
    
    
    if (nzchar(xlab)) {
      
      mtext(
        xlab,
        side = 1,
        line = 2.4,
        cex = 0.90
      )
    }
    
    
    # ------------------------------------------------------------------------
    # Automatic interval annotation
    # ------------------------------------------------------------------------
    
    interval_label <- sprintf(
      "%d\u2013%d ms",
      round(
        1000 * peak_left
      ),
      round(
        1000 * peak_right
      )
    )
    
    
    # Position annotation adaptively
    if (annotation_side == "left") {
      
      text_x <- max(
        0.04,
        peak_mid - 0.12
      )
      
      adj_text <- 1
      
    } else {
      
      text_x <- min(
        1.94,
        peak_mid + 0.14
      )
      
      adj_text <- 0
    }
    
    
    text_y <- min(
      ylim[2] * 0.86,
      peak_y + 0.10 * diff(ylim)
    )
    
    
    # Arrow
    arrows(
      x0 = text_x,
      y0 = text_y - 0.02 * diff(ylim),
      x1 = peak_mid,
      y1 = peak_y,
      length = 0.07,
      angle = 22,
      lwd = 1.15,
      col = col
    )
    
    
    # Label
    text(
      x = text_x,
      y = text_y,
      labels = interval_label,
      col = col,
      cex = annotation_cex,
      adj = c(
        adj_text,
        0
      )
    )
    
    
    invisible(pk)
  }
  
  

  # ==========================================================================
  # 9. PUBLICATION FIGURE -- SINGLE PASS
  # ==========================================================================

  # Common y-range for the five idiosyncratic profiles.
  ylim_I <- c(
    0,
    max(
      5,
      100 * max(profI$tv_share) * 1.15
    )
  )

  grDevices::pdf(
    file = outfile,
    width = 10,
    height = 6.5,
    useDingbats = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # Layout:
  #
  #   1 | 2
  #   1 | 3
  #   1 | 4
  #   1 | 5
  #   1 | 6
  #   7 | 8   <- legends
  #
  # --------------------------------------------------------------------------
  
  layout(
    matrix(
      c(
        1, 2,
        1, 3,
        1, 4,
        1, 5,
        1, 6,
        7, 8
      ),
      nrow = 6,
      byrow = TRUE
    ),
    widths = c(
      0.78,
      1.22
    ),
    heights = c(
      1,
      1,
      1,
      1,
      1,
      0.72
    )
  )
  
  
  # ==========================================================================
# 10. PANEL A -- SHARED PRINCIPAL MEASURE
  # ==========================================================================
  
  par(
    mar = c(
      4.0,
      4.5,
      2.5,
      1.6
    )
  )
  
  
  draw_profile(
    df = profS,
    col = col_shared,
    ylim = c(
      0,
      max(
        6,
        100 * max(
          profS$tv_share
        ) * 1.15
      )
    ),
    show_x_axis = TRUE,
    show_y_axis = TRUE,
    ylab = "Total-variation share (%)",
    xlab = "Time after stimulus onset (s)",
    annotation_cex = 0.85,
    annotation_side = "right"
  )
  
  
  mtext(
    "(a)  Dominant shared principal measure",
    side = 3,
    line = 0.55,
    adj = 0,
    font = 2,
    cex = 0.98
  )
  # ==========================================================================
# 11. PANEL B -- FIVE IDIOSYNCRATIC PRINCIPAL MEASURES
  # ==========================================================================
  
  for (j in seq_along(areas)) {
    
    a <- areas[j]
    
    g <- profI[
      profI$area == a,
      ,
      drop = FALSE
    ]
    
    show_x <- (
      j == length(areas)
    )
    
    
    par(
      mar = c(
        if (show_x) 3.6 else 0.7,
        3.8,
        if (j == 1L) 2.5 else 0.25,
        0.8
      )
    )
    
    
    draw_profile(
      df = g,
      col = col_area[a],
      ylim = ylim_I,
      show_x_axis = show_x,
      show_y_axis = TRUE,
      ylab = "",
      xlab = if (show_x) {
        "Time after stimulus onset (s)"
      } else {
        ""
      },
      annotation_cex = 0.78,
      annotation_side = "right"
    )
    
    
    mtext(
      a,
      side = 2,
      line = 2.35,
      las = 1,
      font = 2,
      cex = 0.88
    )
    
    
    if (j == 1L) {
      
      mtext(
        "(b)  First idiosyncratic principal measures",
        side = 3,
        line = 0.55,
        adj = 0,
        font = 2,
        cex = 0.98
      )
    }
  }
  
  
  # ==========================================================================
# 12. LEGEND UNDER PANEL A
  # ==========================================================================
  
  par(
    mar = c(
      0,
      0,
      0,
      0
    )
  )
  
  plot.new()
  
  legend(
    "center",
    legend = "Shared component",
    col = col_shared,
    lwd = 2.5,
    bty = "n",
    horiz = TRUE,
    cex = 0.85
  )
  
  
  # ==========================================================================
# 13. LEGEND UNDER PANEL B
  # ==========================================================================
  
  par(
    mar = c(
      0,
      0,
      0,
      0
    )
  )
  
  plot.new()
  
  legend(
    "center",
    legend = areas,
    col = col_area[areas],
    lwd = 2.5,
    bty = "n",
    horiz = TRUE,
    cex = 0.78,
    xpd = NA
  )
  
  
  # ==========================================================================
# 14. CLOSE PDF
  # ==========================================================================
  
  dev.off()
  
  
  # ==========================================================================
# 15. CONSOLE INFORMATION
  # ==========================================================================
  
  message(
    "Figure 4 written to: ",
    outfile
  )
  
  message(
    "Figure 4 peak intervals written to: ",
    peakfile
  )
  
  
  # ==========================================================================
# 16. RETURN RESULTS
  # ==========================================================================
  
  invisible(
    list(
      figure = outfile,
      peaks = peak_table,
      shared_profile = profS,
      idiosyncratic_profiles = profI
    )
  )
}
