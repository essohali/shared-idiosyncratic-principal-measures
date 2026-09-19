# ================================================================
# R/07_figures.R
#
# Signature figures for:
# Shared and Idiosyncratic Principal Measures
#
# Figure 1: Shared-idiosyncratic covariance decomposition
# Figure 2: First-order convergence fingerprint
# Figure 3: Principal-measure recovery
# ================================================================


library(dplyr)
library(tidyr)
library(ggplot2)


# ================================================================
# GLOBAL VISUAL IDENTITY
# ================================================================

COL_SHARED <- "#1F5A94"
COL_IDIO   <- "#D97732"
COL_TOTAL  <- "#6F4C9B"
COL_TRUTH  <- "#222222"
COL_LIGHT  <- "#D9D9D9"
COL_GRID   <- "#E8E8E8"
COL_BG     <- "#FAFAFA"


paper_theme <- function(base_size = 12) {
  
  theme_minimal(base_size = base_size) +
    
    theme(
      plot.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      panel.background =
        element_rect(
          fill = COL_BG,
          colour = NA
        ),
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major =
        element_line(
          colour = COL_GRID,
          linewidth = 0.25
        ),
      
      strip.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      strip.text =
        element_text(
          face = "bold",
          colour = "#222222",
          size = base_size
        ),
      
      plot.title =
        element_text(
          face = "bold",
          size = base_size + 3,
          colour = "#111111",
          margin = margin(b = 5)
        ),
      
      plot.subtitle =
        element_text(
          size = base_size,
          colour = "#444444",
          margin = margin(b = 10)
        ),
      
      axis.title =
        element_text(
          face = "bold",
          colour = "#333333"
        ),
      
      axis.text =
        element_text(
          colour = "#444444"
        ),
      
      legend.position =
        "bottom",
      
      legend.title =
        element_blank(),
      
      legend.text =
        element_text(
          size = base_size - 1
        ),
      
      plot.margin =
        margin(10, 14, 10, 10)
    )
}

# ================================================================
# FIGURE 1
# POPULATION VS ESTIMATED COVARIANCE DECOMPOSITION
#
# Population:
#   C_22 = a_2^2 C_S + C_2^I
#
# Estimate:
#   C_22^rec = a_hat_2^2 C_S_hat + C_2_hat^I
#
# Balanced regime, n = 500,
# representative replication = 88
# ================================================================


# ================================================================
# Helper: reconstruct exactly one Monte Carlo replication
# ================================================================

reconstruct_replication <- function(
    scenario_use = "balanced",
    n_use = 500,
    rep_use = 88,
    final = TRUE
) {
  
  design <- make_design(
    final = final
  )
  
  grid <- seq(
    0,
    1,
    length.out = design$grid_n
  )
  
  truth <- scenario_truth(
    design,
    scenario_use,
    grid = grid
  )
  
  truth <- prepare_truth_spectra(
    truth,
    design
  )
  
  scenario_index <- match(
    scenario_use,
    design$scenarios$scenario
  )
  
  if (is.na(scenario_index)) {
    stop("Scenario not found in design$scenarios.")
  }
  
  # Exact seed used in run_one_replication()
  seed <-
    design$base_seed +
    scenario_index * 10000000L +
    n_use * 1000L +
    rep_use
  
  dat <- simulate_dataset(
    n = n_use,
    truth = truth,
    grid = truth$grid,
    seed = seed
  )
  
  fit <- fit_model(
    dat = dat,
    design = design,
    seed = seed + 77L
  )
  
  list(
    design = design,
    truth = truth,
    dat = dat,
    fit = fit,
    seed = seed
  )
}


# ================================================================
# Figure 1
# ================================================================

figure1_truth_vs_estimated <- function(
    scenario_use = "balanced",
    n_use = 500,
    rep_use = 88,
    margin_use = 2,
    final = TRUE
) {
  
  # ==============================================================
  # 1. Reconstruct exact representative replication
  # ==============================================================
  
  obj <- reconstruct_replication(
    scenario_use = scenario_use,
    n_use = n_use,
    rep_use = rep_use,
    final = final
  )
  
  design <- obj$design
  truth  <- obj$truth
  dat    <- obj$dat
  fit    <- obj$fit
  
  r <- margin_use
  tt <- dat$grid
  
  
  # ==============================================================
  # 2. Scenario parameters
  # ==============================================================
  
  scen <- design$scenarios %>%
    dplyr::filter(
      scenario == scenario_use
    )
  
  if (nrow(scen) != 1) {
    stop("Scenario not found uniquely.")
  }
  
  c_scale <- scen$c_scale[[1]]
  
  a_true <- design$a
  
  
  # ==============================================================
  # 3. Population shared covariance C_S
  # ==============================================================
  
  eta1 <- sqrt(2) *
    sin(pi * tt / 2)
  
  eta2 <- sqrt(2) *
    sin(3 * pi * tt / 2)
  
  CS_true <-
    25 * c_scale *
    tcrossprod(eta1) +
    8.75 * c_scale *
    tcrossprod(eta2)
  
  
  # ==============================================================
  # 4. Population idiosyncratic covariance C_r^I
  # ==============================================================
  
  phi1 <- sqrt(2) *
    sin(pi * tt)
  
  phi2 <- sqrt(2) *
    sin(2 * pi * tt)
  
  CV <-
    4 * tcrossprod(phi1) +
    tcrossprod(phi2)
  
  
  # --------------------------------------------------------------
  # Same baseline construction as the DGP
  # --------------------------------------------------------------
  
  B_S <-
    5 * pi / sqrt(2) +
    sqrt(8.75) *
    3 * pi / sqrt(2)
  
  B_I <-
    4 * sqrt(2) * pi
  
  lambda0 <-
    a_true *
    sqrt(c_scale) *
    B_S +
    B_I +
    1
  
  
  # --------------------------------------------------------------
  # Poisson-martingale covariance
  # --------------------------------------------------------------
  
  Kmin <- outer(
    tt,
    tt,
    pmin
  )
  
  CI_true <-
    CV +
    lambda0[r] *
    Kmin
  
  
  # ==============================================================
  # 5. Population decomposition
  # ==============================================================
  
  Cshared_true <-
    a_true[r]^2 *
    CS_true
  
  Ctotal_true <-
    Cshared_true +
    CI_true
  
  
  # ==============================================================
  # 6. Estimated decomposition
  # ==============================================================
  
  CS_hat <- fit$CS_hat
  
  CI_hat <- fit$CI_hat[[r]]
  
  a_hat <- fit$a_hat
  
  Cshared_hat <-
    a_hat[r]^2 *
    CS_hat
  
  # IMPORTANT:
  # This is the reconstructed marginal covariance after PSD
  # projection of the idiosyncratic component.
  Ctotal_hat <-
    Cshared_hat +
    CI_hat
  
  
  # ==============================================================
  # 7. Empirical marginal covariance
  #    Diagnostic only
  # ==============================================================
  
  Crr_emp <-
    empirical_kernel(
      dat$X[[r]]
    )
  
  w <- fit$weights
  
  
  # ==============================================================
  # 8. Weighted Hilbert-Schmidt norm
  # ==============================================================
  
  weighted_HS_norm <- function(A, w) {
    
    sqrt(
      sum(
        A^2 *
          outer(w, w)
      )
    )
  }
  
  
  reconstruction_error <-
    weighted_HS_norm(
      Crr_emp - Ctotal_hat,
      w
    ) /
    weighted_HS_norm(
      Crr_emp,
      w
    )
  
  
  # ==============================================================
  # 9. Population shared-variation proportion
  # ==============================================================
  
  tau_S <-
    33.75 *
    c_scale
  
  rho_true <-
    (
      a_true[r]^2 *
        tau_S
    ) /
    (
      a_true[r]^2 *
        tau_S +
        5 +
        lambda0[r] / 2
    )
  
  
  # ==============================================================
  # 10. Helper: kernel to long-format data
  # ==============================================================
  
  kernel_to_df <- function(
    K,
    row_name,
    component_name
  ) {
    
    expand.grid(
      t = tt,
      u = tt
    ) %>%
      dplyr::mutate(
        value = as.vector(K),
        row = row_name,
        component = component_name
      )
  }
  
  
  # ==============================================================
  # 11. Combine population and estimated kernels
  # ==============================================================
  
  dd <- dplyr::bind_rows(
    
    # ------------------------------------------------------------
    # Population
    # ------------------------------------------------------------
    
    kernel_to_df(
      Ctotal_true,
      "Population",
      "C22"
    ),
    
    kernel_to_df(
      Cshared_true,
      "Population",
      "Shared2"
    ),
    
    kernel_to_df(
      CI_true,
      "Population",
      "Idio2"
    ),
    
    # ------------------------------------------------------------
    # Estimated
    # ------------------------------------------------------------
    
    kernel_to_df(
      Ctotal_hat,
      "Estimated",
      "C22"
    ),
    
    kernel_to_df(
      Cshared_hat,
      "Estimated",
      "Shared2"
    ),
    
    kernel_to_df(
      CI_hat,
      "Estimated",
      "Idio2"
    )
  )
  
  
  # ==============================================================
  # 12. Factor ordering
  # ==============================================================
  
  dd$row <- factor(
    dd$row,
    levels = c(
      "Population",
      "Estimated"
    )
  )
  
  dd$component <- factor(
    dd$component,
    levels = c(
      "C22",
      "Shared2",
      "Idio2"
    )
  )
  
  
  # ==============================================================
  # 13. Mathematical facet labels
  # ==============================================================
  
  component_labeller <- ggplot2::labeller(
    component = ggplot2::as_labeller(
      c(
        C22 =
          "C[22]*group('(',list(t,u),')')",
        
        Shared2 =
          "a[2]^2*C[S]*group('(',list(t,u),')')",
        
        Idio2 =
          "C[2]^I*group('(',list(t,u),')')"
      ),
      default = ggplot2::label_parsed
    )
  )
  
  # ==============================================================
  # 14. Common covariance scale
  #
  # Observed range for this representative replication is roughly
  #
  #   [-0.424, 120.679]
  #
  # Negative values are genuine but very small relative to the
  # positive range. We therefore use an asymmetric diverging scale:
  # blue is reserved for negative covariance, while most colour
  # resolution is allocated to the positive range.
  # ==============================================================
  
  zmin <- min(
    dd$value,
    na.rm = TRUE
  )
  
  zmax <- max(
    dd$value,
    na.rm = TRUE
  )
  
  
  colour_breaks <- c(
    zmin,
    0,
    25,
    50,
    80,
    zmax
  )
  
  colour_breaks <- sort(
    unique(
      pmax(
        zmin,
        pmin(
          zmax,
          colour_breaks
        )
      )
    )
  )
  
  
  # ==============================================================
  # 15. Main plot
  # ==============================================================
  
  p <- ggplot(
    dd,
    aes(
      x = t,
      y = u,
      fill = value
    )
  ) +
    
    geom_raster() +
    
    facet_grid(
      rows = vars(row),
      cols = vars(component),
      labeller = component_labeller,
      switch = "y"
    )+
    
    scale_fill_gradientn(
      colours = c(
        "#3B6FA5",  # negative covariance
        "#F7F7F4",  # zero
        "#FCE2D3",
        "#F6B38D",
        "#E47745",
        "#A93A1F"
      ),
      
      values = scales::rescale(
        colour_breaks,
        from = c(
          zmin,
          zmax
        )
      ),
      
      limits = c(
        zmin,
        zmax
      ),
      
      oob = scales::squish,
      
      name = "Covariance"
    ) +
    
    coord_equal(
      expand = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(
        0,
        0.5,
        1
      ),
      labels = c(
        "0",
        "0.5",
        "1"
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    scale_y_continuous(
      breaks = c(
        0,
        0.5,
        1
      ),
      labels = c(
        "0",
        "0.5",
        "1"
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    labs(
      x = expression(t),
      y = expression(u)
    ) +
    
    theme_minimal(
      base_size = 11
    ) +
    
    theme(
      
      # ----------------------------------------------------------
      # Panel background
      # ----------------------------------------------------------
      
      panel.background =
        element_rect(
          fill = "#FAFAFA",
          colour = NA
        ),
      
      panel.grid =
        element_blank(),
      
      # ----------------------------------------------------------
      # Fine border: same visual grammar as Figures 2--3
      # ----------------------------------------------------------
      
      panel.border =
        element_rect(
          colour = "#555555",
          fill = NA,
          linewidth = 0.45
        ),
      
      # ----------------------------------------------------------
      # Facet strips
      # ----------------------------------------------------------
      
      strip.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      strip.text.x =
        element_text(
          face = "bold",
          size = 12,
          colour = "#111111",
          margin = margin(
            b = 6
          )
        ),
      
      strip.text.y.left =
        element_text(
          face = "bold",
          size = 10.5,
          colour = "#111111",
          angle = 90,
          margin = margin(
            r = 5
          )
        ),
      
      # ----------------------------------------------------------
      # Axes
      # ----------------------------------------------------------
      
      axis.title =
        element_text(
          size = 10.5,
          colour = "#222222"
        ),
      
      axis.text =
        element_text(
          size = 8.5,
          colour = "#333333"
        ),
      
      # ----------------------------------------------------------
      # Panel spacing
      # ----------------------------------------------------------
      
      panel.spacing =
        grid::unit(
          0.6,
          "lines"
        ),
      
      # ----------------------------------------------------------
      # Common covariance colour bar
      # ----------------------------------------------------------
      
      legend.position =
        "right",
      
      legend.title =
        element_text(
          face = "bold",
          size = 9.5
        ),
      
      legend.text =
        element_text(
          size = 8.5
        ),
      
      legend.key.height =
        grid::unit(
          3.2,
          "cm"
        ),
      
      # ----------------------------------------------------------
      # Margins
      # ----------------------------------------------------------
      
      plot.margin =
        margin(
          6,
          8,
          6,
          6
        )
    )
  
  
  # ==============================================================
  # 16. Diagnostics stored in the ggplot object
  # ==============================================================
  
  attr(
    p,
    "figure_info"
  ) <- list(
    
    scenario =
      scenario_use,
    
    n =
      n_use,
    
    replication =
      rep_use,
    
    margin =
      r,
    
    seed =
      obj$seed,
    
    a_true =
      a_true[r],
    
    a_hat =
      a_hat[r],
    
    tau_true =
      tau_S,
    
    tau_hat =
      fit$tau_hat,
    
    rho_true =
      rho_true,
    
    reconstruction_error =
      reconstruction_error,
    
    covariance_range =
      range(
        dd$value,
        na.rm = TRUE
      ),
    
    covariance_min =
      min(
        dd$value,
        na.rm = TRUE
      ),
    
    covariance_max =
      max(
        dd$value,
        na.rm = TRUE
      ),
    
    component_ranges =
      dd %>%
      dplyr::group_by(
        row,
        component
      ) %>%
      dplyr::summarise(
        min = min(
          value,
          na.rm = TRUE
        ),
        max = max(
          value,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
  )
  
  
  return(p)
}
# ================================================================
# FIGURE 2
# FIRST-ORDER CONVERGENCE FINGERPRINT
#
# Numerical validation of Theorems 3--4
# ================================================================
figure2_first_order_convergence <- function(res) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork'.")
  }
  
  # ==============================================================
  # 1. Monte Carlo median errors
  # ==============================================================
  
  dd <- res$metrics %>%
    group_by(
      scenario,
      n
    ) %>%
    summarise(
      Ea    = median(Ea, na.rm = TRUE),
      Etau  = median(Etau, na.rm = TRUE),
      ES    = median(ES, na.rm = TRUE),
      Eeta1 = median(EetaS1, na.rm = TRUE),
      Emu1  = median(EmuS1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(
        Ea,
        Etau,
        ES,
        Eeta1,
        Emu1
      ),
      names_to = "criterion",
      values_to = "error"
    )
  
  
  # ==============================================================
  # 2. Criterion ordering
  # ==============================================================
  
  dd$criterion <- factor(
    dd$criterion,
    levels = c(
      "Ea",
      "Etau",
      "ES",
      "Eeta1",
      "Emu1"
    )
  )
  
  
  # ==============================================================
  # 3. Visual identity
  # ==============================================================
  
  criterion_colors <- c(
    Ea    = "#1F5A94",
    Etau  = "#7A4AA5",
    ES    = "#16836A",
    Eeta1 = "#D96B20",
    Emu1  = "#C63D32"
  )
  
  criterion_shapes <- c(
    Ea    = 16,
    Etau  = 17,
    ES    = 15,
    Eeta1 = 18,
    Emu1  = 8
  )
  
  criterion_labels <- c(
    Ea    = expression(E[a]),
    Etau  = expression(E[tau]),
    ES    = expression(E[S]),
    Eeta1 = expression(E[eta[1]]),
    Emu1  = expression(E[mu[1]])
  )
  
  
  # ==============================================================
  # 4. Helper function for one scenario
  # ==============================================================
  
  make_panel <- function(
    data,
    scenario_use,
    panel_label,
    show_y_title = FALSE
  ) {
    
    ds <- data %>%
      filter(
        scenario == scenario_use
      )
    
    # ------------------------------------------------------------
    # n^{-1/2} reference
    #
    # Anchor slightly above E_S at n = 100
    # ------------------------------------------------------------
    
    anchor <- ds %>%
      filter(
        criterion == "ES",
        n == min(n)
      )
    
    n0 <- anchor$n[1]
    
    y0 <- 1.08 * anchor$error[1]
    
    ref <- tibble(
      n = sort(unique(ds$n))
    ) %>%
      mutate(
        error =
          y0 *
          (n / n0)^(-1 / 2)
      )
    
    
    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    
    ggplot(
      ds,
      aes(
        x = n,
        y = error,
        colour = criterion,
        shape = criterion,
        group = criterion
      )
    ) +
      
      # Reference slope
      geom_line(
        data = ref,
        aes(
          x = n,
          y = error
        ),
        inherit.aes = FALSE,
        colour = "#555555",
        linetype = "22",
        linewidth = 0.85
      ) +
      
      # Estimation curves
      geom_line(
        linewidth = 0.85
      ) +
      
      geom_point(
        size = 2.6,
        stroke = 0.75
      ) +
      
      # Log scales
      scale_x_log10(
        breaks = c(
          100,
          250,
          500,
          1000
        ),
        labels = c(
          "100",
          "250",
          "500",
          "1000"
        )
      ) +
      
      scale_y_log10() +
      
      # Colours
      scale_colour_manual(
        values = criterion_colors,
        labels = criterion_labels
      ) +
      
      # Shapes
      scale_shape_manual(
        values = criterion_shapes,
        labels = criterion_labels
      ) +
      
      labs(
        title = paste0(
          panel_label,
          "  ",
          tools::toTitleCase(
            scenario_use
          )
        ),
        
        x = expression(n),
        
        y =
          if (show_y_title) {
            "Median estimation error"
          } else {
            NULL
          },
        
        colour = NULL,
        shape = NULL
      ) +
      
      theme_minimal(
        base_size = 11
      ) +
      
      theme(
        
        # --------------------------------------------
        # Background
        # --------------------------------------------
        
        panel.background =
          element_rect(
            fill = "#FAFAFA",
            colour = NA
          ),
        
        # --------------------------------------------
        # Border around each panel
        # --------------------------------------------
        
        panel.border =
          element_rect(
            colour = "#555555",
            fill = NA,
            linewidth = 0.45
          ),
        
        # --------------------------------------------
        # Grid
        # --------------------------------------------
        
        panel.grid.minor =
          element_blank(),
        
        panel.grid.major =
          element_line(
            colour = "#E3E3E3",
            linewidth = 0.25
          ),
        
        # --------------------------------------------
        # Panel title
        # --------------------------------------------
        
        plot.title =
          element_text(
            face = "bold",
            size = 11.5,
            hjust = 0.5,
            colour = "#111111",
            margin = margin(
              b = 4
            )
          ),
        
        # --------------------------------------------
        # Axes
        # --------------------------------------------
        
        axis.title =
          element_text(
            size = 10.5,
            colour = "#222222"
          ),
        
        axis.text =
          element_text(
            size = 9,
            colour = "#333333"
          ),
        
        # --------------------------------------------
        # Legend INSIDE each panel
        # --------------------------------------------
        
        legend.position =
          c(
            0.97,
            0.97
          ),
        
        legend.justification =
          c(
            1,
            1
          ),
        
        legend.direction =
          "vertical",
        
        legend.background =
          element_rect(
            fill = scales::alpha(
              "white",
              0.90
            ),
            colour = "#666666",
            linewidth = 0.35
          ),
        
        legend.key =
          element_blank(),
        
        legend.key.height =
          grid::unit(
            0.34,
            "cm"
          ),
        
        legend.key.width =
          grid::unit(
            0.45,
            "cm"
          ),
        
        legend.text =
          element_text(
            size = 8
          ),
        
        legend.spacing.y =
          grid::unit(
            0.01,
            "cm"
          ),
        
        # --------------------------------------------
        # Margins
        # --------------------------------------------
        
        plot.margin =
          margin(
            5,
            5,
            5,
            5
          )
      )
  }
  
  
  # ==============================================================
  # 5. Three panels
  # ==============================================================
  
  p1 <- make_panel(
    data = dd,
    scenario_use = "weak",
    panel_label = "(a)",
    show_y_title = TRUE
  )
  
  p2 <- make_panel(
    data = dd,
    scenario_use = "balanced",
    panel_label = "(b)",
    show_y_title = FALSE
  )
  
  p3 <- make_panel(
    data = dd,
    scenario_use = "strong",
    panel_label = "(c)",
    show_y_title = FALSE
  )
  
  
  # ==============================================================
  # 6. Final publication layout
  # ==============================================================
  
  final_plot <-
    p1 + p2 + p3 +
    
    patchwork::plot_layout(
      nrow = 1,
      widths = c(
        1,
        1,
        1
      )
    )
  
  return(final_plot)
}

# ================================================================
# FIGURE 3
# PRINCIPAL-MEASURE RECOVERY
#
# Signature measure-valued figure
# ================================================================
prepare_principal_measure_data <- function(
    res,
    scenario_use = "balanced",
    n_use = 500,
    k_use = 1
) {
  
  # ==============================================================
  # 1. Check curves
  # ==============================================================
  
  if (is.null(res$curves)) {
    stop("res$curves does not exist.")
  }
  
  if (nrow(res$curves) == 0) {
    stop("res$curves is empty.")
  }
  
  curves <- as.data.frame(res$curves)
  
  required_cols <- c(
    "scenario",
    "n",
    "rep",
    "k",
    "t",
    "estimate",
    "truth"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(curves)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing columns in res$curves: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  
  # ==============================================================
  # 2. Standardise types
  # ==============================================================
  
  curves <- curves %>%
    mutate(
      
      scenario_std =
        tolower(
          trimws(
            as.character(scenario)
          )
        ),
      
      n_std =
        suppressWarnings(
          as.numeric(
            as.character(n)
          )
        ),
      
      k_std =
        suppressWarnings(
          as.numeric(
            as.character(k)
          )
        ),
      
      rep_std =
        suppressWarnings(
          as.numeric(
            as.character(rep)
          )
        ),
      
      t =
        as.numeric(t),
      
      estimate =
        as.numeric(estimate),
      
      truth =
        as.numeric(truth)
    )
  
  
  scenario_target <-
    tolower(
      trimws(
        as.character(scenario_use)
      )
    )
  
  n_target <-
    as.numeric(n_use)
  
  k_target <-
    as.numeric(k_use)
  
  
  # ==============================================================
  # 3. Configuration requested
  # ==============================================================
  
  dd <- curves %>%
    filter(
      scenario_std == scenario_target,
      n_std == n_target,
      k_std == k_target
    )
  
  
  # ==============================================================
  # 4. Diagnostic if configuration is absent
  # ==============================================================
  
  if (nrow(dd) == 0) {
    
    available <- curves %>%
      distinct(
        scenario_std,
        n_std,
        k_std
      ) %>%
      arrange(
        scenario_std,
        n_std,
        k_std
      )
    
    cat(
      "\nAvailable configurations in res$curves:\n\n"
    )
    
    print(
      as.data.frame(available),
      row.names = FALSE
    )
    
    stop(
      paste0(
        "\nRequested configuration not found: scenario = ",
        scenario_use,
        ", n = ",
        n_use,
        ", k = ",
        k_use,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  # ==============================================================
  # 5. Remove incomplete numerical rows
  # ==============================================================
  
  dd <- dd %>%
    filter(
      is.finite(t),
      is.finite(estimate),
      is.finite(truth),
      is.finite(rep_std)
    )
  
  
  if (nrow(dd) == 0) {
    stop(
      "The selected configuration contains no finite curves."
    )
  }
  
  
  # ==============================================================
  # 6. Sign alignment
  # ==============================================================
  
  signs <- dd %>%
    group_by(rep_std) %>%
    summarise(
      
      inner =
        sum(
          estimate * truth,
          na.rm = TRUE
        ),
      
      sign_align =
        ifelse(
          inner >= 0,
          1,
          -1
        ),
      
      .groups = "drop"
    )
  
  
  dd <- dd %>%
    left_join(
      signs,
      by = "rep_std"
    ) %>%
    mutate(
      estimate_aligned =
        sign_align * estimate
    )
  
  
  # ==============================================================
  # 7. Cumulative spectral direction
  # ==============================================================
  
  cumulative <- dd %>%
    group_by(t) %>%
    summarise(
      
      truth =
        first(truth),
      
      median =
        median(
          estimate_aligned,
          na.rm = TRUE
        ),
      
      q05 =
        as.numeric(
          quantile(
            estimate_aligned,
            probs = 0.05,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q25 =
        as.numeric(
          quantile(
            estimate_aligned,
            probs = 0.25,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q75 =
        as.numeric(
          quantile(
            estimate_aligned,
            probs = 0.75,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q95 =
        as.numeric(
          quantile(
            estimate_aligned,
            probs = 0.95,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      .groups = "drop"
    )
  
  
  # ==============================================================
  # 8. Local signed principal mass
  # ==============================================================
  
  increments <- dd %>%
    arrange(
      rep_std,
      t
    ) %>%
    group_by(rep_std) %>%
    mutate(
      
      t_prev =
        lag(t),
      
      eta_prev =
        lag(estimate_aligned),
      
      truth_prev =
        lag(truth),
      
      dt =
        t - t_prev,
      
      t_mid =
        0.5 * (t + t_prev),
      
      mass_est =
        estimate_aligned - eta_prev,
      
      mass_true =
        truth - truth_prev,
      
      local_est =
        mass_est / dt,
      
      local_true =
        mass_true / dt
    ) %>%
    ungroup() %>%
    filter(
      !is.na(t_mid),
      is.finite(dt),
      dt > 0,
      is.finite(local_est),
      is.finite(local_true)
    )
  
  
  if (nrow(increments) == 0) {
    stop(
      "Unable to construct local signed masses."
    )
  }
  
  
  local_mass <- increments %>%
    group_by(t_mid) %>%
    summarise(
      
      truth =
        first(local_true),
      
      median =
        median(
          local_est,
          na.rm = TRUE
        ),
      
      q05 =
        as.numeric(
          quantile(
            local_est,
            probs = 0.05,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q25 =
        as.numeric(
          quantile(
            local_est,
            probs = 0.25,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q75 =
        as.numeric(
          quantile(
            local_est,
            probs = 0.75,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q95 =
        as.numeric(
          quantile(
            local_est,
            probs = 0.95,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      .groups = "drop"
    )
  
  
  # ==============================================================
  # 9. Return
  # ==============================================================
  
  list(
    cumulative = cumulative,
    local_mass = local_mass,
    raw = dd
  )
}

# ================================================================
# FIGURE 3
# PRINCIPAL-MEASURE RECOVERY
#
# Publication version:
#   (a) Cumulative representation, n = 100
#   (b) Cumulative representation, n = 1000
#   (c) Interval-average signed mass, n = 100
#   (d) Interval-average signed mass, n = 1000
#
# Important:
# - common y-scale within each row;
# - no density interpretation is imposed on the principal measure;
# - uncertainty bands are based on retained Monte Carlo curves.
# ================================================================

figure3_principal_measure_recovery <- function(
    res,
    scenario_use = "balanced",
    n_values = c(100, 1000),
    k_use = 1
) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork'.")
  }
  
  if (length(n_values) != 2) {
    stop("n_values must contain exactly two sample sizes.")
  }
  
  
  # ==============================================================
  # 1. Prepare data
  # ==============================================================
  
  dat_small <- prepare_principal_measure_data(
    res = res,
    scenario_use = scenario_use,
    n_use = n_values[1],
    k_use = k_use
  )
  
  dat_large <- prepare_principal_measure_data(
    res = res,
    scenario_use = scenario_use,
    n_use = n_values[2],
    k_use = k_use
  )
  
  
  # ==============================================================
  # 2. Common y-ranges
  #
  # Same scale within each row is essential for an honest
  # comparison of uncertainty contraction.
  # ==============================================================
  
  y_cum <- range(
    c(
      dat_small$cumulative$q05,
      dat_small$cumulative$q95,
      dat_large$cumulative$q05,
      dat_large$cumulative$q95,
      dat_small$cumulative$truth,
      dat_large$cumulative$truth
    ),
    finite = TRUE
  )
  
  y_mass <- range(
    c(
      dat_small$local_mass$q05,
      dat_small$local_mass$q95,
      dat_large$local_mass$q05,
      dat_large$local_mass$q95,
      dat_small$local_mass$truth,
      dat_large$local_mass$truth
    ),
    finite = TRUE
  )
  
  
  # Small padding
  pad_cum <- 0.04 * diff(y_cum)
  pad_mass <- 0.04 * diff(y_mass)
  
  y_cum <- c(
    y_cum[1] - pad_cum,
    y_cum[2] + pad_cum
  )
  
  y_mass <- c(
    y_mass[1] - pad_mass,
    y_mass[2] + pad_mass
  )
  
  
  # ==============================================================
  # 3. Helper: cumulative representation
  # ==============================================================
  
  make_cumulative_panel <- function(
    dat,
    n_use,
    panel_label
  ) {
    
    ggplot(
      dat$cumulative,
      aes(x = t)
    ) +
      
      # ----------------------------------------------------------
    # 90% Monte Carlo envelope
    # ----------------------------------------------------------
    
    geom_ribbon(
      aes(
        ymin = q05,
        ymax = q95
      ),
      fill = COL_SHARED,
      alpha = 0.10
    ) +
      
      # ----------------------------------------------------------
    # 50% Monte Carlo envelope
    # ----------------------------------------------------------
    
    geom_ribbon(
      aes(
        ymin = q25,
        ymax = q75
      ),
      fill = COL_SHARED,
      alpha = 0.25
    ) +
      
      # ----------------------------------------------------------
    # Monte Carlo median
    # ----------------------------------------------------------
    
    geom_line(
      aes(y = median),
      colour = COL_SHARED,
      linewidth = 1.10
    ) +
      
      # ----------------------------------------------------------
    # Population truth
    # ----------------------------------------------------------
    
    geom_line(
      aes(y = truth),
      colour = "#111111",
      linewidth = 0.85
    ) +
      
      # ----------------------------------------------------------
    # Axes
    # ----------------------------------------------------------
    
    scale_x_continuous(
      breaks = c(
        0,
        0.25,
        0.50,
        0.75,
        1
      ),
      labels = c(
        "0",
        "0.25",
        "0.50",
        "0.75",
        "1"
      ),
      limits = c(0, 1),
      expand = expansion(
        mult = c(0.01, 0.01)
      )
    ) +
      
      coord_cartesian(
        ylim = y_cum
      ) +
      
      labs(
        title = paste0(
          panel_label,
          "  Cumulative representation"
        ),
        
        subtitle = paste0(
          "n = ",
          n_use
        ),
        
        x = expression(t),
        
        y = expression(
          F[mu[1]](t) == eta[1](t)
        )
      ) +
      
      paper_theme(
        base_size = 11
      ) +
      
      theme(
        legend.position = "none",
        
        plot.title = element_text(
          face = "bold",
          size = 12,
          colour = "#111111",
          margin = margin(b = 2)
        ),
        
        plot.subtitle = element_text(
          face = "bold",
          size = 10.5,
          colour = COL_SHARED,
          margin = margin(b = 6)
        ),
        
        axis.title = element_text(
          size = 10.5
        ),
        
        axis.text = element_text(
          size = 9
        ),
        
        panel.grid.minor = element_blank(),
        
        panel.grid.major = element_line(
          colour = "#E6E6E6",
          linewidth = 0.25
        ),
        
        # Contour du graphique
        panel.border = element_rect(
          colour = "#555555",
          fill = NA,
          linewidth = 0.7
        ),
        
        plot.margin = margin(
          5, 8, 5, 5
        )
      )
  }
  
  
  # ==============================================================
  # 4. Helper: interval-average signed mass
  #
  # Quantity displayed:
  #
  # mu_1((t_{j-1},t_j]) / (t_j-t_{j-1})
  #
  # This is an interval-average representation of the signed
  # measure. It is NOT interpreted as a density estimator.
  # ==============================================================
  
  make_mass_panel <- function(
    dat,
    n_use,
    panel_label
  ) {
    
    ggplot(
      dat$local_mass,
      aes(x = t_mid)
    ) +
      
      # ----------------------------------------------------------
    # Zero reference
    # ----------------------------------------------------------
    
    geom_hline(
      yintercept = 0,
      colour = "#777777",
      linewidth = 0.35
    ) +
      
      # ----------------------------------------------------------
    # 90% Monte Carlo envelope
    #
    # Deliberately light because finite differencing produces
    # substantially more local variability.
    # ----------------------------------------------------------
    
    geom_ribbon(
      aes(
        ymin = q05,
        ymax = q95
      ),
      fill = COL_IDIO,
      alpha = 0.06
    ) +
      
      # ----------------------------------------------------------
    # 50% Monte Carlo envelope
    # ----------------------------------------------------------
    
    geom_ribbon(
      aes(
        ymin = q25,
        ymax = q75
      ),
      fill = COL_IDIO,
      alpha = 0.21
    ) +
      
      # ----------------------------------------------------------
    # Monte Carlo median
    # ----------------------------------------------------------
    
    geom_line(
      aes(y = median),
      colour = COL_IDIO,
      linewidth = 1.05
    ) +
      
      # ----------------------------------------------------------
    # Population target
    # ----------------------------------------------------------
    
    geom_line(
      aes(y = truth),
      colour = "#111111",
      linewidth = 0.85
    ) +
      
      # ----------------------------------------------------------
    # Axes
    # ----------------------------------------------------------
    
    scale_x_continuous(
      breaks = c(
        0,
        0.25,
        0.50,
        0.75,
        1
      ),
      labels = c(
        "0",
        "0.25",
        "0.50",
        "0.75",
        "1"
      ),
      limits = c(0, 1),
      expand = expansion(
        mult = c(0.01, 0.01)
      )
    ) +
      
      coord_cartesian(
        ylim = y_mass
      ) +
      
      labs(
        title = paste0(
          panel_label,
          "  Interval-average signed mass"
        ),
        
        subtitle = paste0(
          "n = ",
          n_use
        ),
        
        x = expression(t),
        
        y = expression(
          Delta*mu[1] / Delta*t
        )
      ) +
      
      paper_theme(
        base_size = 11
      ) +
      
      theme(
        legend.position = "none",
        
        plot.title = element_text(
          face = "bold",
          size = 12,
          colour = "#111111",
          margin = margin(b = 2)
        ),
        
        plot.subtitle = element_text(
          face = "bold",
          size = 10.5,
          colour = COL_IDIO,
          margin = margin(b = 6)
        ),
        
        axis.title = element_text(
          size = 10.5
        ),
        
        axis.text = element_text(
          size = 9
        ),
        
        panel.grid.minor = element_blank(),
        
        panel.grid.major = element_line(
          colour = "#E6E6E6",
          linewidth = 0.25
        ),
        
        # Contour du graphique
        panel.border = element_rect(
          colour = "#555555",
          fill = NA,
          linewidth = 0.7
        ),
        
        plot.margin = margin(
          5, 8, 5, 5
        )
      )
  }
  
  
  # ==============================================================
  # 5. Four panels
  # ==============================================================
  
  p1 <- make_cumulative_panel(
    dat = dat_small,
    n_use = n_values[1],
    panel_label = "(a)"
  )
  
  p2 <- make_cumulative_panel(
    dat = dat_large,
    n_use = n_values[2],
    panel_label = "(b)"
  )
  
  p3 <- make_mass_panel(
    dat = dat_small,
    n_use = n_values[1],
    panel_label = "(c)"
  )
  
  p4 <- make_mass_panel(
    dat = dat_large,
    n_use = n_values[2],
    panel_label = "(d)"
  )
  
  
  # ==============================================================
  # 6. Publication layout
  #
  # No global title, subtitle, or caption inside the figure.
  # These belong in the manuscript caption.
  # ==============================================================
  
  final_plot <-
    (p1 + p2) /
    (p3 + p4) +
    patchwork::plot_layout(
      widths = c(1, 1),
      heights = c(1, 1)
    ) &
    theme(
      plot.margin = margin(4, 6, 4, 6)
    )
}


# ================================================================
# GENERATE ALL PAPER FIGURES
# ================================================================

make_paper_figures <- function(res) {

  if (!is.list(res) || !all(c("metrics", "curves") %in% names(res))) {
    stop("`res` must contain at least `metrics` and `curves`.")
  }

  # Figure 1 is reconstructed deterministically from the exact Monte Carlo
  # replication used in the manuscript. It does not use res$metrics/curves.
  fig1 <- figure1_truth_vs_estimated(
    scenario_use = "balanced",
    n_use = 500,
    rep_use = 88,
    margin_use = 2,
    final = TRUE
  )

  # Figure 2 uses the final Monte Carlo metrics.
  fig2 <- figure2_first_order_convergence(res)

  # Figure 3 requires 500 retained curve replications for n=100 and n=1000.
  curve_check <- res$curves |>
    dplyr::filter(
      tolower(trimws(as.character(scenario))) == "balanced",
      as.numeric(as.character(n)) %in% c(100, 1000),
      as.numeric(as.character(k)) == 1
    ) |>
    dplyr::mutate(
      n_check = as.numeric(as.character(n)),
      rep_check = as.numeric(as.character(rep))
    ) |>
    dplyr::group_by(n_check) |>
    dplyr::summarise(
      B = dplyr::n_distinct(rep_check),
      min_rep = min(rep_check),
      max_rep = max(rep_check),
      .groups = "drop"
    )

  ok_curve_check <-
    nrow(curve_check) == 2L &&
    all(sort(curve_check$n_check) == c(100, 1000)) &&
    all(curve_check$B == 500L) &&
    all(curve_check$min_rep == 1L) &&
    all(curve_check$max_rep == 500L)

  if (!ok_curve_check) {
    stop(
      paste0(
        "Figure 3 requires balanced, k=1 curves with exactly replications ",
        "1:500 for both n=100 and n=1000."
      )
    )
  }

  fig3 <- figure3_principal_measure_recovery(
    res = res,
    scenario_use = "balanced",
    n_values = c(100, 1000),
    k_use = 1
  )

  list(
    figure1 = fig1,
    figure2 = fig2,
    figure3 = fig3
  )
}


# ================================================================
# EXPORT PAPER FIGURES
# ================================================================

save_paper_figures <- function(res, root = ".") {

  fig_dir <- file.path(root, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  figs <- make_paper_figures(res)

  ggplot2::ggsave(
    filename = file.path(fig_dir, "Figure1_structural_decomposition_recovery.pdf"),
    plot = figs$figure1,
    width = 9.0,
    height = 3.5,
    units = "in"
  )

  ggplot2::ggsave(
    filename = file.path(fig_dir, "Figure2_first_order_convergence.pdf"),
    plot = figs$figure2,
    width = 9.0,
    height = 3.5,
    units = "in"
  )

  ggplot2::ggsave(
    filename = file.path(fig_dir, "Figure3_principal_measure_recovery.pdf"),
    plot = figs$figure3,
    width = 9.0,
    height = 4.6,
    units = "in"
  )

  invisible(figs)
}
