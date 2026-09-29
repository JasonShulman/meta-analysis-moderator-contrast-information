# Simulation Study 4 — frozen-design simulation and validation summaries.
#
# The published run reads data/study4_design/. It does not construct or freeze
# scenarios. Outcome seeds come from make_sim4_config(). Estimators come from
# estimators.R. Contrast metrics come from contrast_information.R.
# Common random numbers draw stochastic components once per replication.
# Beta1 datasets differ only by beta1 * X. Class E families that share
# outcome_scenario_id reuse those draws and scale the study effect by tau.

sim4_het_grid <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
sim4_tau_L <- 0.05
sim4_tau_U <- 0.30
sim4_fragility_omega <- 0.15
sim4_fragility_rho <- 0.50
sim4_rho_plan_primary <- 0.50
sim4_alpha_default <- 0.05

# Round x to the nearest element of {0, .05, ..., .30}, with ties rounded upward.
round_to_nearest_grid_s4 <- function(x, grid = sim4_het_grid) {
  x <- as.numeric(x)
  grid <- sort(as.numeric(grid))
  if (!length(grid)) {
    stop("grid must be non-empty.", call. = FALSE)
  }
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  for (i in which(ok)) {
    d <- abs(grid - x[[i]])
    m <- min(d)
    ties <- which(abs(d - m) <= 1e-12)
    # Ties round upward: take the largest tied grid value.
    out[[i]] <- max(grid[ties])
  }
  out
}

# Near-correct planning values (primary Study 4 regime).
# tau_p = round_0.05(tau_true), omega_p = round_0.05(omega_true), rho_p = .50.
near_correct_plan_s4 <- function(tau_true, omega_true, rho_plan = sim4_rho_plan_primary) {
  list(
    tau_plan = as.numeric(round_to_nearest_grid_s4(tau_true)),
    omega_plan = as.numeric(round_to_nearest_grid_s4(omega_true)),
    rho_plan = as.numeric(rho_plan)
  )
}

relative_difference_s4 <- function(a, b) {
  relative_difference_s3(a, b)
}

# Contrast total T_C and balance B_C from group weights.
# I_C* = T_C * B_C when both groups have positive finite weight.
# Frozen edge cases: B_C = 0 if either group has zero weight.
compute_contrast_balance_s4 <- function(W0, W1) {
  W0 <- as.numeric(W0)
  W1 <- as.numeric(W1)
  T_C <- W0 + W1
  B_C <- rep(NA_real_, length(T_C))
  status <- rep("ok", length(T_C))
  for (i in seq_along(T_C)) {
    w0 <- W0[[i]]
    w1 <- W1[[i]]
    tc <- T_C[[i]]
    if (!is.finite(w0) || !is.finite(w1) || w0 <= 0 || w1 <= 0) {
      B_C[[i]] <- 0
      if (!is.finite(tc) || !is.finite(w0) || !is.finite(w1)) {
        T_C[[i]] <- if (is.finite(w0) && is.finite(w1)) w0 + w1 else NA_real_
      }
      status[[i]] <- "insufficient_group_information"
    } else {
      B_C[[i]] <- (4 * w0 * w1) / (tc^2)
    }
  }
  data.frame(
    T_C = as.numeric(T_C),
    B_C = as.numeric(B_C),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}

# Bottleneck support J_min* and support balance B_J.
# Frozen edge cases: B_J = 0 if either side has zero effective support.
compute_support_diagnostics_s4 <- function(J0_star, J1_star) {
  J0_star <- as.numeric(J0_star)
  J1_star <- as.numeric(J1_star)
  n <- length(J0_star)
  J_min_star <- pmin(J0_star, J1_star)
  B_J <- rep(NA_real_, n)
  status <- rep("ok", n)
  for (i in seq_len(n)) {
    a <- J0_star[[i]]
    b <- J1_star[[i]]
    if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
      B_J[[i]] <- 0
      if (!is.finite(a) || !is.finite(b)) {
        J_min_star[[i]] <- NA_real_
      } else {
        J_min_star[[i]] <- min(a, b)
      }
      status[[i]] <- "zero_group_support"
    } else {
      B_J[[i]] <- (4 * a * b) / ((a + b)^2)
    }
  }
  data.frame(
    J0_star = J0_star,
    J1_star = J1_star,
    J_min_star = as.numeric(J_min_star),
    B_J = as.numeric(B_J),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}

# Information-retention ratio R_H = I_C,H* / I_C,S*.
# If I_C,S* = 0: R_H = NA, status = undefined_sampling_information.
# Values R_H > 1 are retained and not truncated.
compute_information_retention_s4 <- function(I_C_H_star, I_C_S_star) {
  I_C_H_star <- as.numeric(I_C_H_star)
  I_C_S_star <- as.numeric(I_C_S_star)
  n <- length(I_C_S_star)
  if (length(I_C_H_star) == 1L && n > 1L) {
    I_C_H_star <- rep(I_C_H_star, n)
  }
  R_H <- rep(NA_real_, n)
  status <- rep("ok", n)
  for (i in seq_len(n)) {
    den <- I_C_S_star[[i]]
    num <- I_C_H_star[[i]]
    if (!is.finite(den) || den <= 0) {
      R_H[[i]] <- NA_real_
      status[[i]] <- "undefined_sampling_information"
    } else {
      R_H[[i]] <- num / den
    }
  }
  data.frame(
    R_H = as.numeric(R_H),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
}

compute_sampling_information_s4 <- function(study_structure,
                                            rho_plan = sim4_rho_plan_primary,
                                            return_weights = FALSE) {
  compute_contrast_information_s3(
    study_structure = study_structure,
    rho_plan = rho_plan,
    weight_type = "sampling",
    return_weights = return_weights
  )
}

compute_heterogeneity_information_s4 <- function(study_structure,
                                                 tau_plan,
                                                 omega_plan,
                                                 rho_plan = sim4_rho_plan_primary,
                                                 return_weights = FALSE) {
  compute_contrast_information_s3(
    study_structure = study_structure,
    rho_plan = rho_plan,
    tau_plan = tau_plan,
    omega_plan = omega_plan,
    weight_type = "heterogeneity_adjusted",
    return_weights = return_weights
  )
}

enrich_contrast_row_s4 <- function(info) {
  if (is.data.frame(info)) {
    info <- as.list(info[1, , drop = FALSE])
  }
  bal <- compute_contrast_balance_s4(info$W0, info$W1)
  sup <- compute_support_diagnostics_s4(info$J0_star, info$J1_star)
  I_check <- bal$T_C[[1]] * bal$B_C[[1]]
  if (is.finite(info$I_C_star) && is.finite(I_check) &&
      abs(info$I_C_star - I_check) > 1e-8 * max(1, abs(info$I_C_star))) {
    stop("I_C* disagrees with T_C * B_C.", call. = FALSE)
  }
  data.frame(
    weight_type = as.character(info$weight_type),
    rho_plan = as.numeric(info$rho_plan),
    tau_plan = as.numeric(info$tau_plan),
    omega_plan = as.numeric(info$omega_plan),
    J0 = as.integer(info$J0),
    J1 = as.integer(info$J1),
    W0 = as.numeric(info$W0),
    W1 = as.numeric(info$W1),
    T_C = as.numeric(bal$T_C[[1]]),
    B_C = as.numeric(bal$B_C[[1]]),
    I_raw = as.numeric(info$I_raw),
    I_C_star = as.numeric(info$I_C_star),
    predicted_variance = as.numeric(info$predicted_variance),
    predicted_se = as.numeric(info$predicted_se),
    J0_star = as.numeric(sup$J0_star[[1]]),
    J1_star = as.numeric(sup$J1_star[[1]]),
    J_min_star = as.numeric(sup$J_min_star[[1]]),
    B_J = as.numeric(sup$B_J[[1]]),
    nu_C_star = as.numeric(info$nu_C_star),
    info_status = as.character(info$status),
    balance_status = as.character(bal$status[[1]]),
    support_status = as.character(sup$status[[1]]),
    stringsAsFactors = FALSE
  )
}

predict_power_s4 <- function(beta1,
                             predicted_se,
                             nu_C_star,
                             alpha = sim4_alpha_default) {
  predict_contrast_power_s3(
    beta1 = beta1,
    predicted_se = predicted_se,
    nu_C_star = nu_C_star,
    alpha = alpha
  )
}

# Fragility diagnostics F_I and Delta P_H under frozen (tau_L, tau_U, omega, rho).
# F_I = 1 - I_C,H*(tau_U) / I_C,H*(tau_L) with omega_p=.15, rho_p=.50.
# If the low-heterogeneity denominator is zero or non-finite: F_I = NA,
# status = undefined_fragility. No epsilon is added.
compute_fragility_diagnostics_s4 <- function(study_structure,
                                             beta1_levels = c(0.05, 0.10, 0.20),
                                             alpha = sim4_alpha_default) {
  low <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure,
      tau_plan = sim4_tau_L,
      omega_plan = sim4_fragility_omega,
      rho_plan = sim4_fragility_rho
    )
  )
  high <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure,
      tau_plan = sim4_tau_U,
      omega_plan = sim4_fragility_omega,
      rho_plan = sim4_fragility_rho
    )
  )
  den <- low$I_C_star[[1]]
  num <- high$I_C_star[[1]]
  if (!is.finite(den) || den <= 0) {
    F_I <- NA_real_
    fi_status <- "undefined_fragility"
  } else {
    F_I <- 1 - num / den
    fi_status <- "ok"
  }

  pred_low <- predict_power_s4(
    beta1 = beta1_levels,
    predicted_se = low$predicted_se[[1]],
    nu_C_star = low$nu_C_star[[1]],
    alpha = alpha
  )
  pred_high <- predict_power_s4(
    beta1 = beta1_levels,
    predicted_se = high$predicted_se[[1]],
    nu_C_star = high$nu_C_star[[1]],
    alpha = alpha
  )
  delta <- pred_low - pred_high
  names(pred_low) <- NULL
  names(pred_high) <- NULL
  names(delta) <- NULL

  idx10 <- which(abs(as.numeric(beta1_levels) - 0.10) < 1e-12)
  delta_p_h_10 <- if (length(idx10) == 1L) delta[[idx10]] else NA_real_
  p_low_10 <- if (length(idx10) == 1L) pred_low[[idx10]] else NA_real_
  p_high_10 <- if (length(idx10) == 1L) pred_high[[idx10]] else NA_real_

  idx05 <- which(abs(as.numeric(beta1_levels) - 0.05) < 1e-12)
  idx20 <- which(abs(as.numeric(beta1_levels) - 0.20) < 1e-12)

  data.frame(
    F_I = as.numeric(F_I),
    F_I_status = as.character(fi_status),
    I_C_H_star_tau_L = as.numeric(low$I_C_star[[1]]),
    I_C_H_star_tau_U = as.numeric(high$I_C_star[[1]]),
    predicted_se_tau_L = as.numeric(low$predicted_se[[1]]),
    predicted_se_tau_U = as.numeric(high$predicted_se[[1]]),
    nu_C_star_tau_L = as.numeric(low$nu_C_star[[1]]),
    nu_C_star_tau_U = as.numeric(high$nu_C_star[[1]]),
    delta_P_H_0.05 = if (length(idx05) == 1L) as.numeric(delta[[idx05]]) else NA_real_,
    delta_P_H_0.10 = as.numeric(delta_p_h_10),
    delta_P_H_0.20 = if (length(idx20) == 1L) as.numeric(delta[[idx20]]) else NA_real_,
    predicted_power_tau_L_0.10 = as.numeric(p_low_10),
    predicted_power_tau_U_0.10 = as.numeric(p_high_10),
    fragility_tau_L = sim4_tau_L,
    fragility_tau_U = sim4_tau_U,
    fragility_omega_plan = sim4_fragility_omega,
    fragility_rho_plan = sim4_fragility_rho,
    stringsAsFactors = FALSE
  )
}

combine_status_s4 <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- as.character(vals)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) {
    return("ok")
  }
  bad <- vals[vals != "ok"]
  if (!length(bad)) {
    return("ok")
  }
  paste(unique(bad), collapse = "|")
}

# Sampling-only, near-correct, and oracle quantities are stored separately.
compute_metric_profile_s4 <- function(study_structure,
                                      tau_true,
                                      omega_true,
                                      rho_true,
                                      rho_plan_primary = sim4_rho_plan_primary,
                                      alpha = sim4_alpha_default,
                                      scenario_id = NA_integer_) {
  samp <- enrich_contrast_row_s4(
    compute_sampling_information_s4(study_structure, rho_plan = rho_plan_primary)
  )
  nc <- near_correct_plan_s4(tau_true, omega_true, rho_plan_primary)
  het_nc <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure,
      tau_plan = nc$tau_plan,
      omega_plan = nc$omega_plan,
      rho_plan = nc$rho_plan
    )
  )
  het_or <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure,
      tau_plan = tau_true,
      omega_plan = omega_true,
      rho_plan = rho_true
    )
  )
  ret_nc <- compute_information_retention_s4(het_nc$I_C_star, samp$I_C_star)
  ret_or <- compute_information_retention_s4(het_or$I_C_star, samp$I_C_star)
  frag <- compute_fragility_diagnostics_s4(study_structure, alpha = alpha)
  comp <- comparator_quantities_s3(study_structure)

  pwr_s <- predict_power_s4(
    c(0, 0.05, 0.10, 0.20), samp$predicted_se, samp$nu_C_star, alpha
  )
  pwr_nc <- predict_power_s4(
    c(0, 0.05, 0.10, 0.20), het_nc$predicted_se, het_nc$nu_C_star, alpha
  )
  pwr_or <- predict_power_s4(
    c(0, 0.05, 0.10, 0.20), het_or$predicted_se, het_or$nu_C_star, alpha
  )

  primary_status <- combine_status_s4(
    het_nc$info_status, het_nc$balance_status, het_nc$support_status,
    ret_nc$status, frag$F_I_status
  )

  cbind(
    data.frame(
      scenario_id = as.integer(scenario_id),
      tau_true = as.numeric(tau_true),
      omega_true = as.numeric(omega_true),
      rho_true = as.numeric(rho_true),
      tau_plan_near_correct = as.numeric(nc$tau_plan),
      omega_plan_near_correct = as.numeric(nc$omega_plan),
      rho_plan_near_correct = as.numeric(nc$rho_plan),
      # Sampling-only comparator
      I_C_S_star = as.numeric(samp$I_C_star),
      T_C_S = as.numeric(samp$T_C),
      B_C_S = as.numeric(samp$B_C),
      J0_star_S = as.numeric(samp$J0_star),
      J1_star_S = as.numeric(samp$J1_star),
      J_min_star_S = as.numeric(samp$J_min_star),
      B_J_S = as.numeric(samp$B_J),
      nu_C_star_S = as.numeric(samp$nu_C_star),
      predicted_variance_S = as.numeric(samp$predicted_variance),
      predicted_se_S = as.numeric(samp$predicted_se),
      predicted_power_S_0.05 = as.numeric(pwr_s[[2]]),
      predicted_power_S_0.10 = as.numeric(pwr_s[[3]]),
      predicted_power_S_0.20 = as.numeric(pwr_s[[4]]),
      sampling_status = as.character(samp$info_status),
      # Near-correct heterogeneity-adjusted — PRIMARY
      I_C_H_star = as.numeric(het_nc$I_C_star),
      T_C = as.numeric(het_nc$T_C),
      B_C = as.numeric(het_nc$B_C),
      R_H = as.numeric(ret_nc$R_H),
      J0_star = as.numeric(het_nc$J0_star),
      J1_star = as.numeric(het_nc$J1_star),
      J_min_star = as.numeric(het_nc$J_min_star),
      B_J = as.numeric(het_nc$B_J),
      nu_C_star = as.numeric(het_nc$nu_C_star),
      predicted_variance_H = as.numeric(het_nc$predicted_variance),
      predicted_se_H = as.numeric(het_nc$predicted_se),
      predicted_power_H_0.05 = as.numeric(pwr_nc[[2]]),
      predicted_power_H_0.10 = as.numeric(pwr_nc[[3]]),
      predicted_power_H_0.20 = as.numeric(pwr_nc[[4]]),
      R_H_status = as.character(ret_nc$status),
      near_correct_info_status = as.character(het_nc$info_status),
      # Oracle heterogeneity-adjusted — SECONDARY benchmark
      I_C_H_star_oracle = as.numeric(het_or$I_C_star),
      T_C_oracle = as.numeric(het_or$T_C),
      B_C_oracle = as.numeric(het_or$B_C),
      R_H_oracle = as.numeric(ret_or$R_H),
      J0_star_oracle = as.numeric(het_or$J0_star),
      J1_star_oracle = as.numeric(het_or$J1_star),
      J_min_star_oracle = as.numeric(het_or$J_min_star),
      B_J_oracle = as.numeric(het_or$B_J),
      nu_C_star_oracle = as.numeric(het_or$nu_C_star),
      predicted_variance_oracle = as.numeric(het_or$predicted_variance),
      predicted_se_oracle = as.numeric(het_or$predicted_se),
      predicted_power_oracle_0.05 = as.numeric(pwr_or[[2]]),
      predicted_power_oracle_0.10 = as.numeric(pwr_or[[3]]),
      predicted_power_oracle_0.20 = as.numeric(pwr_or[[4]]),
      oracle_info_status = as.character(het_or$info_status),
      primary_status = as.character(primary_status),
      stringsAsFactors = FALSE
    ),
    frag,
    comp[, setdiff(names(comp), c("J0", "J1")), drop = FALSE]
  )
}

# Lightweight diagnostics used only for candidate-pool pair/family selection.
# Computes sampling-only and near-correct quantities. Does not compute the
# practical sensitivity grid (that is stored only for frozen scenarios).
compute_selection_diagnostics_s4 <- function(study_structure,
                                             tau_true,
                                             omega_true,
                                             rho_true,
                                             rho_plan_primary = sim4_rho_plan_primary) {
  samp <- enrich_contrast_row_s4(
    compute_sampling_information_s4(study_structure, rho_plan = rho_plan_primary)
  )
  nc <- near_correct_plan_s4(tau_true, omega_true, rho_plan_primary)
  het <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure,
      tau_plan = nc$tau_plan,
      omega_plan = nc$omega_plan,
      rho_plan = nc$rho_plan
    )
  )
  ret <- compute_information_retention_s4(het$I_C_star, samp$I_C_star)
  comp <- comparator_quantities_s3(study_structure)
  data.frame(
    I_C_S_star = as.numeric(samp$I_C_star),
    T_C_S = as.numeric(samp$T_C),
    B_C_S = as.numeric(samp$B_C),
    J0_star_S = as.numeric(samp$J0_star),
    J1_star_S = as.numeric(samp$J1_star),
    J_min_star_S = as.numeric(samp$J_min_star),
    B_J_S = as.numeric(samp$B_J),
    nu_C_star_S = as.numeric(samp$nu_C_star),
    sampling_status = as.character(samp$info_status),
    I_C_H_star = as.numeric(het$I_C_star),
    T_C = as.numeric(het$T_C),
    B_C = as.numeric(het$B_C),
    R_H = as.numeric(ret$R_H),
    J0_star = as.numeric(het$J0_star),
    J1_star = as.numeric(het$J1_star),
    J_min_star = as.numeric(het$J_min_star),
    B_J = as.numeric(het$B_J),
    nu_C_star = as.numeric(het$nu_C_star),
    R_H_status = as.character(ret$status),
    near_correct_status = as.character(het$info_status),
    tau_plan_near_correct = as.numeric(nc$tau_plan),
    omega_plan_near_correct = as.numeric(nc$omega_plan),
    N = as.numeric(comp$N),
    K = as.numeric(comp$K),
    J = as.integer(comp$J),
    J0 = as.integer(comp$J0),
    J1 = as.integer(comp$J1),
    minority_study_count = as.integer(comp$minority_study_count),
    mean_n = as.numeric(comp$mean_n),
    J_prec = as.numeric(comp$J_prec),
    J_k = as.numeric(comp$J_k),
    stringsAsFactors = FALSE
  )
}

# Frozen practical sensitivity-grid specification.
# Primary grid: tau_p in {0,.05,...,.30} × omega_p in {0,.10,.20,.30} × rho_p=.50.
# Secondary rho sensitivity: rho_p in {.20,.80} at omega_p=.15 across seven tau_p.
sim4_practical_grid_spec <- function() {
  primary <- expand.grid(
    tau_plan = sim4_het_grid,
    omega_plan = c(0, 0.10, 0.20, 0.30),
    rho_plan = sim4_rho_plan_primary,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  primary$grid_family <- "practical_tau_omega"
  primary$variant_id <- sprintf(
    "practical_tau%.2f_omega%.2f_rho%.2f",
    primary$tau_plan, primary$omega_plan, primary$rho_plan
  )
  rho_sens <- expand.grid(
    tau_plan = sim4_het_grid,
    omega_plan = sim4_fragility_omega,
    rho_plan = c(0.20, 0.80),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rho_sens$grid_family <- "rho_sensitivity"
  rho_sens$variant_id <- sprintf(
    "rho_sens_tau%.2f_omega%.2f_rho%.2f",
    rho_sens$tau_plan, rho_sens$omega_plan, rho_sens$rho_plan
  )
  rbind(primary, rho_sens)
}

# Sampling, near-correct, oracle, practical-grid, and rho-sensitivity
# regimes are kept as separate rows.
compute_all_metric_variants_s4 <- function(study_structure,
                                           tau_true,
                                           omega_true,
                                           rho_true,
                                           scenario_id = NA_integer_,
                                           alpha = sim4_alpha_default) {
  samp <- enrich_contrast_row_s4(
    compute_sampling_information_s4(study_structure, rho_plan = sim4_rho_plan_primary)
  )
  nc <- near_correct_plan_s4(tau_true, omega_true)
  het_nc <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure, nc$tau_plan, nc$omega_plan, nc$rho_plan
    )
  )
  het_or <- enrich_contrast_row_s4(
    compute_heterogeneity_information_s4(
      study_structure, tau_true, omega_true, rho_true
    )
  )
  ret_nc <- compute_information_retention_s4(het_nc$I_C_star, samp$I_C_star)
  ret_or <- compute_information_retention_s4(het_or$I_C_star, samp$I_C_star)

  pack <- function(variant_id, regime, is_primary, is_oracle, row, R_H, R_H_status) {
    pwr <- predict_power_s4(
      c(0.05, 0.10, 0.20), row$predicted_se, row$nu_C_star, alpha
    )
    data.frame(
      scenario_id = as.integer(scenario_id),
      variant_id = as.character(variant_id),
      regime = as.character(regime),
      is_primary = isTRUE(is_primary),
      is_oracle = isTRUE(is_oracle),
      weight_type = as.character(row$weight_type),
      rho_plan = as.numeric(row$rho_plan),
      tau_plan = as.numeric(row$tau_plan),
      omega_plan = as.numeric(row$omega_plan),
      W0 = as.numeric(row$W0),
      W1 = as.numeric(row$W1),
      T_C = as.numeric(row$T_C),
      B_C = as.numeric(row$B_C),
      I_C_star = as.numeric(row$I_C_star),
      R_H = as.numeric(R_H),
      J0_star = as.numeric(row$J0_star),
      J1_star = as.numeric(row$J1_star),
      J_min_star = as.numeric(row$J_min_star),
      B_J = as.numeric(row$B_J),
      nu_C_star = as.numeric(row$nu_C_star),
      predicted_variance = as.numeric(row$predicted_variance),
      predicted_se = as.numeric(row$predicted_se),
      predicted_power_0.05 = as.numeric(pwr[[1]]),
      predicted_power_0.10 = as.numeric(pwr[[2]]),
      predicted_power_0.20 = as.numeric(pwr[[3]]),
      info_status = as.character(row$info_status),
      R_H_status = as.character(R_H_status),
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    pack("sampling_rho_0.50", "sampling_only", FALSE, FALSE, samp, NA_real_, NA_character_),
    pack("near_correct", "near_correct", TRUE, FALSE, het_nc, ret_nc$R_H, ret_nc$status),
    pack("oracle", "oracle", FALSE, TRUE, het_or, ret_or$R_H, ret_or$status)
  )

  grid <- sim4_practical_grid_spec()
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, , drop = FALSE]
    het <- enrich_contrast_row_s4(
      compute_heterogeneity_information_s4(
        study_structure,
        tau_plan = g$tau_plan[[1]],
        omega_plan = g$omega_plan[[1]],
        rho_plan = g$rho_plan[[1]]
      )
    )
    ret <- compute_information_retention_s4(het$I_C_star, samp$I_C_star)
    rows[[length(rows) + 1L]] <- pack(
      g$variant_id[[1]], g$grid_family[[1]], FALSE, FALSE, het, ret$R_H, ret$status
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The frozen library uses the Study 3 n_j/k_j bounds and integerization.
# Targeted scenarios are 65-96 (Classes A-D); Class E occupies 97-120.
# Within a Class E family, outcome_scenario_id is the lowest-tau scenario.

# =============================================================================
# Paths
# =============================================================================

sim4_project_root <- function(root = NULL) {
  if (!is.null(root)) {
    return(normalizePath(root, winslash = "/", mustWork = TRUE))
  }
  candidates <- c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../.."))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "study4_functions.R"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot locate reproducibility root (expected R/study4_functions.R).", call. = FALSE)
}

sim4_fixture_dir <- function(root = NULL) {
  file.path(sim4_project_root(root), "data", "study4_design")
}

eq_na_false_s4 <- function(x, value) {
  !is.na(x) & x == value
}

# =============================================================================
# Configuration and seeds
# =============================================================================

sim4_stream_codes <- function() {
  c(u_j = 1L, w_ij = 2L, e_ij = 3L, validation = 4L)
}

is_legal_seed_s4 <- function(seed) {
  seed <- suppressWarnings(as.numeric(seed))
  length(seed) == 1L &&
    is.finite(seed) &&
    seed >= 0 &&
    seed <= .Machine$integer.max &&
    abs(seed - round(seed)) < 1e-9
}

assert_legal_seed_s4 <- function(seed, label = "seed") {
  if (!is_legal_seed_s4(seed)) {
    stop(label, " is not a legal R integer seed: ", seed, call. = FALSE)
  }
  as.integer(round(as.numeric(seed)))
}

make_sim4_config <- function(run_mode = c("smoke", "pilot", "production"),
                             workers = NULL,
                             overwrite = FALSE,
                             resume = TRUE) {
  run_mode <- match.arg(run_mode)
  mode_defaults <- switch(
    run_mode,
    smoke = list(
      n_rep_total = 5L,
      chunk_size = 5L,
      workers = 1L
    ),
    pilot = list(
      n_rep_total = 20L,
      chunk_size = 10L,
      workers = min(3L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
    ),
    production = list(
      n_rep_total = 5000L,
      chunk_size = 500L,
      workers = NA_integer_
    )
  )
  workers_out <- if (is.null(workers)) mode_defaults$workers else as.integer(workers)
  list(
    study_id = 4L,
    run_mode = run_mode,
    beta0 = 0.20,
    beta1_levels = c(0.00, 0.05, 0.10, 0.20),
    confirmatory_beta1 = 0.10,
    alpha = 0.05,
    ci_level = 0.95,
    rho_plan_primary = 0.50,
    n_min = 20L,
    n_max = 3000L,
    k_min = 1L,
    k_max = 30L,
    J_levels = c(24L, 36L, 48L, 72L, 96L, 120L, 144L),
    p1_range = c(0.10, 0.50),
    mean_n_range = c(40, 500),
    n_cv_range = c(0, 0.90),
    mean_k_range = c(2, 10),
    k_cv_range = c(0, 0.90),
    lambda_range = c(0, 1),
    rho_true_range = c(0.10, 0.90),
    tau_range = c(0.02, 0.30),
    omega_range = c(0.02, 0.30),
    nk_relationships = c("independent", "positive", "negative"),
    alignment_directions = c("independent", "minority_favored", "majority_favored"),
    n_broad = 64L,
    n_targeted = 56L,
    n_scenarios = 120L,
    n_candidate_pool = 100000L,
    n_candidate_pool_expanded = 250000L,
    # Construction seeds for the frozen library, not outcome-simulation seeds.
    pair_candidate_seed = 8253026L,
    pair_candidate_seed_expanded = 8253027L,
    scenario_seed_root = 8252026L,
    n_rep_total = as.integer(mode_defaults$n_rep_total),
    chunk_size = as.integer(mode_defaults$chunk_size),
    workers = workers_out,
    seed_root_outcomes = 1100000000L,
    scenario_stride = 1000000L,
    replication_stride = 100L,
    methods = c("mlma", "mlma_cr2"),
    n_methods = 2L,
    type1_band = c(0.025, 0.075),
    failure_flag = 0.01,
    warning_separator = " || ",
    boundary_tolerance = 1e-8,
    # Frozen pair tolerances
    pair_a_ic_s_tol = 0.05,
    pair_a_n_tol = 0.02,
    pair_a_k_tol = 0.02,
    pair_a_rh_sep = 0.25,
    pair_b_tc_tol = 0.05,
    pair_b_jmin_tol = 0.10,
    pair_b_bj_tol = 0.10,
    pair_b_bc_sep = 0.25,
    pair_c_ich_tol = 0.05,
    pair_c_tc_tol = 0.10,
    pair_c_jmin_ratio = 2,
    pair_c_bj_sep = 0.20,
    pair_d_jmin_tol = 0.05,
    pair_d_bj_tol = 0.05,
    pair_d_ich_sep = 0.50,
    class_e_tau = c(0.05, 0.12, 0.20, 0.30),
    class_e_omega = 0.15,
    n_class_e_families = 6L,
    n_pairs_per_class = 4L,
    max_broad_seed_attempts = 100L,
    overwrite = isTRUE(overwrite),
    resume = isTRUE(resume),
    method_row_multiplier = 2L
  )
}

# Deterministic Study 4 outcome-simulation seed.

# seed = seed_root_outcomes + scenario_id * scenario_stride
#        + rep_id * replication_stride + stream_id
sim4_seed <- function(scenario_id,
                      rep_id,
                      stream_id,
                      config = make_sim4_config("smoke")) {
  scenario_id <- as.integer(scenario_id)
  rep_id <- as.integer(rep_id)
  if (is.character(stream_id)) {
    codes <- sim4_stream_codes()
    if (!stream_id %in% names(codes)) {
      stop("Unknown stream '", stream_id, "'.", call. = FALSE)
    }
    stream_id <- codes[[stream_id]]
  }
  stream_id <- as.integer(stream_id)
  if (length(scenario_id) != 1L || is.na(scenario_id) || scenario_id < 1L) {
    stop("scenario_id must be a positive integer.", call. = FALSE)
  }
  if (length(rep_id) != 1L || is.na(rep_id) || rep_id < 1L) {
    stop("rep_id must be a positive integer.", call. = FALSE)
  }
  raw <- as.numeric(config$seed_root_outcomes) +
    as.numeric(scenario_id) * as.numeric(config$scenario_stride) +
    as.numeric(rep_id) * as.numeric(config$replication_stride) +
    as.numeric(stream_id)
  assert_legal_seed_s4(raw, "derived Study 4 outcome seed")
}

run_with_sim4_seed <- function(seed, code) {
  seed <- assert_legal_seed_s4(seed, "run_with_sim4_seed")
  if (exists("run_with_sim3_seed", mode = "function")) {
    return(run_with_sim3_seed(seed, code))
  }
  if (exists("run_with_sim2_seed", mode = "function")) {
    return(run_with_sim2_seed(seed, code))
  }
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop("Package 'withr' is required for Study 4 seeding.", call. = FALSE)
  }
  withr::with_seed(
    seed,
    code,
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
}

# =============================================================================
# Targeted pair selection
# =============================================================================

normalize_pair_row_s4 <- function(row) {
  if (is.data.frame(row)) {
    row <- as.list(row[1, , drop = FALSE])
  }
  if ((is.null(row$N) || (length(row$N) == 1L && is.na(row$N))) &&
      !is.null(row$N_total)) {
    row$N <- row$N_total
  }
  if ((is.null(row$K) || (length(row$K) == 1L && is.na(row$K))) &&
      !is.null(row$K_total)) {
    row$K <- row$K_total
  }
  row
}

pair_qualifies_class_s4 <- function(a, b, class_label, config) {
  a <- normalize_pair_row_s4(a)
  b <- normalize_pair_row_s4(b)
  class_label <- as.character(class_label)
  finite_ok <- function(x) all(is.finite(as.numeric(x)))
  if (identical(class_label, "A")) {
    if (!finite_ok(c(a$I_C_S_star, b$I_C_S_star, a$R_H, b$R_H, a$N, b$N, a$K, b$K))) {
      return(FALSE)
    }
    if (!identical(as.integer(a$J), as.integer(b$J))) return(FALSE)
    if (abs(as.integer(a$J0) - as.integer(b$J0)) > 1L) return(FALSE)
    if (abs(as.integer(a$J1) - as.integer(b$J1)) > 1L) return(FALSE)
    if (!isTRUE(relative_difference_s4(a$I_C_S_star, b$I_C_S_star) <= config$pair_a_ic_s_tol)) {
      return(FALSE)
    }
    if (!isTRUE(relative_difference_s4(a$N, b$N) <= config$pair_a_n_tol)) return(FALSE)
    if (!isTRUE(relative_difference_s4(a$K, b$K) <= config$pair_a_k_tol)) return(FALSE)
    return(isTRUE(abs(a$R_H - b$R_H) >= config$pair_a_rh_sep))
  }
  if (identical(class_label, "B")) {
    if (!finite_ok(c(a$T_C, b$T_C, a$J_min_star, b$J_min_star, a$B_J, b$B_J, a$B_C, b$B_C))) {
      return(FALSE)
    }
    if (!isTRUE(relative_difference_s4(a$T_C, b$T_C) <= config$pair_b_tc_tol)) return(FALSE)
    if (!isTRUE(relative_difference_s4(a$J_min_star, b$J_min_star) <= config$pair_b_jmin_tol)) {
      return(FALSE)
    }
    if (!isTRUE(abs(a$B_J - b$B_J) <= config$pair_b_bj_tol)) return(FALSE)
    return(isTRUE(abs(a$B_C - b$B_C) >= config$pair_b_bc_sep))
  }
  if (identical(class_label, "C")) {
    if (!finite_ok(c(a$I_C_H_star, b$I_C_H_star, a$T_C, b$T_C, a$J_min_star, b$J_min_star, a$B_J, b$B_J))) {
      return(FALSE)
    }
    jmin <- c(a$J_min_star, b$J_min_star)
    if (min(jmin) <= 0) return(FALSE)
    if (!isTRUE(relative_difference_s4(a$I_C_H_star, b$I_C_H_star) <= config$pair_c_ich_tol)) {
      return(FALSE)
    }
    if (!isTRUE(relative_difference_s4(a$T_C, b$T_C) <= config$pair_c_tc_tol)) return(FALSE)
    ratio <- max(jmin) / min(jmin)
    if (!isTRUE(ratio >= config$pair_c_jmin_ratio)) return(FALSE)
    return(isTRUE(abs(a$B_J - b$B_J) >= config$pair_c_bj_sep))
  }
  if (identical(class_label, "D")) {
    if (!finite_ok(c(a$J_min_star, b$J_min_star, a$B_J, b$B_J, a$I_C_H_star, b$I_C_H_star))) {
      return(FALSE)
    }
    if (!isTRUE(relative_difference_s4(a$J_min_star, b$J_min_star) <= config$pair_d_jmin_tol)) {
      return(FALSE)
    }
    if (!isTRUE(abs(a$B_J - b$B_J) <= config$pair_d_bj_tol)) return(FALSE)
    return(isTRUE(relative_difference_s4(a$I_C_H_star, b$I_C_H_star) >= config$pair_d_ich_sep))
  }
  FALSE
}

# =============================================================================
# Class E families
# =============================================================================

class_e_profile_eligibility_s4 <- function(d, profile, q_ic, q_jmin) {
  profile <- as.integer(profile)
  ic <- as.numeric(d$I_C_S_star)
  bc <- as.numeric(d$B_C_S)
  bj <- as.numeric(d$B_J_S)
  jmin <- as.numeric(d$J_min_star_S)
  lambda <- as.numeric(d$alignment_strength)
  align <- as.character(d$alignment_direction)
  low_ic <- is.finite(ic) & ic <= q_ic[[1]]
  high_ic <- is.finite(ic) & ic >= q_ic[[2]]
  mid_ic <- is.finite(ic) & ic > q_ic[[1]] & ic < q_ic[[2]]
  balanced_c <- is.finite(bc) & bc >= 0.80
  poor_c <- is.finite(bc) & bc <= 0.50
  balanced_s <- is.finite(bj) & bj >= 0.80
  bottleneck <- is.finite(jmin) & jmin <= q_jmin & is.finite(bj) & bj <= 0.60
  if (profile == 1L) {
    return(low_ic & balanced_c & balanced_s)
  }
  if (profile == 2L) {
    return(mid_ic & balanced_c & balanced_s)
  }
  if (profile == 3L) {
    return(high_ic & balanced_c & balanced_s)
  }
  if (profile == 4L) {
    return(high_ic & poor_c)
  }
  if (profile == 5L) {
    return(high_ic & bottleneck)
  }
  if (profile == 6L) {
    return(high_ic & is.finite(lambda) & lambda >= 0.75 & align != "independent")
  }
  stop("Unknown Class E profile: ", profile, call. = FALSE)
}

# =============================================================================
# Freeze / validate
# =============================================================================

# TRUE for Study 4 text/CSV design artifacts that must be hashed canonically.
is_s4_text_design_file <- function(path) {
  grepl("\\.(csv|txt|md|tsv)$", path, ignore.case = TRUE, perl = FALSE)
}

# Canonical UTF-8/LF bytes of a text file (CRLF and CR become LF; BOM stripped).

# Frozen Study 4 design hashes are SHA-256 of these bytes, not of raw on-disk
# bytes. Git `core.autocrlf=true` on Windows must not fail validation.
canonical_text_bytes_s4 <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }
  n <- file.info(path, extra_cols = FALSE)$size
  if (!is.finite(n) || n < 0) {
    stop("Cannot read file size: ", path, call. = FALSE)
  }
  raw <- if (n == 0) {
    raw(0)
  } else {
    readBin(path, what = "raw", n = as.integer(n))
  }
  if (length(raw) >= 3L &&
      raw[[1]] == as.raw(0xEF) &&
      raw[[2]] == as.raw(0xBB) &&
      raw[[3]] == as.raw(0xBF)) {
    raw <- raw[-c(1L, 2L, 3L)]
  }
  txt <- if (!length(raw)) {
    ""
  } else {
    rawToChar(raw)
  }
  Encoding(txt) <- "UTF-8"
  txt <- gsub("\r\n", "\n", txt, fixed = TRUE)
  txt <- gsub("\r", "\n", txt, fixed = TRUE)
  bytes <- iconv(txt, from = "UTF-8", to = "UTF-8", toRaw = TRUE)[[1]]
  if (is.null(bytes)) {
    raw(0)
  } else {
    bytes
  }
}

# SHA-256 of canonical UTF-8/LF text bytes. Used for frozen design CSVs.
hash_canonical_text_file_s4 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  bytes <- canonical_text_bytes_s4(path)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for Study 4 design hashing.", call. = FALSE)
  }
  digest::digest(bytes, algo = "sha256", serialize = FALSE)
}

# SHA-256 of a file. Text/CSV design artifacts are canonicalized; RDS stay raw.
hash_file_s4 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  if (is_s4_text_design_file(path)) {
    return(hash_canonical_text_file_s4(path))
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(file = path, algo = "sha256")
  } else {
    paste(unname(tools::md5sum(path)))
  }
}

read_frozen_scenario_library_s4 <- function(root = NULL) {
  dir <- sim4_fixture_dir(root)
  required <- c(
    "scenario_library.csv",
    "broad_scenarios.csv",
    "targeted_scenarios.csv",
    "study_structure.csv",
    "diagnostic_library.csv",
    "targeted_pair_diagnostics.csv",
    "fragility_family_diagnostics.csv",
    "coverage_report.csv",
    "design_hashes.csv",
    "study4_config.csv"
  )
  missing <- required[!file.exists(file.path(dir, required))]
  if (length(missing)) {
    stop(
      "Frozen Study 4 scenario library is incomplete. Missing: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  list(
    dir = dir,
    grid = utils::read.csv(file.path(dir, "scenario_library.csv"), stringsAsFactors = FALSE),
    broad = utils::read.csv(file.path(dir, "broad_scenarios.csv"), stringsAsFactors = FALSE),
    targeted = utils::read.csv(file.path(dir, "targeted_scenarios.csv"), stringsAsFactors = FALSE),
    study_structure = utils::read.csv(
      file.path(dir, "study_structure.csv"), stringsAsFactors = FALSE
    ),
    metrics = utils::read.csv(
      file.path(dir, "diagnostic_library.csv"), stringsAsFactors = FALSE
    ),
    pair_diagnostics = utils::read.csv(
      file.path(dir, "targeted_pair_diagnostics.csv"), stringsAsFactors = FALSE
    ),
    family_diagnostics = utils::read.csv(
      file.path(dir, "fragility_family_diagnostics.csv"), stringsAsFactors = FALSE
    ),
    coverage = utils::read.csv(file.path(dir, "coverage_report.csv"), stringsAsFactors = FALSE),
    hashes = utils::read.csv(file.path(dir, "design_hashes.csv"), stringsAsFactors = FALSE),
    config_df = utils::read.csv(file.path(dir, "study4_config.csv"), stringsAsFactors = FALSE)
  )
}

verify_design_hashes_s4 <- function(root = NULL) {
  lib <- read_frozen_scenario_library_s4(root)
  dir <- lib$dir
  issues <- character(0)
  files <- lib$hashes$file[grepl("\\.csv$", lib$hashes$file)]
  for (f in files) {
    stored <- lib$hashes$sha256[lib$hashes$file == f]
    live <- hash_canonical_text_file_s4(file.path(dir, f))
    if (!identical(as.character(stored), as.character(live))) {
      issues <- c(issues, sprintf("Hash mismatch for %s.", f))
    }
  }
  list(ok = length(issues) == 0L, issues = issues)
}

validate_frozen_scenario_library_s4 <- function(root = NULL,
                                                config = make_sim4_config("smoke")) {
  lib <- read_frozen_scenario_library_s4(root)
  issues <- character(0)
  add <- function(msg) issues <<- c(issues, msg)

  g <- lib$grid
  ss <- lib$study_structure
  if (nrow(g) != 120L) add(sprintf("Expected 120 scenarios, found %d.", nrow(g)))
  n_broad <- sum(g$scenario_set == "broad")
  n_targ <- sum(g$scenario_set == "targeted")
  if (n_broad != 64L) add(sprintf("Expected 64 broad scenarios, found %d.", n_broad))
  if (n_targ != 56L) add(sprintf("Expected 56 targeted scenarios, found %d.", n_targ))
  if (length(unique(g$scenario_id)) != nrow(g)) add("Duplicate scenario_id values.")

  for (cl in c("A", "B", "C", "D")) {
    n_cl <- sum(g$targeted_class == cl, na.rm = TRUE)
    if (n_cl != 8L) add(sprintf("Class %s expected 8 scenarios, found %d.", cl, n_cl))
  }
  n_e <- sum(g$targeted_class == "E", na.rm = TRUE)
  if (n_e != 24L) add(sprintf("Class E expected 24 scenarios, found %d.", n_e))

  for (i in seq_len(nrow(g))) {
    row <- g[i, ]
    ssi <- ss[ss$scenario_id == row$scenario_id, , drop = FALSE]
    if (nrow(ssi) != row$J) {
      add(sprintf("scenario %d: structure rows != J.", row$scenario_id))
      next
    }
    if (sum(ssi$n_j) != row$N_total) {
      add(sprintf("scenario %d: sum(n_j) != N_total.", row$scenario_id))
    }
    if (sum(ssi$k_j) != row$K_total) {
      add(sprintf("scenario %d: sum(k_j) != K_total.", row$scenario_id))
    }
    if (sum(ssi$X == 1L) != row$J1) {
      add(sprintf("scenario %d: minority count mismatch.", row$scenario_id))
    }
    if (sum(ssi$X == 0L) < 1L || sum(ssi$X == 1L) < 1L) {
      add(sprintf("scenario %d: empty moderator group.", row$scenario_id))
    }
    if (any(ssi$n_j < config$n_min | ssi$n_j > config$n_max)) {
      add(sprintf("scenario %d: n_j outside integerization bounds.", row$scenario_id))
    }
    if (any(ssi$k_j < config$k_min | ssi$k_j > config$k_max)) {
      add(sprintf("scenario %d: k_j outside integerization bounds.", row$scenario_id))
    }
    if (identical(row$scenario_set, "broad")) {
      if (!row$J %in% config$J_levels) {
        add(sprintf("broad scenario %d: J not in allowed set.", row$scenario_id))
      }
    }
    recomputed <- compute_selection_diagnostics_s4(
      ssi, row$tau, row$omega, row$rho_true
    )
    if (is.finite(row$I_C_H_star) && is.finite(recomputed$I_C_H_star) &&
        abs(row$I_C_H_star - recomputed$I_C_H_star) > 1e-8 * max(1, abs(row$I_C_H_star))) {
      add(sprintf("scenario %d: stored I_C,H* disagrees with recomputation.", row$scenario_id))
    }
  }

  pd <- lib$pair_diagnostics
  if (nrow(pd) != 16L) add(sprintf("Expected 16 targeted pairs, found %d.", nrow(pd)))
  if (anyDuplicated(c(pd$candidate_id_a, pd$candidate_id_b))) {
    add("Candidate reuse across Classes A-D.")
  }
  for (i in seq_len(nrow(pd))) {
    if (!isTRUE(pd$pass_all_tolerances[[i]])) {
      add(sprintf("Pair class %s id %s fails frozen tolerances.",
                  pd$pair_class[[i]], pd$pair_id[[i]]))
    }
    a <- g[g$scenario_id == pd$scenario_id_a[[i]], , drop = FALSE]
    b <- g[g$scenario_id == pd$scenario_id_b[[i]], , drop = FALSE]
    if (!nrow(a) || !nrow(b)) {
      add(sprintf("Pair %s missing scenario rows.", i))
      next
    }
    aa <- as.list(a[1, ])
    bb <- as.list(b[1, ])
    if (!pair_qualifies_class_s4(aa, bb, pd$pair_class[[i]], config)) {
      add(sprintf("Recomputed pair class %s id %s fails tolerances.",
                  pd$pair_class[[i]], pd$pair_id[[i]]))
    }
  }

  fd <- lib$family_diagnostics
  if (nrow(fd) != 6L) add(sprintf("Expected 6 Class E families, found %d.", nrow(fd)))
  if (any(fd$base_candidate_id %in% c(pd$candidate_id_a, pd$candidate_id_b))) {
    add("Class E base candidate reused in Classes A-D.")
  }
  if (anyDuplicated(fd$base_candidate_id)) add("Duplicate Class E base candidates.")
  for (f in seq_len(nrow(fd))) {
    d_el <- data.frame(
      I_C_S_star = fd$I_C_S_star[[f]],
      B_C_S = fd$B_C_S[[f]],
      B_J_S = fd$B_J_S[[f]],
      J_min_star_S = fd$J_min_star_S[[f]],
      alignment_strength = fd$alignment_strength[[f]],
      alignment_direction = fd$alignment_direction[[f]],
      stringsAsFactors = FALSE
    )
    q_ic <- c(fd$I_C_S_quartile_low[[f]], fd$I_C_S_quartile_high[[f]])
    if (!isTRUE(class_e_profile_eligibility_s4(d_el, fd$profile[[f]], q_ic, fd$J_min_S_q25[[f]]))) {
      add(sprintf("Family %s fails its frozen profile eligibility rule.", fd$family_id[[f]]))
    }
  }
  e <- g[eq_na_false_s4(g$targeted_class, "E"), , drop = FALSE]
  for (f in seq_len(nrow(fd))) {
    fam <- e[e$family_id == fd$family_id[[f]], , drop = FALSE]
    if (nrow(fam) != 4L) add(sprintf("Family %s does not have 4 tau levels.", fd$family_id[[f]]))
    taus <- sort(unique(fam$tau))
    if (!isTRUE(all.equal(taus, sort(config$class_e_tau)))) {
      add(sprintf("Family %s tau levels are not the frozen set.", fd$family_id[[f]]))
    }
    if (any(abs(fam$omega - config$class_e_omega) > 1e-12)) {
      add(sprintf("Family %s omega is not 0.15.", fd$family_id[[f]]))
    }
    if (length(unique(fam$outcome_scenario_id)) != 1L) {
      add(sprintf("Family %s CRN outcome_scenario_id is not shared.", fd$family_id[[f]]))
    }
  }

  cov <- lib$coverage
  if ("pass" %in% names(cov) && any(!as.logical(cov$pass))) {
    add("Stored coverage report contains failing checks.")
  }

  hash_chk <- verify_design_hashes_s4(root)
  if (!isTRUE(hash_chk$ok)) issues <- c(issues, hash_chk$issues)

  list(ok = length(issues) == 0L, issues = issues, library = lib)
}

# Simulation Study 4 — DGM with frozen study structure and common random numbers.
# Structure is never regenerated inside the replication loop.
#
# CRN across beta1: one stochastic draw; only the fixed moderator shift changes.
# CRN across Class E tau: standardized u_j ~ N(0,1) drawn under outcome_scenario_id
# (the family's first scenario) then scaled by true tau. w_ij and e_ij are reused.

# Streams: u_j=1, w_ij=2, e_ij=3. Moderator-effect value is not in the seed.
# Study-level heterogeneity is drawn as standard normals and scaled by tau so
# Class E tau levels that share outcome_scenario_id are exact CRN rescalings.
simulate_base_components_s4 <- function(study_structure,
                                        scenario,
                                        rep_id,
                                        config = make_sim4_config("smoke")) {
  ss <- study_structure[order(study_structure$study), , drop = FALSE]
  J <- nrow(ss)
  k_j <- as.integer(ss$k_j)
  n_j <- as.numeric(ss$n_j)
  if ("outcome_scenario_id" %in% names(scenario) &&
      length(scenario$outcome_scenario_id) &&
      !is.na(scenario$outcome_scenario_id[[1]])) {
    seed_sid <- as.integer(scenario$outcome_scenario_id[[1]])
  } else {
    seed_sid <- as.integer(scenario$scenario_id)
  }
  tau <- as.numeric(scenario$tau)
  omega <- as.numeric(scenario$omega)
  rho <- as.numeric(scenario$rho_true)
  rep_id <- as.integer(rep_id)

  z_u <- run_with_sim4_seed(
    sim4_seed(seed_sid, rep_id, "u_j", config),
    stats::rnorm(J, mean = 0, sd = 1)
  )
  u_j <- as.numeric(z_u) * tau

  K <- sum(k_j)
  w_all <- run_with_sim4_seed(
    sim4_seed(seed_sid, rep_id, "w_ij", config),
    stats::rnorm(K, mean = 0, sd = omega)
  )

  e_list <- run_with_sim4_seed(
    sim4_seed(seed_sid, rep_id, "e_ij", config),
    {
      if (!exists("make_sampling_Sigma", mode = "function")) {
        stop("make_sampling_Sigma() is unavailable. Source R/study1_functions.R.", call. = FALSE)
      }
      out <- vector("list", J)
      for (j in seq_len(J)) {
        kj_j <- k_j[[j]]
        vi_j <- rep(1 / (n_j[[j]] - 3), kj_j)
        Sigma_j <- make_sampling_Sigma(vi_j, rho)
        if (kj_j == 1L) {
          out[[j]] <- stats::rnorm(1L, mean = 0, sd = sqrt(vi_j))
        } else {
          out[[j]] <- as.numeric(
            MASS::mvrnorm(n = 1L, mu = rep(0, kj_j), Sigma = Sigma_j)
          )
        }
      }
      out
    }
  )

  list(
    z_u = as.numeric(z_u),
    u_j = as.numeric(u_j),
    w_all = as.numeric(w_all),
    e_list = e_list,
    study_structure = ss,
    outcome_scenario_id = seed_sid,
    tau = tau,
    omega = omega,
    rho = rho
  )
}

assemble_dataset_s4 <- function(components,
                                scenario,
                                beta1,
                                rep_id,
                                config = make_sim4_config("smoke")) {
  ss <- components$study_structure
  J <- nrow(ss)
  beta0 <- as.numeric(config$beta0)
  beta1 <- as.numeric(beta1)
  rows <- vector("list", J)
  w_offset <- 0L
  for (j in seq_len(J)) {
    kj_j <- as.integer(ss$k_j[[j]])
    n_jj <- as.numeric(ss$n_j[[j]])
    X_j <- as.integer(ss$X[[j]])
    vi_j <- rep(1 / (n_jj - 3), kj_j)
    w_ij <- components$w_all[seq.int(w_offset + 1L, w_offset + kj_j)]
    w_offset <- w_offset + kj_j
    e_ij <- components$e_list[[j]]
    y0 <- beta0 + components$u_j[[j]] + w_ij + e_ij
    yi <- y0 + beta1 * X_j
    rows[[j]] <- data.frame(
      scenario_id = as.integer(scenario$scenario_id),
      scenario_set = as.character(scenario$scenario_set),
      targeted_class = if ("targeted_class" %in% names(scenario)) {
        as.character(scenario$targeted_class)
      } else {
        NA_character_
      },
      pair_id = if ("pair_id" %in% names(scenario)) as.integer(scenario$pair_id) else NA_integer_,
      pair_arm = if ("pair_arm" %in% names(scenario)) as.character(scenario$pair_arm) else NA_character_,
      family_id = if ("family_id" %in% names(scenario)) as.integer(scenario$family_id) else NA_integer_,
      tau_family_level = if ("tau_family_level" %in% names(scenario)) {
        as.numeric(scenario$tau_family_level)
      } else {
        NA_real_
      },
      rep_id = as.integer(rep_id),
      beta1_true = beta1,
      study = as.integer(ss$study[[j]]),
      es_id = seq_len(kj_j),
      X = X_j,
      k_j = kj_j,
      n_j = n_jj,
      vi = vi_j,
      yi = as.numeric(yi),
      y0 = as.numeric(y0),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

simulate_beta_crn_s4 <- function(scenario,
                                 study_structure,
                                 rep_id,
                                 config = make_sim4_config("smoke"),
                                 beta1_levels = NULL) {
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  if (is.null(beta1_levels)) {
    beta1_levels <- config$beta1_levels
  }
  ss <- study_structure[study_structure$scenario_id == scenario$scenario_id, , drop = FALSE]
  if (!nrow(ss)) {
    stop("No study structure found for scenario_id ", scenario$scenario_id, call. = FALSE)
  }
  components <- simulate_base_components_s4(ss, scenario, rep_id, config)
  datasets <- lapply(beta1_levels, function(b) {
    assemble_dataset_s4(components, scenario, b, rep_id, config)
  })
  names(datasets) <- sprintf("beta1_%s", format(beta1_levels, trim = TRUE, scientific = FALSE))
  list(
    components = components,
    datasets = datasets,
    beta1_levels = as.numeric(beta1_levels)
  )
}

replication_diagnostics_s4 <- function(scenario, study_structure, rep_id, datasets, components = NULL) {
  ss <- study_structure[study_structure$scenario_id == scenario$scenario_id, , drop = FALSE]
  d0 <- datasets[[1]]
  data.frame(
    scenario_id = as.integer(scenario$scenario_id),
    scenario_set = as.character(scenario$scenario_set),
    outcome_scenario_id = if (!is.null(components)) {
      as.integer(components$outcome_scenario_id)
    } else if ("outcome_scenario_id" %in% names(scenario)) {
      as.integer(scenario$outcome_scenario_id)
    } else {
      as.integer(scenario$scenario_id)
    },
    rep_id = as.integer(rep_id),
    J = nrow(ss),
    N_total = as.integer(sum(ss$n_j)),
    K_total = as.integer(sum(ss$k_j)),
    rho_true = as.numeric(scenario$rho_true),
    tau = as.numeric(scenario$tau),
    omega = as.numeric(scenario$omega),
    n_effects = nrow(d0),
    y0_mean = mean(d0$y0),
    y0_sd = stats::sd(d0$y0),
    stringsAsFactors = FALSE
  )
}

# Model fitting (MLMA + CR2 only)

condition_adapter_s4 <- function(scenario, beta1_true) {
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  data.frame(
    condition_id = as.integer(scenario$scenario_id),
    condition_set = as.character(scenario$scenario_set),
    J = as.integer(scenario$J),
    k_imbalance = "study4",
    moderator_split = "study4",
    n_imbalance = "study4",
    beta1_true = as.numeric(beta1_true),
    allocation_pattern = if ("alignment_direction" %in% names(scenario)) {
      as.character(scenario$alignment_direction)
    } else {
      "independent"
    },
    stringsAsFactors = FALSE
  )
}

fitting_config_s4 <- function(scenario, config_s4) {
  cfg <- make_sim2_config("smoke")
  cfg$beta0 <- config_s4$beta0
  cfg$alpha <- config_s4$alpha
  cfg$ci_level <- config_s4$ci_level
  cfg$rho_assumed_aggregation <- 0.50
  cfg$rho_assumed_che <- as.numeric(scenario$rho_true)
  cfg$rho_true <- as.numeric(scenario$rho_true)
  cfg$tau <- as.numeric(scenario$tau)
  cfg$omega <- as.numeric(scenario$omega)
  cfg$warning_separator <- config_s4$warning_separator
  cfg$boundary_tolerance <- config_s4$boundary_tolerance
  cfg
}

remap_method_result_s4 <- function(row, scenario, beta1_true) {
  data.frame(
    scenario_id = as.integer(scenario$scenario_id),
    scenario_set = as.character(scenario$scenario_set),
    targeted_class = if ("targeted_class" %in% names(scenario)) {
      as.character(scenario$targeted_class)
    } else {
      NA_character_
    },
    pair_id = if ("pair_id" %in% names(scenario)) as.integer(scenario$pair_id) else NA_integer_,
    pair_arm = if ("pair_arm" %in% names(scenario)) as.character(scenario$pair_arm) else NA_character_,
    family_id = if ("family_id" %in% names(scenario)) as.integer(scenario$family_id) else NA_integer_,
    tau_family_level = if ("tau_family_level" %in% names(scenario)) {
      as.numeric(scenario$tau_family_level)
    } else {
      NA_real_
    },
    rep_id = as.integer(row$rep_id),
    beta1_true = as.numeric(beta1_true),
    method = as.character(row$method),
    estimate = as.numeric(row$estimate),
    se = as.numeric(row$se),
    ci_lower = as.numeric(row$ci_lower),
    ci_upper = as.numeric(row$ci_upper),
    p_value = as.numeric(row$p_value),
    df = as.numeric(row$df),
    test_statistic = if (is.finite(row$se) && isTRUE(row$se > 0) && is.finite(row$estimate)) {
      as.numeric(row$estimate / row$se)
    } else {
      NA_real_
    },
    reject = if (is.finite(row$p_value)) as.integer(row$p_value < 0.05) else NA_integer_,
    converged = as.logical(row$converged),
    status = as.character(row$status),
    error = as.character(row$error),
    warning_count = as.integer(row$warning_count),
    warnings = as.character(row$warnings),
    fit_runtime_sec = as.numeric(row$fit_runtime_sec),
    variance_study = as.numeric(row$variance_study),
    variance_effect = as.numeric(row$variance_effect),
    variance_residual = as.numeric(row$variance_residual),
    boundary_study = as.logical(row$boundary_study),
    boundary_effect = as.logical(row$boundary_effect),
    boundary_residual = as.logical(row$boundary_residual),
    stringsAsFactors = FALSE
  )
}

fit_mlma_s4 <- function(data) {
  fit_mlma_base_s2(data)
}

extract_cr2_s4 <- function(fit_capture, data, condition, rep_id, config) {
  extract_cr2_s2(
    fit_capture = fit_capture,
    data = data,
    method_label = "mlma_cr2",
    condition = condition,
    rep_id = rep_id,
    config = config
  )
}

empty_method_pair_s4 <- function(scenario, rep_id, beta1_true, error, config_s4) {
  cond <- condition_adapter_s4(scenario, beta1_true)
  cfg2 <- fitting_config_s4(scenario, config_s4)
  mlma <- finalize_method_result_s2(
    empty_method_result_s2(
      "mlma", cond, rep_id, error = error, status = "error", config = cfg2
    ),
    cfg2
  )
  cr2 <- finalize_method_result_s2(
    empty_method_result_s2(
      "mlma_cr2", cond, rep_id, error = error, status = "error", config = cfg2
    ),
    cfg2
  )
  rbind(
    remap_method_result_s4(mlma, scenario, beta1_true),
    remap_method_result_s4(cr2, scenario, beta1_true)
  )
}

# Fit MLMA once and extract model-based MLMA plus MLMA+CR2.

# Always returns exactly two rows. Failures are captured per method.
fit_primary_method_s4 <- function(data,
                                  scenario,
                                  rep_id,
                                  beta1_true,
                                  config_s4 = make_sim4_config("smoke")) {
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  cond <- condition_adapter_s4(scenario, beta1_true)
  cfg2 <- fitting_config_s4(scenario, config_s4)

  mlma_cap <- tryCatch(
    fit_mlma_s4(data),
    error = function(e) {
      list(
        value = NULL,
        error = conditionMessage(e),
        warnings = character(0),
        elapsed_sec = NA_real_
      )
    }
  )

  row_mlma <- tryCatch(
    extract_mlma_model_s2(mlma_cap, cond, as.integer(rep_id), cfg2),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "mlma", cond, as.integer(rep_id),
          error = conditionMessage(e),
          warnings = mlma_cap$warnings,
          elapsed_sec = mlma_cap$elapsed_sec,
          status = "error",
          config = cfg2
        ),
        cfg2
      )
    }
  )

  row_cr2 <- tryCatch(
    extract_cr2_s4(mlma_cap, data, cond, as.integer(rep_id), cfg2),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "mlma_cr2", cond, as.integer(rep_id),
          error = conditionMessage(e),
          warnings = mlma_cap$warnings,
          elapsed_sec = mlma_cap$elapsed_sec,
          status = "error",
          config = cfg2
        ),
        cfg2
      )
    }
  )

  out <- rbind(
    remap_method_result_s4(row_mlma, scenario, beta1_true),
    remap_method_result_s4(row_cr2, scenario, beta1_true)
  )
  rownames(out) <- NULL
  if (nrow(out) != 2L) {
    stop("fit_primary_method_s4() must return exactly two rows.", call. = FALSE)
  }
  if (!identical(as.character(out$method), c("mlma", "mlma_cr2"))) {
    stop("fit_primary_method_s4() returned unexpected method labels.", call. = FALSE)
  }
  out
}

# Simulation Study 4 — replication runner and confirmatory summaries.

run_one_replication_s4 <- function(scenario,
                                   study_structure,
                                   rep_id,
                                   config = make_sim4_config("smoke"),
                                   return_data = FALSE) {
  if (is.data.frame(scenario)) {
    if (nrow(scenario) != 1L) {
      stop("scenario must be a one-row data frame or list.", call. = FALSE)
    }
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  rep_id <- as.integer(rep_id)
  sim <- tryCatch(
    simulate_beta_crn_s4(
      scenario = scenario,
      study_structure = study_structure,
      rep_id = rep_id,
      config = config
    ),
    error = function(e) e
  )

  beta1_levels <- config$beta1_levels
  if (inherits(sim, "error")) {
    err <- conditionMessage(sim)
    method_results <- do.call(
      rbind,
      lapply(beta1_levels, function(b) {
        empty_method_pair_s4(scenario, rep_id, b, err, config)
      })
    )
    diagnostics <- data.frame(
      scenario_id = as.integer(scenario$scenario_id),
      scenario_set = as.character(scenario$scenario_set),
      rep_id = rep_id,
      dgm_error = err,
      stringsAsFactors = FALSE
    )
    out <- list(method_results = method_results, replication_diagnostics = diagnostics)
    if (isTRUE(return_data)) out$datasets <- NULL
    return(out)
  }

  method_list <- vector("list", length(sim$datasets))
  for (i in seq_along(sim$datasets)) {
    method_list[[i]] <- fit_primary_method_s4(
      data = sim$datasets[[i]],
      scenario = scenario,
      rep_id = rep_id,
      beta1_true = sim$beta1_levels[[i]],
      config_s4 = config
    )
  }
  method_results <- do.call(rbind, method_list)
  rownames(method_results) <- NULL
  diagnostics <- replication_diagnostics_s4(
    scenario, study_structure, rep_id, sim$datasets, sim$components
  )
  out <- list(
    method_results = method_results,
    replication_diagnostics = diagnostics
  )
  if (isTRUE(return_data)) {
    out$datasets <- sim$datasets
    out$components <- sim$components
  }
  out
}

summarize_condition_performance_s4 <- function(results, config) {
  if (!nrow(results)) {
    stop("results is empty.", call. = FALSE)
  }
  group_cols <- c("scenario_id", "scenario_set", "beta1_true", "method")
  missing <- setdiff(group_cols, names(results))
  if (length(missing)) {
    stop("results missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  split_keys <- interaction(results[, group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(seq_len(nrow(results)), split_keys), function(idx) {
    d <- results[idx, , drop = FALSE]
    design <- d[1, group_cols, drop = FALSE]
    extra_cols <- intersect(
      c("targeted_class", "pair_id", "pair_arm", "family_id", "tau_family_level", "J"),
      names(d)
    )
    extra <- d[1, extra_cols, drop = FALSE]
    beta1_true <- as.numeric(design$beta1_true)
    n_attempted <- nrow(d)
    n_success <- sum(d$status == "success", na.rm = TRUE)
    n_error <- sum(d$status %in% c("error", "upstream_fit_failure"), na.rm = TRUE)
    n_invalid <- sum(d$status == "invalid_result", na.rm = TRUE)
    d_ok <- d[d$status == "success", , drop = FALSE]
    if (n_success > 0L) {
      mean_estimate <- mean(d_ok$estimate)
      bias <- mean_estimate - beta1_true
      empirical_var <- stats::var(d_ok$estimate)
      empirical_sd <- stats::sd(d_ok$estimate)
      mean_se <- mean(d_ok$se)
      se_ratio <- mean_se / empirical_sd
      rmse <- sqrt(mean((d_ok$estimate - beta1_true)^2))
      covered <- (d_ok$ci_lower <= beta1_true) & (d_ok$ci_upper >= beta1_true)
      coverage <- mean(covered)
      mean_ci_width <- mean(d_ok$ci_upper - d_ok$ci_lower, na.rm = TRUE)
      reject_rate <- mean(d_ok$p_value < config$alpha)
      mcse_reject <- sqrt(reject_rate * (1 - reject_rate) / n_success)
      dfs <- d_ok$df
      if (all(is.na(dfs))) {
        mean_df <- median_df <- sd_df <- min_df <- max_df <- NA_real_
      } else {
        mean_df <- mean(dfs, na.rm = TRUE)
        median_df <- stats::median(dfs, na.rm = TRUE)
        sd_df <- stats::sd(dfs, na.rm = TRUE)
        min_df <- min(dfs, na.rm = TRUE)
        max_df <- max(dfs, na.rm = TRUE)
      }
      boundary_study_rate <- mean(d_ok$boundary_study %in% TRUE, na.rm = TRUE)
      boundary_effect_rate <- mean(d_ok$boundary_effect %in% TRUE, na.rm = TRUE)
      mean_fit_runtime_sec <- mean(d_ok$fit_runtime_sec, na.rm = TRUE)
    } else {
      mean_estimate <- bias <- empirical_var <- empirical_sd <- mean_se <- se_ratio <-
        rmse <- coverage <- mean_ci_width <- reject_rate <- mcse_reject <- mean_df <-
        median_df <- sd_df <- min_df <- max_df <- boundary_study_rate <-
        boundary_effect_rate <- mean_fit_runtime_sec <- NA_real_
    }
    type1_error <- if (isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    power <- if (!isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    elevated_failure <- isTRUE(((n_attempted - n_success) / n_attempted) > config$failure_flag)
    cbind(
      design,
      extra,
      data.frame(
        n_attempted = n_attempted,
        n_success = n_success,
        n_error = n_error,
        n_invalid = n_invalid,
        error_rate = (n_attempted - n_success) / n_attempted,
        elevated_failure_flag = elevated_failure,
        mean_estimate = mean_estimate,
        bias = bias,
        empirical_variance = empirical_var,
        empirical_sd = empirical_sd,
        mean_se = mean_se,
        se_ratio = se_ratio,
        rmse = rmse,
        coverage = coverage,
        mean_ci_width = mean_ci_width,
        reject_rate = reject_rate,
        mcse_reject_rate = mcse_reject,
        type1_error = type1_error,
        power = power,
        mean_df = mean_df,
        median_df = median_df,
        sd_df = sd_df,
        min_df = min_df,
        max_df = max_df,
        boundary_study_rate = boundary_study_rate,
        boundary_effect_rate = boundary_effect_rate,
        mean_fit_runtime_sec = mean_fit_runtime_sec,
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(out$scenario_id, out$beta1_true, out$method), , drop = FALSE]
}

ols_calibration_s4 <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (sum(ok) < 3L) {
    return(list(intercept = NA_real_, slope = NA_real_, r2 = NA_real_, pearson = NA_real_))
  }
  fit <- stats::lm(observed[ok] ~ predicted[ok])
  sm <- summary(fit)
  pearson <- if (stats::sd(observed[ok]) > 0 && stats::sd(predicted[ok]) > 0) {
    as.numeric(stats::cor(observed[ok], predicted[ok], method = "pearson"))
  } else {
    NA_real_
  }
  list(
    intercept = unname(stats::coef(fit)[[1]]),
    slope = unname(stats::coef(fit)[[2]]),
    r2 = sm$r.squared,
    pearson = pearson
  )
}

prediction_error_metrics_s4 <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) {
    return(data.frame(
      n_cells = 0L,
      rmse = NA_real_,
      mae = NA_real_,
      mean_signed_error = NA_real_,
      abs_mean_signed_error = NA_real_,
      spearman = NA_real_,
      pearson = NA_real_,
      calibration_intercept = NA_real_,
      calibration_slope = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  e <- predicted[ok] - observed[ok]
  cal <- ols_calibration_s4(observed, predicted)
  spearman <- if (sum(ok) >= 3L && stats::sd(observed[ok]) > 0 && stats::sd(predicted[ok]) > 0) {
    as.numeric(stats::cor(observed[ok], predicted[ok], method = "spearman"))
  } else {
    NA_real_
  }
  data.frame(
    n_cells = as.integer(sum(ok)),
    rmse = sqrt(mean(e^2)),
    mae = mean(abs(e)),
    mean_signed_error = mean(e),
    abs_mean_signed_error = abs(mean(e)),
    spearman = spearman,
    pearson = cal$pearson,
    calibration_intercept = cal$intercept,
    calibration_slope = cal$slope,
    stringsAsFactors = FALSE
  )
}

merge_primary_diagnostics_s4 <- function(method_summary, scenario_grid) {
  keep <- unique(scenario_grid[, intersect(
    c(
      "scenario_id", "scenario_set", "targeted_class", "pair_id", "pair_arm",
      "family_id", "tau_family_level", "J", "J0", "J1", "N_total", "K_total",
      "mean_n_target", "alignment_direction", "alignment_strength",
      "nk_relationship", "rho_true", "tau", "omega",
      "I_C_S_star", "I_C_H_star", "R_H", "T_C", "B_C",
      "J0_star", "J1_star", "J_min_star", "B_J", "nu_C_star",
      "predicted_variance_S", "predicted_variance_H", "predicted_variance_oracle",
      "predicted_se_S", "predicted_se_H", "predicted_se_oracle",
      "predicted_power_S_0.05", "predicted_power_S_0.10", "predicted_power_S_0.20",
      "predicted_power_H_0.05", "predicted_power_H_0.10", "predicted_power_H_0.20",
      "predicted_power_oracle_0.05", "predicted_power_oracle_0.10",
      "predicted_power_oracle_0.20",
      "F_I", "delta_P_H_0.10", "J_prec", "J_k", "minority_study_count",
      "primary_status"
    ),
    names(scenario_grid)
  ), drop = FALSE])
  merge(method_summary, keep, by = "scenario_id", all.x = TRUE, suffixes = c("", "_design"))
}

predicted_power_for_row_s4 <- function(row, regime) {
  b <- as.numeric(row$beta1_true)
  pick <- function(p05, p10, p20) {
    if (isTRUE(all.equal(b, 0.05))) return(as.numeric(p05))
    if (isTRUE(all.equal(b, 0.10))) return(as.numeric(p10))
    if (isTRUE(all.equal(b, 0.20))) return(as.numeric(p20))
    NA_real_
  }
  if (identical(regime, "sampling_only")) {
    return(pick(row$predicted_power_S_0.05, row$predicted_power_S_0.10, row$predicted_power_S_0.20))
  }
  if (identical(regime, "near_correct")) {
    return(pick(row$predicted_power_H_0.05, row$predicted_power_H_0.10, row$predicted_power_H_0.20))
  }
  if (identical(regime, "oracle")) {
    return(pick(
      row$predicted_power_oracle_0.05, row$predicted_power_oracle_0.10,
      row$predicted_power_oracle_0.20
    ))
  }
  NA_real_
}

summarize_power_calibration_s4 <- function(method_summary,
                                           scenario_grid,
                                           method = "mlma_cr2") {
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  ms <- ms[ms$method == method, , drop = FALSE]
  ms <- ms[!isTRUE(all.equal(ms$beta1_true, 0)) & ms$beta1_true > 0, , drop = FALSE]
  if (!nrow(ms)) {
    return(data.frame())
  }
  nulls <- method_summary[
    method_summary$method == method & abs(method_summary$beta1_true) < 1e-12,
    c("scenario_id", "type1_error"),
    drop = FALSE
  ]
  names(nulls)[2] <- "type1_at_null"
  ms <- merge(ms, nulls, by = "scenario_id", all.x = TRUE)
  ms$in_type1_band <- is.finite(ms$type1_at_null) &
    ms$type1_at_null >= 0.025 & ms$type1_at_null <= 0.075

  strata <- list(
    list(name = "overall_nonnull", idx = rep(TRUE, nrow(ms))),
    list(name = "beta1_0.05", idx = abs(ms$beta1_true - 0.05) < 1e-12),
    list(name = "beta1_0.10", idx = abs(ms$beta1_true - 0.10) < 1e-12),
    list(name = "beta1_0.20", idx = abs(ms$beta1_true - 0.20) < 1e-12),
    list(name = "broad_primary", idx = ms$scenario_set == "broad" | ms$scenario_set_design == "broad"),
    list(name = "targeted_secondary", idx = ms$scenario_set == "targeted" | ms$scenario_set_design == "targeted"),
    list(name = "type1_robustness_band", idx = ms$in_type1_band %in% TRUE)
  )
  # scenario_set may have been suffixed
  if (!"scenario_set" %in% names(ms) && "scenario_set_design" %in% names(ms)) {
    strata[[5]]$idx <- ms$scenario_set_design == "broad"
    strata[[6]]$idx <- ms$scenario_set_design == "targeted"
  }

  regimes <- c("sampling_only", "near_correct", "oracle")
  out <- list()
  k <- 0L
  for (rg in regimes) {
    pred <- vapply(seq_len(nrow(ms)), function(i) predicted_power_for_row_s4(ms[i, ], rg), numeric(1))
    for (st in strata) {
      k <- k + 1L
      idx <- st$idx & is.finite(ms$power)
      met <- prediction_error_metrics_s4(ms$power[idx], pred[idx])
      out[[k]] <- cbind(
        data.frame(
          regime = rg,
          stratum = st$name,
          method = method,
          stringsAsFactors = FALSE
        ),
        met
      )
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

summarize_variance_calibration_s4 <- function(method_summary,
                                              scenario_grid,
                                              method = "mlma_cr2") {
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  ms <- ms[ms$method == method & abs(ms$beta1_true) < 1e-12, , drop = FALSE]
  if (!nrow(ms)) return(data.frame())
  regimes <- data.frame(
    regime = c("sampling_only", "near_correct", "oracle"),
    pred_col = c("predicted_variance_S", "predicted_variance_H", "predicted_variance_oracle"),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(regimes))
  for (i in seq_len(nrow(regimes))) {
    pc <- regimes$pred_col[[i]]
    if (!pc %in% names(ms)) {
      rows[[i]] <- data.frame(regime = regimes$regime[[i]], n_cells = 0L, stringsAsFactors = FALSE)
      next
    }
    met <- prediction_error_metrics_s4(ms$empirical_variance, ms[[pc]])
    ratio <- ms[[pc]] / ms$empirical_variance
    rows[[i]] <- cbind(
      data.frame(
        regime = regimes$regime[[i]],
        method = method,
        mean_predicted_empirical_ratio = mean(ratio[is.finite(ratio)]),
        stringsAsFactors = FALSE
      ),
      met
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_retention_validation_s4 <- function(method_summary,
                                              scenario_grid,
                                              method = "mlma_cr2") {
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  cr2_0 <- ms[ms$method == method & abs(ms$beta1_true) < 1e-12, , drop = FALSE]
  cr2_0 <- cr2_0[cr2_0$scenario_set == "broad" | cr2_0$scenario_set_design == "broad", , drop = FALSE]
  vif <- cr2_0$empirical_variance / cr2_0$predicted_variance_S
  log_vif <- log(vif)
  ok <- is.finite(cr2_0$R_H) & is.finite(log_vif)
  spear <- if (sum(ok) >= 3L) {
    as.numeric(stats::cor(cr2_0$R_H[ok], log_vif[ok], method = "spearman"))
  } else {
    NA_real_
  }
  pear <- if (sum(ok) >= 3L && stats::sd(cr2_0$R_H[ok]) > 0 && stats::sd(log_vif[ok]) > 0) {
    as.numeric(stats::cor(cr2_0$R_H[ok], log_vif[ok], method = "pearson"))
  } else {
    NA_real_
  }
  cr2_10 <- ms[ms$method == method & abs(ms$beta1_true - 0.10) < 1e-12, , drop = FALSE]
  cr2_10 <- cr2_10[cr2_10$scenario_id %in% cr2_0$scenario_id, , drop = FALSE]
  m10 <- merge(
    cr2_0[, c("scenario_id", "R_H"), drop = FALSE],
    cr2_10[, c("scenario_id", "power", "predicted_power_S_0.10"), drop = FALSE],
    by = "scenario_id"
  )
  optimism <- m10$predicted_power_S_0.10 - m10$power
  ok2 <- is.finite(m10$R_H) & is.finite(optimism)
  spear_opt <- if (sum(ok2) >= 3L) {
    as.numeric(stats::cor(m10$R_H[ok2], optimism[ok2], method = "spearman"))
  } else {
    NA_real_
  }
  data.frame(
    n_broad_null = as.integer(sum(ok)),
    spearman_R_H_logVIF_S = spear,
    pearson_R_H_logVIF_S = pear,
    expected_direction = "negative",
    n_broad_beta10 = as.integer(sum(ok2)),
    spearman_R_H_sampling_power_optimism = spear_opt,
    stringsAsFactors = FALSE
  )
}

summarize_support_validation_s4 <- function(method_summary,
                                            scenario_grid,
                                            pair_diagnostics = NULL,
                                            method = "mlma_cr2") {
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  cr2_0 <- ms[ms$method == method & abs(ms$beta1_true) < 1e-12, , drop = FALSE]
  assoc <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    data.frame(
      n = as.integer(sum(ok)),
      spearman = if (sum(ok) >= 3L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
        as.numeric(stats::cor(x[ok], y[ok], method = "spearman"))
      } else {
        NA_real_
      },
      pearson = if (sum(ok) >= 3L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
        as.numeric(stats::cor(x[ok], y[ok], method = "pearson"))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
  rows <- rbind(
    cbind(data.frame(diagnostic = "J_min_star", target = "mean_df", stringsAsFactors = FALSE),
          assoc(cr2_0$J_min_star, cr2_0$mean_df)),
    cbind(data.frame(diagnostic = "B_J", target = "mean_df", stringsAsFactors = FALSE),
          assoc(cr2_0$B_J, cr2_0$mean_df)),
    cbind(data.frame(diagnostic = "nu_C_star", target = "mean_df", stringsAsFactors = FALSE),
          assoc(cr2_0$nu_C_star, cr2_0$mean_df)),
    cbind(data.frame(diagnostic = "J_min_star", target = "mean_ci_width", stringsAsFactors = FALSE),
          assoc(cr2_0$J_min_star, cr2_0$mean_ci_width)),
    cbind(data.frame(diagnostic = "B_J", target = "mean_ci_width", stringsAsFactors = FALSE),
          assoc(cr2_0$B_J, cr2_0$mean_ci_width)),
    cbind(data.frame(diagnostic = "nu_C_star", target = "mean_ci_width", stringsAsFactors = FALSE),
          assoc(cr2_0$nu_C_star, cr2_0$mean_ci_width))
  )
  rownames(rows) <- NULL
  rows
}

summarize_targeted_pairs_s4 <- function(method_summary,
                                        scenario_grid,
                                        pair_diagnostics,
                                        method = "mlma_cr2") {
  if (is.null(pair_diagnostics) || !nrow(pair_diagnostics)) {
    return(NULL)
  }
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  ms <- ms[ms$method == method, , drop = FALSE]
  pick <- function(sid, beta) {
    ms[ms$scenario_id == sid & abs(ms$beta1_true - beta) < 1e-12, , drop = FALSE]
  }
  rows <- vector("list", nrow(pair_diagnostics))
  for (i in seq_len(nrow(pair_diagnostics))) {
    pd <- pair_diagnostics[i, ]
    a0 <- pick(pd$scenario_id_a, 0)
    b0 <- pick(pd$scenario_id_b, 0)
    a10 <- pick(pd$scenario_id_a, 0.10)
    b10 <- pick(pd$scenario_id_b, 0.10)
    a05 <- pick(pd$scenario_id_a, 0.05)
    b05 <- pick(pd$scenario_id_b, 0.05)
    a20 <- pick(pd$scenario_id_a, 0.20)
    b20 <- pick(pd$scenario_id_b, 0.20)
    cl <- as.character(pd$pair_class)
    # Identify the "higher diagnostic" arm per frozen directional hypothesis.
    if (identical(cl, "A")) {
      expected <- "lower_R_H_has_larger_VIF_S"
      higher_is_a <- isTRUE(pd$R_H_a < pd$R_H_b)
    } else if (identical(cl, "B")) {
      expected <- "higher_B_C_has_lower_empirical_variance"
      higher_is_a <- isTRUE(pd$B_C_a > pd$B_C_b)
    } else if (identical(cl, "C")) {
      expected <- "greater_J_min_and_B_J_has_greater_mean_df"
      higher_is_a <- isTRUE(pd$J_min_star_a > pd$J_min_star_b)
    } else {
      expected <- "higher_I_C_H_has_lower_empirical_variance"
      higher_is_a <- isTRUE(pd$I_C_H_star_a > pd$I_C_H_star_b)
    }
    var_a <- if (nrow(a0)) a0$empirical_variance[[1]] else NA_real_
    var_b <- if (nrow(b0)) b0$empirical_variance[[1]] else NA_real_
    df_a <- if (nrow(a0)) a0$mean_df[[1]] else NA_real_
    df_b <- if (nrow(b0)) b0$mean_df[[1]] else NA_real_
    p10_a <- if (nrow(a10)) a10$power[[1]] else NA_real_
    p10_b <- if (nrow(b10)) b10$power[[1]] else NA_real_
    vif_a <- if (nrow(a0)) a0$empirical_variance[[1]] / a0$predicted_variance_S[[1]] else NA_real_
    vif_b <- if (nrow(b0)) b0$empirical_variance[[1]] / b0$predicted_variance_S[[1]] else NA_real_
    direction_ok <- NA
    if (identical(cl, "A")) {
      if (isTRUE(pd$R_H_a < pd$R_H_b)) {
        direction_ok <- isTRUE(vif_a > vif_b)
      } else {
        direction_ok <- isTRUE(vif_b > vif_a)
      }
    } else if (identical(cl, "B")) {
      if (isTRUE(pd$B_C_a > pd$B_C_b)) {
        direction_ok <- isTRUE(var_a < var_b)
      } else {
        direction_ok <- isTRUE(var_b < var_a)
      }
    } else if (identical(cl, "C")) {
      # Arm with greater J_min* should have greater mean df.
      if (isTRUE(pd$J_min_star_a > pd$J_min_star_b)) {
        direction_ok <- isTRUE(df_a > df_b)
      } else {
        direction_ok <- isTRUE(df_b > df_a)
      }
    } else if (identical(cl, "D")) {
      if (isTRUE(pd$I_C_H_star_a > pd$I_C_H_star_b)) {
        direction_ok <- isTRUE(var_a < var_b)
      } else {
        direction_ok <- isTRUE(var_b < var_a)
      }
    }
    rows[[i]] <- data.frame(
      pair_class = cl,
      pair_id = as.integer(pd$pair_id),
      scenario_id_a = as.integer(pd$scenario_id_a),
      scenario_id_b = as.integer(pd$scenario_id_b),
      expected_direction = expected,
      empirical_variance_a = var_a,
      empirical_variance_b = var_b,
      empirical_variance_diff = var_a - var_b,
      mean_df_a = df_a,
      mean_df_b = df_b,
      mean_df_diff = df_a - df_b,
      power_0.10_a = p10_a,
      power_0.10_b = p10_b,
      power_0.10_diff = p10_a - p10_b,
      power_0.05_diff = if (nrow(a05) && nrow(b05)) a05$power[[1]] - b05$power[[1]] else NA_real_,
      power_0.20_diff = if (nrow(a20) && nrow(b20)) a20$power[[1]] - b20$power[[1]] else NA_real_,
      VIF_S_a = vif_a,
      VIF_S_b = vif_b,
      direction_supported = as.logical(direction_ok),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_fragility_families_s4 <- function(method_summary,
                                            scenario_grid,
                                            family_diagnostics,
                                            method = "mlma_cr2") {
  if (is.null(family_diagnostics) || !nrow(family_diagnostics)) {
    return(NULL)
  }
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  ms <- ms[ms$method == method, , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(family_diagnostics))) {
    fd <- family_diagnostics[i, ]
    fam <- ms[ms$family_id == fd$family_id | ms$family_id_design == fd$family_id, , drop = FALSE]
    if (!nrow(fam) && "family_id" %in% names(ms)) {
      fam <- ms[ms$family_id == fd$family_id, , drop = FALSE]
    }
    p10 <- fam[abs(fam$beta1_true - 0.10) < 1e-12, , drop = FALSE]
    getp <- function(tau) {
      r <- p10[abs(p10$tau - tau) < 1e-12 | abs(p10$tau_family_level - tau) < 1e-12, ]
      if (!nrow(r)) return(list(power = NA_real_, var = NA_real_, pred = NA_real_, ich = NA_real_))
      list(
        power = r$power[[1]],
        var = r$empirical_variance[[1]],
        pred = r$predicted_power_H_0.10[[1]],
        ich = r$I_C_H_star[[1]]
      )
    }
    t05 <- getp(0.05)
    t12 <- getp(0.12)
    t20 <- getp(0.20)
    t30 <- getp(0.30)
    delta_emp <- t05$power - t30$power
    rows[[i]] <- data.frame(
      family_id = as.integer(fd$family_id),
      profile = if ("profile" %in% names(fd)) as.integer(fd$profile) else NA_integer_,
      profile_label = if ("profile_label" %in% names(fd)) as.character(fd$profile_label) else NA_character_,
      F_I = as.numeric(fd$F_I),
      delta_P_H_0.10 = as.numeric(fd$delta_P_H_0.10),
      power_0.10_tau_0.05 = t05$power,
      power_0.10_tau_0.12 = t12$power,
      power_0.10_tau_0.20 = t20$power,
      power_0.10_tau_0.30 = t30$power,
      var_tau_0.05 = t05$var,
      var_tau_0.12 = t12$var,
      var_tau_0.20 = t20$var,
      var_tau_0.30 = t30$var,
      pred_power_tau_0.05 = t05$pred,
      pred_power_tau_0.12 = t12$pred,
      pred_power_tau_0.20 = t20$pred,
      pred_power_tau_0.30 = t30$pred,
      I_C_H_star_tau_0.05 = t05$ich,
      I_C_H_star_tau_0.30 = t30$ich,
      delta_P_emp_0.10 = delta_emp,
      signed_sensitivity_error = as.numeric(fd$delta_P_H_0.10) - delta_emp,
      abs_sensitivity_error = abs(as.numeric(fd$delta_P_H_0.10) - delta_emp),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_comparator_performance_s4 <- function(method_summary,
                                                scenario_grid,
                                                method = "mlma_cr2") {
  ms <- merge_primary_diagnostics_s4(method_summary, scenario_grid)
  cr2_0 <- ms[ms$method == method & abs(ms$beta1_true) < 1e-12, , drop = FALSE]
  cr2_10 <- ms[ms$method == method & abs(ms$beta1_true - 0.10) < 1e-12, , drop = FALSE]
  is_broad <- function(d) {
    if ("scenario_set" %in% names(d)) {
      return(d$scenario_set == "broad")
    }
    d$scenario_set_design == "broad"
  }
  cr2_0 <- cr2_0[is_broad(cr2_0), , drop = FALSE]
  cr2_10 <- cr2_10[cr2_10$scenario_id %in% cr2_0$scenario_id, , drop = FALSE]
  m <- merge(
    cr2_0[, intersect(c(
      "scenario_id", "empirical_variance", "N_total", "J", "K_total", "mean_n_target",
      "minority_study_count", "J_prec", "J_k", "I_C_S_star", "I_C_H_star", "R_H",
      "B_C", "J_min_star", "B_J"
    ), names(cr2_0)), drop = FALSE],
    cr2_10[, c("scenario_id", "power"), drop = FALSE],
    by = "scenario_id"
  )
  comps <- c(
    "N_total", "J", "K_total", "mean_n_target", "minority_study_count",
    "J_prec", "J_k", "I_C_S_star", "I_C_H_star", "R_H", "B_C", "J_min_star", "B_J"
  )
  comps <- intersect(comps, names(m))
  assoc <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
      return(NA_real_)
    }
    as.numeric(stats::cor(x[ok], y[ok], method = "spearman"))
  }
  rows <- lapply(comps, function(nm) {
    data.frame(
      comparator = nm,
      spearman_vs_empirical_variance = assoc(m[[nm]], m$empirical_variance),
      spearman_vs_power_0.10 = assoc(m[[nm]], m$power),
      n = nrow(m),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_type1_robustness_s4 <- function(method_summary, method = "mlma_cr2") {
  d <- method_summary[method_summary$method == method & abs(method_summary$beta1_true) < 1e-12, ]
  in_band <- is.finite(d$type1_error) & d$type1_error >= 0.025 & d$type1_error <= 0.075
  data.frame(
    n_scenarios = nrow(d),
    n_in_band = sum(in_band, na.rm = TRUE),
    n_below_band = sum(is.finite(d$type1_error) & d$type1_error < 0.025),
    n_above_band = sum(is.finite(d$type1_error) & d$type1_error > 0.075),
    mean_type1 = mean(d$type1_error, na.rm = TRUE),
    min_type1 = min(d$type1_error, na.rm = TRUE),
    max_type1 = max(d$type1_error, na.rm = TRUE),
    band_low = 0.025,
    band_high = 0.075,
    stringsAsFactors = FALSE
  )
}

summarize_heterogeneity_sensitivity_s4 <- function(metric_library) {
  if (is.null(metric_library) || !nrow(metric_library)) {
    return(data.frame())
  }
  g <- metric_library
  fam <- unique(g$regime[g$regime %in% c("practical_tau_omega", "rho_sensitivity")])
  rows <- lapply(fam, function(rg) {
    d <- g[g$regime == rg, , drop = FALSE]
    data.frame(
      regime = rg,
      n_rows = nrow(d),
      n_scenarios = length(unique(d$scenario_id)),
      n_variants = length(unique(d$variant_id)),
      mean_I_C_star = mean(d$I_C_star, na.rm = TRUE),
      mean_predicted_power_0.10 = mean(d$predicted_power_0.10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Simulation Study 4 — chunking, checkpointing, parallel execution
# Chunks are scenario-level (not beta1-level) so CRN components are generated once.

source_sim4_functions <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  source(file.path(root, "R", "study1_functions.R"), local = FALSE)
  source(file.path(root, "R", "study2_design.R"), local = FALSE)
  source(file.path(root, "R", "study2_functions.R"), local = FALSE)
  source(file.path(root, "R", "estimators.R"), local = FALSE)
  source(file.path(root, "R", "contrast_information.R"), local = FALSE)
  source(file.path(root, "R", "study3_functions.R"), local = FALSE)
  source(file.path(root, "R", "study4_functions.R"), local = FALSE)
  invisible(TRUE)
}

sim4_design_hash <- function(root = NULL) {
  dir <- sim4_fixture_dir(root)
  path <- file.path(dir, "design_hashes.csv")
  if (!file.exists(path)) {
    stop("No Study 4 design_hashes.csv found under ", dir, call. = FALSE)
  }
  hashes <- utils::read.csv(path, stringsAsFactors = FALSE)
  payload <- paste(hashes$file, hashes$sha256, sep = "=", collapse = "\n")
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(payload, algo = "sha256")
  } else {
    sprintf("raw-%d", abs(as.integer(sum(utf8ToInt(payload) %% 1000003L))))
  }
}

chunk_rds_path_s4 <- function(out_dir, scenario_set, scenario_id, chunk_id) {
  file.path(
    out_dir,
    "raw_chunks",
    as.character(scenario_set),
    sprintf("scenario_%03d_chunk_%03d.rds", as.integer(scenario_id), as.integer(chunk_id))
  )
}

make_chunk_grid_s4 <- function(scenario_grid, n_rep_total, chunk_size) {
  n_rep_total <- as.integer(n_rep_total)
  chunk_size <- as.integer(chunk_size)
  if (n_rep_total < 1L || chunk_size < 1L) {
    stop("n_rep_total and chunk_size must be >= 1.", call. = FALSE)
  }
  if (n_rep_total %% chunk_size != 0L) {
    stop("n_rep_total must be divisible by chunk_size.", call. = FALSE)
  }
  required <- c("scenario_id", "scenario_set", "J")
  missing <- setdiff(required, names(scenario_grid))
  if (length(missing)) {
    stop("scenario_grid missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  n_chunks <- n_rep_total / chunk_size
  rows <- vector("list", nrow(scenario_grid) * n_chunks)
  k <- 0L
  for (i in seq_len(nrow(scenario_grid))) {
    sc <- scenario_grid[i, , drop = FALSE]
    for (ch in seq_len(n_chunks)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        scenario_id = as.integer(sc$scenario_id),
        scenario_set = as.character(sc$scenario_set),
        targeted_class = if ("targeted_class" %in% names(sc)) as.character(sc$targeted_class) else NA_character_,
        pair_id = if ("pair_id" %in% names(sc)) as.integer(sc$pair_id) else NA_integer_,
        family_id = if ("family_id" %in% names(sc)) as.integer(sc$family_id) else NA_integer_,
        chunk_id = as.integer(ch),
        J = as.integer(sc$J),
        rho_true = if ("rho_true" %in% names(sc)) as.numeric(sc$rho_true) else NA_real_,
        tau = if ("tau" %in% names(sc)) as.numeric(sc$tau) else NA_real_,
        omega = if ("omega" %in% names(sc)) as.numeric(sc$omega) else NA_real_,
        rep_start = as.integer((ch - 1L) * chunk_size + 1L),
        rep_end = as.integer(ch * chunk_size),
        n_reps = as.integer(chunk_size),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (anyDuplicated(out[, c("scenario_id", "chunk_id")])) {
    stop("Chunk grid has duplicate scenario/chunk IDs.", call. = FALSE)
  }
  out
}

assert_chunk_coverage_s4 <- function(chunk_grid, n_rep_total, n_scenarios = NULL) {
  issues <- character(0)
  if (!is.null(n_scenarios) && length(unique(chunk_grid$scenario_id)) != n_scenarios) {
    issues <- c(issues, sprintf(
      "Expected %d scenarios in chunk grid, found %d.",
      n_scenarios, length(unique(chunk_grid$scenario_id))
    ))
  }
  split_ids <- split(chunk_grid, chunk_grid$scenario_id)
  for (nm in names(split_ids)) {
    ch <- split_ids[[nm]]
    ch <- ch[order(ch$rep_start), , drop = FALSE]
    reps <- unlist(mapply(seq, ch$rep_start, ch$rep_end, SIMPLIFY = FALSE))
    if (anyDuplicated(reps)) {
      issues <- c(issues, sprintf("scenario %s: overlapping replication IDs.", nm))
    }
    if (!identical(as.integer(sort(reps)), seq_len(as.integer(n_rep_total)))) {
      issues <- c(issues, sprintf("scenario %s: replication IDs do not cover 1:%s exactly once.",
                                  nm, n_rep_total))
    }
    if (any(ch$n_reps != (ch$rep_end - ch$rep_start + 1L))) {
      issues <- c(issues, sprintf("scenario %s: n_reps does not match rep interval.", nm))
    }
  }
  list(ok = length(issues) == 0L, issues = issues)
}

validate_chunk_object_s4 <- function(obj, chunk, n_beta = 4L, n_methods = 2L) {
  if (!is.list(obj) || !all(c("metadata", "method_results", "replication_diagnostics") %in% names(obj))) {
    return(list(ok = FALSE, reason = "Chunk object lacks required top-level elements."))
  }
  meta <- obj$metadata
  if (is.data.frame(meta)) meta <- as.list(meta[1, , drop = FALSE])
  if (!identical(as.integer(meta$scenario_id), as.integer(chunk$scenario_id)) ||
      !identical(as.integer(meta$chunk_id), as.integer(chunk$chunk_id))) {
    return(list(ok = FALSE, reason = "metadata identifiers do not match chunk job."))
  }
  n_reps <- as.integer(chunk$n_reps)
  mr <- obj$method_results
  dd <- obj$replication_diagnostics
  expected_method_rows <- n_reps * n_beta * n_methods
  if (!is.data.frame(mr) || nrow(mr) != expected_method_rows) {
    return(list(ok = FALSE, reason = sprintf(
      "method_results has %s rows; expected %d.",
      if (is.data.frame(mr)) nrow(mr) else "NA", expected_method_rows
    )))
  }
  if (!is.data.frame(dd) || nrow(dd) != n_reps) {
    return(list(ok = FALSE, reason = sprintf(
      "replication_diagnostics has %s rows; expected %d.",
      if (is.data.frame(dd)) nrow(dd) else "NA", n_reps
    )))
  }
  list(ok = TRUE, reason = NA_character_)
}

scenario_from_chunk_row_s4 <- function(chunk, scenario_grid) {
  sid <- as.integer(chunk$scenario_id[[1]])
  rows <- scenario_grid[scenario_grid$scenario_id == sid, , drop = FALSE]
  if (nrow(rows) != 1L) {
    stop("Could not uniquely recover scenario_id ", sid, call. = FALSE)
  }
  rows
}

run_one_chunk_s4 <- function(chunk,
                             scenario_grid,
                             study_structure,
                             root,
                             out_dir,
                             config,
                             design_hash = NULL) {
  chunk <- as.list(chunk)
  out_file <- chunk_rds_path_s4(
    out_dir, chunk$scenario_set, chunk$scenario_id, chunk$chunk_id
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (!exists("run_one_replication_s4", mode = "function")) {
    source_sim4_functions(root)
  }
  if (is.null(design_hash)) {
    design_hash <- tryCatch(sim4_design_hash(root), error = function(e) NA_character_)
  }

  if (file.exists(out_file)) {
    if (!isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      stop("Chunk file exists and resume/overwrite are FALSE: ", out_file, call. = FALSE)
    }
    if (isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      obj <- tryCatch(readRDS(out_file), error = function(e) e)
      if (!inherits(obj, "error")) {
        chk <- validate_chunk_object_s4(
          obj, chunk,
          n_beta = length(config$beta1_levels),
          n_methods = config$n_methods
        )
        if (isTRUE(chk$ok)) {
          log_row <- data.frame(
            scenario_id = as.integer(chunk$scenario_id),
            chunk_id = as.integer(chunk$chunk_id),
            scenario_set = as.character(chunk$scenario_set),
            status = "skipped",
            started_at = NA_character_,
            completed_at = NA_character_,
            elapsed_sec = NA_real_,
            n_method_rows = nrow(obj$method_results),
            n_diagnostic_rows = nrow(obj$replication_diagnostics),
            n_errors = sum(obj$method_results$status != "success", na.rm = TRUE),
            file_path = out_file,
            file_sha256 = hash_file_s4(out_file),
            error_message = NA_character_,
            stringsAsFactors = FALSE
          )
          return(list(status = "skipped", log = log_row, object = obj))
        }
      }
    }
  }

  started_at <- Sys.time()
  err_msg <- NA_character_
  scenario <- scenario_from_chunk_row_s4(chunk, scenario_grid)
  n_beta <- length(config$beta1_levels)
  n_methods <- config$n_methods

  result <- tryCatch(
    {
      reps <- seq.int(as.integer(chunk$rep_start), as.integer(chunk$rep_end))
      method_list <- vector("list", length(reps))
      diag_list <- vector("list", length(reps))
      for (i in seq_along(reps)) {
        one <- run_one_replication_s4(
          scenario = scenario,
          study_structure = study_structure,
          rep_id = reps[[i]],
          config = config,
          return_data = FALSE
        )
        method_list[[i]] <- one$method_results
        diag_list[[i]] <- one$replication_diagnostics
      }
      list(
        method_results = do.call(rbind, method_list),
        replication_diagnostics = do.call(rbind, diag_list)
      )
    },
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )

  completed_at <- Sys.time()
  elapsed_sec <- as.numeric(difftime(completed_at, started_at, units = "secs"))

  if (is.null(result)) {
    expand <- expand.grid(
      rep_id = seq.int(chunk$rep_start, chunk$rep_end),
      beta1_true = config$beta1_levels,
      method = config$methods,
      stringsAsFactors = FALSE
    )
    method_results <- data.frame(
      scenario_id = as.integer(chunk$scenario_id),
      scenario_set = as.character(chunk$scenario_set),
      rep_id = expand$rep_id,
      beta1_true = expand$beta1_true,
      method = expand$method,
      estimate = NA_real_,
      se = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      status = "error",
      error = err_msg,
      stringsAsFactors = FALSE
    )
    replication_diagnostics <- data.frame(
      scenario_id = as.integer(chunk$scenario_id),
      rep_id = seq.int(chunk$rep_start, chunk$rep_end),
      chunk_error = err_msg,
      stringsAsFactors = FALSE
    )
    chunk_status <- "failed"
  } else {
    method_results <- result$method_results
    replication_diagnostics <- result$replication_diagnostics
    rownames(method_results) <- NULL
    rownames(replication_diagnostics) <- NULL
    chunk_status <- "completed"
  }

  metadata <- list(
    scenario_id = as.integer(chunk$scenario_id),
    scenario_set = as.character(chunk$scenario_set),
    chunk_id = as.integer(chunk$chunk_id),
    rep_start = as.integer(chunk$rep_start),
    rep_end = as.integer(chunk$rep_end),
    n_reps = as.integer(chunk$n_reps),
    n_beta = n_beta,
    n_methods = n_methods,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S"),
    completed_at = format(completed_at, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = elapsed_sec,
    host = as.character(Sys.info()[["nodename"]]),
    git_commit = sim2_git_commit(root),
    design_hash = design_hash
  )
  obj <- list(
    metadata = metadata,
    method_results = method_results,
    replication_diagnostics = replication_diagnostics
  )
  save_ok <- tryCatch(
    {
      atomic_save_rds_s2(obj, out_file)
      TRUE
    },
    error = function(e) {
      err_msg <<- paste(na.omit(c(err_msg, conditionMessage(e))), collapse = "; ")
      FALSE
    }
  )
  if (!save_ok && identical(chunk_status, "completed")) {
    chunk_status <- "failed"
  }
  log_row <- data.frame(
    scenario_id = as.integer(chunk$scenario_id),
    chunk_id = as.integer(chunk$chunk_id),
    scenario_set = as.character(chunk$scenario_set),
    status = chunk_status,
    started_at = metadata$started_at,
    completed_at = metadata$completed_at,
    elapsed_sec = elapsed_sec,
    n_method_rows = nrow(method_results),
    n_diagnostic_rows = nrow(replication_diagnostics),
    n_errors = sum(method_results$status != "success", na.rm = TRUE),
    file_path = out_file,
    file_sha256 = if (file.exists(out_file)) hash_file_s4(out_file) else NA_character_,
    error_message = err_msg,
    stringsAsFactors = FALSE
  )
  list(status = chunk_status, log = log_row, object = obj)
}

run_chunks_s4 <- function(chunk_grid,
                          scenario_grid,
                          study_structure,
                          root,
                          out_dir,
                          config) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(out_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)
  design_hash <- tryCatch(sim4_design_hash(root), error = function(e) NA_character_)
  workers <- as.integer(config$workers)
  if (is.na(workers) || workers < 1L) {
    stop("config$workers must be a positive integer.", call. = FALSE)
  }
  wall_start <- Sys.time()
  if (workers == 1L) {
    results <- lapply(seq_len(nrow(chunk_grid)), function(i) {
      run_one_chunk_s4(
        chunk = chunk_grid[i, , drop = FALSE],
        scenario_grid = scenario_grid,
        study_structure = study_structure,
        root = root,
        out_dir = out_dir,
        config = config,
        design_hash = design_hash
      )
    })
  } else {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required for parallel Study 4 runs.",
           call. = FALSE)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    results <- future.apply::future_lapply(
      seq_len(nrow(chunk_grid)),
      function(i) {
        run_one_chunk_s4(
          chunk = chunk_grid[i, , drop = FALSE],
          scenario_grid = scenario_grid,
          study_structure = study_structure,
          root = root,
          out_dir = out_dir,
          config = config,
          design_hash = design_hash
        )
      },
      # Seeds are set inside simulate_base_components_s4() via withr::with_seed().
      # Do not add a parallel RNG stream.
      future.seed = NULL
    )
  }
  wall_elapsed_sec <- as.numeric(difftime(Sys.time(), wall_start, units = "secs"))
  chunk_log <- do.call(rbind, lapply(results, function(x) x$log))
  rownames(chunk_log) <- NULL
  list(chunk_log = chunk_log, wall_elapsed_sec = wall_elapsed_sec, results = results)
}

collect_sim4_chunk_results <- function(out_dir, chunk_grid) {
  method_list <- vector("list", nrow(chunk_grid))
  diag_list <- vector("list", nrow(chunk_grid))
  for (i in seq_len(nrow(chunk_grid))) {
    ch <- chunk_grid[i, , drop = FALSE]
    path <- chunk_rds_path_s4(out_dir, ch$scenario_set[[1]], ch$scenario_id[[1]], ch$chunk_id[[1]])
    obj <- readRDS(path)
    method_list[[i]] <- obj$method_results
    diag_list[[i]] <- obj$replication_diagnostics
  }
  mr <- do.call(rbind, method_list)
  dd <- do.call(rbind, diag_list)
  rownames(mr) <- NULL
  rownames(dd) <- NULL
  list(method_results = mr, replication_diagnostics = dd)
}

write_sim4_summaries <- function(out_dir,
                                 config,
                                 scenario_grid,
                                 metric_library,
                                 method_results,
                                 pair_diagnostics = NULL,
                                 family_diagnostics = NULL,
                                 scientific = TRUE) {
  method_summary <- summarize_condition_performance_s4(method_results, config)
  if (isTRUE(scientific)) {
    utils::write.csv(
      method_summary,
      file.path(out_dir, "summaries", "condition_method_summary.csv"),
      row.names = FALSE
    )
    mc <- method_summary[, intersect(
      c(
        "scenario_id", "beta1_true", "method", "n_success", "n_attempted",
        "reject_rate", "mcse_reject_rate", "empirical_sd", "mean_df",
        "elevated_failure_flag", "error_rate"
      ),
      names(method_summary)
    ), drop = FALSE]
    utils::write.csv(mc, file.path(out_dir, "summaries", "monte_carlo_uncertainty.csv"), row.names = FALSE)
  }

  if (!isTRUE(scientific)) {
    op_cols <- intersect(
      c(
        "scenario_id", "scenario_set", "beta1_true", "method",
        "n_attempted", "n_success", "n_error", "n_invalid",
        "error_rate", "elevated_failure_flag", "mean_fit_runtime_sec"
      ),
      names(method_summary)
    )
    utils::write.csv(
      method_summary[, op_cols, drop = FALSE],
      file.path(out_dir, "summaries", "condition_method_summary.csv"),
      row.names = FALSE
    )
    mc <- method_summary[, intersect(
      c(
        "scenario_id", "beta1_true", "method", "n_success", "n_attempted",
        "elevated_failure_flag", "error_rate", "mean_fit_runtime_sec"
      ),
      names(method_summary)
    ), drop = FALSE]
    utils::write.csv(mc, file.path(out_dir, "summaries", "monte_carlo_uncertainty.csv"), row.names = FALSE)
    return(list(method_summary = method_summary, scientific = FALSE))
  }

  pwr <- summarize_power_calibration_s4(method_summary, scenario_grid)
  vcal <- summarize_variance_calibration_s4(method_summary, scenario_grid)
  ret <- summarize_retention_validation_s4(method_summary, scenario_grid)
  sup <- summarize_support_validation_s4(method_summary, scenario_grid, pair_diagnostics)
  pairs <- summarize_targeted_pairs_s4(method_summary, scenario_grid, pair_diagnostics)
  fam <- summarize_fragility_families_s4(method_summary, scenario_grid, family_diagnostics)
  comps <- summarize_comparator_performance_s4(method_summary, scenario_grid)
  het <- summarize_heterogeneity_sensitivity_s4(metric_library)
  t1 <- summarize_type1_robustness_s4(method_summary)

  write_one <- function(df, name) {
    if (!is.null(df) && is.data.frame(df) && nrow(df)) {
      utils::write.csv(df, file.path(out_dir, "summaries", name), row.names = FALSE)
    }
  }
  write_one(pwr, "power_calibration_summary.csv")
  write_one(vcal, "variance_calibration_summary.csv")
  write_one(comps, "comparator_performance_summary.csv")
  write_one(ret, "retention_validation_summary.csv")
  write_one(sup, "support_validation_summary.csv")
  write_one(pairs, "targeted_pair_summary.csv")
  write_one(fam, "fragility_family_summary.csv")
  write_one(het, "heterogeneity_sensitivity_summary.csv")
  write_one(t1, "type1_robustness_summary.csv")

  list(
    method_summary = method_summary,
    power_calibration = pwr,
    variance_calibration = vcal,
    retention = ret,
    support = sup,
    pairs = pairs,
    families = fam,
    comparators = comps,
    heterogeneity = het,
    type1 = t1,
    scientific = TRUE
  )
}

# Class E fragility from the frozen library and condition-method summary.

derive_study4_class_e_fragility <- function(scenario_library, condition_method_summary) {
  num <- function(x) suppressWarnings(as.numeric(x))
  lib <- scenario_library
  cms <- condition_method_summary
  required_lib <- c(
    "scenario_id", "targeted_class", "family_id", "tau", "omega",
    "delta_P_H_0.10", "predicted_power_tau_L_0.10", "predicted_power_tau_U_0.10",
    "predicted_power_H_0.10"
  )
  required_cms <- c("scenario_id", "method", "beta1_true", "power", "family_id", "tau_family_level")
  missing_lib <- setdiff(required_lib, names(lib))
  missing_cms <- setdiff(required_cms, names(cms))
  if (length(missing_lib) || length(missing_cms)) {
    stop(
      "Class E derivation is missing columns: ",
      paste(c(missing_lib, missing_cms), collapse = ", "),
      call. = FALSE
    )
  }

  e <- lib[!is.na(lib$targeted_class) & lib$targeted_class == "E", , drop = FALSE]
  if (nrow(e) != 24L) {
    stop("Expected 24 Class E scenarios in the scenario library, found ", nrow(e), call. = FALSE)
  }
  e$family_id <- as.integer(e$family_id)
  e$scenario_id <- as.integer(e$scenario_id)
  e$tau <- num(e$tau)
  e$omega <- num(e$omega)
  e$delta_P_H_0.10 <- num(e$delta_P_H_0.10)
  e$predicted_power_tau_L_0.10 <- num(e$predicted_power_tau_L_0.10)
  e$predicted_power_tau_U_0.10 <- num(e$predicted_power_tau_U_0.10)
  e$predicted_power_H_0.10 <- num(e$predicted_power_H_0.10)
  families <- sort(unique(e$family_id))
  if (!identical(families, 1:6)) {
    stop(
      "Expected Class E family_id values 1 through 6, found ",
      paste(families, collapse = ","),
      call. = FALSE
    )
  }

  cr <- cms[cms$method == "mlma_cr2" & abs(num(cms$beta1_true) - 0.10) < 1e-12, , drop = FALSE]
  if (nrow(cr) != 120L) {
    stop("Expected 120 mlma_cr2 beta1=0.10 rows, found ", nrow(cr), call. = FALSE)
  }
  cr$scenario_id <- as.integer(cr$scenario_id)
  cr$power <- num(cr$power)
  cr$family_id <- as.integer(cr$family_id)
  cr$tau_family_level <- num(cr$tau_family_level)

  tau_levels <- c(`0.05` = 0.05, `0.12` = 0.12, `0.20` = 0.20, `0.30` = 0.30)
  rows <- vector("list", length(families))
  for (i in seq_along(families)) {
    fam <- e[e$family_id == families[[i]], , drop = FALSE]
    if (nrow(fam) != 4L) {
      stop("Family ", families[[i]], " has ", nrow(fam), " scenarios.", call. = FALSE)
    }
    endpoint_delta <- fam$predicted_power_tau_L_0.10 - fam$predicted_power_tau_U_0.10
    stored_delta <- fam$delta_P_H_0.10
    if (length(unique(round(endpoint_delta, 12))) != 1L ||
        length(unique(round(stored_delta, 12))) != 1L ||
        max(abs(endpoint_delta - stored_delta)) > 1e-10) {
      stop(
        "Family ", families[[i]], " Delta P_H(.10) disagrees with the endpoint difference.",
        call. = FALSE
      )
    }
    stored_delta <- stored_delta[[1L]]
    omega_u <- unique(round(fam$omega, 12))
    if (length(omega_u) != 1L) {
      stop("Family ", families[[i]], " does not have one omega.", call. = FALSE)
    }
    picked <- vector("list", length(tau_levels))
    for (k in seq_along(tau_levels)) {
      hit <- fam[abs(fam$tau - unname(tau_levels[[k]])) < 1e-8, , drop = FALSE]
      if (nrow(hit) != 1L) {
        stop(
          "Family ", families[[i]], " does not have one scenario at tau=",
          names(tau_levels)[[k]], call. = FALSE
        )
      }
      emp <- cr[cr$scenario_id == hit$scenario_id, , drop = FALSE]
      if (nrow(emp) != 1L || !is.finite(emp$power)) {
        stop(
          "Expected one finite mlma_cr2 beta1=0.10 power for scenario ",
          hit$scenario_id, call. = FALSE
        )
      }
      if (!is.na(emp$family_id) && emp$family_id != families[[i]]) {
        stop(
          "condition_method_summary family_id disagrees with the library for scenario ",
          hit$scenario_id, call. = FALSE
        )
      }
      if (is.finite(emp$tau_family_level) &&
          abs(emp$tau_family_level - unname(tau_levels[[k]])) > 1e-8) {
        stop(
          "condition_method_summary tau_family_level disagrees with library tau for scenario ",
          hit$scenario_id, call. = FALSE
        )
      }
      picked[[k]] <- list(
        scenario_id = hit$scenario_id,
        tau = hit$tau,
        power = emp$power,
        predicted_power_H = hit$predicted_power_H_0.10
      )
    }
    delta_emp <- picked[[1L]]$power - picked[[4L]]$power
    rows[[i]] <- data.frame(
      family_id = families[[i]],
      n_scenarios = nrow(fam),
      scenario_id_tau_0.05 = picked[[1L]]$scenario_id,
      scenario_id_tau_0.12 = picked[[2L]]$scenario_id,
      scenario_id_tau_0.20 = picked[[3L]]$scenario_id,
      scenario_id_tau_0.30 = picked[[4L]]$scenario_id,
      tau_0.05 = picked[[1L]]$tau,
      tau_0.12 = picked[[2L]]$tau,
      tau_0.20 = picked[[3L]]$tau,
      tau_0.30 = picked[[4L]]$tau,
      omega = omega_u,
      delta_P_H_0.10 = stored_delta,
      power_0.10_tau_0.05 = picked[[1L]]$power,
      power_0.10_tau_0.12 = picked[[2L]]$power,
      power_0.10_tau_0.20 = picked[[3L]]$power,
      power_0.10_tau_0.30 = picked[[4L]]$power,
      predicted_power_H_0.10_tau_0.05 = picked[[1L]]$predicted_power_H,
      predicted_power_H_0.10_tau_0.12 = picked[[2L]]$predicted_power_H,
      predicted_power_H_0.10_tau_0.20 = picked[[3L]]$predicted_power_H,
      predicted_power_H_0.10_tau_0.30 = picked[[4L]]$predicted_power_H,
      delta_P_emp_0.10 = delta_emp,
      signed_error = stored_delta - delta_emp,
      abs_error = abs(stored_delta - delta_emp),
      stringsAsFactors = FALSE
    )
  }
  families_df <- do.call(rbind, rows)
  rownames(families_df) <- NULL
  families_df$n_families <- nrow(families_df)
  families_df$family_mae <- mean(families_df$abs_error)
  families_df$family_rmse <- sqrt(mean(families_df$signed_error^2))
  families_df$family_mean_signed_error <- mean(families_df$signed_error)
  families_df$spearman_delta_P_H_vs_emp <- as.numeric(stats::cor(
    families_df$delta_P_H_0.10, families_df$delta_P_emp_0.10, method = "spearman"
  ))
  families_df$method <- "mlma_cr2"
  families_df$beta1_true <- 0.10

  provenance <- data.frame(
    input_scenario_library = "docs/Simulation_Study_4/data/scenario_library.csv",
    input_condition_method_summary = "output/sim_study_4_production/summaries/condition_method_summary.csv",
    not_used_as_inputs = paste(
      "fragility_family_summary_reconstructed.csv",
      "manuscript prose",
      "fragility_family_summary.csv empirical columns",
      sep = "; "
    ),
    join_keys = paste(
      "Class E: scenario_library.targeted_class == E",
      "family: scenario_library.family_id",
      "tau endpoint: scenario_library.tau in {0.05, 0.12, 0.20, 0.30}",
      "empirical power: condition_method_summary.scenario_id + method == mlma_cr2 + beta1_true == 0.10",
      sep = "; "
    ),
    formula_delta_P_H = "prespecified scenario_library.delta_P_H_0.10, required to equal predicted_power_tau_L_0.10 - predicted_power_tau_U_0.10 within each family",
    formula_delta_P_emp = "power(tau=0.05) - power(tau=0.30) for mlma_cr2 at beta1_true=0.10",
    formula_signed_error = "delta_P_H_0.10 - delta_P_emp_0.10",
    formula_mae = "mean(abs(signed_error)) over the 6 families",
    formula_rmse = "sqrt(mean(signed_error^2)) over the 6 families",
    formula_spearman = "spearman correlation of delta_P_H_0.10 and delta_P_emp_0.10 over the 6 families",
    n_library_rows = nrow(lib),
    n_class_e_rows = nrow(e),
    n_families = nrow(families_df),
    n_scenarios_per_family = 4L,
    n_condition_rows_mlma_cr2_beta10 = nrow(cr),
    n_condition_rows_joined = 24L,
    n_endpoint_rows_used_for_loss = 12L,
    stringsAsFactors = FALSE
  )
  list(families = families_df, provenance = provenance)
}

compare_class_e_summaries <- function(derived, reconstruction) {
  num <- function(x) suppressWarnings(as.numeric(x))
  pairs <- c(
    delta_P_H_0.10 = "delta_P_H_0.10",
    power_0.10_tau_0.05 = "power_0.10_tau_0.05",
    power_0.10_tau_0.12 = "power_0.10_tau_0.12",
    power_0.10_tau_0.20 = "power_0.10_tau_0.20",
    power_0.10_tau_0.30 = "power_0.10_tau_0.30",
    predicted_power_H_0.10_tau_0.05 = "pred_power_tau_0.05",
    predicted_power_H_0.10_tau_0.12 = "pred_power_tau_0.12",
    predicted_power_H_0.10_tau_0.20 = "pred_power_tau_0.20",
    predicted_power_H_0.10_tau_0.30 = "pred_power_tau_0.30",
    delta_P_emp_0.10 = "delta_P_emp_0.10",
    signed_error = "signed_sensitivity_error",
    abs_error = "abs_sensitivity_error",
    scenario_id_tau_0.05 = "scenario_id_tau_0.05",
    scenario_id_tau_0.12 = "scenario_id_tau_0.12",
    scenario_id_tau_0.20 = "scenario_id_tau_0.20",
    scenario_id_tau_0.30 = "scenario_id_tau_0.30",
    family_mae = "family_mae_delta_P",
    family_rmse = "family_rmse_delta_P",
    spearman_delta_P_H_vs_emp = "spearman_delta_P_H_vs_emp"
  )
  rows <- list()
  for (fid in derived$family_id) {
    d <- derived[derived$family_id == fid, , drop = FALSE]
    r <- reconstruction[as.integer(reconstruction$family_id) == as.integer(fid), , drop = FALSE]
    if (nrow(d) != 1L || nrow(r) != 1L) {
      stop("Audit expected one derived and one reconstructed row for family ", fid, call. = FALSE)
    }
    for (nm in names(pairs)) {
      rows[[length(rows) + 1L]] <- data.frame(
        family_id = as.integer(fid),
        quantity = nm,
        derived = num(d[[nm]]),
        reconstruction = num(r[[pairs[[nm]]]]),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$abs_diff <- abs(out$derived - out$reconstruction)
  out$match <- is.finite(out$abs_diff) & out$abs_diff <= 1e-8
  list(comparisons = out, ok = all(out$match), max_abs_diff = max(out$abs_diff))
}

# Class E and planning-range post-processing uses the frozen scenario,
# family, and diagnostic libraries plus the condition-method summary.

fragility_family_summary_from_library <- function(scenario_library,
                                                    condition_method_summary,
                                                    family_diagnostics) {
  num <- function(x) suppressWarnings(as.numeric(x))
  cms <- condition_method_summary
  scen <- scenario_library
  fam <- family_diagnostics
  cr2_10 <- cms[
    cms$method == "mlma_cr2" & abs(num(cms$beta1_true) - 0.10) < 1e-12,
    ,
    drop = FALSE
  ]
  if (!nrow(cr2_10)) {
    stop("No mlma_cr2 beta1=0.10 rows in the condition-method summary.", call. = FALSE)
  }
  lookup_row <- function(sid) {
    hit <- cr2_10[as.integer(cr2_10$scenario_id) == as.integer(sid), , drop = FALSE]
    if (nrow(hit) != 1L) {
      stop("Expected one mlma_cr2 beta1=0.10 row for scenario ", sid, call. = FALSE)
    }
    hit
  }
  lookup_scen <- function(sid) {
    hit <- scen[as.integer(scen$scenario_id) == as.integer(sid), , drop = FALSE]
    if (nrow(hit) != 1L) {
      stop("Expected one scenario_library row for scenario ", sid, call. = FALSE)
    }
    hit
  }
  recon_rows <- vector("list", nrow(fam))
  for (i in seq_len(nrow(fam))) {
    fd <- fam[i, , drop = FALSE]
    sids <- as.integer(c(
      fd$scenario_id_tau_0.05, fd$scenario_id_tau_0.12,
      fd$scenario_id_tau_0.20, fd$scenario_id_tau_0.30
    ))
    powers <- numeric(4L)
    vars <- numeric(4L)
    preds <- numeric(4L)
    ich <- numeric(4L)
    n_success <- integer(4L)
    n_attempted <- integer(4L)
    true_tau <- numeric(4L)
    for (k in seq_len(4L)) {
      emp <- lookup_row(sids[[k]])
      sl <- lookup_scen(sids[[k]])
      powers[[k]] <- num(emp$power)
      vars[[k]] <- num(emp$empirical_variance)
      preds[[k]] <- num(sl$predicted_power_H_0.10)
      ich[[k]] <- num(sl$I_C_H_star)
      n_success[[k]] <- as.integer(emp$n_success)
      n_attempted[[k]] <- as.integer(emp$n_attempted)
      true_tau[[k]] <- num(sl$tau)
    }
    delta_emp <- powers[[1]] - powers[[4]]
    delta_h <- num(fd$delta_P_H_0.10)
    signed <- delta_h - delta_emp
    recon_rows[[i]] <- data.frame(
      family_id = as.integer(fd$family_id),
      profile = as.integer(fd$profile),
      profile_label = as.character(fd$profile_label),
      base_candidate_id = as.integer(fd$base_candidate_id),
      scenario_id_tau_0.05 = sids[[1]],
      scenario_id_tau_0.12 = sids[[2]],
      scenario_id_tau_0.20 = sids[[3]],
      scenario_id_tau_0.30 = sids[[4]],
      true_tau_0.05 = true_tau[[1]],
      true_tau_0.12 = true_tau[[2]],
      true_tau_0.20 = true_tau[[3]],
      true_tau_0.30 = true_tau[[4]],
      omega_family = num(fd$omega_family),
      F_I = num(fd$F_I),
      delta_P_H_0.10 = delta_h,
      power_0.10_tau_0.05 = powers[[1]],
      power_0.10_tau_0.12 = powers[[2]],
      power_0.10_tau_0.20 = powers[[3]],
      power_0.10_tau_0.30 = powers[[4]],
      var_tau_0.05 = vars[[1]],
      var_tau_0.12 = vars[[2]],
      var_tau_0.20 = vars[[3]],
      var_tau_0.30 = vars[[4]],
      pred_power_tau_0.05 = preds[[1]],
      pred_power_tau_0.12 = preds[[2]],
      pred_power_tau_0.20 = preds[[3]],
      pred_power_tau_0.30 = preds[[4]],
      I_C_H_star_tau_0.05 = ich[[1]],
      I_C_H_star_tau_0.30 = ich[[4]],
      n_success_tau_0.05 = n_success[[1]],
      n_success_tau_0.12 = n_success[[2]],
      n_success_tau_0.20 = n_success[[3]],
      n_success_tau_0.30 = n_success[[4]],
      n_attempted_tau_0.05 = n_attempted[[1]],
      n_attempted_tau_0.12 = n_attempted[[2]],
      n_attempted_tau_0.20 = n_attempted[[3]],
      n_attempted_tau_0.30 = n_attempted[[4]],
      delta_P_emp_0.10 = delta_emp,
      signed_sensitivity_error = signed,
      abs_sensitivity_error = abs(signed),
      method = "mlma_cr2",
      beta1_true = 0.10,
      reconstruction_status = "documentation_reconstruction",
      official_empirical_columns = "NA_in_fragility_family_summary.csv",
      reconstruction_inputs = paste(
        "condition_method_summary.csv (mlma_cr2, beta1=0.10);",
        "scenario_library.csv (predicted_power_H_0.10, I_C_H_star, tau);",
        "fragility_family_diagnostics.csv (profile, scenario IDs, F_I, delta_P_H_0.10)"
      ),
      stringsAsFactors = FALSE
    )
  }
  recon <- do.call(rbind, recon_rows)
  rownames(recon) <- NULL
  mae <- mean(recon$abs_sensitivity_error)
  rmse <- sqrt(mean(recon$signed_sensitivity_error^2))
  mean_signed <- mean(recon$signed_sensitivity_error)
  spearman_pearson <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) {
      return(c(spearman = NA_real_, pearson = NA_real_))
    }
    c(
      spearman = as.numeric(stats::cor(x, y, method = "spearman")),
      pearson = as.numeric(stats::cor(x, y, method = "pearson"))
    )
  }
  cor_dp <- spearman_pearson(recon$delta_P_H_0.10, recon$delta_P_emp_0.10)
  cor_fi <- spearman_pearson(recon$F_I, recon$delta_P_emp_0.10)
  recon$family_mae_delta_P <- mae
  recon$family_rmse_delta_P <- rmse
  recon$family_mean_signed_error <- mean_signed
  recon$spearman_delta_P_H_vs_emp <- unname(cor_dp[["spearman"]])
  recon$pearson_delta_P_H_vs_emp <- unname(cor_dp[["pearson"]])
  recon$spearman_F_I_vs_emp <- unname(cor_fi[["spearman"]])
  recon$pearson_F_I_vs_emp <- unname(cor_fi[["pearson"]])
  recon
}

heterogeneity_planning_grid_s4 <- function(diagnostic_library) {
  num <- function(x) suppressWarnings(as.numeric(x))
  diag <- diagnostic_library
  prac <- diag[diag$regime == "practical_tau_omega", , drop = FALSE]
  rho <- diag[diag$regime == "rho_sensitivity", , drop = FALSE]
  agg_cells <- function(d, regime_label) {
    keys <- unique(d[, c("tau_plan", "omega_plan", "rho_plan"), drop = FALSE])
    keys <- keys[order(keys$tau_plan, keys$omega_plan, keys$rho_plan), , drop = FALSE]
    rows <- vector("list", nrow(keys))
    for (j in seq_len(nrow(keys))) {
      hit <- d[
        abs(d$tau_plan - keys$tau_plan[[j]]) < 1e-12 &
          abs(d$omega_plan - keys$omega_plan[[j]]) < 1e-12 &
          abs(d$rho_plan - keys$rho_plan[[j]]) < 1e-12,
        ,
        drop = FALSE
      ]
      rows[[j]] <- data.frame(
        table_kind = "grid_cell",
        regime = regime_label,
        tau_plan = keys$tau_plan[[j]],
        omega_plan = keys$omega_plan[[j]],
        rho_plan = keys$rho_plan[[j]],
        n_scenarios = length(unique(hit$scenario_id)),
        n_rows = nrow(hit),
        mean_I_C_star = mean(num(hit$I_C_star), na.rm = TRUE),
        mean_predicted_power_0.10 = mean(num(hit$predicted_power_0.10), na.rm = TRUE),
        aggregation = "mean over 120 scenarios at this planning cell",
        quantity_type = "design_only_predicted",
        not_a_confidence_or_prediction_interval = TRUE,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  marginal <- function(d, regime_label, vary, label) {
    vals <- sort(unique(d[[vary]]))
    rows <- lapply(vals, function(v) {
      hit <- d[abs(d[[vary]] - v) < 1e-12, , drop = FALSE]
      data.frame(
        table_kind = label,
        regime = regime_label,
        tau_plan = if (identical(vary, "tau_plan")) v else NA_real_,
        omega_plan = if (identical(vary, "omega_plan")) v else NA_real_,
        rho_plan = if (identical(vary, "rho_plan")) v else NA_real_,
        n_scenarios = length(unique(hit$scenario_id)),
        n_rows = nrow(hit),
        mean_I_C_star = mean(num(hit$I_C_star), na.rm = TRUE),
        mean_predicted_power_0.10 = mean(num(hit$predicted_power_0.10), na.rm = TRUE),
        aggregation = paste0(
          "mean over scenarios and the other planning factors in ", regime_label
        ),
        quantity_type = "design_only_predicted",
        not_a_confidence_or_prediction_interval = TRUE,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
  grid_doc <- rbind(
    agg_cells(prac, "practical_tau_omega"),
    agg_cells(rho, "rho_sensitivity"),
    marginal(prac, "practical_tau_omega", "tau_plan", "marginal_tau"),
    marginal(prac, "practical_tau_omega", "omega_plan", "marginal_omega"),
    marginal(rho, "rho_sensitivity", "rho_plan", "marginal_rho")
  )
  rownames(grid_doc) <- NULL
  grid_doc
}

heterogeneity_planning_range_s4 <- function(diagnostic_library) {
  num <- function(x) suppressWarnings(as.numeric(x))
  diag <- diagnostic_library
  prac <- diag[diag$regime == "practical_tau_omega", , drop = FALSE]
  rho <- diag[diag$regime == "rho_sensitivity", , drop = FALSE]
  range_one <- function(d, vary, holds, grid_label) {
    d$group_key <- apply(
      d[, c("scenario_id", holds), drop = FALSE],
      1L,
      function(r) paste(r, collapse = "|")
    )
    rng <- tapply(num(d$predicted_power_0.10), d$group_key, function(z) {
      z <- z[is.finite(z)]
      if (length(z) < 2L) return(NA_real_)
      max(z) - min(z)
    })
    rng <- as.numeric(rng)
    rng <- rng[is.finite(rng)]
    data.frame(
      varied_component = vary,
      held_components = paste(holds, collapse = ","),
      grid = grid_label,
      n_groups = length(rng),
      mean_range_predicted_power_0.10 = mean(rng),
      median_range_predicted_power_0.10 = stats::median(rng),
      max_range_predicted_power_0.10 = max(rng),
      quantity_type = "design_only_predicted_power_range",
      not_a_confidence_or_prediction_interval = TRUE,
      stringsAsFactors = FALSE
    )
  }
  range_doc <- rbind(
    range_one(prac, "tau_plan", c("omega_plan", "rho_plan"), "practical_tau_omega"),
    range_one(prac, "omega_plan", c("tau_plan", "rho_plan"), "practical_tau_omega"),
    range_one(rho, "rho_plan", c("tau_plan", "omega_plan"), "rho_sensitivity")
  )
  rownames(range_doc) <- NULL
  range_doc
}

write_study4_derived_tables <- function(out_dir,
                                                   scenario_library,
                                                   condition_method_summary,
                                                   family_diagnostics,
                                                   diagnostic_library) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  derived <- derive_study4_class_e_fragility(
    scenario_library, condition_method_summary
  )
  recon <- fragility_family_summary_from_library(
    scenario_library, condition_method_summary, family_diagnostics
  )
  audit <- compare_class_e_summaries(derived$families, recon)
  if (!isTRUE(audit$ok)) {
    stop("Class E derivation does not match the reconstruction.", call. = FALSE)
  }
  grid_doc <- heterogeneity_planning_grid_s4(diagnostic_library)
  range_doc <- heterogeneity_planning_range_s4(diagnostic_library)
  utils::write.csv(
    derived$families,
    file.path(out_dir, "fragility_family_summary_derived.csv"),
    row.names = FALSE, na = "NA"
  )
  utils::write.csv(
    derived$provenance,
    file.path(out_dir, "fragility_family_summary_derived_provenance.csv"),
    row.names = FALSE, na = "NA"
  )
  utils::write.csv(
    audit$comparisons,
    file.path(out_dir, "fragility_family_summary_derived_vs_reconstruction_audit.csv"),
    row.names = FALSE, na = "NA"
  )
  utils::write.csv(
    recon,
    file.path(out_dir, "fragility_family_summary_reconstructed.csv"),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
  )
  utils::write.csv(
    grid_doc,
    file.path(out_dir, "heterogeneity_sensitivity_grid_documentation.csv"),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
  )
  utils::write.csv(
    range_doc,
    file.path(out_dir, "heterogeneity_sensitivity_range_documentation.csv"),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
  )
  invisible(list(
    derived = derived,
    reconstruction = recon,
    audit = audit,
    planning_grid = grid_doc,
    planning_range = range_doc
  ))
}
