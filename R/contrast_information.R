# Contrast-specific information metrics used in Simulation Studies 3 and 4.
# These quantities depend only on the study design, not on simulated outcomes
# or fitted models.

sim3_df_epsilon <- 1e-10

# Relative difference: |a-b| / ((a+b)/2).
relative_difference_s3 <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  denom <- 0.5 * (a + b)
  out <- abs(a - b) / denom
  out[!is.finite(denom) | denom <= 0] <- NA_real_
  out
}

# Design-only sampling-information study weight.
#
# w_j^(S) = k_j (n_j - 3) / (1 + (k_j - 1) rho_plan)
study_sampling_weight_s3 <- function(n_j, k_j, rho_plan) {
  n_j <- as.numeric(n_j)
  k_j <- as.numeric(k_j)
  rho_plan <- as.numeric(rho_plan)
  if (length(rho_plan) != 1L || !is.finite(rho_plan)) {
    stop("rho_plan must be a single finite number.", call. = FALSE)
  }
  if (length(n_j) != length(k_j)) {
    stop("n_j and k_j must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(n_j)) || any(!is.finite(k_j))) {
    stop("n_j and k_j must be finite.", call. = FALSE)
  }
  if (any(n_j <= 3)) {
    stop("All n_j must be greater than 3.", call. = FALSE)
  }
  if (any(k_j < 1)) {
    stop("All k_j must be at least 1.", call. = FALSE)
  }
  denom <- 1 + (k_j - 1) * rho_plan
  if (any(denom <= 0)) {
    stop("1 + (k_j - 1) * rho_plan must be positive.", call. = FALSE)
  }
  k_j * (n_j - 3) / denom
}

# Heterogeneity-adjusted planning weight (secondary framework).
#
# w_j^(H) = [tau_plan^2 + omega_plan^2 / k_j + (1+(k_j-1)rho_plan)/(k_j(n_j-3))]^{-1}
study_heterogeneity_weight_s3 <- function(n_j, k_j, rho_plan, tau_plan, omega_plan) {
  n_j <- as.numeric(n_j)
  k_j <- as.numeric(k_j)
  rho_plan <- as.numeric(rho_plan)
  tau_plan <- as.numeric(tau_plan)
  omega_plan <- as.numeric(omega_plan)
  if (length(rho_plan) != 1L || !is.finite(rho_plan)) {
    stop("rho_plan must be a single finite number.", call. = FALSE)
  }
  if (length(tau_plan) != 1L || !is.finite(tau_plan) || tau_plan < 0) {
    stop("tau_plan must be a single finite non-negative number.", call. = FALSE)
  }
  if (length(omega_plan) != 1L || !is.finite(omega_plan) || omega_plan < 0) {
    stop("omega_plan must be a single finite non-negative number.", call. = FALSE)
  }
  w_s_inv <- (1 + (k_j - 1) * rho_plan) / (k_j * (n_j - 3))
  v <- tau_plan^2 + (omega_plan^2) / k_j + w_s_inv
  if (any(!is.finite(v) | v <= 0)) {
    stop("Heterogeneity-adjusted variance components must be positive.", call. = FALSE)
  }
  1 / v
}

# Kish-style effective cluster support for a vector of weights.
#
# J* = (sum w)^2 / sum(w^2). Empty or all-zero weight vectors return 0.
effective_cluster_support_s3 <- function(w) {
  w <- as.numeric(w)
  if (!length(w)) {
    return(0)
  }
  w <- w[is.finite(w)]
  if (!length(w)) {
    return(0)
  }
  sw <- sum(w)
  sw2 <- sum(w^2)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) {
    return(0)
  }
  (sw^2) / sw2
}

normalize_study_structure_s3 <- function(study_structure) {
  if (is.null(study_structure) || !nrow(study_structure)) {
    stop("study_structure is empty.", call. = FALSE)
  }
  required <- c("study", "X", "n_j", "k_j")
  missing <- setdiff(required, names(study_structure))
  if (length(missing)) {
    stop(
      "study_structure missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out <- data.frame(
    study = as.integer(study_structure$study),
    X = as.integer(study_structure$X),
    n_j = as.numeric(study_structure$n_j),
    k_j = as.numeric(study_structure$k_j),
    stringsAsFactors = FALSE
  )
  if (any(is.na(out$X)) || any(!out$X %in% c(0L, 1L))) {
    stop("X must be a binary study-level moderator in {0, 1}.", call. = FALSE)
  }
  out
}

# Contrast-specific information, effective cluster support, and anticipated df.
# weight_type = "sampling" uses design-only weights; "heterogeneity_adjusted"
# additionally uses tau_plan and omega_plan.
compute_contrast_information_s3 <- function(study_structure,
                                            rho_plan = 0.50,
                                            tau_plan = NULL,
                                            omega_plan = NULL,
                                            weight_type = c("sampling", "heterogeneity_adjusted"),
                                            return_weights = FALSE) {
  weight_type <- match.arg(weight_type)
  ss <- normalize_study_structure_s3(study_structure)

  if (identical(weight_type, "sampling")) {
    w <- study_sampling_weight_s3(ss$n_j, ss$k_j, rho_plan)
    tau_out <- NA_real_
    omega_out <- NA_real_
  } else {
    if (is.null(tau_plan) || is.null(omega_plan)) {
      stop(
        "tau_plan and omega_plan are required for heterogeneity_adjusted weights.",
        call. = FALSE
      )
    }
    w <- study_heterogeneity_weight_s3(
      ss$n_j, ss$k_j, rho_plan, tau_plan, omega_plan
    )
    tau_out <- as.numeric(tau_plan)
    omega_out <- as.numeric(omega_plan)
  }

  g0 <- ss$X == 0L
  g1 <- ss$X == 1L
  w0 <- w[g0]
  w1 <- w[g1]
  W0 <- sum(w0)
  W1 <- sum(w1)
  J0_star <- effective_cluster_support_s3(w0)
  J1_star <- effective_cluster_support_s3(w1)
  J0 <- as.integer(sum(g0))
  J1 <- as.integer(sum(g1))

  empty_pred <- function(I_C_star, I_raw, nu, status) {
    data.frame(
      weight_type = weight_type,
      rho_plan = as.numeric(rho_plan),
      tau_plan = tau_out,
      omega_plan = omega_out,
      J0 = J0,
      J1 = J1,
      W0 = as.numeric(W0),
      W1 = as.numeric(W1),
      I_raw = as.numeric(I_raw),
      I_C_star = as.numeric(I_C_star),
      predicted_variance = NA_real_,
      predicted_se = NA_real_,
      J0_star = as.numeric(J0_star),
      J1_star = as.numeric(J1_star),
      nu_C_star = as.numeric(nu),
      status = as.character(status),
      stringsAsFactors = FALSE
    )
  }

  if (!is.finite(W0) || !is.finite(W1) || W0 <= 0 || W1 <= 0) {
    row <- empty_pred(
      I_C_star = 0,
      I_raw = 0,
      nu = NA_real_,
      status = "insufficient_group_information"
    )
    if (isTRUE(return_weights)) {
      attr(row, "weights") <- w
    }
    return(row)
  }

  I_raw <- (W0 * W1) / (W0 + W1)
  I_C_star <- 4 * I_raw
  predicted_variance <- 4 / I_C_star
  predicted_se <- 2 / sqrt(I_C_star)
  predicted_se_w <- sqrt(1 / W0 + 1 / W1)
  if (is.finite(predicted_se) && is.finite(predicted_se_w) &&
      abs(predicted_se - predicted_se_w) > 1e-10 * max(1, predicted_se_w)) {
    stop("Contrast standard-error identity does not hold.", call. = FALSE)
  }

  V0 <- 1 / W0
  V1 <- 1 / W1
  d0 <- J0_star - 1
  d1 <- J1_star - 1
  eps <- sim3_df_epsilon

  # Effective support near one.
  # Mathematical limit is zero when either group's (J_g* - 1) is non-positive.
  if (d0 <= eps || d1 <= eps) {
    nu_C_star <- 0
    status <- "zero_cluster_support"
  } else {
    denom <- (V0^2) / d0 + (V1^2) / d1
    nu_C_star <- ((V0 + V1)^2) / denom
    if (!is.finite(nu_C_star)) {
      nu_C_star <- NA_real_
      status <- "numerical_failure"
    } else if (nu_C_star < 0) {
      # Negative df is not a defined Satterthwaite support quantity.
      nu_C_star <- NA_real_
      status <- "numerical_failure"
    } else if (nu_C_star == 0) {
      status <- "zero_cluster_support"
    } else if (nu_C_star < 1) {
      # Retain very small positive df; do not floor.
      status <- "very_low_df"
    } else {
      status <- "ok"
    }
  }

  row <- data.frame(
    weight_type = weight_type,
    rho_plan = as.numeric(rho_plan),
    tau_plan = tau_out,
    omega_plan = omega_out,
    J0 = J0,
    J1 = J1,
    W0 = as.numeric(W0),
    W1 = as.numeric(W1),
    I_raw = as.numeric(I_raw),
    I_C_star = as.numeric(I_C_star),
    predicted_variance = as.numeric(predicted_variance),
    predicted_se = as.numeric(predicted_se),
    J0_star = as.numeric(J0_star),
    J1_star = as.numeric(J1_star),
    nu_C_star = as.numeric(nu_C_star),
    status = as.character(status),
    stringsAsFactors = FALSE
  )
  if (isTRUE(return_weights)) {
    attr(row, "weights") <- w
  }
  row
}

# Theoretical two-sided noncentral-t power (primary power prediction).
#
# Does not floor nu_C_star. Returns NA when df is non-positive or non-finite.
predict_contrast_power_s3 <- function(beta1,
                                      predicted_se,
                                      nu_C_star,
                                      alpha = 0.05) {
  beta1 <- as.numeric(beta1)
  predicted_se <- as.numeric(predicted_se)
  nu_C_star <- as.numeric(nu_C_star)
  alpha <- as.numeric(alpha)
  n <- max(length(beta1), length(predicted_se), length(nu_C_star))
  if (length(beta1) == 1L) beta1 <- rep(beta1, n)
  if (length(predicted_se) == 1L) predicted_se <- rep(predicted_se, n)
  if (length(nu_C_star) == 1L) nu_C_star <- rep(nu_C_star, n)
  if (length(alpha) == 1L) alpha <- rep(alpha, n)

  out <- rep(NA_real_, n)
  ok <- is.finite(beta1) &
    is.finite(predicted_se) & predicted_se > 0 &
    is.finite(nu_C_star) & nu_C_star > 0 &
    is.finite(alpha) & alpha > 0 & alpha < 1
  if (!any(ok)) {
    return(out)
  }
  delta <- beta1[ok] / predicted_se[ok]
  crit <- stats::qt(1 - alpha[ok] / 2, df = nu_C_star[ok])
  # Two-sided: P(T < -c) + P(T > c) with noncentrality delta.
  lower <- stats::pt(-crit, df = nu_C_star[ok], ncp = delta)
  upper <- stats::pt(crit, df = nu_C_star[ok], ncp = delta, lower.tail = FALSE)
  pow <- lower + upper
  pow[!is.finite(pow)] <- NA_real_
  out[ok] <- pow
  out
}

comparator_quantities_s3 <- function(study_structure) {
  ss <- normalize_study_structure_s3(study_structure)
  n_j <- ss$n_j
  k_j <- ss$k_j
  X <- ss$X
  J <- nrow(ss)
  J0 <- as.integer(sum(X == 0L))
  J1 <- as.integer(sum(X == 1L))
  N <- sum(n_j)
  K <- sum(k_j)
  precision_j <- n_j - 3
  j_prec <- (sum(precision_j)^2) / sum(precision_j^2)
  j_k <- (sum(k_j)^2) / sum(k_j^2)
  min_count <- min(J0, J1)
  min_prop <- min_count / J
  minority_is_1 <- J1 <= J0
  min_prec_share <- if (minority_is_1) {
    sum(precision_j[X == 1L]) / sum(precision_j)
  } else {
    sum(precision_j[X == 0L]) / sum(precision_j)
  }
  min_k_share <- if (minority_is_1) {
    sum(k_j[X == 1L]) / sum(k_j)
  } else {
    sum(k_j[X == 0L]) / sum(k_j)
  }
  data.frame(
    J = as.integer(J),
    J0 = J0,
    J1 = J1,
    minority_study_count = as.integer(min_count),
    minority_prop = as.numeric(min_prop),
    N = as.numeric(N),
    K = as.numeric(K),
    mean_n = mean(n_j),
    mean_k = mean(k_j),
    n_cv = stats::sd(n_j) / mean(n_j),
    k_cv = stats::sd(k_j) / mean(k_j),
    J_prec = as.numeric(j_prec),
    J_prec_over_J = as.numeric(j_prec / J),
    J_k = as.numeric(j_k),
    J_k_over_J = as.numeric(j_k / J),
    minority_precision_share = as.numeric(min_prec_share),
    minority_effect_share = as.numeric(min_k_share),
    largest_precision_share = max(precision_j) / sum(precision_j),
    largest_effect_share = max(k_j) / sum(k_j),
    spearman_n_k = if (exists("safe_cor_s2", mode = "function")) {
      safe_cor_s2(n_j, k_j, "spearman")
    } else if (stats::sd(n_j) > 0 && stats::sd(k_j) > 0) {
      as.numeric(stats::cor(n_j, k_j, method = "spearman"))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

# Metric-variant catalog for every frozen scenario (long format).
#
# Sampling: rho_plan in {0.20, 0.50, 0.80} plus oracle rho_true.
# Heterogeneity-adjusted: tau_plan = omega_plan in {0.05, 0.10, 0.20} crossed
# with rho_plan in {0.20, 0.50, 0.80}, plus oracle (true tau, omega, rho).
sim3_metric_variant_spec <- function(rho_true = NA_real_,
                                     tau_true = NA_real_,
                                     omega_true = NA_real_) {
  sampling <- data.frame(
    variant_id = c(
      "sampling_rho_0.20",
      "sampling_rho_0.50",
      "sampling_rho_0.80",
      "sampling_rho_oracle"
    ),
    weight_type = "sampling",
    rho_plan = c(0.20, 0.50, 0.80, as.numeric(rho_true)),
    tau_plan = NA_real_,
    omega_plan = NA_real_,
    is_primary = c(FALSE, TRUE, FALSE, FALSE),
    is_oracle_sampling = c(FALSE, FALSE, FALSE, TRUE),
    is_default_het = FALSE,
    is_oracle_het = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- expand.grid(
    het = c(0.05, 0.10, 0.20),
    rho = c(0.20, 0.50, 0.80),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  het_sens <- data.frame(
    variant_id = sprintf(
      "het_tau%.2f_omega%.2f_rho%.2f",
      grid$het, grid$het, grid$rho
    ),
    weight_type = "heterogeneity_adjusted",
    rho_plan = grid$rho,
    tau_plan = grid$het,
    omega_plan = grid$het,
    is_primary = FALSE,
    is_oracle_sampling = FALSE,
    is_default_het = abs(grid$het - 0.10) < 1e-12 & abs(grid$rho - 0.50) < 1e-12,
    is_oracle_het = FALSE,
    stringsAsFactors = FALSE
  )
  het_oracle <- data.frame(
    variant_id = "het_oracle",
    weight_type = "heterogeneity_adjusted",
    rho_plan = as.numeric(rho_true),
    tau_plan = as.numeric(tau_true),
    omega_plan = as.numeric(omega_true),
    is_primary = FALSE,
    is_oracle_sampling = FALSE,
    is_default_het = FALSE,
    is_oracle_het = TRUE,
    stringsAsFactors = FALSE
  )
  rbind(sampling, het_sens, het_oracle)
}

compute_all_metric_variants_s3 <- function(study_structure,
                                           rho_true,
                                           tau_true,
                                           omega_true,
                                           scenario_id = NA_integer_) {
  spec <- sim3_metric_variant_spec(rho_true, tau_true, omega_true)
  rows <- vector("list", nrow(spec))
  for (i in seq_len(nrow(spec))) {
    v <- spec[i, , drop = FALSE]
    info <- compute_contrast_information_s3(
      study_structure = study_structure,
      rho_plan = v$rho_plan[[1]],
      tau_plan = if (is.na(v$tau_plan[[1]])) NULL else v$tau_plan[[1]],
      omega_plan = if (is.na(v$omega_plan[[1]])) NULL else v$omega_plan[[1]],
      weight_type = v$weight_type[[1]]
    )
    rows[[i]] <- cbind(
      data.frame(
        scenario_id = as.integer(scenario_id),
        variant_id = v$variant_id[[1]],
        is_primary = v$is_primary[[1]],
        is_oracle_sampling = v$is_oracle_sampling[[1]],
        is_default_het = v$is_default_het[[1]],
        is_oracle_het = v$is_oracle_het[[1]],
        stringsAsFactors = FALSE
      ),
      info
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

predict_power_grid_s3 <- function(metric_row,
                                  beta1_levels = c(0, 0.05, 0.10, 0.20),
                                  alpha = 0.05) {
  if (is.data.frame(metric_row)) {
    if (nrow(metric_row) != 1L) {
      stop("metric_row must have exactly one row.", call. = FALSE)
    }
    metric_row <- as.list(metric_row[1, , drop = FALSE])
  }
  data.frame(
    beta1_true = as.numeric(beta1_levels),
    predicted_se = as.numeric(metric_row$predicted_se),
    nu_C_star = as.numeric(metric_row$nu_C_star),
    predicted_ncp = as.numeric(beta1_levels) / as.numeric(metric_row$predicted_se),
    predicted_power = predict_contrast_power_s3(
      beta1 = beta1_levels,
      predicted_se = metric_row$predicted_se,
      nu_C_star = metric_row$nu_C_star,
      alpha = alpha
    ),
    metric_status = as.character(metric_row$status),
    stringsAsFactors = FALSE
  )
}

safe_cor_s2 <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    return(NA_real_)
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::cor(x, y, method = method))
}
