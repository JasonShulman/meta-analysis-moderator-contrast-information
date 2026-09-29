# Shared estimator functions for Simulation Studies 2-4.
# Study 2 and Study 3 use aggregate HKSJ, MLMA, MLMA+CR2, and CHE-RVE.
# Study 4 uses MLMA and MLMA+CR2.

# =============================================================================
# Warning / error capture
# =============================================================================

capture_conditions <- function(expr) {
  warnings <- character(0)
  started <- proc.time()[[3]]
  value <- NULL
  error <- NA_character_

  value <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }
  )

  list(
    value = value,
    error = error,
    warnings = warnings,
    elapsed_sec = as.numeric(proc.time()[[3]] - started)
  )
}

collapse_warnings_s2 <- function(warnings, separator = " || ") {
  if (!length(warnings) || all(is.na(warnings))) {
    return(NA_character_)
  }
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (!length(warnings)) {
    return(NA_character_)
  }
  paste(warnings, collapse = separator)
}

empty_method_result_s2 <- function(method,
                                   condition,
                                   rep_id,
                                   error = NA_character_,
                                   warnings = character(0),
                                   elapsed_sec = NA_real_,
                                   status = "error",
                                   config = make_sim2_config("smoke")) {
  condition <- normalize_condition_s2(condition)
  data.frame(
    condition_id = as.integer(condition$condition_id),
    condition_set = as.character(condition$condition_set),
    rep_id = as.integer(rep_id),
    method = as.character(method),
    estimate = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    p_value = NA_real_,
    df = NA_real_,
    converged = FALSE,
    status = as.character(status),
    error = if (length(error)) as.character(error[[1]]) else NA_character_,
    warning_count = as.integer(length(warnings[!is.na(warnings) & nzchar(warnings)])),
    warnings = collapse_warnings_s2(warnings, config$warning_separator),
    fit_runtime_sec = as.numeric(elapsed_sec),
    variance_study = NA_real_,
    variance_effect = NA_real_,
    variance_residual = NA_real_,
    boundary_study = NA,
    boundary_effect = NA,
    boundary_residual = NA,
    stringsAsFactors = FALSE
  )
}

boundary_flag_s2 <- function(variance, tolerance) {
  if (length(variance) != 1L || is.na(variance)) {
    return(NA)
  }
  isTRUE(as.numeric(variance) <= as.numeric(tolerance))
}

# =============================================================================
# clubSandwich interface check
# =============================================================================

.sim2_clubsandwich_api_ok <- FALSE

# Verify clubSandwich::coef_test formals once per session.
#
# `cluster` is accepted via `...` in clubSandwich 0.7.0, so it is not required
# as a named formal. Live output columns are checked in assert_coef_test_columns_s2().
verify_clubsandwich_api_s2 <- function() {
  if (isTRUE(.sim2_clubsandwich_api_ok)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    stop("Package 'clubSandwich' is required for Study 2 CR2 inference.", call. = FALSE)
  }
  fml <- names(formals(clubSandwich::coef_test))
  required_args <- c("obj", "vcov", "test", "coefs")
  missing_args <- setdiff(required_args, fml)
  if (length(missing_args)) {
    stop(
      "clubSandwich::coef_test() API changed; missing formals: ",
      paste(missing_args, collapse = ", "),
      call. = FALSE
    )
  }
  if (!("..." %in% fml)) {
    stop(
      "clubSandwich::coef_test() API changed; expected `...` for cluster=.",
      call. = FALSE
    )
  }
  .sim2_clubsandwich_api_ok <<- TRUE
  invisible(TRUE)
}

assert_coef_test_columns_s2 <- function(ct_df) {
  required <- c("beta", "SE", "df_Satt", "p_Satt")
  missing <- setdiff(required, names(ct_df))
  if (length(missing)) {
    stop(
      "clubSandwich::coef_test() output is missing expected columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# =============================================================================
# Result validation
# =============================================================================

validate_method_result_s2 <- function(row, config) {
  if (nrow(row) != 1L) {
    stop("validate_method_result_s2() expects a one-row data frame.", call. = FALSE)
  }
  # Upstream / hard errors are left unchanged.
  if (row$status %in% c("error", "upstream_fit_failure")) {
    return(row)
  }

  requires_df <- row$method %in% c("aggregate_hksj", "mlma_cr2", "che_rve")
  ok <- is.finite(row$estimate) &&
    is.finite(row$se) && row$se > 0 &&
    is.finite(row$ci_lower) &&
    is.finite(row$ci_upper) &&
    row$ci_lower <= row$ci_upper &&
    is.finite(row$p_value) &&
    row$p_value >= 0 && row$p_value <= 1

  if (requires_df) {
    ok <- ok && is.finite(row$df) && row$df > 0
  }

  if (!isTRUE(ok)) {
    row$status <- "invalid_result"
    row$converged <- FALSE
    detail <- paste0(
      "Invalid inferential quantities for method '", row$method, "'",
      if (requires_df) " (including df checks)" else "",
      "."
    )
    row$error <- if (is.na(row$error) || !nzchar(row$error)) {
      detail
    } else {
      paste(row$error, detail, sep = config$warning_separator)
    }
  } else if (identical(row$status, "success") || identical(row$status, "")) {
    row$status <- "success"
  }
  row
}

finalize_method_result_s2 <- function(row, config) {
  row$boundary_study <- boundary_flag_s2(row$variance_study, config$boundary_tolerance)
  row$boundary_effect <- boundary_flag_s2(row$variance_effect, config$boundary_tolerance)
  row$boundary_residual <- boundary_flag_s2(
    row$variance_residual,
    config$boundary_tolerance
  )
  validate_method_result_s2(row, config)
}

# =============================================================================
# Aggregation + Knapp-Hartung
# =============================================================================

# Dependence-aware study-level aggregation for HKSJ.
aggregate_study_effects_s2 <- function(data, rho_assumed = 0.50) {
  required <- c("study", "yi", "vi", "X", "k_j", "n_j")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "aggregate_study_effects_s2() missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(data)) {
    stop("data is empty.", call. = FALSE)
  }

  studies <- unique(as.integer(data$study))
  rows <- lapply(studies, function(s) {
    d <- data[as.integer(data$study) == s, , drop = FALSE]
    if (length(unique(d$X)) != 1L ||
        length(unique(d$n_j)) != 1L ||
        length(unique(d$k_j)) != 1L ||
        length(unique(d$vi)) != 1L) {
      stop(
        "X, n_j, k_j, and vi must be constant within study ", s, ".",
        call. = FALSE
      )
    }
    k_j <- as.integer(unique(d$k_j))
    if (nrow(d) != k_j) {
      stop(
        "Study ", s, " has ", nrow(d), " rows but k_j = ", k_j, ".",
        call. = FALSE
      )
    }
    vi <- as.numeric(unique(d$vi))
    data.frame(
      study = as.integer(s),
      X = as.integer(unique(d$X)),
      k_j = k_j,
      n_j = as.numeric(unique(d$n_j)),
      ybar = mean(as.numeric(d$yi)),
      vbar = vi * (rho_assumed + (1 - rho_assumed) / k_j),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

fit_aggregate_hksj_s2 <- function(data, condition, rep_id, config) {
  method <- "aggregate_hksj"
  condition <- normalize_condition_s2(condition)

  cap_agg <- capture_conditions({
    aggregate_study_effects_s2(data, rho_assumed = config$rho_assumed_aggregation)
  })
  if (!is.na(cap_agg$error)) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method, condition, rep_id,
        error = cap_agg$error,
        warnings = cap_agg$warnings,
        elapsed_sec = cap_agg$elapsed_sec,
        status = "error",
        config = config
      ),
      config
    ))
  }

  agg <- cap_agg$value
  cap <- capture_conditions({
    metafor::rma.uni(
      yi = ybar,
      vi = vbar,
      mods = ~ X,
      method = "REML",
      test = "knha",
      data = agg
    )
  })
  elapsed <- cap_agg$elapsed_sec + cap$elapsed_sec
  warns <- c(cap_agg$warnings, cap$warnings)

  if (!is.na(cap$error) || is.null(cap$value)) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method, condition, rep_id,
        error = cap$error,
        warnings = warns,
        elapsed_sec = elapsed,
        status = "error",
        config = config
      ),
      config
    ))
  }

  fit <- cap$value
  cn <- names(stats::coef(fit))
  if (!"X" %in% cn) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method, condition, rep_id,
        error = "Could not locate coefficient 'X' in rma.uni() fit.",
        warnings = warns,
        elapsed_sec = elapsed,
        status = "error",
        config = config
      ),
      config
    ))
  }
  idx <- which(cn == "X")
  converged <- if (is.null(fit$converged)) TRUE else isTRUE(fit$converged)

  row <- data.frame(
    condition_id = as.integer(condition$condition_id),
    condition_set = as.character(condition$condition_set),
    rep_id = as.integer(rep_id),
    method = method,
    estimate = as.numeric(stats::coef(fit)[idx]),
    se = as.numeric(fit$se[idx]),
    ci_lower = as.numeric(fit$ci.lb[idx]),
    ci_upper = as.numeric(fit$ci.ub[idx]),
    p_value = as.numeric(fit$pval[idx]),
    df = as.numeric(fit$k - fit$p),
    converged = converged,
    status = "success",
    error = NA_character_,
    warning_count = as.integer(length(warns)),
    warnings = collapse_warnings_s2(warns, config$warning_separator),
    fit_runtime_sec = elapsed,
    variance_study = NA_real_,
    variance_effect = NA_real_,
    variance_residual = as.numeric(fit$tau2),
    boundary_study = NA,
    boundary_effect = NA,
    boundary_residual = NA,
    stringsAsFactors = FALSE
  )
  finalize_method_result_s2(row, config)
}

# =============================================================================
# MLMA
# =============================================================================

# Fit the Study 2 MLMA working model once (diagonal V).
fit_mlma_base_s2 <- function(data) {
  capture_conditions({
    dat <- data
    dat$study <- factor(dat$study)
    dat$es_id <- factor(dat$es_id)
    metafor::rma.mv(
      yi = yi,
      V = vi,
      mods = ~ X,
      random = ~ 1 | study / es_id,
      method = "REML",
      sparse = TRUE,
      data = dat
    )
  })
}

extract_mlma_model_s2 <- function(fit_capture, condition, rep_id, config) {
  method <- "mlma"
  condition <- normalize_condition_s2(condition)

  if (!is.na(fit_capture$error) || is.null(fit_capture$value)) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method, condition, rep_id,
        error = fit_capture$error,
        warnings = fit_capture$warnings,
        elapsed_sec = fit_capture$elapsed_sec,
        status = "upstream_fit_failure",
        config = config
      ),
      config
    ))
  }

  fit <- fit_capture$value
  cn <- names(stats::coef(fit))
  if (!"X" %in% cn) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method, condition, rep_id,
        error = "Could not locate coefficient 'X' in rma.mv() fit.",
        warnings = fit_capture$warnings,
        elapsed_sec = fit_capture$elapsed_sec,
        status = "error",
        config = config
      ),
      config
    ))
  }
  idx <- which(cn == "X")
  converged <- if (is.null(fit$converged)) TRUE else isTRUE(fit$converged)
  sigma2 <- as.numeric(fit$sigma2)

  row <- data.frame(
    condition_id = as.integer(condition$condition_id),
    condition_set = as.character(condition$condition_set),
    rep_id = as.integer(rep_id),
    method = method,
    estimate = as.numeric(stats::coef(fit)[idx]),
    se = as.numeric(fit$se[idx]),
    ci_lower = as.numeric(fit$ci.lb[idx]),
    ci_upper = as.numeric(fit$ci.ub[idx]),
    p_value = as.numeric(fit$pval[idx]),
    df = NA_real_,
    converged = converged,
    status = "success",
    error = NA_character_,
    warning_count = as.integer(length(fit_capture$warnings)),
    warnings = collapse_warnings_s2(fit_capture$warnings, config$warning_separator),
    fit_runtime_sec = as.numeric(fit_capture$elapsed_sec),
    variance_study = sigma2[[1]],
    variance_effect = sigma2[[2]],
    variance_residual = NA_real_,
    boundary_study = NA,
    boundary_effect = NA,
    boundary_residual = NA,
    stringsAsFactors = FALSE
  )
  finalize_method_result_s2(row, config)
}

# =============================================================================
# CR2/Satterthwaite extraction
# =============================================================================

extract_cr2_s2 <- function(fit_capture,
                           data,
                           method_label,
                           condition,
                           rep_id,
                           config) {
  condition <- normalize_condition_s2(condition)

  if (!is.na(fit_capture$error) || is.null(fit_capture$value)) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method_label, condition, rep_id,
        error = fit_capture$error,
        warnings = fit_capture$warnings,
        elapsed_sec = fit_capture$elapsed_sec,
        status = "upstream_fit_failure",
        config = config
      ),
      config
    ))
  }

  verify_clubsandwich_api_s2()
  fit <- fit_capture$value
  sigma2 <- as.numeric(fit$sigma2)
  converged <- if (is.null(fit$converged)) TRUE else isTRUE(fit$converged)

  cap <- capture_conditions({
    clubSandwich::coef_test(
      fit,
      vcov = "CR2",
      cluster = data$study,
      test = "Satterthwaite",
      coefs = "X"
    )
  })
  elapsed <- as.numeric(fit_capture$elapsed_sec) + as.numeric(cap$elapsed_sec)
  warns <- c(fit_capture$warnings, cap$warnings)

  if (!is.na(cap$error) || is.null(cap$value)) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method_label, condition, rep_id,
        error = cap$error,
        warnings = warns,
        elapsed_sec = elapsed,
        status = "error",
        config = config
      ),
      config
    ))
  }

  ct_df <- as.data.frame(cap$value)
  assert_coef_test_columns_s2(ct_df)
  if (nrow(ct_df) != 1L) {
    # coefs="X" should return one row; fall back to name matching if needed.
    rn <- rownames(ct_df)
    if ("X" %in% rn) {
      ct_df <- ct_df["X", , drop = FALSE]
    } else if ("Coef" %in% names(ct_df)) {
      ct_df <- ct_df[ct_df$Coef == "X", , drop = FALSE]
    }
  }
  if (nrow(ct_df) != 1L) {
    return(finalize_method_result_s2(
      empty_method_result_s2(
        method_label, condition, rep_id,
        error = "Could not isolate coefficient 'X' from coef_test() output.",
        warnings = warns,
        elapsed_sec = elapsed,
        status = "error",
        config = config
      ),
      config
    ))
  }

  estimate <- as.numeric(ct_df$beta[[1]])
  se <- as.numeric(ct_df$SE[[1]])
  df <- as.numeric(ct_df$df_Satt[[1]])
  p_value <- as.numeric(ct_df$p_Satt[[1]])

  if (is.finite(df) && df > 0 && is.finite(se)) {
    critical <- stats::qt(1 - (1 - config$ci_level) / 2, df)
    ci_lower <- estimate - critical * se
    ci_upper <- estimate + critical * se
  } else {
    ci_lower <- NA_real_
    ci_upper <- NA_real_
  }

  row <- data.frame(
    condition_id = as.integer(condition$condition_id),
    condition_set = as.character(condition$condition_set),
    rep_id = as.integer(rep_id),
    method = as.character(method_label),
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    df = df,
    converged = converged,
    status = "success",
    error = NA_character_,
    warning_count = as.integer(length(warns)),
    warnings = collapse_warnings_s2(warns, config$warning_separator),
    fit_runtime_sec = elapsed,
    variance_study = sigma2[[1]],
    variance_effect = sigma2[[2]],
    variance_residual = NA_real_,
    boundary_study = NA,
    boundary_effect = NA,
    boundary_residual = NA,
    stringsAsFactors = FALSE
  )
  finalize_method_result_s2(row, config)
}

# =============================================================================
# CHE-RVE working model
# =============================================================================

# vcalc supplies the working covariance and rma.mv fits the nested model.
fit_che_base_s2 <- function(data, config, return_V = FALSE) {
  capture_conditions({
    dat <- data
    dat$study <- factor(dat$study)
    dat$es_id <- factor(dat$es_id)

    V_che <- metafor::vcalc(
      vi = vi,
      cluster = study,
      obs = es_id,
      rho = config$rho_assumed_che,
      data = dat,
      sparse = TRUE,
      checkpd = TRUE,
      nearpd = FALSE
    )

    fit <- metafor::rma.mv(
      yi = yi,
      V = V_che,
      mods = ~ X,
      random = ~ 1 | study / es_id,
      method = "REML",
      sparse = TRUE,
      data = dat
    )

    if (isTRUE(return_V)) {
      list(fit = fit, V_che = V_che)
    } else {
      fit
    }
  })
}

# =============================================================================
# Fit all Study 2 methods
# =============================================================================

# Aggregate HKSJ; one MLMA fit feeds MLMA/CR2; one CHE fit feeds CHE-RVE.
fit_all_methods_s2 <- function(data, condition, rep_id, config) {
  condition <- normalize_condition_s2(condition)
  rep_id <- as.integer(rep_id)

  row_hksj <- tryCatch(
    fit_aggregate_hksj_s2(data, condition, rep_id, config),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "aggregate_hksj", condition, rep_id,
          error = conditionMessage(e),
          status = "error",
          config = config
        ),
        config
      )
    }
  )

  mlma_cap <- tryCatch(
    fit_mlma_base_s2(data),
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
    extract_mlma_model_s2(mlma_cap, condition, rep_id, config),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "mlma", condition, rep_id,
          error = conditionMessage(e),
          warnings = mlma_cap$warnings,
          elapsed_sec = mlma_cap$elapsed_sec,
          status = "error",
          config = config
        ),
        config
      )
    }
  )

  row_mlma_cr2 <- tryCatch(
    extract_cr2_s2(
      fit_capture = mlma_cap,
      data = data,
      method_label = "mlma_cr2",
      condition = condition,
      rep_id = rep_id,
      config = config
    ),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "mlma_cr2", condition, rep_id,
          error = conditionMessage(e),
          warnings = mlma_cap$warnings,
          elapsed_sec = mlma_cap$elapsed_sec,
          status = "error",
          config = config
        ),
        config
      )
    }
  )

  che_cap <- tryCatch(
    fit_che_base_s2(data, config, return_V = FALSE),
    error = function(e) {
      list(
        value = NULL,
        error = conditionMessage(e),
        warnings = character(0),
        elapsed_sec = NA_real_
      )
    }
  )

  if (is.list(che_cap$value) && !inherits(che_cap$value, "rma") &&
      !is.null(che_cap$value$fit)) {
    che_cap$value <- che_cap$value$fit
  }

  row_che <- tryCatch(
    extract_cr2_s2(
      fit_capture = che_cap,
      data = data,
      method_label = "che_rve",
      condition = condition,
      rep_id = rep_id,
      config = config
    ),
    error = function(e) {
      finalize_method_result_s2(
        empty_method_result_s2(
          "che_rve", condition, rep_id,
          error = conditionMessage(e),
          warnings = che_cap$warnings,
          elapsed_sec = che_cap$elapsed_sec,
          status = "error",
          config = config
        ),
        config
      )
    }
  )

  out <- rbind(row_hksj, row_mlma, row_mlma_cr2, row_che)
  rownames(out) <- NULL
  if (nrow(out) != 4L) {
    stop("fit_all_methods_s2() must return exactly four rows.", call. = FALSE)
  }
  if (!identical(as.character(out$method),
                 c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve"))) {
    stop("fit_all_methods_s2() returned unexpected method labels.", call. = FALSE)
  }
  out
}
