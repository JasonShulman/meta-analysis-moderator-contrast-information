# Simulation Study 1 data generation, fitting, chunk execution, and summaries.

# Data-generating mechanism for Simulation Study 1
# Fisher z: y_ij = beta0 + beta1 * X_j + u_j + w_ij + e_ij
#
# Level 3: u_j  ~ N(0, tau^2)              between-study heterogeneity
# Level 2: w_ij ~ N(0, omega^2)            within-study true-effect heterogeneity
# Level 1: e_j  ~ MVN(0, Sigma_j)          sampling error with rho dependence

# Fixed k_j templates with sum_j k_j == 4 * J exactly.
make_kj_template <- function(J, imbalance = c("balanced", "moderate", "severe")) {
  imbalance <- match.arg(imbalance)
  if (!J %in% c(20L, 40L, 80L)) {
    stop("J must be 20, 40, or 80 for fixed k_j templates.", call. = FALSE)
  }

  kj <- switch(
    imbalance,
    balanced = rep(4L, J),
    moderate = {
      # Evenly distributed across {1, 3, 5, 7}; J/4 studies of each size.
      n_each <- J / 4L
      rep(c(1L, 3L, 5L, 7L), each = n_each)
    },
    severe = {
      if (J == 20L) {
        c(rep(1L, 10), rep(2L, 5), rep(9L, 3), 16L, 17L)
      } else if (J == 40L) {
        c(rep(1L, 20), rep(2L, 10), rep(9L, 6), rep(16L, 3), 18L)
      } else {
        # J == 80
        c(rep(1L, 40), rep(2L, 20), rep(9L, 12), rep(16L, 6), rep(18L, 2))
      }
    }
  )

  if (length(kj) != J) {
    stop("k_j template length does not equal J.", call. = FALSE)
  }
  if (sum(kj) != 4L * J) {
    stop("k_j template does not sum to 4J.", call. = FALSE)
  }

  as.integer(kj)
}

# Assign a binary study-level moderator with exact subgroup counts.
# X is assigned independently of k_j. The majority group receives X = 0.
assign_moderator <- function(J, split = c("50/50", "70/30", "85/15")) {
  split <- match.arg(split)

  n1 <- switch(
    split,
    "50/50" = as.integer(round(0.50 * J)),
    "70/30" = as.integer(round(0.30 * J)),
    "85/15" = as.integer(round(0.15 * J))
  )
  n0 <- J - n1

  expected <- switch(
    split,
    "50/50" = list(`20` = c(10L, 10L), `40` = c(20L, 20L), `80` = c(40L, 40L)),
    "70/30" = list(`20` = c(14L, 6L),  `40` = c(28L, 12L), `80` = c(56L, 24L)),
    "85/15" = list(`20` = c(17L, 3L),  `40` = c(34L, 6L),  `80` = c(68L, 12L))
  )
  key <- as.character(as.integer(J))
  if (!is.null(expected[[key]])) {
    # Compare as integers so numeric J=20 matches integer templates.
    if (!identical(as.integer(c(n0, n1)), as.integer(expected[[key]]))) {
      stop("Moderator counts do not match the split template.", call. = FALSE)
    }
  }

  sample(c(rep(0L, as.integer(n0)), rep(1L, as.integer(n1))))
}

# Build a within-study sampling-error covariance matrix.
# Diagonal: sampling variances v_i.
# Off-diagonal: rho * sqrt(v_i * v_i').
make_sampling_Sigma <- function(v, rho) {
  k <- length(v)
  if (k == 0L) {
    stop("v must have length at least 1.", call. = FALSE)
  }
  if (k == 1L) {
    return(matrix(v, nrow = 1L, ncol = 1L))
  }

  s <- sqrt(v)
  Sigma <- rho * tcrossprod(s)
  diag(Sigma) <- v
  Sigma
}

simulate_meta_dataset <- function(J,
                                  imbalance = c("balanced", "moderate", "severe"),
                                  split = c("50/50", "70/30", "85/15"),
                                  beta0 = 0.20,
                                  beta1 = 0,
                                  tau = 0.10,
                                  omega = 0.10,
                                  rho = 0.50,
                                  n_j = 100,
                                  shuffle_kj = TRUE,
                                  return_components = FALSE) {
  imbalance <- match.arg(imbalance)
  split <- match.arg(split)

  kj <- make_kj_template(J, imbalance)
  if (shuffle_kj) {
    kj <- sample(kj)
  }

  Xj <- assign_moderator(J, split)

  if (length(n_j) == 1L) {
    n_j <- rep(as.numeric(n_j), J)
  }
  if (length(n_j) != J) {
    stop("n_j must be a scalar or a vector of length J.", call. = FALSE)
  }
  if (any(n_j <= 3)) {
    stop("All n_j must be greater than 3 so that vi = 1/(n_j - 3) is defined.",
         call. = FALSE)
  }

  # Level 3: one between-study random effect per cluster.
  uj <- stats::rnorm(J, mean = 0, sd = tau)

  rows <- vector("list", J)
  for (j in seq_len(J)) {
    kj_j <- kj[j]
    vi_j <- rep(1 / (n_j[j] - 3), kj_j)

    # Level 2: within-study true-effect heterogeneity.
    wij <- stats::rnorm(kj_j, mean = 0, sd = omega)

    # Level 1: correlated sampling errors.
    Sigma_j <- make_sampling_Sigma(vi_j, rho)
    if (kj_j == 1L) {
      eij <- stats::rnorm(1L, mean = 0, sd = sqrt(vi_j))
    } else {
      eij <- as.numeric(MASS::mvrnorm(n = 1L, mu = rep(0, kj_j), Sigma = Sigma_j))
    }

    yij <- beta0 + beta1 * Xj[j] + uj[j] + wij + eij

    df_j <- data.frame(
      study = j,
      es_id = seq_len(kj_j),
      X = Xj[j],
      k_j = kj_j,
      n_j = n_j[j],
      vi = vi_j,
      yi = yij,
      stringsAsFactors = FALSE
    )
    if (return_components) {
      df_j$u_j <- uj[j]
      df_j$w_ij <- wij
      df_j$e_ij <- eij
    }
    rows[[j]] <- df_j
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Study 1 fits independent-effects RE, three-level MLMA, and MLMA with CR2.
# POMADE, CHE, and cluster wild bootstrap are not fit.

empty_method_row <- function(method, error = NA_character_) {
  data.frame(
    method = method,
    estimate = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    p_value = NA_real_,
    df = NA_real_,
    converged = FALSE,
    error = error,
    stringsAsFactors = FALSE
  )
}

extract_beta1_model_based <- function(fit, method) {
  cn <- names(stats::coef(fit))
  x_name <- if ("X" %in% cn) {
    "X"
  } else {
    # Fallback if metafor names the coefficient differently.
    grep("^X$|X$", cn, value = TRUE)[1]
  }
  if (is.na(x_name) || !length(x_name)) {
    stop("Could not locate the X coefficient in the fitted model.", call. = FALSE)
  }

  idx <- which(cn == x_name)
  est <- as.numeric(stats::coef(fit)[idx])
  se <- as.numeric(fit$se[idx])
  # Model-based Wald CIs are stored on the fitted object.
  ci_lower <- as.numeric(fit$ci.lb[idx])
  ci_upper <- as.numeric(fit$ci.ub[idx])
  p_val <- as.numeric(fit$pval[idx])
  # Some metafor fits may leave $converged NULL; treat successful fit as converged.
  converged <- if (is.null(fit$converged)) TRUE else isTRUE(fit$converged)

  data.frame(
    method = method,
    estimate = est,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_val,
    df = NA_real_,
    converged = converged,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# Intentionally treats every effect size as independent.
fit_naive_re <- function(data) {
  method <- "naive_re"
  tryCatch(
    {
      fit <- metafor::rma.uni(
        yi = yi,
        vi = vi,
        mods = ~ X,
        data = data,
        method = "REML"
      )
      extract_beta1_model_based(fit, method)
    },
    error = function(e) empty_method_row(method, conditionMessage(e))
  )
}

# Three-level MLMA: ~ 1 | study / es_id.
fit_mlma_model <- function(data) {
  dat <- data
  dat$study <- factor(dat$study)
  dat$es_id <- factor(dat$es_id)

  metafor::rma.mv(
    yi = yi,
    V = vi,
    mods = ~ X,
    random = ~ 1 | study / es_id,
    data = dat,
    method = "REML",
    sparse = TRUE
  )
}

# Three-level multilevel meta-analysis with model-based inference.
fit_mlma <- function(data) {
  method <- "mlma"
  tryCatch(
    {
      fit <- fit_mlma_model(data)
      extract_beta1_model_based(fit, method)
    },
    error = function(e) empty_method_row(method, conditionMessage(e))
  )
}

# Three-level MLMA with clubSandwich CR2 robust inference for X.
# Fits the same rma.mv() model as fit_mlma(), then applies CR2 cluster-robust
# SEs / Satterthwaite tests with cluster = study.
fit_mlma_cr2 <- function(data) {
  method <- "mlma_cr2"
  tryCatch(
    {
      fit <- fit_mlma_model(data)
      cluster <- data$study

      ct <- clubSandwich::coef_test(
        fit,
        vcov = "CR2",
        cluster = cluster,
        test = "Satterthwaite"
      )
      ct_df <- as.data.frame(ct)
      row_names <- rownames(ct_df)
      x_idx <- which(row_names == "X")
      if (length(x_idx) != 1L && "Coef" %in% names(ct_df)) {
        x_idx <- which(ct_df$Coef == "X")
      }
      if (length(x_idx) != 1L) {
        stop("Could not locate the X row in clubSandwich::coef_test() output.",
             call. = FALSE)
      }

      est <- ct_df$beta[x_idx]
      se <- ct_df$SE[x_idx]
      df_satt <- ct_df$df_Satt[x_idx]
      # Prefer p_Satt; fall back to any p-value column present.
      p_val <- if ("p_Satt" %in% names(ct_df)) {
        ct_df$p_Satt[x_idx]
      } else if ("p_val" %in% names(ct_df)) {
        ct_df$p_val[x_idx]
      } else {
        NA_real_
      }

      if (is.finite(df_satt) && df_satt > 0) {
        tcrit <- stats::qt(0.975, df = df_satt)
        ci_lower <- est - tcrit * se
        ci_upper <- est + tcrit * se
      } else {
        ci_lower <- NA_real_
        ci_upper <- NA_real_
      }

      converged <- if (is.null(fit$converged)) TRUE else isTRUE(fit$converged)

      data.frame(
        method = method,
        estimate = as.numeric(est),
        se = as.numeric(se),
        ci_lower = as.numeric(ci_lower),
        ci_upper = as.numeric(ci_upper),
        p_value = as.numeric(p_val),
        df = as.numeric(df_satt),
        converged = converged,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) empty_method_row(method, conditionMessage(e))
  )
}

fit_all_models <- function(data) {
  rbind(
    fit_naive_re(data),
    fit_mlma(data),
    fit_mlma_cr2(data)
  )
}

run_one_replication <- function(rep_id,
                                J = 20L,
                                imbalance = "balanced",
                                split = "50/50",
                                beta0 = 0.20,
                                beta1 = 0,
                                tau = 0.10,
                                omega = 0.10,
                                rho = 0.50,
                                n_j = 100,
                                seed_base = 20260000L) {
  rep_id <- as.integer(rep_id)

  tryCatch(
    {
      set.seed(as.integer(seed_base) + rep_id)

      dat <- simulate_meta_dataset(
        J = J,
        imbalance = imbalance,
        split = split,
        beta0 = beta0,
        beta1 = beta1,
        tau = tau,
        omega = omega,
        rho = rho,
        n_j = n_j,
        shuffle_kj = TRUE
      )

      fits <- fit_all_models(dat)

      data.frame(
        rep_id = rep_id,
        method = fits$method,
        beta1_true = beta1,
        estimate = fits$estimate,
        se = fits$se,
        ci_lower = fits$ci_lower,
        ci_upper = fits$ci_upper,
        p_value = fits$p_value,
        df = fits$df,
        converged = fits$converged,
        error = fits$error,
        status = ifelse(
          is.na(fits$error) & is.finite(fits$estimate) & !is.na(fits$converged) &
            fits$converged,
          "success",
          "fail"
        ),
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      # Entire replication failed (e.g., data generation); one fail row per method.
      methods <- c("naive_re", "mlma", "mlma_cr2")
      data.frame(
        rep_id = rep_id,
        method = methods,
        beta1_true = beta1,
        estimate = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        p_value = NA_real_,
        df = NA_real_,
        converged = FALSE,
        error = conditionMessage(e),
        status = "fail",
        stringsAsFactors = FALSE
      )
    }
  )
}

mcse_proportion <- function(p, n) {
  if (is.na(p) || is.na(n) || n <= 0) {
    return(NA_real_)
  }
  as.numeric(sqrt(p * (1 - p) / n))
}

# Operating characteristics by method.
# Successful fits are those with status == "success". Coverage and rejection
# rates are computed only over successful fits. Type I error is reported when
# beta1_true == 0; power when beta1_true != 0 (NA otherwise).
summarize_simulation_performance <- function(results, alpha = 0.05) {
  if (!nrow(results)) {
    stop("results is empty.", call. = FALSE)
  }

  beta1_true <- unique(results$beta1_true)
  if (length(beta1_true) != 1L) {
    stop("summarize_simulation_performance() expects a single beta1_true value.",
         call. = FALSE)
  }
  beta1_true <- beta1_true[[1]]

  methods <- unique(results$method)
  rows <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    m <- methods[[i]]
    d <- results[results$method == m, , drop = FALSE]
    n_attempted <- nrow(d)
    ok <- d$status == "success"
    n_success <- sum(ok)
    n_fail <- n_attempted - n_success
    error_rate <- n_fail / n_attempted

    d_ok <- d[ok, , drop = FALSE]

    if (n_success > 0L) {
      mean_est <- mean(d_ok$estimate)
      bias <- mean_est - beta1_true
      mean_se <- mean(d_ok$se)
      ci_width <- d_ok$ci_upper - d_ok$ci_lower
      mean_ci_width <- mean(ci_width)
      covered <- (d_ok$ci_lower <= beta1_true) & (d_ok$ci_upper >= beta1_true)
      coverage <- mean(covered)
      rejected <- d_ok$p_value < alpha
      reject_rate <- mean(rejected)
      if (all(is.na(d_ok$df))) {
        mean_df <- NA_real_
        median_df <- NA_real_
        min_df <- NA_real_
        max_df <- NA_real_
      } else {
        mean_df <- mean(d_ok$df, na.rm = TRUE)
        median_df <- stats::median(d_ok$df, na.rm = TRUE)
        min_df <- min(d_ok$df, na.rm = TRUE)
        max_df <- max(d_ok$df, na.rm = TRUE)
      }
      mcse_coverage <- mcse_proportion(coverage, n_success)
      mcse_reject <- mcse_proportion(reject_rate, n_success)
    } else {
      mean_est <- NA_real_
      bias <- NA_real_
      mean_se <- NA_real_
      mean_ci_width <- NA_real_
      coverage <- NA_real_
      reject_rate <- NA_real_
      mean_df <- NA_real_
      median_df <- NA_real_
      min_df <- NA_real_
      max_df <- NA_real_
      mcse_coverage <- NA_real_
      mcse_reject <- NA_real_
    }

    type_I_error <- if (isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    power <- if (!isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_

    rows[[i]] <- data.frame(
      method = m,
      beta1_true = beta1_true,
      n_attempted = n_attempted,
      n_success = n_success,
      n_fail = n_fail,
      error_rate = error_rate,
      mean_estimate = mean_est,
      bias = bias,
      mean_se = mean_se,
      mean_ci_width = mean_ci_width,
      coverage = coverage,
      mcse_coverage = mcse_coverage,
      type_I_error = type_I_error,
      power = power,
      mcse_reject_rate = mcse_reject,
      mean_df = mean_df,
      median_df = median_df,
      min_df = min_df,
      max_df = max_df,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

make_study1_condition_grid <- function(J = c(20L, 40L),
                                      imbalance = c("balanced", "severe"),
                                      split = c("50/50", "85/15"),
                                      beta1 = c(0, 0.30),
                                      beta0 = 0.20,
                                      rho = 0.50,
                                      tau = 0.10,
                                      omega = 0.10,
                                      n_j = 100,
                                      n_rep = 250L) {
  grid <- expand.grid(
    J = as.integer(J),
    imbalance = imbalance,
    split = split,
    beta1 = beta1,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$J, grid$imbalance, grid$split, grid$beta1), , drop = FALSE]
  rownames(grid) <- NULL

  grid$condition_id <- seq_len(nrow(grid))
  grid$beta0 <- beta0
  grid$rho <- rho
  grid$tau <- tau
  grid$omega <- omega
  grid$n_j <- n_j
  grid$n_rep <- as.integer(n_rep)

  grid[, c(
    "condition_id", "J", "imbalance", "split", "beta1",
    "beta0", "rho", "tau", "omega", "n_j", "n_rep"
  )]
}

summarize_study1_performance <- function(results, alpha = 0.05) {
  ids <- sort(unique(results$condition_id))
  pieces <- vector("list", length(ids))

  for (i in seq_along(ids)) {
    id <- ids[[i]]
    d <- results[results$condition_id == id, , drop = FALSE]
    perf <- summarize_simulation_performance(d, alpha = alpha)
    design <- d[1, c("condition_id", "J", "imbalance", "split"), drop = FALSE]
    pieces[[i]] <- cbind(
      design[rep(1L, nrow(perf)), , drop = FALSE],
      perf,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL

  out[, c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "n_attempted", "n_success", "n_fail", "error_rate",
    "mean_estimate", "bias", "mean_se", "mean_ci_width",
    "coverage", "mcse_coverage",
    "type_I_error", "power", "mcse_reject_rate",
    "mean_df", "median_df", "min_df", "max_df"
  )]
}

# Deterministic seed_base for a condition (chunk-invariant).
# set.seed(seed_base + rep_id) inside run_one_replication(), so the same
# (condition_id, rep_id) always draws the same data regardless of chunk_size.
# chunk_id is stored on the job for checkpointing/logging only.
condition_seed_base <- function(condition_id, seed_offset = 202700000L) {
  as.integer(seed_offset + as.integer(condition_id) * 100000L)
}

make_chunk_grid <- function(condition_grid,
                            n_rep_total = 100L,
                            chunk_size = 25L,
                            seed_offset = 202700000L) {
  n_rep_total <- as.integer(n_rep_total)
  chunk_size <- as.integer(chunk_size)
  if (n_rep_total < 1L || chunk_size < 1L) {
    stop("n_rep_total and chunk_size must be >= 1.", call. = FALSE)
  }
  if (n_rep_total %% chunk_size != 0L) {
    stop("n_rep_total must be divisible by chunk_size.", call. = FALSE)
  }

  n_chunks <- n_rep_total / chunk_size
  rows <- vector("list", nrow(condition_grid) * n_chunks)
  k <- 0L

  for (i in seq_len(nrow(condition_grid))) {
    cond <- condition_grid[i, , drop = FALSE]
    for (ch in seq_len(n_chunks)) {
      k <- k + 1L
      rep_start <- as.integer((ch - 1L) * chunk_size + 1L)
      rep_end <- as.integer(ch * chunk_size)
      rows[[k]] <- data.frame(
        condition_id = cond$condition_id,
        chunk_id = as.integer(ch),
        J = cond$J,
        imbalance = cond$imbalance,
        split = cond$split,
        beta1_true = cond$beta1,
        beta0 = cond$beta0,
        rho = cond$rho,
        tau = cond$tau,
        omega = cond$omega,
        n_j = cond$n_j,
        rep_start = rep_start,
        rep_end = rep_end,
        n_rep_chunk = as.integer(chunk_size),
        seed_base = condition_seed_base(cond$condition_id, seed_offset),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

chunk_rds_path <- function(out_dir, condition_id, chunk_id) {
  file.path(
    out_dir,
    "raw_chunks",
    sprintf(
      "condition_%02d_chunk_%02d.rds",
      as.integer(condition_id),
      as.integer(chunk_id)
    )
  )
}

source_simulation_functions <- function(root) {
  source(file.path(root, "R", "study1_functions.R"), local = FALSE)
}

# Run one chunk, reusing an existing checkpoint when allowed.
run_one_chunk <- function(chunk, root, out_dir, overwrite = FALSE) {
  chunk <- as.list(chunk)
  out_file <- chunk_rds_path(out_dir, chunk$condition_id, chunk$chunk_id)

  # Checkpoint: reuse existing chunk output.
  if (!isTRUE(overwrite) && file.exists(out_file)) {
    raw <- tryCatch(
      readRDS(out_file),
      error = function(e) NULL
    )
    file_mb <- tryCatch(
      file.info(out_file)$size / (1024^2),
      error = function(e) NA_real_
    )
    n_success <- if (!is.null(raw)) sum(raw$status == "success", na.rm = TRUE) else NA_integer_
    n_fail <- if (!is.null(raw)) sum(raw$status != "success", na.rm = TRUE) else NA_integer_
    n_tot <- if (is.finite(n_success) && is.finite(n_fail)) n_success + n_fail else NA_real_

    log_row <- data.frame(
      condition_id = chunk$condition_id,
      chunk_id = chunk$chunk_id,
      J = chunk$J,
      imbalance = chunk$imbalance,
      split = chunk$split,
      beta1_true = chunk$beta1_true,
      n_rep_chunk = chunk$n_rep_chunk,
      start_time = NA_character_,
      end_time = NA_character_,
      elapsed_sec = NA_real_,
      sec_per_rep = NA_real_,
      n_success = as.integer(n_success),
      n_fail = as.integer(n_fail),
      error_rate = if (is.finite(n_tot) && n_tot > 0) n_fail / n_tot else NA_real_,
      output_file = out_file,
      output_file_size_mb = as.numeric(file_mb),
      status = "skipped",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    return(list(status = "skipped", log = log_row, raw = raw))
  }

  start_time <- Sys.time()
  err_msg <- NA_character_

  raw <- tryCatch(
    {
      # Fresh workers need the DGM / fitting / performance functions.
      source_simulation_functions(root)

      reps <- seq.int(chunk$rep_start, chunk$rep_end)
      rows <- vector("list", length(reps))
      for (j in seq_along(reps)) {
        r <- reps[[j]]
        rows[[j]] <- run_one_replication(
          rep_id = r,
          J = chunk$J,
          imbalance = chunk$imbalance,
          split = chunk$split,
          beta0 = chunk$beta0,
          beta1 = chunk$beta1_true,
          tau = chunk$tau,
          omega = chunk$omega,
          rho = chunk$rho,
          n_j = chunk$n_j,
          seed_base = chunk$seed_base
        )
      }
      out <- do.call(rbind, rows)
      out$condition_id <- chunk$condition_id
      out$chunk_id <- chunk$chunk_id
      out$J <- chunk$J
      out$imbalance <- chunk$imbalance
      out$split <- chunk$split
      out$beta0 <- chunk$beta0
      out$rho <- chunk$rho
      out$tau <- chunk$tau
      out$omega <- chunk$omega
      out$n_j <- chunk$n_j
      rownames(out) <- NULL
      out
    },
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )

  end_time <- Sys.time()
  elapsed_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (is.null(raw)) {
    # Record a failed chunk placeholder so combine/summary can proceed.
    methods <- c("naive_re", "mlma", "mlma_cr2")
    expand <- expand.grid(
      rep_id = seq.int(chunk$rep_start, chunk$rep_end),
      method = methods,
      stringsAsFactors = FALSE
    )
    raw <- data.frame(
      condition_id = chunk$condition_id,
      chunk_id = chunk$chunk_id,
      J = chunk$J,
      imbalance = chunk$imbalance,
      split = chunk$split,
      beta0 = chunk$beta0,
      beta1_true = chunk$beta1_true,
      rho = chunk$rho,
      tau = chunk$tau,
      omega = chunk$omega,
      n_j = chunk$n_j,
      rep_id = expand$rep_id,
      method = expand$method,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      converged = FALSE,
      error = err_msg,
      status = "fail",
      stringsAsFactors = FALSE
    )
    chunk_status <- "failed"
  } else {
    chunk_status <- "completed"
  }

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  save_ok <- tryCatch(
    {
      saveRDS(raw, out_file)
      TRUE
    },
    error = function(e) {
      detail <- sprintf(
        "saveRDS failed for output_file='%s': %s",
        out_file,
        conditionMessage(e)
      )
      err_msg <<- paste(na.omit(c(err_msg, detail)), collapse = "; ")
      FALSE
    }
  )
  if (!save_ok && identical(chunk_status, "completed")) {
    chunk_status <- "failed"
  }
  if (!save_ok && identical(chunk_status, "failed") && is.na(err_msg)) {
    err_msg <- sprintf("saveRDS failed for output_file='%s'", out_file)
  }

  file_mb <- if (file.exists(out_file)) {
    file.info(out_file)$size / (1024^2)
  } else {
    NA_real_
  }

  n_success <- sum(raw$status == "success", na.rm = TRUE)
  n_fail <- sum(raw$status != "success", na.rm = TRUE)
  n_tot <- n_success + n_fail

  log_row <- data.frame(
    condition_id = chunk$condition_id,
    chunk_id = chunk$chunk_id,
    J = chunk$J,
    imbalance = chunk$imbalance,
    split = chunk$split,
    beta1_true = chunk$beta1_true,
    n_rep_chunk = chunk$n_rep_chunk,
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = elapsed_sec,
    sec_per_rep = elapsed_sec / as.numeric(chunk$n_rep_chunk),
    n_success = as.integer(n_success),
    n_fail = as.integer(n_fail),
    error_rate = if (n_tot > 0) n_fail / n_tot else NA_real_,
    output_file = out_file,
    output_file_size_mb = as.numeric(file_mb),
    status = chunk_status,
    error_message = err_msg,
    stringsAsFactors = FALSE
  )

  list(status = chunk_status, log = log_row, raw = raw)
}

run_chunks_parallel <- function(chunk_grid,
                                root,
                                out_dir,
                                workers = 4L,
                                overwrite = FALSE) {
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Packages future and future.apply are required.", call. = FALSE)
  }

  dir.create(file.path(out_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = as.integer(workers))

  message(sprintf(
    "Running %d chunks with future::multisession (workers = %d)...",
    nrow(chunk_grid), as.integer(workers)
  ))

  wall_start <- Sys.time()
  results <- future.apply::future_lapply(
    X = seq_len(nrow(chunk_grid)),
    FUN = function(i) {
      run_one_chunk(
        chunk = chunk_grid[i, , drop = FALSE],
        root = root,
        out_dir = out_dir,
        overwrite = overwrite
      )
    },
    # Seeds are set inside run_one_replication(). Do not add a parallel RNG stream.
    future.seed = NULL,
    future.packages = c("metafor", "clubSandwich", "MASS")
  )
  wall_end <- Sys.time()
  wall_elapsed <- as.numeric(difftime(wall_end, wall_start, units = "secs"))

  chunk_log <- do.call(rbind, lapply(results, function(x) x$log))
  rownames(chunk_log) <- NULL

  list(
    chunk_log = chunk_log,
    wall_elapsed_sec = wall_elapsed,
    wall_start = wall_start,
    wall_end = wall_end,
    results = results
  )
}

combine_chunk_rds <- function(out_dir, chunk_grid = NULL) {
  chunk_dir <- file.path(out_dir, "raw_chunks")
  files <- if (is.null(chunk_grid)) {
    list.files(chunk_dir, pattern = "^condition_.*_chunk_.*\\.rds$", full.names = TRUE)
  } else {
    vapply(
      seq_len(nrow(chunk_grid)),
      function(i) {
        chunk_rds_path(out_dir, chunk_grid$condition_id[i], chunk_grid$chunk_id[i])
      },
      character(1)
    )
  }
  files <- files[file.exists(files)]
  if (!length(files)) {
    stop("No chunk RDS files found to combine.", call. = FALSE)
  }

  pieces <- lapply(files, function(f) {
    tryCatch(
      readRDS(f),
      error = function(e) {
        warning("Could not read ", f, ": ", conditionMessage(e))
        NULL
      }
    )
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

aggregate_condition_runtime <- function(chunk_log) {
  ids <- sort(unique(chunk_log$condition_id))
  rows <- lapply(ids, function(id) {
    d <- chunk_log[chunk_log$condition_id == id, , drop = FALSE]
    timed <- d[!is.na(d$elapsed_sec), , drop = FALSE]
    data.frame(
      condition_id = id,
      J = d$J[1],
      imbalance = d$imbalance[1],
      split = d$split[1],
      beta1_true = d$beta1_true[1],
      n_chunks = nrow(d),
      n_chunks_completed = sum(d$status == "completed"),
      n_chunks_skipped = sum(d$status == "skipped"),
      n_chunks_failed = sum(d$status == "failed"),
      n_rep_total = sum(d$n_rep_chunk, na.rm = TRUE),
      total_chunk_elapsed_sec = sum(timed$elapsed_sec, na.rm = TRUE),
      mean_sec_per_rep = {
        if (nrow(timed) && sum(timed$n_rep_chunk, na.rm = TRUE) > 0) {
          sum(timed$elapsed_sec, na.rm = TRUE) / sum(timed$n_rep_chunk, na.rm = TRUE)
        } else {
          NA_real_
        }
      },
      n_success = sum(d$n_success, na.rm = TRUE),
      n_fail = sum(d$n_fail, na.rm = TRUE),
      error_rate = {
        ns <- sum(d$n_success, na.rm = TRUE)
        nf <- sum(d$n_fail, na.rm = TRUE)
        if ((ns + nf) > 0) nf / (ns + nf) else NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ensure_sim_study_dirs <- function(out_dir) {
  subdirs <- c(
    "raw_chunks", "figure_data", "tables", "figures", "logs"
  )
  for (sd in subdirs) {
    dir.create(file.path(out_dir, sd), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(out_dir)
}

write_sim_study_tables <- function(summary_tab,
                                   condition_grid,
                                   condition_log,
                                   tables_dir,
                                   figure_data_dir = NULL) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(figure_data_dir)) {
    dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (nrow(summary_tab)) {
    summary_tab$J <- factor(summary_tab$J, levels = c(20, 40, 80))
    summary_tab$imbalance <- factor(
      summary_tab$imbalance,
      levels = c("balanced", "moderate", "severe")
    )
    summary_tab$split <- factor(
      summary_tab$split,
      levels = c("50/50", "70/30", "85/15")
    )
    summary_tab$method <- factor(
      summary_tab$method,
      levels = c("naive_re", "mlma", "mlma_cr2")
    )
  }

  design_grid <- condition_grid
  utils::write.csv(design_grid, file.path(tables_dir, "design_grid.csv"),
                   row.names = FALSE)

  type1 <- summary_tab[as.numeric(as.character(summary_tab$beta1_true)) == 0, , drop = FALSE]
  type1_out <- type1[, intersect(c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "type_I_error", "mcse_reject_rate", "coverage", "mcse_coverage",
    "mean_df", "median_df", "min_df", "max_df", "n_success", "error_rate"
  ), names(type1)), drop = FALSE]
  utils::write.csv(type1_out, file.path(tables_dir, "table_type1_error.csv"),
                   row.names = FALSE)

  p10 <- summary_tab[abs(as.numeric(as.character(summary_tab$beta1_true)) - 0.10) < 1e-8, , drop = FALSE]
  p10_out <- p10[, intersect(c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "power", "mcse_reject_rate", "coverage", "mcse_coverage",
    "mean_df", "median_df", "min_df", "max_df", "n_success", "error_rate"
  ), names(p10)), drop = FALSE]
  utils::write.csv(p10_out, file.path(tables_dir, "table_power_beta10.csv"),
                   row.names = FALSE)

  p30 <- summary_tab[abs(as.numeric(as.character(summary_tab$beta1_true)) - 0.30) < 1e-8, , drop = FALSE]
  p30_out <- p30[, intersect(c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "power", "mcse_reject_rate", "coverage", "mcse_coverage",
    "mean_df", "median_df", "min_df", "max_df", "n_success", "error_rate"
  ), names(p30)), drop = FALSE]
  utils::write.csv(p30_out, file.path(tables_dir, "table_power_beta30.csv"),
                   row.names = FALSE)

  cov_out <- summary_tab[, intersect(c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "coverage", "mcse_coverage", "mean_ci_width", "n_success"
  ), names(summary_tab)), drop = FALSE]
  utils::write.csv(cov_out, file.path(tables_dir, "table_coverage.csv"),
                   row.names = FALSE)

  cr2 <- summary_tab[as.character(summary_tab$method) == "mlma_cr2", , drop = FALSE]
  cr2_out <- cr2[, intersect(c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "mean_df", "median_df", "min_df", "max_df",
    "power", "type_I_error", "coverage", "n_success"
  ), names(cr2)), drop = FALSE]
  utils::write.csv(cr2_out, file.path(tables_dir, "table_cr2_df.csv"),
                   row.names = FALSE)

  utils::write.csv(
    condition_log,
    file.path(tables_dir, "table_runtime.csv"),
    row.names = FALSE
  )

  cmp_cols <- c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "coverage", "type_I_error", "power", "mean_ci_width", "mean_df",
    "bias", "error_rate"
  )
  cmp <- summary_tab[, intersect(cmp_cols, names(summary_tab)), drop = FALSE]
  utils::write.csv(cmp, file.path(tables_dir, "table_method_comparison.csv"),
                   row.names = FALSE)

  tables <- list(
    design_grid = design_grid,
    type1 = type1_out,
    power_beta10 = p10_out,
    power_beta30 = p30_out,
    coverage = cov_out,
    cr2_df = cr2_out,
    runtime = condition_log,
    method_comparison = cmp
  )

  if (!is.null(figure_data_dir)) {
    for (nm in names(tables)) {
      utils::write.csv(
        tables[[nm]],
        file.path(figure_data_dir, paste0(nm, ".csv")),
        row.names = FALSE
      )
    }
    utils::write.csv(
      summary_tab,
      file.path(figure_data_dir, "full_summary.csv"),
      row.names = FALSE
    )
  }

  tables
}

ALPHA <- 0.05
J_LEVELS <- c(20L, 40L, 80L)
IMBALANCE_LEVELS <- c("balanced", "moderate", "severe")
SPLIT_LEVELS <- c("50/50", "70/30", "85/15")
METHOD_LEVELS <- c("naive_re", "mlma", "mlma_cr2")

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

mc_interval <- function(est, mcse) {
  lower <- clamp01(est - 1.96 * mcse)
  upper <- clamp01(est + 1.96 * mcse)
  list(lower = lower, upper = upper)
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7))
}

apply_factor_order <- function(df) {
  df$J <- factor(as.integer(df$J), levels = J_LEVELS)
  df$imbalance <- factor(as.character(df$imbalance), levels = IMBALANCE_LEVELS)
  df$split <- factor(as.character(df$split), levels = SPLIT_LEVELS)
  df$beta1_true <- as.numeric(df$beta1_true)
  df$method <- factor(as.character(df$method), levels = METHOD_LEVELS)
  df
}

order_summary_rows <- function(df) {
  df <- apply_factor_order(df)
  df[order(df$condition_id, df$method), , drop = FALSE]
}

summarize_condition_methods <- function(raw, alpha = ALPHA) {
  required <- c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "rep_id", "estimate", "se", "ci_lower", "ci_upper", "p_value", "df",
    "status", "converged", "error"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop(
      "Raw chunk missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cond_ids <- unique(raw$condition_id)
  if (length(cond_ids) != 1L) {
    stop("summarize_condition_methods() expects a single condition_id.", call. = FALSE)
  }
  beta1_vals <- unique(raw$beta1_true)
  if (length(beta1_vals) != 1L) {
    stop("summarize_condition_methods() expects a single beta1_true.", call. = FALSE)
  }
  true_beta1 <- as.numeric(beta1_vals[[1L]])

  methods <- METHOD_LEVELS
  rows <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    m <- methods[[i]]
    d <- raw[raw$method == m, , drop = FALSE]
    n_attempted <- nrow(d)
    ok <- !is.na(d$status) & d$status == "success"
    n_success <- sum(ok)
    n_failed <- n_attempted - n_success
    failure_rate <- if (n_attempted > 0L) n_failed / n_attempted else NA_real_

    d_ok <- d[ok, , drop = FALSE]

    if (n_success > 0L) {
      est <- d_ok$estimate
      mean_estimate <- mean(est)
      bias <- mean(est - true_beta1)
      absolute_bias <- abs(bias)
      empirical_sd <- if (n_success >= 2L) stats::sd(est) else NA_real_
      rmse <- sqrt(mean((est - true_beta1)^2))
      mean_model_se <- mean(d_ok$se)
      median_model_se <- stats::median(d_ok$se)
      se_calibration_ratio <- if (is.finite(empirical_sd) && empirical_sd > 0) {
        mean_model_se / empirical_sd
      } else {
        NA_real_
      }

      covered <- (d_ok$ci_lower <= true_beta1) & (d_ok$ci_upper >= true_beta1)
      coverage <- mean(covered)
      ci_width <- d_ok$ci_upper - d_ok$ci_lower
      mean_ci_width <- mean(ci_width)
      median_ci_width <- stats::median(ci_width)
      ci_q <- safe_quantile(ci_width, c(0.05, 0.25, 0.75, 0.95))

      rejected <- d_ok$p_value < alpha
      rejection_rate <- mean(rejected)

      df_ok <- d_ok$df
      if (all(!is.finite(df_ok))) {
        mean_df <- median_df <- min_df <- max_df <- NA_real_
        df_p05 <- df_p25 <- df_p75 <- df_p95 <- NA_real_
      } else {
        mean_df <- mean(df_ok, na.rm = TRUE)
        median_df <- stats::median(df_ok, na.rm = TRUE)
        min_df <- min(df_ok, na.rm = TRUE)
        max_df <- max(df_ok, na.rm = TRUE)
        df_q <- safe_quantile(df_ok, c(0.05, 0.25, 0.75, 0.95))
        df_p05 <- df_q[1]
        df_p25 <- df_q[2]
        df_p75 <- df_q[3]
        df_p95 <- df_q[4]
      }

      mcse_bias <- if (is.finite(empirical_sd) && n_success > 0L) {
        empirical_sd / sqrt(n_success)
      } else {
        NA_real_
      }
      mcse_empirical_sd <- if (is.finite(empirical_sd) && n_success > 1L) {
        empirical_sd / sqrt(2 * (n_success - 1))
      } else {
        NA_real_
      }
      mcse_coverage <- sqrt(coverage * (1 - coverage) / n_success)
      mcse_rejection_rate <- sqrt(rejection_rate * (1 - rejection_rate) / n_success)
      cov_ci <- mc_interval(coverage, mcse_coverage)
      rej_ci <- mc_interval(rejection_rate, mcse_rejection_rate)
    } else {
      mean_estimate <- bias <- absolute_bias <- empirical_sd <- rmse <- NA_real_
      mean_model_se <- median_model_se <- se_calibration_ratio <- NA_real_
      coverage <- mean_ci_width <- median_ci_width <- NA_real_
      ci_q <- rep(NA_real_, 4L)
      rejection_rate <- NA_real_
      mean_df <- median_df <- min_df <- max_df <- NA_real_
      df_p05 <- df_p25 <- df_p75 <- df_p95 <- NA_real_
      mcse_bias <- mcse_empirical_sd <- mcse_coverage <- mcse_rejection_rate <- NA_real_
      cov_ci <- list(lower = NA_real_, upper = NA_real_)
      rej_ci <- list(lower = NA_real_, upper = NA_real_)
    }

    relative_bias <- if (isTRUE(all.equal(true_beta1, 0))) {
      NA_real_
    } else if (is.finite(bias)) {
      bias / true_beta1
    } else {
      NA_real_
    }

    reject_interp <- if (isTRUE(all.equal(true_beta1, 0))) {
      "type1_error"
    } else {
      "power"
    }

    type_I_error <- if (identical(reject_interp, "type1_error")) rejection_rate else NA_real_
    power <- if (identical(reject_interp, "power")) rejection_rate else NA_real_

    rows[[i]] <- data.frame(
      condition_id = as.integer(cond_ids[[1L]]),
      J = as.integer(raw$J[1L]),
      imbalance = as.character(raw$imbalance[1L]),
      split = as.character(raw$split[1L]),
      beta1_true = true_beta1,
      true_beta1 = true_beta1,
      method = m,
      n_attempted = as.integer(n_attempted),
      n_success = as.integer(n_success),
      n_fail = as.integer(n_failed),
      n_failed = as.integer(n_failed),
      error_rate = failure_rate,
      failure_rate = failure_rate,
      n_with_warning = NA_integer_,
      warning_rate = NA_real_,
      mean_estimate = mean_estimate,
      bias = bias,
      absolute_bias = absolute_bias,
      relative_bias = relative_bias,
      empirical_sd = empirical_sd,
      rmse = rmse,
      mean_se = mean_model_se,
      mean_model_se = mean_model_se,
      median_model_se = median_model_se,
      se_calibration_ratio = se_calibration_ratio,
      coverage = coverage,
      mean_ci_width = mean_ci_width,
      median_ci_width = median_ci_width,
      ci_width_p05 = ci_q[1],
      ci_width_p25 = ci_q[2],
      ci_width_p75 = ci_q[3],
      ci_width_p95 = ci_q[4],
      rejection_rate = rejection_rate,
      rejection_interpretation = reject_interp,
      type_I_error = type_I_error,
      power = power,
      mean_df = mean_df,
      median_df = median_df,
      min_df = min_df,
      df_p05 = df_p05,
      df_p25 = df_p25,
      df_p75 = df_p75,
      df_p95 = df_p95,
      max_df = max_df,
      mcse_bias = mcse_bias,
      mcse_empirical_sd = mcse_empirical_sd,
      mcse_coverage = mcse_coverage,
      mcse_rejection_rate = mcse_rejection_rate,
      mcse_reject_rate = mcse_rejection_rate,
      coverage_mc_lower = cov_ci$lower,
      coverage_mc_upper = cov_ci$upper,
      rejection_mc_lower = rej_ci$lower,
      rejection_mc_upper = rej_ci$upper,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

extract_failures <- function(raw, seed_lookup) {
  fails <- raw[!is.na(raw$status) & raw$status != "success", , drop = FALSE]
  if (!nrow(fails)) {
    return(fails[0, , drop = FALSE])
  }

  mlma_fail_reps <- unique(fails$rep_id[fails$method == "mlma"])
  cr2_fail_reps <- unique(fails$rep_id[fails$method == "mlma_cr2"])
  both_fail <- intersect(mlma_fail_reps, cr2_fail_reps)

  seed_base <- seed_lookup[as.character(fails$condition_id[1L])]
  seed_vals <- if (is.finite(seed_base)) {
    as.integer(seed_base) + as.integer(fails$rep_id)
  } else {
    rep(NA_integer_, nrow(fails))
  }

  data.frame(
    condition_id = as.integer(fails$condition_id),
    replication_id = as.integer(fails$rep_id),
    rep_id = as.integer(fails$rep_id),
    seed = seed_vals,
    seed_base = if (is.finite(seed_base)) as.integer(seed_base) else NA_integer_,
    chunk_id = if ("chunk_id" %in% names(fails)) as.integer(fails$chunk_id) else NA_integer_,
    chunk_filename = if ("chunk_id" %in% names(fails)) {
      sprintf(
        "condition_%02d_chunk_%02d.rds",
        as.integer(fails$condition_id),
        as.integer(fails$chunk_id)
      )
    } else {
      NA_character_
    },
    J = as.integer(fails$J),
    imbalance = as.character(fails$imbalance),
    split = as.character(fails$split),
    true_beta1 = as.numeric(fails$beta1_true),
    beta1_true = as.numeric(fails$beta1_true),
    method = as.character(fails$method),
    error_class = NA_character_,
    error_message = if ("error" %in% names(fails)) as.character(fails$error) else NA_character_,
    warning_message = NA_character_,
    converged = if ("converged" %in% names(fails)) as.logical(fails$converged) else NA,
    status = as.character(fails$status),
    cr2_downstream_of_mlma_failure = ifelse(
      fails$method == "mlma_cr2" & fails$rep_id %in% both_fail,
      TRUE,
      ifelse(fails$method == "mlma_cr2", FALSE, NA)
    ),
    note = "warning_message unavailable: production raw output did not retain warnings.",
    stringsAsFactors = FALSE
  )
}

# Extended summaries use a named vector of condition seed bases.
write_study1_extended_summaries <- function(raw, out_dir, seed_lookup) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ids <- sort(unique(raw$condition_id))
  summary_pieces <- lapply(ids, function(id) {
    summarize_condition_methods(raw[raw$condition_id == id, , drop = FALSE])
  })
  extended <- do.call(rbind, summary_pieces)
  rownames(extended) <- NULL
  extended <- order_summary_rows(extended)

  failure_pieces <- lapply(ids, function(id) {
    extract_failures(raw[raw$condition_id == id, , drop = FALSE], seed_lookup)
  })
  nonempty <- failure_pieces[vapply(failure_pieces, nrow, integer(1)) > 0L]
  if (!length(nonempty)) {
    failures <- data.frame(
      condition_id = integer(0),
      replication_id = integer(0),
      rep_id = integer(0),
      seed = integer(0),
      seed_base = integer(0),
      chunk_id = integer(0),
      chunk_filename = character(0),
      J = integer(0),
      imbalance = character(0),
      split = character(0),
      true_beta1 = numeric(0),
      beta1_true = numeric(0),
      method = character(0),
      error_class = character(0),
      error_message = character(0),
      warning_message = character(0),
      converged = logical(0),
      status = character(0),
      cr2_downstream_of_mlma_failure = logical(0),
      note = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    failures <- do.call(rbind, nonempty)
    rownames(failures) <- NULL
  }

  for (nm in c("J", "imbalance", "split", "method")) {
    extended[[nm]] <- if (nm == "J") {
      as.integer(as.character(extended[[nm]]))
    } else {
      as.character(extended[[nm]])
    }
  }

  table_se <- extended[, c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "empirical_sd", "mean_model_se", "se_calibration_ratio", "rmse", "n_success"
  ), drop = FALSE]
  table_mc <- extended[, c(
    "condition_id", "J", "imbalance", "split", "beta1_true", "method",
    "coverage", "mcse_coverage", "coverage_mc_lower", "coverage_mc_upper",
    "rejection_rate", "rejection_interpretation",
    "mcse_rejection_rate", "rejection_mc_lower", "rejection_mc_upper",
    "n_success"
  ), drop = FALSE]

  utils::write.csv(
    extended,
    file.path(out_dir, "sim_study_1_summary_extended.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    failures,
    file.path(out_dir, "sim_study_1_failure_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    table_se,
    file.path(out_dir, "table_se_calibration.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    table_mc,
    file.path(out_dir, "table_monte_carlo_uncertainty.csv"),
    row.names = FALSE
  )
  invisible(list(
    extended = extended,
    failures = failures,
    table_se = table_se,
    table_mc = table_mc
  ))
}
