# Simulation Study 2 design constants, seeds, templates, and condition grids.

sim2_project_root <- function(root = NULL) {
  if (!is.null(root)) {
    return(normalizePath(root, winslash = "/", mustWork = TRUE))
  }
  candidates <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "../..")
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "study2_design.R"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot locate reproducibility root (expected R/study2_design.R).", call. = FALSE)
}

sim2_fixture_dir <- function(root = NULL) {
  file.path(sim2_project_root(root), "data", "study2_design")
}

read_sim2_fixture <- function(filename, root = NULL) {
  path <- file.path(sim2_fixture_dir(root), filename)
  if (!file.exists(path)) {
    stop("Missing Study 2 fixture: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

make_sim2_config <- function(run_mode = c("smoke", "pilot", "production"),
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
      n_rep_total = 100L,
      chunk_size = 25L,
      workers = min(10L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
    ),
    production = list(
      n_rep_total = 5000L,
      chunk_size = 500L,
      workers = NA_integer_
    )
  )

  workers_out <- if (is.null(workers)) mode_defaults$workers else as.integer(workers)

  list(
    study_id = 2L,
    run_mode = run_mode,
    beta0 = 0.20,
    beta1_levels = c(0.00, 0.10),
    rho_true = 0.50,
    rho_assumed_aggregation = 0.50,
    rho_assumed_che = 0.50,
    tau = 0.10,
    omega = 0.10,
    alpha = 0.05,
    ci_level = 0.95,
    J_levels = c(20L, 40L, 80L),
    k_imbalance_levels = c("balanced", "moderate", "severe"),
    split_levels = c("50_50", "70_30", "85_15"),
    n_imbalance_levels = c("balanced", "moderate", "severe", "extreme"),
    allocation_patterns_core = "independent",
    allocation_patterns_alignment = c(
      "n_majority", "n_minority", "nk_majority", "nk_minority"
    ),
    mean_n = 300L,
    n_min = 20L,
    n_max = 2000L,
    sigma_log_n = c(
      balanced = 0.0,
      moderate = 0.5,
      severe = 0.9,
      extreme = 1.4
    ),
    n_rep_total = as.integer(mode_defaults$n_rep_total),
    chunk_size = as.integer(mode_defaults$chunk_size),
    workers = workers_out,
    seed_root = 600000000L,
    condition_stride = 1000000L,
    replication_stride = 100L,
    warning_separator = " || ",
    boundary_tolerance = 1e-8,
    overwrite = isTRUE(overwrite),
    resume = isTRUE(resume),
    core_conditions = 216L,
    alignment_new_conditions = 16L,
    unique_run_conditions = 232L,
    stress_analysis_cells = 20L
  )
}

# Named component codes for Study 2 deterministic seeds.
# Codes occupy the trailing offset in:
#   seed_root + condition_id * condition_stride + rep_id * replication_stride + code
# with replication_stride = 100, so codes must be unique integers in 1..99.
# Structural components (n, k, X, alignment ties) are separate so that changing
# one allocation implementation does not alter unrelated random draws.
sim2_component_codes <- function() {
  c(
    n_allocation = 1L,
    k_allocation = 2L,
    x_assignment = 3L,
    alignment_ties = 4L,
    u_j = 5L,
    w_ij = 6L,
    e_ij = 7L,
    validation = 8L
  )
}

sim2_seed <- function(condition_id,
                      rep_id,
                      component,
                      config = make_sim2_config("smoke")) {
  condition_id <- as.integer(condition_id)
  rep_id <- as.integer(rep_id)
  if (length(condition_id) != 1L || is.na(condition_id) || condition_id < 1L) {
    stop("condition_id must be a positive integer.", call. = FALSE)
  }
  if (length(rep_id) != 1L || is.na(rep_id) || rep_id < 1L) {
    stop("rep_id must be a positive integer.", call. = FALSE)
  }

  codes <- sim2_component_codes()
  if (is.character(component)) {
    if (length(component) != 1L || !component %in% names(codes)) {
      stop("Unknown component: ", component, call. = FALSE)
    }
    code <- unname(codes[[component]])
  } else {
    code <- as.integer(component)
    if (length(code) != 1L || is.na(code) || code < 1L || code >= config$replication_stride) {
      stop("component code is outside the legal range.", call. = FALSE)
    }
  }

  seed <- as.numeric(config$seed_root) +
    as.numeric(condition_id) * as.numeric(config$condition_stride) +
    as.numeric(rep_id) * as.numeric(config$replication_stride) +
    as.numeric(code)

  if (!is.finite(seed) || seed > .Machine$integer.max || seed < -.Machine$integer.max) {
    stop("Computed seed is outside the valid R integer range: ", seed, call. = FALSE)
  }
  as.integer(seed)
}

sim2_rep_seed_base <- function(condition_id,
                               rep_id,
                               config = make_sim2_config("smoke")) {
  as.integer(
    as.numeric(config$seed_root) +
      as.numeric(as.integer(condition_id)) * as.numeric(config$condition_stride) +
      as.numeric(as.integer(rep_id)) * as.numeric(config$replication_stride)
  )
}

# Largest-remainder integerization to an exact total under lower/upper
# bounds. Fractional ties are resolved by original index.
integerize_bounded_largest_remainder <- function(x,
                                                 target_total,
                                                 lower,
                                                 upper,
                                                 tol = 1e-8) {
  if (!length(x)) {
    stop("x must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x must contain only finite values.", call. = FALSE)
  }
  target_total <- as.integer(target_total)
  lower <- as.integer(lower)
  upper <- as.integer(upper)
  if (lower > upper) {
    stop("lower must be <= upper.", call. = FALSE)
  }
  n <- length(x)
  if (target_total < n * lower || target_total > n * upper) {
    stop("target_total is infeasible under the stated bounds.", call. = FALSE)
  }
  if (any(x < lower - tol | x > upper + tol)) {
    stop("Continuous values lie outside [lower, upper].", call. = FALSE)
  }

  x_clamp <- pmin(as.numeric(upper), pmax(as.numeric(lower), as.numeric(x)))
  floored <- floor(x_clamp + tol)
  floored <- pmax(as.numeric(lower), pmin(as.numeric(upper), floored))
  remaining <- as.integer(target_total - sum(floored))
  if (remaining < 0L) {
    stop(
      "Floored sum exceeds target_total; cannot integerize under bounds.",
      call. = FALSE
    )
  }

  frac <- x_clamp - floored
  # Stable descending fractional order: secondary key is original index.
  ord <- order(-frac, seq_along(frac), method = "radix")

  for (idx in ord) {
    if (remaining <= 0L) {
      break
    }
    if (floored[[idx]] < upper) {
      floored[[idx]] <- floored[[idx]] + 1
      remaining <- remaining - 1L
    }
  }

  if (remaining > 0L) {
    guard <- 0L
    max_guard <- remaining * n + 1L
    while (remaining > 0L) {
      guard <- guard + 1L
      if (guard > max_guard) {
        stop("Cannot reach target_total without violating the upper bound.", call. = FALSE)
      }
      eligible <- which(floored < upper)
      if (!length(eligible)) {
        stop("Cannot reach target_total without violating the upper bound.", call. = FALSE)
      }
      # Preserve the original fractional ranking among remaining eligible cells.
      el_ord <- eligible[order(-frac[eligible], eligible, method = "radix")]
      for (idx in el_ord) {
        if (remaining <= 0L) {
          break
        }
        floored[[idx]] <- floored[[idx]] + 1
        remaining <- remaining - 1L
        # Classic largest-remainder adds at most one per pass over the ranked list.
      }
    }
  }

  out <- as.integer(floored)
  if (sum(out) != target_total) {
    stop("Integerized values do not sum to target_total.", call. = FALSE)
  }
  if (any(out < lower) || any(out > upper)) {
    stop("Integerized values violate bounds.", call. = FALSE)
  }
  out
}

expand_uniroot_bracket <- function(f, lo, hi, expand_factor = 2, max_expands = 60L) {
  f_lo <- f(lo)
  f_hi <- f(hi)
  expands <- 0L
  while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && expands < max_expands) {
    expands <- expands + 1L
    if (f_hi < 0) {
      # Scale still too small at hi; raise the upper bracket.
      hi <- hi * expand_factor
      f_hi <- f(hi)
    } else if (f_lo > 0) {
      # Scale still too large at lo; lower the lower bracket.
      lo <- lo / expand_factor
      f_lo <- f(lo)
    } else {
      hi <- hi * expand_factor
      lo <- lo / expand_factor
      f_lo <- f(lo)
      f_hi <- f(hi)
    }
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= 0) {
      break
    }
  }
  if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
    stop(
      "Failed to find a valid uniroot bracket for the sample-size scale.",
      call. = FALSE
    )
  }
  c(lo = lo, hi = hi)
}

# For non-balanced levels, uses bounded-lognormal quantiles scaled so that
# sum(n_j) = J * mean_n exactly after integerization.
make_nj_template_s2 <- function(J,
                                imbalance = c("balanced", "moderate", "severe", "extreme"),
                                mean_n = 300L,
                                n_min = 20L,
                                n_max = 2000L,
                                sigma_map = c(
                                  balanced = 0.0,
                                  moderate = 0.5,
                                  severe = 0.9,
                                  extreme = 1.4
                                )) {
  imbalance <- match.arg(imbalance)
  J <- as.integer(J)
  mean_n <- as.integer(mean_n)
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)

  if (!J %in% c(20L, 40L, 80L)) {
    stop("J must be 20, 40, or 80 for Study 2 sample-size templates.", call. = FALSE)
  }
  if (n_min > n_max) {
    stop("n_min must be <= n_max.", call. = FALSE)
  }
  if (J * n_min > J * mean_n || J * mean_n > J * n_max) {
    stop(
      "Infeasible mean_n under bounds: require n_min <= mean_n <= n_max.",
      call. = FALSE
    )
  }
  if (!imbalance %in% names(sigma_map)) {
    stop("Missing sigma for imbalance level '", imbalance, "'.", call. = FALSE)
  }

  if (identical(imbalance, "balanced")) {
    out <- rep(mean_n, J)
  } else {
    sigma <- as.numeric(sigma_map[[imbalance]])
    if (!is.finite(sigma) || sigma < 0) {
      stop("sigma_map values must be finite and nonnegative.", call. = FALSE)
    }
    probs <- (seq_len(J) - 0.5) / J
    z <- stats::qnorm(probs)
    q <- exp(sigma * z)
    target <- as.integer(J * mean_n)

    objective <- function(scale) {
      sum(pmin(n_max, pmax(n_min, scale * q))) - target
    }
    bracket <- expand_uniroot_bracket(objective, lo = 1e-6, hi = 1)
    scale <- stats::uniroot(objective, interval = bracket)$root
    x_continuous <- pmin(n_max, pmax(n_min, scale * q))
    out <- integerize_bounded_largest_remainder(
      x = x_continuous,
      target_total = target,
      lower = n_min,
      upper = n_max
    )
  }

  out <- as.integer(sort(out))
  if (any(out < n_min) || any(out > n_max)) {
    stop("n_j template violates bounds.", call. = FALSE)
  }
  if (sum(out) != J * mean_n) {
    stop("n_j template does not sum to J * mean_n.", call. = FALSE)
  }
  out
}

# Precision-concentration diagnostic J_prec = (sum(n-3))^2 / sum((n-3)^2).
j_precision_s2 <- function(n_j) {
  p <- as.numeric(n_j) - 3
  (sum(p)^2) / sum(p^2)
}

summarize_nj_templates_s2 <- function(config = make_sim2_config("smoke")) {
  rows <- lapply(config$J_levels, function(J) {
    do.call(rbind, lapply(config$n_imbalance_levels, function(imb) {
      nj <- make_nj_template_s2(
        J = J,
        imbalance = imb,
        mean_n = config$mean_n,
        n_min = config$n_min,
        n_max = config$n_max,
        sigma_map = config$sigma_log_n
      )
      prec <- as.numeric(nj) - 3
      data.frame(
        J = as.integer(J),
        n_imbalance = imb,
        sigma_log_n = unname(config$sigma_log_n[[imb]]),
        n_min = as.integer(min(nj)),
        n_q1 = as.numeric(stats::quantile(nj, 0.25, names = FALSE, type = 7)),
        n_median = as.numeric(stats::median(nj)),
        n_mean = mean(as.numeric(nj)),
        n_q3 = as.numeric(stats::quantile(nj, 0.75, names = FALSE, type = 7)),
        n_max = as.integer(max(nj)),
        n_total = as.integer(sum(nj)),
        cv_sample = stats::sd(as.numeric(nj)) / mean(as.numeric(nj)),
        J_precision = j_precision_s2(nj),
        J_precision_ratio = j_precision_s2(nj) / as.numeric(J),
        max_precision_share = max(prec) / sum(prec),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

make_kj_template_s2 <- function(J,
                                imbalance = c("balanced", "moderate", "severe")) {
  imbalance <- match.arg(imbalance)
  if (!exists("make_kj_template", mode = "function")) {
    stop("make_kj_template() is unavailable. Source R/study1_functions.R.", call. = FALSE)
  }
  kj <- make_kj_template(J = as.integer(J), imbalance = imbalance)
  as.integer(sort(kj))
}

# Exact moderator-group counts for Study 2 (lookup table; no rounding).
moderator_group_counts_s2 <- function(J,
                                      split = c("50_50", "70_30", "85_15")) {
  split <- match.arg(split)
  J <- as.integer(J)
  key <- paste(J, split, sep = "|")
  table <- c(
    "20|50_50" = "10,10",
    "20|70_30" = "14,6",
    "20|85_15" = "17,3",
    "40|50_50" = "20,20",
    "40|70_30" = "28,12",
    "40|85_15" = "34,6",
    "80|50_50" = "40,40",
    "80|70_30" = "56,24",
    "80|85_15" = "68,12"
  )
  if (!key %in% names(table)) {
    stop(
      "Unsupported J/split combination for moderator counts: J=", J,
      ", split=", split,
      call. = FALSE
    )
  }
  parts <- as.integer(strsplit(table[[key]], ",", fixed = TRUE)[[1]])
  list(n0 = parts[[1]], n1 = parts[[2]])
}

# Core 216-condition factorial grid for Simulation Study 2.
# Nested order (first factor slowest): J, k_imbalance, moderator_split,
# n_imbalance, beta1_true.
make_core_condition_grid_s2 <- function(config = make_sim2_config("smoke")) {
  # expand.grid varies the first argument fastest; reverse scientific order.
  grid <- expand.grid(
    beta1_true = config$beta1_levels,
    n_imbalance = config$n_imbalance_levels,
    moderator_split = config$split_levels,
    k_imbalance = config$k_imbalance_levels,
    J = as.integer(config$J_levels),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$condition_id <- seq_len(nrow(grid))
  grid$condition_set <- "core"
  grid$allocation_pattern <- "independent"
  grid <- grid[, c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern"
  ), drop = FALSE]
  grid$J <- as.integer(grid$J)
  rownames(grid) <- NULL
  if (nrow(grid) != config$core_conditions) {
    stop(
      "Core grid has ", nrow(grid), " rows; expected ",
      config$core_conditions, ".",
      call. = FALSE
    )
  }
  grid
}

# Twenty alignment-analysis cells map onto 16 additional simulation conditions.
make_alignment_condition_grid_s2 <- function(core_grid = NULL,
                                             config = make_sim2_config("smoke")) {
  if (is.null(core_grid)) {
    core_grid <- make_core_condition_grid_s2(config)
  }

  base <- core_grid[
    core_grid$J %in% c(20L, 40L) &
      core_grid$k_imbalance == "severe" &
      core_grid$moderator_split == "85_15" &
      core_grid$n_imbalance == "extreme",
    ,
    drop = FALSE
  ]
  base <- base[order(base$J, base$beta1_true), , drop = FALSE]
  if (nrow(base) != 4L) {
    stop("Expected 4 alignment base cells from the core grid; found ",
         nrow(base), ".", call. = FALSE)
  }

  patterns <- c(
    "independent",
    config$allocation_patterns_alignment
  )
  analysis_rows <- list()
  unique_rows <- list()
  stress_cell_id <- 0L
  next_condition_id <- 217L

  for (i in seq_len(nrow(base))) {
    b <- base[i, , drop = FALSE]
    for (pattern in patterns) {
      stress_cell_id <- stress_cell_id + 1L
      reuses_core <- identical(pattern, "independent")
      if (reuses_core) {
        condition_id <- as.integer(b$condition_id)
      } else {
        condition_id <- next_condition_id
        next_condition_id <- next_condition_id + 1L
        unique_rows[[length(unique_rows) + 1L]] <- data.frame(
          condition_id = condition_id,
          condition_set = "alignment",
          J = as.integer(b$J),
          k_imbalance = b$k_imbalance,
          moderator_split = b$moderator_split,
          n_imbalance = b$n_imbalance,
          beta1_true = as.numeric(b$beta1_true),
          allocation_pattern = pattern,
          stringsAsFactors = FALSE
        )
      }
      analysis_rows[[stress_cell_id]] <- data.frame(
        stress_cell_id = stress_cell_id,
        condition_id = condition_id,
        reuses_core_condition = reuses_core,
        J = as.integer(b$J),
        k_imbalance = b$k_imbalance,
        moderator_split = b$moderator_split,
        n_imbalance = b$n_imbalance,
        beta1_true = as.numeric(b$beta1_true),
        allocation_pattern = pattern,
        stringsAsFactors = FALSE
      )
    }
  }

  analysis_grid <- do.call(rbind, analysis_rows)
  unique_run_grid <- do.call(rbind, unique_rows)
  rownames(analysis_grid) <- NULL
  rownames(unique_run_grid) <- NULL

  if (nrow(analysis_grid) != config$stress_analysis_cells) {
    stop("Alignment analysis grid must have ",
         config$stress_analysis_cells, " rows.", call. = FALSE)
  }
  if (nrow(unique_run_grid) != config$alignment_new_conditions) {
    stop("Alignment unique-run grid must have ",
         config$alignment_new_conditions, " rows.", call. = FALSE)
  }
  if (!identical(as.integer(unique_run_grid$condition_id), 217:232)) {
    stop("Aligned unique condition IDs must be 217:232.", call. = FALSE)
  }

  list(
    analysis_grid = analysis_grid,
    unique_run_grid = unique_run_grid
  )
}

make_unique_run_grid_s2 <- function(config = make_sim2_config("smoke")) {
  core <- make_core_condition_grid_s2(config)
  alignment <- make_alignment_condition_grid_s2(core, config)
  out <- rbind(core, alignment$unique_run_grid)
  rownames(out) <- NULL
  if (nrow(out) != config$unique_run_conditions) {
    stop(
      "Unique run grid has ", nrow(out), " rows; expected ",
      config$unique_run_conditions, ".",
      call. = FALSE
    )
  }
  out
}

normalize_logical_like <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.character(x) || is.factor(x)) {
    xc <- tolower(as.character(x))
    if (all(xc %in% c("true", "false", "t", "f", "1", "0"))) {
      return(xc %in% c("true", "t", "1"))
    }
  }
  x
}
