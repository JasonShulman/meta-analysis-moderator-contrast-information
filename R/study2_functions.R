# Simulation Study 2 allocation, data generation, replication, chunks, and summaries.
#
# Estimators are fit by estimators.R (aggregate + Knapp-Hartung, MLMA, MLMA+CR2, CHE-RVE).
# Component seeds use withr::with_seed() and are not combined across components.

# Evaluate under an explicit component seed.
run_with_sim2_seed <- function(seed, code) {
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop(
      "Package 'withr' is required for Study 2 component seeding. ",
      "Install it before running allocation or DGM code.",
      call. = FALSE
    )
  }
  withr::with_seed(
    seed,
    code,
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
}

same_multiset <- function(x, y) {
  if (is.numeric(x) && is.numeric(y)) {
    return(isTRUE(identical(sort(as.numeric(x)), sort(as.numeric(y)))))
  }
  isTRUE(identical(sort(as.character(x)), sort(as.character(y))))
}

moderator_template_s2 <- function(J, moderator_split) {
  counts <- moderator_group_counts_s2(J, moderator_split)
  c(rep(0L, counts$n0), rep(1L, counts$n1))
}

assign_moderator_labels_s2 <- function(J, moderator_split, seed) {
  template <- moderator_template_s2(J, moderator_split)
  as.integer(run_with_sim2_seed(seed, sample(template)))
}

# When aligning k_j with n_j, shuffle within tied k_j blocks so ties are
# not resolved by their original ordering.
shuffle_within_ties <- function(x) {
  x <- as.vector(x)
  if (!length(x)) {
    return(x)
  }
  x_ord <- x[order(x, method = "radix")]
  out <- x_ord
  blocks <- split(seq_along(x_ord), x_ord)
  for (idx in blocks) {
    if (length(idx) > 1L) {
      out[idx] <- x_ord[sample(idx)]
    }
  }
  out
}

allocate_independent_s2 <- function(J,
                                    nj_template,
                                    kj_template,
                                    moderator_split,
                                    condition_id,
                                    rep_id,
                                    config) {
  seed_n <- sim2_seed(condition_id, rep_id, "n_allocation", config)
  seed_k <- sim2_seed(condition_id, rep_id, "k_allocation", config)
  seed_x <- sim2_seed(condition_id, rep_id, "x_assignment", config)

  list(
    n_j = as.integer(run_with_sim2_seed(seed_n, sample(as.integer(nj_template)))),
    k_j = as.integer(run_with_sim2_seed(seed_k, sample(as.integer(kj_template)))),
    X = assign_moderator_labels_s2(J, moderator_split, seed_x)
  )
}

# n-only alignment (n_majority / n_minority).
# alignment_ties draw order:
# 1. permute n values within the minority group
# 2. permute n values within the majority group
# k_j is permuted independently via k_allocation.
allocate_n_aligned_s2 <- function(J,
                                  nj_template,
                                  kj_template,
                                  moderator_split,
                                  minority_gets = c("smallest", "largest"),
                                  condition_id,
                                  rep_id,
                                  config) {
  minority_gets <- match.arg(minority_gets)
  seed_x <- sim2_seed(condition_id, rep_id, "x_assignment", config)
  seed_k <- sim2_seed(condition_id, rep_id, "k_allocation", config)
  seed_tie <- sim2_seed(condition_id, rep_id, "alignment_ties", config)

  X <- assign_moderator_labels_s2(J, moderator_split, seed_x)
  counts <- moderator_group_counts_s2(J, moderator_split)
  n1 <- counts$n1

  n_sorted <- sort(as.integer(nj_template))
  if (identical(minority_gets, "smallest")) {
    n_min <- n_sorted[seq_len(n1)]
    n_maj <- n_sorted[seq.int(n1 + 1L, J)]
  } else {
    n_min <- n_sorted[seq.int(J - n1 + 1L, J)]
    n_maj <- n_sorted[seq_len(J - n1)]
  }

  min_studies <- which(X == 1L)
  maj_studies <- which(X == 0L)

  n_assigns <- run_with_sim2_seed(seed_tie, {
    list(min = sample(n_min), maj = sample(n_maj))
  })

  n_j <- integer(J)
  n_j[min_studies] <- n_assigns$min
  n_j[maj_studies] <- n_assigns$maj

  list(
    n_j = as.integer(n_j),
    k_j = as.integer(run_with_sim2_seed(seed_k, sample(as.integer(kj_template)))),
    X = X
  )
}

# Joint n/k alignment (nk_majority / nk_minority).
# alignment_ties draw order (single with_seed block):
# 1. within-tie shuffle while building the positively aligned k vector
# 2. permute minority pairs onto minority studies
# 3. permute majority pairs onto majority studies
allocate_nk_aligned_s2 <- function(J,
                                   nj_template,
                                   kj_template,
                                   moderator_split,
                                   minority_gets = c("smallest", "largest"),
                                   condition_id,
                                   rep_id,
                                   config) {
  minority_gets <- match.arg(minority_gets)
  seed_x <- sim2_seed(condition_id, rep_id, "x_assignment", config)
  seed_tie <- sim2_seed(condition_id, rep_id, "alignment_ties", config)

  X <- assign_moderator_labels_s2(J, moderator_split, seed_x)
  counts <- moderator_group_counts_s2(J, moderator_split)
  n1 <- counts$n1
  min_studies <- which(X == 1L)
  maj_studies <- which(X == 0L)

  assigned <- run_with_sim2_seed(seed_tie, {
    n_sorted <- sort(as.integer(nj_template))
    k_aligned <- shuffle_within_ties(sort(as.integer(kj_template)))
    pairs <- data.frame(
      n_j = n_sorted,
      k_j = as.integer(k_aligned),
      stringsAsFactors = FALSE
    )

    if (identical(minority_gets, "smallest")) {
      pairs_min <- pairs[seq_len(n1), , drop = FALSE]
      pairs_maj <- pairs[seq.int(n1 + 1L, J), , drop = FALSE]
    } else {
      pairs_min <- pairs[seq.int(J - n1 + 1L, J), , drop = FALSE]
      pairs_maj <- pairs[seq_len(J - n1), , drop = FALSE]
    }

    pairs_min <- pairs_min[sample.int(nrow(pairs_min)), , drop = FALSE]
    pairs_maj <- pairs_maj[sample.int(nrow(pairs_maj)), , drop = FALSE]

    n_j <- integer(J)
    k_j <- integer(J)
    n_j[min_studies] <- pairs_min$n_j
    k_j[min_studies] <- pairs_min$k_j
    n_j[maj_studies] <- pairs_maj$n_j
    k_j[maj_studies] <- pairs_maj$k_j
    list(n_j = n_j, k_j = k_j)
  })

  list(
    n_j = as.integer(assigned$n_j),
    k_j = as.integer(assigned$k_j),
    X = X
  )
}

# Independent allocations can show chance associations, so they return NA.
alignment_rule_holds_s2 <- function(structure_df, allocation_pattern) {
  X <- structure_df$X
  n_j <- structure_df$n_j
  k_j <- structure_df$k_j

  if (identical(allocation_pattern, "independent")) {
    return(NA)
  }
  if (identical(allocation_pattern, "n_majority")) {
    return(isTRUE(max(n_j[X == 1L]) <= min(n_j[X == 0L])))
  }
  if (identical(allocation_pattern, "n_minority")) {
    return(isTRUE(min(n_j[X == 1L]) >= max(n_j[X == 0L])))
  }
  if (identical(allocation_pattern, "nk_majority")) {
    n_ok <- isTRUE(max(n_j[X == 1L]) <= min(n_j[X == 0L]))
    # Pair membership: the joint (n,k) multiset must equal sort(n) x sort(k).
    pairs_ok <- same_multiset(
      paste(n_j, k_j, sep = ":"),
      paste(sort(n_j), sort(k_j), sep = ":")
    )
    return(isTRUE(n_ok && pairs_ok))
  }
  if (identical(allocation_pattern, "nk_minority")) {
    n_ok <- isTRUE(min(n_j[X == 1L]) >= max(n_j[X == 0L]))
    pairs_ok <- same_multiset(
      paste(n_j, k_j, sep = ":"),
      paste(sort(n_j), sort(k_j), sep = ":")
    )
    return(isTRUE(n_ok && pairs_ok))
  }
  stop("Unknown allocation_pattern: ", allocation_pattern, call. = FALSE)
}

validate_study_structure_s2 <- function(structure_df,
                                        nj_template,
                                        kj_template,
                                        moderator_split,
                                        allocation_pattern) {
  J <- nrow(structure_df)
  if (!identical(as.integer(structure_df$study), seq_len(J))) {
    stop("Study IDs must be fixed as 1:J in order.", call. = FALSE)
  }
  if (!same_multiset(structure_df$n_j, nj_template)) {
    stop("n_j multiset does not match the sample-size template.", call. = FALSE)
  }
  if (!same_multiset(structure_df$k_j, kj_template)) {
    stop("k_j multiset does not match the effect-count template.", call. = FALSE)
  }

  counts <- moderator_group_counts_s2(J, moderator_split)
  if (sum(structure_df$X == 0L) != counts$n0 ||
      sum(structure_df$X == 1L) != counts$n1) {
    stop("Moderator counts are not exact.", call. = FALSE)
  }
  if (sum(structure_df$n_j) != 300L * J) {
    stop("Total N must equal 300J.", call. = FALSE)
  }
  if (sum(structure_df$k_j) != 4L * J) {
    stop("Total K must equal 4J.", call. = FALSE)
  }
  if (!identical(
    as.character(structure_df$allocation_pattern),
    rep(allocation_pattern, J)
  )) {
    stop("allocation_pattern column is inconsistent.", call. = FALSE)
  }

  if (!identical(allocation_pattern, "independent")) {
    if (!isTRUE(alignment_rule_holds_s2(structure_df, allocation_pattern))) {
      stop(
        "Allocation does not satisfy the alignment rule for pattern '",
        allocation_pattern, "'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# Assign study-level sample sizes, effect counts, and moderator labels.
# Study IDs are fixed as 1:J. Randomness uses distinct component seeds derived
# from (condition_id, rep_id) via sim2_seed().
# Independent core conditions permute n_j, k_j, and X independently. Aligned
# conditions first assign X, then allocate templates by the alignment rule, then
# randomize within groups using alignment_ties.
assign_study_structure_s2 <- function(J,
                                      nj_template,
                                      kj_template,
                                      moderator_split = c("50_50", "70_30", "85_15"),
                                      allocation_pattern = c(
                                        "independent",
                                        "n_majority",
                                        "n_minority",
                                        "nk_majority",
                                        "nk_minority"
                                      ),
                                      condition_id,
                                      rep_id,
                                      config = make_sim2_config("smoke")) {
  moderator_split <- match.arg(moderator_split)
  allocation_pattern <- match.arg(allocation_pattern)
  J <- as.integer(J)
  nj_template <- as.integer(nj_template)
  kj_template <- as.integer(kj_template)

  if (length(nj_template) != J || length(kj_template) != J) {
    stop("nj_template and kj_template must have length J.", call. = FALSE)
  }
  if (sum(nj_template) != 300L * J) {
    stop("nj_template must sum to 300J.", call. = FALSE)
  }
  if (sum(kj_template) != 4L * J) {
    stop("kj_template must sum to 4J.", call. = FALSE)
  }

  allocated <- switch(
    allocation_pattern,
    independent = allocate_independent_s2(
      J, nj_template, kj_template, moderator_split,
      condition_id, rep_id, config
    ),
    n_majority = allocate_n_aligned_s2(
      J, nj_template, kj_template, moderator_split,
      minority_gets = "smallest",
      condition_id, rep_id, config
    ),
    n_minority = allocate_n_aligned_s2(
      J, nj_template, kj_template, moderator_split,
      minority_gets = "largest",
      condition_id, rep_id, config
    ),
    nk_majority = allocate_nk_aligned_s2(
      J, nj_template, kj_template, moderator_split,
      minority_gets = "smallest",
      condition_id, rep_id, config
    ),
    nk_minority = allocate_nk_aligned_s2(
      J, nj_template, kj_template, moderator_split,
      minority_gets = "largest",
      condition_id, rep_id, config
    )
  )

  structure_df <- data.frame(
    study = seq_len(J),
    n_j = as.integer(allocated$n_j),
    k_j = as.integer(allocated$k_j),
    X = as.integer(allocated$X),
    n_rank = as.numeric(rank(allocated$n_j, ties.method = "average")),
    k_rank = as.numeric(rank(allocated$k_j, ties.method = "average")),
    allocation_pattern = allocation_pattern,
    stringsAsFactors = FALSE
  )

  validate_study_structure_s2(
    structure_df = structure_df,
    nj_template = nj_template,
    kj_template = kj_template,
    moderator_split = moderator_split,
    allocation_pattern = allocation_pattern
  )

  structure_df
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

summarize_design_diagnostics_s2 <- function(study_structure,
                                            condition_id = NULL,
                                            rep_id = NULL) {
  if (!nrow(study_structure)) {
    stop("study_structure is empty.", call. = FALSE)
  }

  if (is.null(condition_id)) {
    condition_id <- if ("condition_id" %in% names(study_structure)) {
      study_structure$condition_id[[1]]
    } else {
      NA_integer_
    }
  }
  if (is.null(rep_id)) {
    rep_id <- if ("rep_id" %in% names(study_structure)) {
      study_structure$rep_id[[1]]
    } else {
      NA_integer_
    }
  }

  n_j <- as.numeric(study_structure$n_j)
  k_j <- as.numeric(study_structure$k_j)
  X <- as.integer(study_structure$X)
  J <- length(n_j)
  precision_j <- n_j - 3
  q_j <- k_j * precision_j

  precision_X0 <- sum(precision_j[X == 0L])
  precision_X1 <- sum(precision_j[X == 1L])
  effects_X0 <- sum(k_j[X == 0L])
  effects_X1 <- sum(k_j[X == 1L])
  q_X0 <- sum(q_j[X == 0L])
  q_X1 <- sum(q_j[X == 1L])
  j_prec <- (sum(precision_j)^2) / sum(precision_j^2)
  j_k <- (sum(k_j)^2) / sum(k_j^2)

  data.frame(
    condition_id = as.integer(condition_id),
    rep_id = as.integer(rep_id),
    J = as.integer(J),
    K_total = as.integer(sum(k_j)),
    N_total = as.integer(sum(n_j)),
    n_min = as.numeric(min(n_j)),
    n_median = as.numeric(stats::median(n_j)),
    n_max = as.numeric(max(n_j)),
    n_mean = mean(n_j),
    n_cv = stats::sd(n_j) / mean(n_j),
    max_precision_share = max(precision_j) / sum(precision_j),
    J_precision = j_prec,
    J_precision_ratio = j_prec / as.numeric(J),
    precision_X0 = precision_X0,
    precision_X1 = precision_X1,
    precision_X1_prop = precision_X1 / (precision_X0 + precision_X1),
    precision_X1_X0_ratio = precision_X1 / precision_X0,
    k_min = as.numeric(min(k_j)),
    k_median = as.numeric(stats::median(k_j)),
    k_max = as.numeric(max(k_j)),
    k_mean = mean(k_j),
    max_k_share = max(k_j) / sum(k_j),
    J_k = j_k,
    J_k_ratio = j_k / as.numeric(J),
    effects_X0 = as.numeric(effects_X0),
    effects_X1 = as.numeric(effects_X1),
    effects_X1_prop = effects_X1 / (effects_X0 + effects_X1),
    pearson_n_k = safe_cor_s2(n_j, k_j, "pearson"),
    spearman_n_k = safe_cor_s2(n_j, k_j, "spearman"),
    max_q_share = max(q_j) / sum(q_j),
    q_X0 = q_X0,
    q_X1 = q_X1,
    q_X1_prop = q_X1 / (q_X0 + q_X1),
    allocation_pattern = as.character(study_structure$allocation_pattern[[1]]),
    stringsAsFactors = FALSE
  )
}

normalize_condition_s2 <- function(condition) {
  if (is.data.frame(condition)) {
    if (nrow(condition) != 1L) {
      stop("condition must be a one-row data frame or list.", call. = FALSE)
    }
    condition <- as.list(condition[1, , drop = FALSE])
  }
  required <- c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern"
  )
  missing <- setdiff(required, names(condition))
  if (length(missing)) {
    stop(
      "condition is missing required field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  condition$condition_id <- as.integer(condition$condition_id)
  condition$J <- as.integer(condition$J)
  condition$beta1_true <- as.numeric(condition$beta1_true)
  condition$condition_set <- as.character(condition$condition_set)
  condition$k_imbalance <- as.character(condition$k_imbalance)
  condition$moderator_split <- as.character(condition$moderator_split)
  condition$n_imbalance <- as.character(condition$n_imbalance)
  condition$allocation_pattern <- as.character(condition$allocation_pattern)
  condition
}

draw_u_j_s2 <- function(J, tau, condition_id, rep_id, config) {
  seed <- sim2_seed(condition_id, rep_id, "u_j", config)
  run_with_sim2_seed(seed, stats::rnorm(J, mean = 0, sd = tau))
}

# Documented order: studies 1:J, and within study es_id 1:k_j.
draw_w_ij_s2 <- function(k_j, omega, condition_id, rep_id, config) {
  seed <- sim2_seed(condition_id, rep_id, "w_ij", config)
  K <- sum(as.integer(k_j))
  as.numeric(run_with_sim2_seed(seed, stats::rnorm(K, mean = 0, sd = omega)))
}

# Draw sampling errors under the e_ij component seed.
# Uses Study 1's make_sampling_Sigma() + MASS::mvrnorm() (or rnorm for k=1).
# Documented order: studies 1:J.
draw_e_ij_s2 <- function(n_j, k_j, rho, condition_id, rep_id, config) {
  if (!exists("make_sampling_Sigma", mode = "function")) {
    stop("make_sampling_Sigma() is unavailable. Source R/study1_functions.R.", call. = FALSE)
  }
  seed <- sim2_seed(condition_id, rep_id, "e_ij", config)
  n_j <- as.numeric(n_j)
  k_j <- as.integer(k_j)
  J <- length(k_j)

  run_with_sim2_seed(seed, {
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
  })
}

simulate_meta_dataset_s2 <- function(condition,
                                     rep_id,
                                     config = make_sim2_config("smoke"),
                                     return_components = FALSE) {
  condition <- normalize_condition_s2(condition)
  rep_id <- as.integer(rep_id)

  J <- condition$J
  nj_template <- make_nj_template_s2(
    J = J,
    imbalance = condition$n_imbalance,
    mean_n = config$mean_n,
    n_min = config$n_min,
    n_max = config$n_max,
    sigma_map = config$sigma_log_n
  )
  kj_template <- make_kj_template_s2(
    J = J,
    imbalance = condition$k_imbalance
  )

  study_structure <- assign_study_structure_s2(
    J = J,
    nj_template = nj_template,
    kj_template = kj_template,
    moderator_split = condition$moderator_split,
    allocation_pattern = condition$allocation_pattern,
    condition_id = condition$condition_id,
    rep_id = rep_id,
    config = config
  )
  study_structure$condition_id <- condition$condition_id
  study_structure$condition_set <- condition$condition_set
  study_structure$rep_id <- rep_id
  study_structure$beta1_true <- condition$beta1_true

  u_j <- draw_u_j_s2(
    J = J,
    tau = config$tau,
    condition_id = condition$condition_id,
    rep_id = rep_id,
    config = config
  )
  w_all <- draw_w_ij_s2(
    k_j = study_structure$k_j,
    omega = config$omega,
    condition_id = condition$condition_id,
    rep_id = rep_id,
    config = config
  )
  e_list <- draw_e_ij_s2(
    n_j = study_structure$n_j,
    k_j = study_structure$k_j,
    rho = config$rho_true,
    condition_id = condition$condition_id,
    rep_id = rep_id,
    config = config
  )

  rows <- vector("list", J)
  w_offset <- 0L
  for (j in seq_len(J)) {
    kj_j <- as.integer(study_structure$k_j[[j]])
    n_jj <- as.numeric(study_structure$n_j[[j]])
    X_j <- as.integer(study_structure$X[[j]])
    vi_j <- rep(1 / (n_jj - 3), kj_j)
    w_ij <- w_all[seq.int(w_offset + 1L, w_offset + kj_j)]
    w_offset <- w_offset + kj_j
    e_ij <- e_list[[j]]
    true_theta_ij <- config$beta0 + condition$beta1_true * X_j + u_j[[j]] + w_ij
    yi <- true_theta_ij + e_ij

    df_j <- data.frame(
      condition_id = condition$condition_id,
      condition_set = condition$condition_set,
      rep_id = rep_id,
      study = as.integer(j),
      es_id = seq_len(kj_j),
      X = X_j,
      k_j = kj_j,
      n_j = n_jj,
      vi = vi_j,
      yi = as.numeric(yi),
      allocation_pattern = condition$allocation_pattern,
      stringsAsFactors = FALSE
    )
    if (isTRUE(return_components)) {
      df_j$u_j <- u_j[[j]]
      df_j$w_ij <- as.numeric(w_ij)
      df_j$e_ij <- as.numeric(e_ij)
      df_j$true_theta_ij <- as.numeric(true_theta_ij)
    }
    rows[[j]] <- df_j
  }

  data <- do.call(rbind, rows)
  rownames(data) <- NULL

  diagnostics <- summarize_design_diagnostics_s2(
    study_structure = study_structure,
    condition_id = condition$condition_id,
    rep_id = rep_id
  )

  list(
    data = data,
    study_structure = study_structure,
    diagnostics = diagnostics
  )
}

# Replication execution and performance summaries.

run_one_replication_s2 <- function(condition,
                                   rep_id,
                                   config = make_sim2_config("smoke"),
                                   return_data = FALSE) {
  condition <- normalize_condition_s2(condition)
  rep_id <- as.integer(rep_id)
  if (length(rep_id) != 1L || is.na(rep_id) || rep_id < 1L) {
    stop("rep_id must be a positive integer.", call. = FALSE)
  }

  rep_seed_base <- sim2_rep_seed_base(
    condition_id = condition$condition_id,
    rep_id = rep_id,
    config = config
  )

  sim <- tryCatch(
    simulate_meta_dataset_s2(
      condition = condition,
      rep_id = rep_id,
      config = config,
      return_components = FALSE
    ),
    error = function(e) e
  )

  if (inherits(sim, "error")) {
    err <- conditionMessage(sim)
    methods <- c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve")
    method_results <- do.call(
      rbind,
      lapply(methods, function(m) {
        empty_method_result_s2(
          method = m,
          condition = condition,
          rep_id = rep_id,
          error = err,
          status = "error",
          config = config
        )
      })
    )
    method_results$rep_seed_base <- rep_seed_base
    diagnostics <- data.frame(
      condition_id = as.integer(condition$condition_id),
      condition_set = as.character(condition$condition_set),
      rep_id = rep_id,
      rep_seed_base = rep_seed_base,
      J = as.integer(condition$J),
      K_total = NA_integer_,
      N_total = NA_integer_,
      allocation_pattern = as.character(condition$allocation_pattern),
      dgm_error = err,
      stringsAsFactors = FALSE
    )
    out <- list(
      method_results = method_results,
      dataset_diagnostics = diagnostics
    )
    if (isTRUE(return_data)) {
      out$simulated_data <- NULL
    }
    return(out)
  }

  method_results <- fit_all_methods_s2(
    data = sim$data,
    condition = condition,
    rep_id = rep_id,
    config = config
  )
  method_results$condition_set <- as.character(condition$condition_set)
  method_results$rep_seed_base <- rep_seed_base

  diagnostics <- sim$diagnostics
  diagnostics$condition_set <- as.character(condition$condition_set)
  diagnostics$rep_seed_base <- rep_seed_base

  # Stable column order: identifiers first.
  id_cols <- c("condition_id", "condition_set", "rep_id", "rep_seed_base")
  method_results <- method_results[, c(
    id_cols,
    setdiff(names(method_results), id_cols)
  ), drop = FALSE]
  diagnostics <- diagnostics[, c(
    id_cols,
    setdiff(names(diagnostics), id_cols)
  ), drop = FALSE]

  out <- list(
    method_results = method_results,
    dataset_diagnostics = diagnostics
  )
  if (isTRUE(return_data)) {
    out$simulated_data <- sim$data
  }
  out
}

attach_condition_factors_s2 <- function(results, condition_grid) {
  keep <- c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern"
  )
  missing <- setdiff(keep, names(condition_grid))
  if (length(missing)) {
    stop(
      "condition_grid missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  design <- unique(condition_grid[, keep, drop = FALSE])
  # Drop overlapping non-key columns from results before join.
  overlap <- setdiff(intersect(names(results), keep), "condition_id")
  if (length(overlap)) {
    results <- results[, setdiff(names(results), overlap), drop = FALSE]
  }
  merged <- merge(results, design, by = "condition_id", all.x = TRUE, sort = FALSE)
  if (any(is.na(merged$J))) {
    stop("Some condition_id values were not found in condition_grid.", call. = FALSE)
  }
  merged
}

summarize_method_performance_s2 <- function(results, config, condition_grid = NULL) {
  if (!nrow(results)) {
    stop("results is empty.", call. = FALSE)
  }
  if (!is.null(condition_grid)) {
    results <- attach_condition_factors_s2(results, condition_grid)
  }
  group_cols <- c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern", "method"
  )
  missing <- setdiff(group_cols, names(results))
  if (length(missing)) {
    stop(
      "results missing grouping columns: ",
      paste(missing, collapse = ", "),
      ". Pass condition_grid or pre-join design factors.",
      call. = FALSE
    )
  }

  split_keys <- interaction(results[, group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(seq_len(nrow(results)), split_keys), function(idx) {
    d <- results[idx, , drop = FALSE]
    design <- d[1, group_cols, drop = FALSE]
    beta1_true <- as.numeric(design$beta1_true)

    n_attempted <- nrow(d)
    n_success <- sum(d$status == "success", na.rm = TRUE)
    n_error <- sum(d$status %in% c("error", "upstream_fit_failure"), na.rm = TRUE)
    n_invalid <- sum(d$status == "invalid_result", na.rm = TRUE)
    error_rate <- (n_attempted - n_success) / n_attempted
    warning_rate <- mean(d$warning_count > 0, na.rm = TRUE)

    d_ok <- d[d$status == "success", , drop = FALSE]
    if (n_success > 0L) {
      mean_estimate <- mean(d_ok$estimate)
      bias <- mean_estimate - beta1_true
      absolute_bias <- abs(bias)
      empirical_sd <- stats::sd(d_ok$estimate)
      mean_se <- mean(d_ok$se)
      se_ratio <- mean_se / empirical_sd
      rmse <- sqrt(mean((d_ok$estimate - beta1_true)^2))
      covered <- (d_ok$ci_lower <= beta1_true) & (d_ok$ci_upper >= beta1_true)
      coverage <- mean(covered)
      mean_ci_width <- mean(d_ok$ci_upper - d_ok$ci_lower)
      reject_rate <- mean(d_ok$p_value < config$alpha)
      mcse_reject_rate <- sqrt(reject_rate * (1 - reject_rate) / n_success)
      mcse_coverage <- sqrt(coverage * (1 - coverage) / n_success)
      mcse_bias <- empirical_sd / sqrt(n_success)
      mcse_empirical_sd <- if (n_success > 1L) {
        empirical_sd / sqrt(2 * (n_success - 1))
      } else {
        NA_real_
      }

      dfs <- d_ok$df
      if (all(is.na(dfs))) {
        mean_df <- median_df <- sd_df <- min_df <- q05_df <- q10_df <- max_df <- NA_real_
        prop_df_lt_2 <- prop_df_lt_4 <- prop_df_lt_10 <- NA_real_
      } else {
        mean_df <- mean(dfs, na.rm = TRUE)
        median_df <- stats::median(dfs, na.rm = TRUE)
        sd_df <- stats::sd(dfs, na.rm = TRUE)
        min_df <- min(dfs, na.rm = TRUE)
        q05_df <- as.numeric(stats::quantile(dfs, 0.05, na.rm = TRUE, names = FALSE))
        q10_df <- as.numeric(stats::quantile(dfs, 0.10, na.rm = TRUE, names = FALSE))
        max_df <- max(dfs, na.rm = TRUE)
        prop_df_lt_2 <- mean(dfs < 2, na.rm = TRUE)
        prop_df_lt_4 <- mean(dfs < 4, na.rm = TRUE)
        prop_df_lt_10 <- mean(dfs < 10, na.rm = TRUE)
      }

      mean_variance_study <- mean(d_ok$variance_study, na.rm = TRUE)
      mean_variance_effect <- mean(d_ok$variance_effect, na.rm = TRUE)
      mean_variance_residual <- mean(d_ok$variance_residual, na.rm = TRUE)
      boundary_study_rate <- mean(d_ok$boundary_study %in% TRUE, na.rm = TRUE)
      boundary_effect_rate <- mean(d_ok$boundary_effect %in% TRUE, na.rm = TRUE)
      boundary_residual_rate <- mean(d_ok$boundary_residual %in% TRUE, na.rm = TRUE)
      mean_fit_runtime_sec <- mean(d_ok$fit_runtime_sec, na.rm = TRUE)
    } else {
      mean_estimate <- bias <- absolute_bias <- empirical_sd <- mean_se <- se_ratio <-
        rmse <- coverage <- mean_ci_width <- reject_rate <- mcse_reject_rate <-
        mcse_coverage <- mcse_bias <- mcse_empirical_sd <- mean_df <- median_df <-
        sd_df <- min_df <- q05_df <- q10_df <- max_df <- prop_df_lt_2 <-
        prop_df_lt_4 <- prop_df_lt_10 <- mean_variance_study <-
        mean_variance_effect <- mean_variance_residual <- boundary_study_rate <-
        boundary_effect_rate <- boundary_residual_rate <- mean_fit_runtime_sec <-
        NA_real_
    }

    type1_error <- if (isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    power <- if (!isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_

    cbind(
      design,
      data.frame(
        n_attempted = n_attempted,
        n_success = n_success,
        n_error = n_error,
        n_invalid = n_invalid,
        error_rate = error_rate,
        warning_rate = warning_rate,
        mean_estimate = mean_estimate,
        bias = bias,
        absolute_bias = absolute_bias,
        empirical_sd = empirical_sd,
        mean_se = mean_se,
        se_ratio = se_ratio,
        rmse = rmse,
        coverage = coverage,
        mcse_coverage = mcse_coverage,
        mean_ci_width = mean_ci_width,
        reject_rate = reject_rate,
        mcse_reject_rate = mcse_reject_rate,
        mcse_bias = mcse_bias,
        mcse_empirical_sd = mcse_empirical_sd,
        type1_error = type1_error,
        power = power,
        mean_df = mean_df,
        median_df = median_df,
        sd_df = sd_df,
        min_df = min_df,
        q05_df = q05_df,
        q10_df = q10_df,
        max_df = max_df,
        prop_df_lt_2 = prop_df_lt_2,
        prop_df_lt_4 = prop_df_lt_4,
        prop_df_lt_10 = prop_df_lt_10,
        mean_variance_study = mean_variance_study,
        mean_variance_effect = mean_variance_effect,
        mean_variance_residual = mean_variance_residual,
        boundary_study_rate = boundary_study_rate,
        boundary_effect_rate = boundary_effect_rate,
        boundary_residual_rate = boundary_residual_rate,
        mean_fit_runtime_sec = mean_fit_runtime_sec,
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(out$condition_id, out$method), , drop = FALSE]
}

summarize_design_diagnostics_s2_by_condition <- function(diagnostics,
                                                        condition_grid = NULL) {
  if (!nrow(diagnostics)) {
    stop("diagnostics is empty.", call. = FALSE)
  }
  if (!is.null(condition_grid)) {
    diagnostics <- attach_condition_factors_s2(diagnostics, condition_grid)
  }
  group_cols <- c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern"
  )
  missing <- setdiff(group_cols, names(diagnostics))
  if (length(missing)) {
    stop("diagnostics missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  split_keys <- interaction(diagnostics[, group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(seq_len(nrow(diagnostics)), split_keys), function(idx) {
    d <- diagnostics[idx, , drop = FALSE]
    design <- d[1, group_cols, drop = FALSE]
    num_cols <- c(
      "N_total", "K_total", "n_min", "n_median", "n_max", "n_mean", "n_cv",
      "max_precision_share", "J_precision", "J_precision_ratio",
      "precision_X0", "precision_X1", "precision_X1_prop", "precision_X1_X0_ratio",
      "k_min", "k_median", "k_max", "k_mean", "max_k_share", "J_k", "J_k_ratio",
      "effects_X0", "effects_X1", "effects_X1_prop",
      "pearson_n_k", "spearman_n_k", "max_q_share", "q_X0", "q_X1", "q_X1_prop"
    )
    num_cols <- intersect(num_cols, names(d))
    means <- lapply(num_cols, function(nm) mean(d[[nm]], na.rm = TRUE))
    names(means) <- paste0("mean_", num_cols)
    cbind(
      design,
      data.frame(
        n_rep = nrow(d),
        N_total_ok = all(d$N_total == 300L * as.integer(design$J), na.rm = TRUE),
        K_total_ok = all(d$K_total == 4L * as.integer(design$J), na.rm = TRUE),
        means,
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(out$condition_id), , drop = FALSE]
}

# Build the 20-cell alignment stress-test summary table.
summarize_alignment_stress_s2 <- function(method_summary,
                                          alignment_analysis = NULL,
                                          root = NULL) {
  if (is.null(alignment_analysis)) {
    alignment_analysis <- read_sim2_fixture(
      "condition_grid_alignment_analysis.csv",
      root = root
    )
    alignment_analysis$reuses_core_condition <- normalize_logical_like(
      alignment_analysis$reuses_core_condition
    )
  }
  keep_methods <- method_summary
  joined <- merge(
    alignment_analysis,
    keep_methods,
    by = c(
      "condition_id", "J", "k_imbalance", "moderator_split",
      "n_imbalance", "beta1_true", "allocation_pattern"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  if (length(unique(joined$stress_cell_id)) != 20L) {
    stop(
      "Alignment stress join did not yield 20 stress_cell_id values.",
      call. = FALSE
    )
  }
  joined[order(joined$stress_cell_id, joined$method), , drop = FALSE]
}

# Chunking and parallel execution.
#
# Component seeds are sim2_seed() in study2_design.R.
# Parallel orchestration uses future.seed = NULL and those explicit component seeds.
source_sim2_functions <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  source(file.path(root, "R", "study1_functions.R"), local = FALSE)
  source(file.path(root, "R", "study2_design.R"), local = FALSE)
  source(file.path(root, "R", "study2_functions.R"), local = FALSE)
  source(file.path(root, "R", "estimators.R"), local = FALSE)
  invisible(TRUE)
}

# Hash authoritative Study 2 design fixtures for chunk metadata.
sim2_design_hash <- function(root = NULL) {
  dir <- sim2_fixture_dir(root)
  files <- sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
  if (!length(files)) {
    stop("No Study 2 design fixtures found under ", dir, call. = FALSE)
  }
  payload <- paste(
    vapply(files, function(f) {
      paste(basename(f), paste(readLines(f, warn = FALSE), collapse = "\n"), sep = "=")
    }, character(1)),
    collapse = "\n"
  )
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(payload, algo = "sha256")
  } else {
    sprintf("raw-%d", abs(as.integer(sum(utf8ToInt(payload) %% 1000003L))))
  }
}

# Best-effort git commit hash (non-fatal).
sim2_git_commit <- function(root = NULL) {
  root <- sim2_project_root(root)
  git_dir <- file.path(root, ".git")
  if (!dir.exists(git_dir)) {
    return(NA_character_)
  }
  tryCatch(
    {
      out <- system2(
        "git",
        c("-C", root, "rev-parse", "HEAD"),
        stdout = TRUE,
        stderr = FALSE
      )
      if (length(out)) out[[1]] else NA_character_
    },
    error = function(e) NA_character_
  )
}

chunk_rds_path_s2 <- function(out_dir, condition_set, condition_id, chunk_id) {
  file.path(
    out_dir,
    "raw_chunks",
    as.character(condition_set),
    sprintf(
      "condition_%03d_chunk_%02d.rds",
      as.integer(condition_id),
      as.integer(chunk_id)
    )
  )
}

make_chunk_grid_s2 <- function(condition_grid,
                               n_rep_total,
                               chunk_size) {
  n_rep_total <- as.integer(n_rep_total)
  chunk_size <- as.integer(chunk_size)
  if (n_rep_total < 1L || chunk_size < 1L) {
    stop("n_rep_total and chunk_size must be >= 1.", call. = FALSE)
  }
  if (n_rep_total %% chunk_size != 0L) {
    stop("n_rep_total must be divisible by chunk_size.", call. = FALSE)
  }
  required <- c(
    "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
    "n_imbalance", "beta1_true", "allocation_pattern"
  )
  missing <- setdiff(required, names(condition_grid))
  if (length(missing)) {
    stop(
      "condition_grid missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
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
        condition_id = as.integer(cond$condition_id),
        condition_set = as.character(cond$condition_set),
        chunk_id = as.integer(ch),
        J = as.integer(cond$J),
        k_imbalance = as.character(cond$k_imbalance),
        moderator_split = as.character(cond$moderator_split),
        n_imbalance = as.character(cond$n_imbalance),
        beta1_true = as.numeric(cond$beta1_true),
        allocation_pattern = as.character(cond$allocation_pattern),
        rep_start = rep_start,
        rep_end = rep_end,
        n_reps = as.integer(chunk_size),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

condition_from_chunk_row_s2 <- function(chunk) {
  chunk <- as.list(chunk)
  data.frame(
    condition_id = as.integer(chunk$condition_id),
    condition_set = as.character(chunk$condition_set),
    J = as.integer(chunk$J),
    k_imbalance = as.character(chunk$k_imbalance),
    moderator_split = as.character(chunk$moderator_split),
    n_imbalance = as.character(chunk$n_imbalance),
    beta1_true = as.numeric(chunk$beta1_true),
    allocation_pattern = as.character(chunk$allocation_pattern),
    stringsAsFactors = FALSE
  )
}

validate_chunk_object_s2 <- function(obj, chunk) {
  chunk <- as.list(chunk)
  if (!is.list(obj) ||
      !all(c("metadata", "method_results", "dataset_diagnostics") %in% names(obj))) {
    return(list(ok = FALSE, reason = "Chunk object lacks required top-level elements."))
  }
  meta <- obj$metadata
  if (is.data.frame(meta)) {
    if (nrow(meta) != 1L) {
      return(list(ok = FALSE, reason = "metadata must be one row."))
    }
    meta <- as.list(meta[1, , drop = FALSE])
  }
  if (!is.list(meta)) {
    return(list(ok = FALSE, reason = "metadata must be a list or one-row data.frame."))
  }

  required_meta <- c(
    "condition_id", "chunk_id", "rep_start", "rep_end", "n_reps"
  )
  missing_meta <- setdiff(required_meta, names(meta))
  if (length(missing_meta)) {
    return(list(
      ok = FALSE,
      reason = paste("metadata missing:", paste(missing_meta, collapse = ", "))
    ))
  }

  if (!identical(as.integer(meta$condition_id), as.integer(chunk$condition_id)) ||
      !identical(as.integer(meta$chunk_id), as.integer(chunk$chunk_id)) ||
      !identical(as.integer(meta$rep_start), as.integer(chunk$rep_start)) ||
      !identical(as.integer(meta$rep_end), as.integer(chunk$rep_end))) {
    return(list(ok = FALSE, reason = "metadata identifiers do not match chunk job."))
  }

  n_reps <- as.integer(chunk$n_reps)
  mr <- obj$method_results
  dd <- obj$dataset_diagnostics
  if (!is.data.frame(mr) || !is.data.frame(dd)) {
    return(list(ok = FALSE, reason = "method_results/diagnostics must be data.frames."))
  }
  if (nrow(mr) != n_reps * 4L) {
    return(list(
      ok = FALSE,
      reason = sprintf(
        "method_results has %d rows; expected %d.",
        nrow(mr), n_reps * 4L
      )
    ))
  }
  if (nrow(dd) != n_reps) {
    return(list(
      ok = FALSE,
      reason = sprintf(
        "dataset_diagnostics has %d rows; expected %d.",
        nrow(dd), n_reps
      )
    ))
  }
  if (anyDuplicated(mr[, c("condition_id", "rep_id", "method")]) > 0L) {
    return(list(ok = FALSE, reason = "duplicate method-result keys."))
  }
  if (anyDuplicated(dd[, c("condition_id", "rep_id")]) > 0L) {
    return(list(ok = FALSE, reason = "duplicate diagnostic keys."))
  }
  list(ok = TRUE, reason = NA_character_)
}

# Atomic chunk write.
atomic_save_rds_s2 <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  saveRDS(object, tmp)
  # On some platforms rename fails if destination exists.
  if (file.exists(path)) {
    unlink(path)
  }
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok) || !file.exists(path)) {
    # Fallback copy+unlink if rename fails across devices.
    if (file.exists(tmp)) {
      file.copy(tmp, path, overwrite = TRUE)
      unlink(tmp)
    }
  }
  if (!file.exists(path)) {
    stop("Failed atomic save for ", path, call. = FALSE)
  }
  invisible(path)
}

# Run one chunk of replications with checkpoint/resume semantics.
run_one_chunk_s2 <- function(chunk,
                             root,
                             out_dir,
                             config,
                             design_hash = NULL) {
  chunk <- as.list(chunk)
  out_file <- chunk_rds_path_s2(
    out_dir = out_dir,
    condition_set = chunk$condition_set,
    condition_id = chunk$condition_id,
    chunk_id = chunk$chunk_id
  )

  if (!exists("run_one_replication_s2", mode = "function") ||
      !exists("make_chunk_grid_s2", mode = "function")) {
    source_sim2_functions(root)
  }
  if (is.null(design_hash)) {
    design_hash <- sim2_design_hash(root)
  }

  if (file.exists(out_file)) {
    if (!isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      stop(
        "Chunk file exists and resume/overwrite are FALSE: ", out_file,
        call. = FALSE
      )
    }
    if (isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      obj <- tryCatch(readRDS(out_file), error = function(e) e)
      if (!inherits(obj, "error")) {
        chk <- validate_chunk_object_s2(obj, chunk)
        if (isTRUE(chk$ok)) {
          log_row <- data.frame(
            condition_id = as.integer(chunk$condition_id),
            chunk_id = as.integer(chunk$chunk_id),
            condition_set = as.character(chunk$condition_set),
            status = "skipped",
            started_at = NA_character_,
            completed_at = NA_character_,
            elapsed_sec = NA_real_,
            n_method_rows = nrow(obj$method_results),
            n_diagnostic_rows = nrow(obj$dataset_diagnostics),
            n_errors = sum(obj$method_results$status != "success", na.rm = TRUE),
            n_warnings = sum(obj$method_results$warning_count > 0, na.rm = TRUE),
            file_path = out_file,
            file_size = file.info(out_file)$size,
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
  condition <- condition_from_chunk_row_s2(chunk)

  result <- tryCatch(
    {
      reps <- seq.int(as.integer(chunk$rep_start), as.integer(chunk$rep_end))
      method_list <- vector("list", length(reps))
      diag_list <- vector("list", length(reps))
      for (i in seq_along(reps)) {
        one <- run_one_replication_s2(
          condition = condition,
          rep_id = reps[[i]],
          config = config,
          return_data = FALSE
        )
        method_list[[i]] <- one$method_results
        diag_list[[i]] <- one$dataset_diagnostics
      }
      list(
        method_results = do.call(rbind, method_list),
        dataset_diagnostics = do.call(rbind, diag_list)
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
    # Record the failed chunk window so the summary can continue.
    methods <- c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve")
    expand <- expand.grid(
      rep_id = seq.int(chunk$rep_start, chunk$rep_end),
      method = methods,
      stringsAsFactors = FALSE
    )
    method_results <- data.frame(
      condition_id = as.integer(chunk$condition_id),
      condition_set = as.character(chunk$condition_set),
      rep_id = as.integer(expand$rep_id),
      rep_seed_base = NA_integer_,
      method = expand$method,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      converged = FALSE,
      status = "error",
      error = err_msg,
      warning_count = 0L,
      warnings = NA_character_,
      fit_runtime_sec = NA_real_,
      variance_study = NA_real_,
      variance_effect = NA_real_,
      variance_residual = NA_real_,
      boundary_study = NA,
      boundary_effect = NA,
      boundary_residual = NA,
      stringsAsFactors = FALSE
    )
    dataset_diagnostics <- data.frame(
      condition_id = as.integer(chunk$condition_id),
      condition_set = as.character(chunk$condition_set),
      rep_id = seq.int(chunk$rep_start, chunk$rep_end),
      rep_seed_base = NA_integer_,
      chunk_error = err_msg,
      stringsAsFactors = FALSE
    )
    chunk_status <- "failed"
  } else {
    method_results <- result$method_results
    dataset_diagnostics <- result$dataset_diagnostics
    rownames(method_results) <- NULL
    rownames(dataset_diagnostics) <- NULL
    chunk_status <- "completed"
  }

  metadata <- list(
    condition_id = as.integer(chunk$condition_id),
    condition_set = as.character(chunk$condition_set),
    chunk_id = as.integer(chunk$chunk_id),
    rep_start = as.integer(chunk$rep_start),
    rep_end = as.integer(chunk$rep_end),
    n_reps = as.integer(chunk$n_reps),
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S"),
    completed_at = format(completed_at, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = elapsed_sec,
    host = as.character(Sys.info()[["nodename"]]),
    pid = as.integer(Sys.getpid()),
    worker = NA_character_,
    git_commit = sim2_git_commit(root),
    design_hash = design_hash
  )

  obj <- list(
    metadata = metadata,
    method_results = method_results,
    dataset_diagnostics = dataset_diagnostics
  )

  save_ok <- tryCatch(
    {
      atomic_save_rds_s2(obj, out_file)
      TRUE
    },
    error = function(e) {
      err_msg <<- paste(
        na.omit(c(err_msg, paste("atomic_save_rds_s2 failed:", conditionMessage(e)))),
        collapse = "; "
      )
      FALSE
    }
  )
  if (!save_ok && identical(chunk_status, "completed")) {
    chunk_status <- "failed"
  }

  log_row <- data.frame(
    condition_id = as.integer(chunk$condition_id),
    chunk_id = as.integer(chunk$chunk_id),
    condition_set = as.character(chunk$condition_set),
    status = chunk_status,
    started_at = metadata$started_at,
    completed_at = metadata$completed_at,
    elapsed_sec = elapsed_sec,
    n_method_rows = nrow(method_results),
    n_diagnostic_rows = nrow(dataset_diagnostics),
    n_errors = sum(method_results$status != "success", na.rm = TRUE),
    n_warnings = sum(method_results$warning_count > 0, na.rm = TRUE),
    file_path = out_file,
    file_size = if (file.exists(out_file)) file.info(out_file)$size else NA_real_,
    error_message = err_msg,
    stringsAsFactors = FALSE
  )

  list(status = chunk_status, log = log_row, object = obj)
}

run_chunks_s2 <- function(chunk_grid,
                          root,
                          out_dir,
                          config) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(out_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)
  design_hash <- sim2_design_hash(root)
  workers <- as.integer(config$workers)
  if (is.na(workers) || workers < 1L) {
    stop("config$workers must be a positive integer for run_chunks_s2().", call. = FALSE)
  }

  wall_start <- Sys.time()

  if (workers == 1L) {
    results <- lapply(seq_len(nrow(chunk_grid)), function(i) {
      run_one_chunk_s2(
        chunk = chunk_grid[i, , drop = FALSE],
        root = root,
        out_dir = out_dir,
        config = config,
        design_hash = design_hash
      )
    })
  } else {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages future and future.apply are required for workers > 1.", call. = FALSE)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)

    results <- future.apply::future_lapply(
      X = seq_len(nrow(chunk_grid)),
      FUN = function(i) {
        # Fresh workers need Study 2 functions.
        source_sim2_functions(root)
        run_one_chunk_s2(
          chunk = chunk_grid[i, , drop = FALSE],
          root = root,
          out_dir = out_dir,
          config = config,
          design_hash = design_hash
        )
      },
      # Seeds are set inside run_one_replication_s2() via withr::with_seed().
      # Do not add a parallel RNG stream.
      future.seed = NULL,
      future.packages = c("metafor", "clubSandwich", "MASS", "withr")
    )
  }

  wall_end <- Sys.time()
  chunk_log <- do.call(rbind, lapply(results, function(x) x$log))
  rownames(chunk_log) <- NULL

  list(
    chunk_log = chunk_log,
    wall_elapsed_sec = as.numeric(difftime(wall_end, wall_start, units = "secs")),
    wall_start = wall_start,
    wall_end = wall_end,
    results = results
  )
}

collect_sim2_chunk_results <- function(out_dir, chunk_grid) {
  method_list <- vector("list", nrow(chunk_grid))
  diag_list <- vector("list", nrow(chunk_grid))
  for (i in seq_len(nrow(chunk_grid))) {
    ch <- chunk_grid[i, , drop = FALSE]
    path <- chunk_rds_path_s2(
      out_dir,
      ch$condition_set[[1]],
      ch$condition_id[[1]],
      ch$chunk_id[[1]]
    )
    if (!file.exists(path)) {
      stop("Missing chunk file: ", path, call. = FALSE)
    }
    obj <- readRDS(path)
    chk <- validate_chunk_object_s2(obj, ch)
    if (!isTRUE(chk$ok)) {
      stop("Invalid chunk file ", path, ": ", chk$reason, call. = FALSE)
    }
    method_list[[i]] <- obj$method_results
    diag_list[[i]] <- obj$dataset_diagnostics
  }
  methods <- do.call(rbind, method_list)
  diags <- do.call(rbind, diag_list)
  rownames(methods) <- NULL
  rownames(diags) <- NULL
  list(method_results = methods, dataset_diagnostics = diags)
}

# Focused tables retain the design factors used in the manuscript.
sim2_table_design_cols <- function() {
  c(
    "condition_id", "J", "k_imbalance", "moderator_split", "n_imbalance",
    "beta1_true", "allocation_pattern", "method"
  )
}

# Select columns that exist on a data frame, preserving order.
sim2_select_cols <- function(dat, cols) {
  dat[, intersect(cols, names(dat)), drop = FALSE]
}

write_sim2_tables <- function(method_summary,
                              design_summary,
                              alignment_summary,
                              condition_grid,
                              tables_dir,
                              config,
                              figure_data_dir = NULL,
                              condition_runtime = NULL) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(figure_data_dir)) {
    dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
  }

  design_cols <- sim2_table_design_cols()

  table_design <- condition_grid
  table_n_templates <- summarize_nj_templates_s2(config)

  cmp <- sim2_select_cols(method_summary, c(
    design_cols,
    "coverage", "type1_error", "power", "mean_ci_width", "mean_df",
    "bias", "se_ratio", "error_rate", "warning_rate"
  ))

  type1 <- method_summary[abs(as.numeric(method_summary$beta1_true)) < 1e-12, , drop = FALSE]
  type1_out <- sim2_select_cols(type1, c(
    design_cols,
    "type1_error", "mcse_reject_rate", "coverage", "mcse_coverage",
    "mean_df", "median_df", "min_df", "max_df", "n_success", "error_rate"
  ))

  power <- method_summary[
    abs(as.numeric(method_summary$beta1_true) - 0.10) < 1e-12, ,
    drop = FALSE
  ]
  power_out <- sim2_select_cols(power, c(
    design_cols,
    "power", "mcse_reject_rate", "coverage", "mcse_coverage",
    "mean_df", "median_df", "min_df", "max_df", "n_success", "error_rate"
  ))

  cov_out <- sim2_select_cols(method_summary, c(
    design_cols, "coverage", "mcse_coverage", "mean_ci_width", "n_success"
  ))

  se_out <- sim2_select_cols(method_summary, c(
    design_cols, "mean_se", "empirical_sd", "se_ratio", "n_success"
  ))

  df_out <- sim2_select_cols(method_summary, c(
    design_cols,
    "mean_df", "median_df", "sd_df", "min_df", "q05_df", "q10_df", "max_df",
    "prop_df_lt_2", "prop_df_lt_4", "prop_df_lt_10",
    "power", "type1_error", "coverage", "n_success"
  ))

  fail_out <- sim2_select_cols(method_summary, c(
    design_cols,
    "n_attempted", "n_success", "n_error", "n_invalid",
    "error_rate", "warning_rate",
    "boundary_study_rate", "boundary_effect_rate", "boundary_residual_rate"
  ))

  align_out <- NULL
  if (!is.null(alignment_summary) && nrow(alignment_summary)) {
    align_out <- sim2_select_cols(alignment_summary, c(
      design_cols, "stress_cell_id", "reuses_core_condition", "condition_set",
      "coverage", "type1_error", "power", "mean_ci_width", "mean_df",
      "se_ratio", "bias", "error_rate", "warning_rate",
      "boundary_study_rate", "boundary_effect_rate", "boundary_residual_rate",
      "n_success"
    ))
  }

  runtime_out <- condition_runtime

  utils::write.csv(table_design, file.path(tables_dir, "table_design.csv"),
                   row.names = FALSE)
  utils::write.csv(table_n_templates, file.path(tables_dir, "table_n_templates.csv"),
                   row.names = FALSE)
  utils::write.csv(cmp, file.path(tables_dir, "table_method_comparison.csv"),
                   row.names = FALSE)
  utils::write.csv(type1_out, file.path(tables_dir, "table_type1_error.csv"),
                   row.names = FALSE)
  utils::write.csv(cov_out, file.path(tables_dir, "table_coverage.csv"),
                   row.names = FALSE)
  utils::write.csv(se_out, file.path(tables_dir, "table_se_calibration.csv"),
                   row.names = FALSE)
  utils::write.csv(power_out, file.path(tables_dir, "table_power.csv"),
                   row.names = FALSE)
  utils::write.csv(df_out, file.path(tables_dir, "table_df.csv"),
                   row.names = FALSE)
  utils::write.csv(fail_out, file.path(tables_dir, "table_failures.csv"),
                   row.names = FALSE)
  if (!is.null(align_out)) {
    utils::write.csv(align_out, file.path(tables_dir, "table_alignment.csv"),
                     row.names = FALSE)
  }
  if (!is.null(runtime_out) && nrow(runtime_out)) {
    utils::write.csv(runtime_out, file.path(tables_dir, "table_runtime.csv"),
                     row.names = FALSE)
  }

  tables <- list(
    design = table_design,
    n_templates = table_n_templates,
    method_comparison = cmp,
    type1 = type1_out,
    coverage = cov_out,
    se_calibration = se_out,
    power = power_out,
    df = df_out,
    failures = fail_out
  )
  if (!is.null(align_out)) tables$alignment <- align_out
  if (!is.null(runtime_out) && nrow(runtime_out)) tables$runtime <- runtime_out

  if (!is.null(figure_data_dir)) {
    for (nm in names(tables)) {
      utils::write.csv(
        tables[[nm]],
        file.path(figure_data_dir, paste0(nm, ".csv")),
        row.names = FALSE
      )
    }
    utils::write.csv(
      method_summary,
      file.path(figure_data_dir, "full_summary.csv"),
      row.names = FALSE
    )
    if (!is.null(design_summary) && nrow(design_summary)) {
      utils::write.csv(
        design_summary,
        file.path(figure_data_dir, "design_summary.csv"),
        row.names = FALSE
      )
    }
  }

  invisible(tables)
}
