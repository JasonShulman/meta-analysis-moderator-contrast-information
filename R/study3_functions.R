# Simulation Study 3 outcome simulation, summaries, and heterogeneity sensitivity.
#
# The published run reads data/study3_design/. It does not construct or freeze scenarios.
# Outcome seeds use seed_root, scenario_stride, and replication_stride from make_sim3_config().
# Estimators are fit by estimators.R. Contrast metrics are computed by contrast_information.R.
# Common random numbers: stochastic components are drawn once per scenario-replication,
# then the beta1 datasets differ only by beta1 * X_j.

# =============================================================================
# Paths
# =============================================================================

sim3_project_root <- function(root = NULL) {
  if (!is.null(root)) {
    return(normalizePath(root, winslash = "/", mustWork = TRUE))
  }
  candidates <- c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../.."))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "study3_functions.R"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot locate reproducibility root (expected R/study3_functions.R).", call. = FALSE)
}

sim3_fixture_dir <- function(root = NULL) {
  file.path(sim3_project_root(root), "data", "study3_design")
}

# =============================================================================
# Configuration and seeds
# =============================================================================

sim3_stream_codes <- function() {
  c(u_j = 1L, w_ij = 2L, e_ij = 3L, validation = 4L)
}

make_sim3_config <- function(run_mode = c("smoke", "pilot", "production"),
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
    study_id = 3L,
    run_mode = run_mode,
    beta0 = 0.20,
    beta1_levels = c(0.00, 0.05, 0.10, 0.20),
    alpha = 0.05,
    ci_level = 0.95,
    rho_plan_primary = 0.50,
    rho_assumed_aggregation = 0.50,
    n_min = 20L,
    n_max = 3000L,
    k_min = 1L,
    k_max = 30L,
    J_levels = c(24L, 30L, 36L, 48L, 60L, 72L, 96L, 120L),
    p1_range = c(0.08, 0.45),
    mean_n_range = c(80, 800),
    sigma_log_n_range = c(0, 1.4),
    mean_k_range = c(1.5, 8),
    sigma_log_k_range = c(0, 1.2),
    lambda_range = c(0.35, 1),
    rho_true_range = c(0.10, 0.90),
    tau_range = c(0.03, 0.25),
    omega_range = c(0.03, 0.25),
    nk_modes = c("negative", "independent", "positive"),
    alignment_patterns = c(
      "independent", "n_majority", "n_minority", "nk_majority", "nk_minority"
    ),
    n_broad = 80L,
    n_targeted = 20L,
    n_candidate_pool = 50000L,
    n_rep_total = as.integer(mode_defaults$n_rep_total),
    chunk_size = as.integer(mode_defaults$chunk_size),
    workers = workers_out,
    seed_root = 900000000L,
    scenario_stride = 1000000L,
    replication_stride = 100L,
    # Construction seeds for the frozen library, not outcome-simulation seeds.
    scenario_seed_root = 81032026L,
    candidate_seed_root = 81032999L,
    warning_separator = " || ",
    boundary_tolerance = 1e-8,
    pair_ic_div_min = 0.30,
    pair_class_c_ic_tol = 0.02,
    pair_class_d_nu_tol = 0.05,
    pair_class_e_conc_tol = 0.01,
    pair_class_b_share_min = 0.10,
    overwrite = isTRUE(overwrite),
    resume = isTRUE(resume)
  )
}

sim3_seed <- function(scenario_id,
                      rep_id,
                      stream_id,
                      config = make_sim3_config("smoke")) {
  scenario_id <- as.integer(scenario_id)
  rep_id <- as.integer(rep_id)
  if (is.character(stream_id)) {
    codes <- sim3_stream_codes()
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
  as.integer(
    as.numeric(config$seed_root) +
      as.numeric(scenario_id) * as.numeric(config$scenario_stride) +
      as.numeric(rep_id) * as.numeric(config$replication_stride) +
      as.numeric(stream_id)
  )
}

run_with_sim3_seed <- function(seed, code) {
  if (exists("run_with_sim2_seed", mode = "function")) {
    return(run_with_sim2_seed(seed, code))
  }
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop("Package 'withr' is required for Study 3 seeding.", call. = FALSE)
  }
  withr::with_seed(
    seed,
    code,
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
}

primary_metrics_s3 <- function(study_structure, rho_plan = 0.50) {
  compute_contrast_information_s3(
    study_structure = study_structure,
    rho_plan = rho_plan,
    weight_type = "sampling"
  )
}

pair_satisfies_class_s3 <- function(diag_row, class_label, config) {
  d <- as.list(diag_row)
  if (any(c(d$status_a, d$status_b) %in% c("insufficient_group_information", "numerical_failure"))) {
    return(FALSE)
  }
  ic <- as.numeric(d$I_C_rel_diff)
  nu <- as.numeric(d$nu_rel_diff)
  j1s <- as.numeric(d$J1_star_rel_diff)
  if (identical(class_label, "A")) {
    same_marg <- isTRUE(d$J_a == d$J_b) &&
      isTRUE(d$N_a == d$N_b) &&
      isTRUE(d$K_a == d$K_b) &&
      isTRUE(d$J1_a == d$J1_b)
    return(isTRUE(same_marg) && isTRUE(ic >= config$pair_ic_div_min))
  }
  if (identical(class_label, "B")) {
    same_marg <- isTRUE(d$J_a == d$J_b) &&
      isTRUE(d$N_a == d$N_b) &&
      isTRUE(d$K_a == d$K_b) &&
      isTRUE(d$J1_a == d$J1_b)
    return(
      isTRUE(same_marg) &&
        isTRUE(ic >= config$pair_ic_div_min) &&
        isTRUE(as.numeric(d$minority_share_abs_diff) >= config$pair_class_b_share_min)
    )
  }
  if (identical(class_label, "C")) {
    return(
      isTRUE(ic <= config$pair_class_c_ic_tol) &&
        (isTRUE(j1s >= config$pair_ic_div_min) || isTRUE(nu >= config$pair_ic_div_min))
    )
  }
  if (identical(class_label, "D")) {
    return(
      isTRUE(nu <= config$pair_class_d_nu_tol) &&
        isTRUE(ic >= config$pair_ic_div_min)
    )
  }
  if (identical(class_label, "E")) {
    same_counts <- isTRUE(d$J_a == d$J_b) &&
      isTRUE(d$N_a == d$N_b) &&
      isTRUE(d$K_a == d$K_b) &&
      isTRUE(d$J1_a == d$J1_b)
    return(
      isTRUE(same_counts) &&
        isTRUE(as.numeric(d$J_prec_rel_diff) <= config$pair_class_e_conc_tol) &&
        isTRUE(as.numeric(d$J_k_rel_diff) <= config$pair_class_e_conc_tol) &&
        isTRUE(ic >= config$pair_ic_div_min)
    )
  }
  FALSE
}

hash_file_s3 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(file = path, algo = "sha256")
  } else {
    paste(unname(tools::md5sum(path)))
  }
}

read_frozen_scenario_library_s3 <- function(root = NULL) {
  dir <- sim3_fixture_dir(root)
  required <- c(
    "scenario_grid_all.csv",
    "scenario_study_structure.csv",
    "scenario_metric_library.csv",
    "targeted_pair_diagnostics.csv"
  )
  missing <- required[!file.exists(file.path(dir, required))]
  if (length(missing)) {
    stop(
      "Frozen Study 3 scenario library is incomplete. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    dir = dir,
    grid = utils::read.csv(file.path(dir, "scenario_grid_all.csv"), stringsAsFactors = FALSE),
    study_structure = utils::read.csv(
      file.path(dir, "scenario_study_structure.csv"),
      stringsAsFactors = FALSE
    ),
    metrics = utils::read.csv(
      file.path(dir, "scenario_metric_library.csv"),
      stringsAsFactors = FALSE
    ),
    diagnostics = utils::read.csv(
      file.path(dir, "targeted_pair_diagnostics.csv"),
      stringsAsFactors = FALSE
    ),
    hashes = if (file.exists(file.path(dir, "scenario_hashes.csv"))) {
      utils::read.csv(file.path(dir, "scenario_hashes.csv"), stringsAsFactors = FALSE)
    } else {
      NULL
    }
  )
}

validate_frozen_scenario_library_s3 <- function(root = NULL,
                                                config = make_sim3_config("smoke")) {
  lib <- read_frozen_scenario_library_s3(root)
  issues <- character(0)
  add <- function(msg) issues <<- c(issues, msg)

  g <- lib$grid
  ss <- lib$study_structure
  if (nrow(g) != 100L) add(sprintf("Expected 100 scenarios, found %d.", nrow(g)))
  n_broad <- sum(g$scenario_set == "broad")
  n_targ <- sum(g$scenario_set == "targeted")
  if (n_broad != 80L) add(sprintf("Expected 80 broad scenarios, found %d.", n_broad))
  if (n_targ != 20L) add(sprintf("Expected 20 targeted scenarios, found %d.", n_targ))
  if (length(unique(g$scenario_id)) != nrow(g)) add("Duplicate scenario_id values.")

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
    if (any(ssi$n_j < config$n_min | ssi$n_j > config$n_max)) {
      add(sprintf("scenario %d: n_j outside bounds.", row$scenario_id))
    }
    if (any(ssi$k_j < config$k_min | ssi$k_j > config$k_max)) {
      add(sprintf("scenario %d: k_j outside bounds.", row$scenario_id))
    }
    if (identical(row$scenario_set, "broad")) {
      if (!row$J %in% config$J_levels) {
        add(sprintf("broad scenario %d: J not in allowed set.", row$scenario_id))
      }
      if (row$J1 < 3L) add(sprintf("broad scenario %d: J1 < 3.", row$scenario_id))
      if (row$J1 > floor(row$J / 2)) {
        add(sprintf("broad scenario %d: J1 > floor(J/2).", row$scenario_id))
      }
    }
    recomputed <- primary_metrics_s3(ssi)
    if (is.finite(row$I_C_star) && is.finite(recomputed$I_C_star) &&
        abs(row$I_C_star - recomputed$I_C_star) > 1e-8 * max(1, abs(row$I_C_star))) {
      add(sprintf("scenario %d: stored I_C* disagrees with recomputation.", row$scenario_id))
    }
  }

  d <- lib$diagnostics
  if (nrow(d) != 10L) add(sprintf("Expected 10 targeted pairs, found %d.", nrow(d)))
  for (i in seq_len(nrow(d))) {
    if (!pair_satisfies_class_s3(d[i, ], d$pair_class[[i]], config)) {
      add(sprintf(
        "Pair %s class %s fails matching/divergence tolerances.",
        d$pair_id[[i]], d$pair_class[[i]]
      ))
    }
  }

  list(ok = length(issues) == 0L, issues = issues, library = lib)
}

# Simulation Study 3 — DGM with frozen study structure and common random numbers.
# Structure is never regenerated inside the replication loop.

# Draw u_j, w_ij, and e_ij once per scenario/replication.
# beta1 is deliberately absent from the seed.
simulate_base_components_s3 <- function(study_structure,
                                        scenario,
                                        rep_id,
                                        config = make_sim3_config("smoke")) {
  ss <- study_structure[order(study_structure$study), , drop = FALSE]
  J <- nrow(ss)
  k_j <- as.integer(ss$k_j)
  n_j <- as.numeric(ss$n_j)
  scenario_id <- as.integer(scenario$scenario_id)
  tau <- as.numeric(scenario$tau)
  omega <- as.numeric(scenario$omega)
  rho <- as.numeric(scenario$rho_true)
  rep_id <- as.integer(rep_id)

  u_j <- run_with_sim3_seed(
    sim3_seed(scenario_id, rep_id, "u_j", config),
    stats::rnorm(J, mean = 0, sd = tau)
  )

  K <- sum(k_j)
  w_all <- run_with_sim3_seed(
    sim3_seed(scenario_id, rep_id, "w_ij", config),
    stats::rnorm(K, mean = 0, sd = omega)
  )

  e_list <- run_with_sim3_seed(
    sim3_seed(scenario_id, rep_id, "e_ij", config),
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
    u_j = as.numeric(u_j),
    w_all = as.numeric(w_all),
    e_list = e_list,
    study_structure = ss
  )
}

assemble_dataset_s3 <- function(components,
                                scenario,
                                beta1,
                                rep_id,
                                config = make_sim3_config("smoke")) {
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
      pair_id = if ("pair_id" %in% names(scenario)) as.integer(scenario$pair_id) else NA_integer_,
      pair_class = if ("pair_class" %in% names(scenario)) as.character(scenario$pair_class) else NA_character_,
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

# Stochastic components are drawn once. Datasets differ only by beta1 * X_j.
simulate_meta_datasets_crn_s3 <- function(scenario,
                                          study_structure,
                                          rep_id,
                                          config = make_sim3_config("smoke"),
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
  components <- simulate_base_components_s3(ss, scenario, rep_id, config)
  datasets <- lapply(beta1_levels, function(b) {
    assemble_dataset_s3(components, scenario, b, rep_id, config)
  })
  names(datasets) <- sprintf("beta1_%s", format(beta1_levels, trim = TRUE, scientific = FALSE))
  list(
    components = components,
    datasets = datasets,
    beta1_levels = as.numeric(beta1_levels)
  )
}

replication_diagnostics_s3 <- function(scenario, study_structure, rep_id, datasets) {
  ss <- study_structure[study_structure$scenario_id == scenario$scenario_id, , drop = FALSE]
  d0 <- datasets[[1]]
  data.frame(
    scenario_id = as.integer(scenario$scenario_id),
    scenario_set = as.character(scenario$scenario_set),
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

# Model fitting
# Reuses the Study 2 fitters. CHE uses scenario-specific rho; aggregation uses .50.

condition_adapter_s3 <- function(scenario, beta1_true) {
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  data.frame(
    condition_id = as.integer(scenario$scenario_id),
    condition_set = as.character(scenario$scenario_set),
    J = as.integer(scenario$J),
    k_imbalance = "study3",
    moderator_split = "study3",
    n_imbalance = "study3",
    beta1_true = as.numeric(beta1_true),
    allocation_pattern = if ("alignment_pattern" %in% names(scenario)) {
      as.character(scenario$alignment_pattern)
    } else {
      "independent"
    },
    stringsAsFactors = FALSE
  )
}

# Reuse the Study 2 fitter configuration with Study 3-specific constants.
fitting_config_s3 <- function(scenario, config_s3) {
  cfg <- make_sim2_config("smoke")
  cfg$beta0 <- config_s3$beta0
  cfg$alpha <- config_s3$alpha
  cfg$ci_level <- config_s3$ci_level
  cfg$rho_assumed_aggregation <- config_s3$rho_assumed_aggregation
  cfg$rho_assumed_che <- as.numeric(scenario$rho_true)
  cfg$rho_true <- as.numeric(scenario$rho_true)
  cfg$tau <- as.numeric(scenario$tau)
  cfg$omega <- as.numeric(scenario$omega)
  cfg$warning_separator <- config_s3$warning_separator
  cfg$boundary_tolerance <- config_s3$boundary_tolerance
  cfg
}

remap_method_result_s3 <- function(row, scenario, beta1_true) {
  data.frame(
    scenario_id = as.integer(scenario$scenario_id),
    scenario_set = as.character(scenario$scenario_set),
    pair_id = if ("pair_id" %in% names(scenario)) as.integer(scenario$pair_id) else NA_integer_,
    pair_class = if ("pair_class" %in% names(scenario)) as.character(scenario$pair_class) else NA_character_,
    rep_id = as.integer(row$rep_id),
    beta1_true = as.numeric(beta1_true),
    method = as.character(row$method),
    estimate = as.numeric(row$estimate),
    se = as.numeric(row$se),
    ci_lower = as.numeric(row$ci_lower),
    ci_upper = as.numeric(row$ci_upper),
    p_value = as.numeric(row$p_value),
    df = as.numeric(row$df),
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

# MLMA and MLMA+CR2 share the identical rma.mv fit via fit_all_methods_s2().
fit_all_methods_s3 <- function(data,
                               scenario,
                               rep_id,
                               beta1_true,
                               config_s3 = make_sim3_config("smoke")) {
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  cond <- condition_adapter_s3(scenario, beta1_true)
  cfg2 <- fitting_config_s3(scenario, config_s3)
  rows <- fit_all_methods_s2(
    data = data,
    condition = cond,
    rep_id = as.integer(rep_id),
    config = cfg2
  )
  out <- do.call(
    rbind,
    lapply(seq_len(nrow(rows)), function(i) {
      remap_method_result_s3(rows[i, , drop = FALSE], scenario, beta1_true)
    })
  )
  rownames(out) <- NULL
  if (nrow(out) != 4L) {
    stop("fit_all_methods_s3() must return exactly four rows.", call. = FALSE)
  }
  out
}

# Simulation Study 3 — replication runner, performance, and prediction error.

# Run one Study 3 replication: one stochastic draw, four beta1 datasets, four methods.
run_one_replication_s3 <- function(scenario,
                                   study_structure,
                                   rep_id,
                                   config = make_sim3_config("smoke"),
                                   return_data = FALSE) {
  if (is.data.frame(scenario)) {
    if (nrow(scenario) != 1L) {
      stop("scenario must be a one-row data frame or list.", call. = FALSE)
    }
    scenario <- as.list(scenario[1, , drop = FALSE])
  }
  rep_id <- as.integer(rep_id)
  sim <- tryCatch(
    simulate_meta_datasets_crn_s3(
      scenario = scenario,
      study_structure = study_structure,
      rep_id = rep_id,
      config = config
    ),
    error = function(e) e
  )

  methods <- c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve")
  beta1_levels <- config$beta1_levels

  if (inherits(sim, "error")) {
    err <- conditionMessage(sim)
    expand <- expand.grid(
      beta1_true = beta1_levels,
      method = methods,
      stringsAsFactors = FALSE
    )
    method_results <- data.frame(
      scenario_id = as.integer(scenario$scenario_id),
      scenario_set = as.character(scenario$scenario_set),
      pair_id = if ("pair_id" %in% names(scenario)) as.integer(scenario$pair_id) else NA_integer_,
      pair_class = if ("pair_class" %in% names(scenario)) as.character(scenario$pair_class) else NA_character_,
      rep_id = rep_id,
      beta1_true = expand$beta1_true,
      method = expand$method,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      converged = FALSE,
      status = "error",
      error = err,
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
    method_list[[i]] <- fit_all_methods_s3(
      data = sim$datasets[[i]],
      scenario = scenario,
      rep_id = rep_id,
      beta1_true = sim$beta1_levels[[i]],
      config_s3 = config
    )
  }
  method_results <- do.call(rbind, method_list)
  rownames(method_results) <- NULL
  diagnostics <- replication_diagnostics_s3(
    scenario, study_structure, rep_id, sim$datasets
  )
  out <- list(
    method_results = method_results,
    replication_diagnostics = diagnostics
  )
  if (isTRUE(return_data)) {
    out$datasets <- sim$datasets
  }
  out
}

summarize_method_performance_s3 <- function(results, config) {
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
    extra <- d[1, intersect(c("pair_id", "pair_class", "J"), names(d)), drop = FALSE]
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
        rmse <- coverage <- reject_rate <- mcse_reject <- mean_df <- median_df <-
        sd_df <- min_df <- max_df <- boundary_study_rate <- boundary_effect_rate <-
        mean_fit_runtime_sec <- NA_real_
    }
    type1_error <- if (isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    power <- if (!isTRUE(all.equal(beta1_true, 0))) reject_rate else NA_real_
    cbind(
      design,
      extra,
      data.frame(
        n_attempted = n_attempted,
        n_success = n_success,
        n_error = n_error,
        n_invalid = n_invalid,
        error_rate = (n_attempted - n_success) / n_attempted,
        mean_estimate = mean_estimate,
        bias = bias,
        empirical_variance = empirical_var,
        empirical_sd = empirical_sd,
        mean_se = mean_se,
        se_ratio = se_ratio,
        rmse = rmse,
        coverage = coverage,
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

ols_calibration_s3 <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (sum(ok) < 3L) {
    return(list(intercept = NA_real_, slope = NA_real_, r2 = NA_real_))
  }
  fit <- stats::lm(observed[ok] ~ predicted[ok])
  sm <- summary(fit)
  list(
    intercept = unname(stats::coef(fit)[[1]]),
    slope = unname(stats::coef(fit)[[2]]),
    r2 = sm$r.squared
  )
}

prediction_error_metrics_s3 <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) {
    return(data.frame(
      n = 0L,
      rmse = NA_real_,
      mae = NA_real_,
      mean_signed_error = NA_real_,
      median_ape = NA_real_,
      spearman = NA_real_,
      max_abs_error = NA_real_,
      calibration_intercept = NA_real_,
      calibration_slope = NA_real_,
      r2 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  e <- predicted[ok] - observed[ok]
  cal <- ols_calibration_s3(observed, predicted)
  data.frame(
    n = as.integer(sum(ok)),
    rmse = sqrt(mean(e^2)),
    mae = mean(abs(e)),
    mean_signed_error = mean(e),
    median_ape = stats::median(abs(e) / pmax(abs(observed[ok]), 1e-12)),
    spearman = if (sum(ok) >= 3L && stats::sd(observed[ok]) > 0 && stats::sd(predicted[ok]) > 0) {
      as.numeric(stats::cor(observed[ok], predicted[ok], method = "spearman"))
    } else {
      NA_real_
    },
    max_abs_error = max(abs(e)),
    calibration_intercept = cal$intercept,
    calibration_slope = cal$slope,
    r2 = cal$r2,
    stringsAsFactors = FALSE
  )
}

drop_overlapping_metric_cols_s3 <- function(left, right, keys = "scenario_id") {
  overlap <- setdiff(intersect(names(left), names(right)), keys)
  if (length(overlap)) {
    left <- left[, setdiff(names(left), overlap), drop = FALSE]
  }
  left
}

summarize_variance_prediction_s3 <- function(method_summary,
                                             metric_library,
                                             variant_id = "sampling_rho_0.50") {
  mlma0 <- method_summary[
    method_summary$method == "mlma" & abs(method_summary$beta1_true) < 1e-12,
    ,
    drop = FALSE
  ]
  met <- metric_library[metric_library$variant_id == variant_id, , drop = FALSE]
  mlma0 <- drop_overlapping_metric_cols_s3(mlma0, met)
  merged <- merge(mlma0, met, by = "scenario_id", suffixes = c("", "_metric"))
  merged$prediction_ratio <- merged$predicted_variance / merged$empirical_variance
  merged$abs_error <- abs(merged$predicted_variance - merged$empirical_variance)
  merged$rel_error <- (merged$predicted_variance - merged$empirical_variance) /
    merged$empirical_variance
  err <- prediction_error_metrics_s3(merged$empirical_variance, merged$predicted_variance)
  list(by_scenario = merged, overall = cbind(
    data.frame(
      variant_id = variant_id,
      method = NA_character_,
      target = "mlma_empirical_variance",
      stringsAsFactors = FALSE
    ),
    err
  ))
}

summarize_df_prediction_s3 <- function(method_summary,
                                       metric_library,
                                       variant_id = "sampling_rho_0.50",
                                       method = "mlma_cr2") {
  cr2_0 <- method_summary[
    method_summary$method == method & abs(method_summary$beta1_true) < 1e-12,
    ,
    drop = FALSE
  ]
  met <- metric_library[metric_library$variant_id == variant_id, , drop = FALSE]
  cr2_0 <- drop_overlapping_metric_cols_s3(cr2_0, met)
  merged <- merge(cr2_0, met, by = "scenario_id", suffixes = c("", "_metric"))
  err <- prediction_error_metrics_s3(merged$mean_df, merged$nu_C_star)
  list(by_scenario = merged, overall = cbind(
    data.frame(variant_id = variant_id, method = method, target = "mean_satterthwaite_df",
               stringsAsFactors = FALSE),
    err
  ))
}

summarize_power_prediction_s3 <- function(method_summary,
                                          metric_library,
                                          variant_id = "sampling_rho_0.50",
                                          method = "mlma_cr2") {
  pwr <- method_summary[
    method_summary$method == method & method_summary$beta1_true > 0,
    ,
    drop = FALSE
  ]
  met <- metric_library[metric_library$variant_id == variant_id, , drop = FALSE]
  pwr <- drop_overlapping_metric_cols_s3(pwr, met)
  merged <- merge(pwr, met, by = "scenario_id", suffixes = c("", "_metric"))
  merged$predicted_power <- predict_contrast_power_s3(
    beta1 = merged$beta1_true,
    predicted_se = merged$predicted_se,
    nu_C_star = merged$nu_C_star
  )
  overall <- prediction_error_metrics_s3(merged$power, merged$predicted_power)
  by_beta <- do.call(rbind, lapply(split(merged, merged$beta1_true), function(d) {
    cbind(
      data.frame(beta1_true = d$beta1_true[[1]], stringsAsFactors = FALSE),
      prediction_error_metrics_s3(d$power, d$predicted_power)
    )
  }))
  list(
    by_scenario = merged,
    overall = cbind(
      data.frame(variant_id = variant_id, method = method, target = "empirical_power",
                 stringsAsFactors = FALSE),
      overall
    ),
    by_beta1 = by_beta
  )
}

summarize_targeted_pairs_s3 <- function(method_summary, diagnostics) {
  if (is.null(diagnostics) || !nrow(diagnostics)) {
    return(diagnostics)
  }
  cr2 <- method_summary[method_summary$method == "mlma_cr2", , drop = FALSE]
  mlma0 <- method_summary[
    method_summary$method == "mlma" & abs(method_summary$beta1_true) < 1e-12,
    c("scenario_id", "empirical_variance", "mean_df"),
    drop = FALSE
  ]
  names(mlma0) <- c("scenario_id", "emp_var", "mean_df_mlma")
  cr2_0 <- cr2[abs(cr2$beta1_true) < 1e-12, c("scenario_id", "mean_df", "type1_error")]
  names(cr2_0) <- c("scenario_id", "mean_df_cr2", "type1_error")
  power <- cr2[cr2$beta1_true > 0, c("scenario_id", "beta1_true", "power")]
  rows <- vector("list", nrow(diagnostics))
  for (i in seq_len(nrow(diagnostics))) {
    d <- diagnostics[i, , drop = FALSE]
    a <- d$scenario_id_a
    b <- d$scenario_id_b
    va <- mlma0[mlma0$scenario_id == a, , drop = FALSE]
    vb <- mlma0[mlma0$scenario_id == b, , drop = FALSE]
    da <- cr2_0[cr2_0$scenario_id == a, , drop = FALSE]
    db <- cr2_0[cr2_0$scenario_id == b, , drop = FALSE]
    pa <- power[power$scenario_id == a, , drop = FALSE]
    pb <- power[power$scenario_id == b, , drop = FALSE]
    rows[[i]] <- cbind(
      d,
      data.frame(
        emp_var_a = if (nrow(va)) va$emp_var[[1]] else NA_real_,
        emp_var_b = if (nrow(vb)) vb$emp_var[[1]] else NA_real_,
        mean_df_cr2_a = if (nrow(da)) da$mean_df_cr2[[1]] else NA_real_,
        mean_df_cr2_b = if (nrow(db)) db$mean_df_cr2[[1]] else NA_real_,
        type1_a = if (nrow(da)) da$type1_error[[1]] else NA_real_,
        type1_b = if (nrow(db)) db$type1_error[[1]] else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Prespecified metric variants for Study 3 theoretical-power reconstruction.

# Includes the primary design-only metric, the default heterogeneity-adjusted
# metric (also present in the planning grid), the full tau=omega x rho grid,
# oracle design-only rho, and the oracle heterogeneity-adjusted benchmark.
sim3_power_metric_variant_ids <- function() {
  grid <- expand.grid(
    het = c(0.05, 0.10, 0.20),
    rho = c(0.20, 0.50, 0.80),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  c(
    "sampling_rho_0.50",
    sprintf("het_tau%.2f_omega%.2f_rho%.2f", grid$het, grid$het, grid$rho),
    "sampling_rho_oracle",
    "het_oracle"
  )
}

label_metric_variant_s3 <- function(variant_id,
                                    is_primary = FALSE,
                                    is_default_het = FALSE,
                                    is_oracle_sampling = FALSE,
                                    is_oracle_het = FALSE) {
  variant_id <- as.character(variant_id)
  is_primary <- as.logical(is_primary)
  is_default_het <- as.logical(is_default_het)
  is_oracle_sampling <- as.logical(is_oracle_sampling)
  is_oracle_het <- as.logical(is_oracle_het)
  out <- variant_id
  out[is_primary %in% TRUE] <- "primary_design_only"
  out[is_default_het %in% TRUE] <- "default_heterogeneity_adjusted"
  out[is_oracle_sampling %in% TRUE] <- "oracle_design_only"
  out[is_oracle_het %in% TRUE] <- "oracle_heterogeneity_adjusted"
  out
}

empirical_power_bin_s3 <- function(power) {
  power <- as.numeric(power)
  out <- rep(NA_character_, length(power))
  ok <- is.finite(power)
  out[ok & power < 0.20] <- "<0.20"
  out[ok & power >= 0.20 & power < 0.40] <- "0.20-0.40"
  out[ok & power >= 0.40 & power < 0.60] <- "0.40-0.60"
  out[ok & power >= 0.60 & power <= 0.80] <- "0.60-0.80"
  out[ok & power > 0.80] <- ">0.80"
  factor(
    out,
    levels = c("<0.20", "0.20-0.40", "0.40-0.60", "0.60-0.80", ">0.80")
  )
}

abs_misspec_bin_s3 <- function(x, breaks, labels) {
  x <- as.numeric(x)
  out <- cut(
    x,
    breaks = breaks,
    labels = labels,
    right = FALSE,
    include.lowest = TRUE
  )
  as.character(out)
}

power_error_summary_row_s3 <- function(observed, predicted, meta) {
  err <- prediction_error_metrics_s3(observed, predicted)
  cbind(
    meta,
    data.frame(
      n_rows = err$n,
      rmse_power = err$rmse,
      mae_power = err$mae,
      mean_signed_error = err$mean_signed_error,
      spearman = err$spearman,
      calibration_intercept = err$calibration_intercept,
      calibration_slope = err$calibration_slope,
      max_absolute_error = err$max_abs_error,
      stringsAsFactors = FALSE
    )
  )
}

variant_planning_meta_s3 <- function(met) {
  oracle_rho <- isTRUE(as.logical(met$is_oracle_sampling[[1]])) ||
    isTRUE(as.logical(met$is_oracle_het[[1]]))
  oracle_het <- isTRUE(as.logical(met$is_oracle_het[[1]]))
  sampling <- identical(as.character(met$weight_type[[1]]), "sampling")
  data.frame(
    variant_id = as.character(met$variant_id[[1]]),
    metric_variant = label_metric_variant_s3(
      met$variant_id[[1]],
      is_primary = met$is_primary[[1]],
      is_default_het = met$is_default_het[[1]],
      is_oracle_sampling = met$is_oracle_sampling[[1]],
      is_oracle_het = met$is_oracle_het[[1]]
    ),
    weight_type = as.character(met$weight_type[[1]]),
    rho_plan = if (oracle_rho) NA_real_ else as.numeric(met$rho_plan[[1]]),
    tau_plan = if (oracle_het || sampling) NA_real_ else as.numeric(met$tau_plan[[1]]),
    omega_plan = if (oracle_het || sampling) NA_real_ else as.numeric(met$omega_plan[[1]]),
    oracle_rho = oracle_rho,
    oracle_heterogeneity = oracle_het,
    stringsAsFactors = FALSE
  )
}

# Merge a frozen metric variant onto empirical MLMA+CR2 power; no calibration is fit.
merge_variant_power_s3 <- function(method_summary,
                                   metric_library,
                                   variant_id,
                                   method = "mlma_cr2") {
  pwr <- method_summary[
    method_summary$method == method & method_summary$beta1_true > 0,
    ,
    drop = FALSE
  ]
  met <- metric_library[metric_library$variant_id == variant_id, , drop = FALSE]
  if (!nrow(met)) {
    stop("Unknown metric variant_id: ", variant_id, call. = FALSE)
  }
  pwr <- drop_overlapping_metric_cols_s3(pwr, met)
  merged <- merge(pwr, met, by = "scenario_id", suffixes = c("", "_metric"))
  merged$predicted_power <- predict_contrast_power_s3(
    beta1 = merged$beta1_true,
    predicted_se = merged$predicted_se,
    nu_C_star = merged$nu_C_star
  )
  merged$power_error <- merged$predicted_power - merged$power
  merged$tau_plan_eval <- ifelse(is.na(merged$tau_plan), 0, merged$tau_plan)
  merged$omega_plan_eval <- ifelse(is.na(merged$omega_plan), 0, merged$omega_plan)
  merged$abs_rho_misspec <- abs(merged$rho_plan - merged$rho_true)
  merged$abs_tau_misspec <- abs(merged$tau_plan_eval - merged$tau)
  merged$abs_omega_misspec <- abs(merged$omega_plan_eval - merged$omega)
  merged$signed_rho_misspec <- merged$rho_plan - merged$rho_true
  merged$signed_tau_misspec <- merged$tau_plan_eval - merged$tau
  merged$signed_omega_misspec <- merged$omega_plan_eval - merged$omega
  merged$empirical_power_bin <- empirical_power_bin_s3(merged$power)
  merged$metric_variant <- label_metric_variant_s3(
    merged$variant_id,
    is_primary = merged$is_primary,
    is_default_het = merged$is_default_het,
    is_oracle_sampling = merged$is_oracle_sampling,
    is_oracle_het = merged$is_oracle_het
  )
  merged
}

# Theoretical noncentral-t power vs empirical MLMA+CR2 power for metric variants.

# Reconstructs predictions from frozen metric-library SE/df and existing
# scenario-level empirical power. No new simulation and no fitted correction.
summarize_metric_variant_power_s3 <- function(method_summary,
                                              metric_library,
                                              variant_ids = sim3_power_metric_variant_ids(),
                                              method = "mlma_cr2") {
  variant_ids <- unique(as.character(variant_ids))
  detail_list <- vector("list", length(variant_ids))
  summary_list <- list()
  range_list <- list()
  miss_list <- list()
  assoc_list <- list()

  beta_groups <- c(0.05, 0.10, 0.20)
  rho_breaks <- c(0, 0.10, 0.25, 0.40, Inf)
  rho_labels <- c("[0,0.10)", "[0.10,0.25)", "[0.25,0.40)", "[0.40,Inf)")
  het_breaks <- c(0, 0.05, 0.10, 0.15, Inf)
  het_labels <- c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,Inf)")

  for (i in seq_along(variant_ids)) {
    vid <- variant_ids[[i]]
    merged <- merge_variant_power_s3(method_summary, metric_library, vid, method)
    detail_list[[i]] <- merged
    meta <- variant_planning_meta_s3(merged)
    n_scen_all <- length(unique(merged$scenario_id[is.finite(merged$power) & is.finite(merged$predicted_power)]))

    add_row <- function(d, beta1_group, n_scen) {
      power_error_summary_row_s3(
        d$power,
        d$predicted_power,
        cbind(meta, data.frame(
          beta1_group = beta1_group,
          n_scenarios = as.integer(n_scen),
          method = method,
          stringsAsFactors = FALSE
        ))
      )
    }

    summary_list[[length(summary_list) + 1L]] <- add_row(merged, "all_non_null", n_scen_all)
    for (b in beta_groups) {
      d <- merged[abs(merged$beta1_true - b) < 1e-12, , drop = FALSE]
      n_scen <- length(unique(d$scenario_id[is.finite(d$power) & is.finite(d$predicted_power)]))
      summary_list[[length(summary_list) + 1L]] <- add_row(d, sprintf("%.2f", b), n_scen)
    }

    for (bn in levels(merged$empirical_power_bin)) {
      d <- merged[as.character(merged$empirical_power_bin) == bn, , drop = FALSE]
      if (!nrow(d)) next
      n_scen <- length(unique(d$scenario_id[is.finite(d$power) & is.finite(d$predicted_power)]))
      range_list[[length(range_list) + 1L]] <- cbind(
        add_row(d, "all_non_null", n_scen),
        data.frame(empirical_power_range = bn, stringsAsFactors = FALSE)
      )
    }

    miss_specs <- list(
      list(col = "abs_rho_misspec", quantity = "abs_rho_plan_minus_true", breaks = rho_breaks, labels = rho_labels),
      list(col = "abs_tau_misspec", quantity = "abs_tau_plan_minus_true", breaks = het_breaks, labels = het_labels),
      list(col = "abs_omega_misspec", quantity = "abs_omega_plan_minus_true", breaks = het_breaks, labels = het_labels)
    )
    for (ms in miss_specs) {
      bins <- abs_misspec_bin_s3(merged[[ms$col]], ms$breaks, ms$labels)
      for (bn in ms$labels) {
        d <- merged[bins == bn, , drop = FALSE]
        if (!nrow(d)) next
        n_scen <- length(unique(d$scenario_id[is.finite(d$power) & is.finite(d$predicted_power)]))
        miss_list[[length(miss_list) + 1L]] <- cbind(
          add_row(d, "all_non_null", n_scen),
          data.frame(
            misspec_quantity = ms$quantity,
            misspec_bin = bn,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    assoc_one <- function(abs_x, signed_x, quantity) {
      ok_abs <- is.finite(merged$power_error) & is.finite(abs_x)
      ok_signed <- is.finite(merged$power_error) & is.finite(signed_x)
      spearman_safe <- function(a, b) {
        if (length(a) < 3L || stats::sd(a) == 0 || stats::sd(b) == 0) {
          return(NA_real_)
        }
        as.numeric(stats::cor(a, b, method = "spearman"))
      }
      cbind(
        meta,
        data.frame(
          misspec_quantity = quantity,
          n_rows = as.integer(sum(ok_abs)),
          spearman_abs_error_vs_abs_misspec = spearman_safe(
            abs(merged$power_error[ok_abs]), abs_x[ok_abs]
          ),
          spearman_signed_error_vs_signed_misspec = spearman_safe(
            merged$power_error[ok_signed], signed_x[ok_signed]
          ),
          stringsAsFactors = FALSE
        )
      )
    }
    assoc_list[[length(assoc_list) + 1L]] <- rbind(
      assoc_one(merged$abs_rho_misspec, merged$signed_rho_misspec, "rho_plan_minus_true"),
      assoc_one(merged$abs_tau_misspec, merged$signed_tau_misspec, "tau_plan_minus_true"),
      assoc_one(merged$abs_omega_misspec, merged$signed_omega_misspec, "omega_plan_minus_true")
    )
  }

  by_scenario <- do.call(rbind, detail_list)
  rownames(by_scenario) <- NULL
  overall <- do.call(rbind, summary_list)
  rownames(overall) <- NULL
  by_range <- if (length(range_list)) do.call(rbind, range_list) else overall[0, ]
  rownames(by_range) <- NULL
  misspec <- if (length(miss_list)) do.call(rbind, miss_list) else overall[0, ]
  rownames(misspec) <- NULL
  assoc <- do.call(rbind, assoc_list)
  rownames(assoc) <- NULL

  col_front <- c(
    "metric_variant", "variant_id", "beta1_group", "n_scenarios", "n_rows",
    "rmse_power", "mae_power", "mean_signed_error", "spearman",
    "calibration_intercept", "calibration_slope", "max_absolute_error",
    "weight_type", "rho_plan", "tau_plan", "omega_plan",
    "oracle_rho", "oracle_heterogeneity", "method"
  )
  overall <- overall[, col_front, drop = FALSE]
  range_cols <- c(col_front, "empirical_power_range")
  by_range <- by_range[, intersect(range_cols, names(by_range)), drop = FALSE]
  miss_cols <- c(col_front, "misspec_quantity", "misspec_bin")
  misspec <- misspec[, intersect(miss_cols, names(misspec)), drop = FALSE]

  list(
    overall = overall,
    by_empirical_range = by_range,
    misspecification = misspec,
    misspecification_association = assoc,
    by_scenario = by_scenario
  )
}

# Simulation Study 3 — chunking, checkpointing, parallel execution
# Chunks are scenario-level (not beta1-level) so CRN components are generated once.

source_sim3_functions <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  source(file.path(root, "R", "study1_functions.R"), local = FALSE)
  source(file.path(root, "R", "study2_design.R"), local = FALSE)
  source(file.path(root, "R", "study2_functions.R"), local = FALSE)
  source(file.path(root, "R", "estimators.R"), local = FALSE)
  source(file.path(root, "R", "contrast_information.R"), local = FALSE)
  source(file.path(root, "R", "study3_functions.R"), local = FALSE)
  invisible(TRUE)
}

sim3_design_hash <- function(root = NULL) {
  dir <- sim3_fixture_dir(root)
  files <- sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
  if (!length(files)) {
    stop("No Study 3 scenario fixtures found under ", dir, call. = FALSE)
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

chunk_rds_path_s3 <- function(out_dir, scenario_set, scenario_id, chunk_id) {
  file.path(
    out_dir,
    "raw_chunks",
    as.character(scenario_set),
    sprintf("scenario_%03d_chunk_%02d.rds", as.integer(scenario_id), as.integer(chunk_id))
  )
}

make_chunk_grid_s3 <- function(scenario_grid, n_rep_total, chunk_size) {
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
        pair_id = if ("pair_id" %in% names(sc)) as.integer(sc$pair_id) else NA_integer_,
        pair_class = if ("pair_class" %in% names(sc)) as.character(sc$pair_class) else NA_character_,
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
  out
}

validate_chunk_object_s3 <- function(obj, chunk, n_beta = 4L, n_methods = 4L) {
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

scenario_from_chunk_row_s3 <- function(chunk, scenario_grid) {
  sid <- as.integer(chunk$scenario_id[[1]])
  rows <- scenario_grid[scenario_grid$scenario_id == sid, , drop = FALSE]
  if (nrow(rows) != 1L) {
    stop("Could not uniquely recover scenario_id ", sid, call. = FALSE)
  }
  rows
}

run_one_chunk_s3 <- function(chunk,
                             scenario_grid,
                             study_structure,
                             root,
                             out_dir,
                             config,
                             design_hash = NULL) {
  chunk <- as.list(chunk)
  out_file <- chunk_rds_path_s3(
    out_dir, chunk$scenario_set, chunk$scenario_id, chunk$chunk_id
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (!exists("run_one_replication_s3", mode = "function")) {
    source_sim3_functions(root)
  }
  if (is.null(design_hash)) {
    design_hash <- sim3_design_hash(root)
  }

  if (file.exists(out_file)) {
    if (!isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      stop("Chunk file exists and resume/overwrite are FALSE: ", out_file, call. = FALSE)
    }
    if (isTRUE(config$resume) && !isTRUE(config$overwrite)) {
      obj <- tryCatch(readRDS(out_file), error = function(e) e)
      if (!inherits(obj, "error")) {
        chk <- validate_chunk_object_s3(obj, chunk, n_beta = length(config$beta1_levels))
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
  scenario <- scenario_from_chunk_row_s3(chunk, scenario_grid)
  n_beta <- length(config$beta1_levels)
  n_methods <- 4L

  result <- tryCatch(
    {
      reps <- seq.int(as.integer(chunk$rep_start), as.integer(chunk$rep_end))
      method_list <- vector("list", length(reps))
      diag_list <- vector("list", length(reps))
      for (i in seq_along(reps)) {
        one <- run_one_replication_s3(
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
      method = c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve"),
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
    error_message = err_msg,
    stringsAsFactors = FALSE
  )
  list(status = chunk_status, log = log_row, object = obj)
}

run_chunks_s3 <- function(chunk_grid,
                          scenario_grid,
                          study_structure,
                          root,
                          out_dir,
                          config) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(out_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)
  design_hash <- sim3_design_hash(root)
  workers <- as.integer(config$workers)
  if (is.na(workers) || workers < 1L) {
    stop("config$workers must be a positive integer.", call. = FALSE)
  }
  wall_start <- Sys.time()
  if (workers == 1L) {
    results <- lapply(seq_len(nrow(chunk_grid)), function(i) {
      run_one_chunk_s3(
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
      stop("Packages 'future' and 'future.apply' are required for parallel Study 3 runs.",
           call. = FALSE)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    results <- future.apply::future_lapply(
      seq_len(nrow(chunk_grid)),
      function(i) {
        run_one_chunk_s3(
          chunk = chunk_grid[i, , drop = FALSE],
          scenario_grid = scenario_grid,
          study_structure = study_structure,
          root = root,
          out_dir = out_dir,
          config = config,
          design_hash = design_hash
        )
      },
      # Seeds are set inside simulate_base_components_s3() via withr::with_seed().
      # Do not add a parallel RNG stream.
      future.seed = NULL
    )
  }
  wall_elapsed_sec <- as.numeric(difftime(Sys.time(), wall_start, units = "secs"))
  chunk_log <- do.call(rbind, lapply(results, function(x) x$log))
  rownames(chunk_log) <- NULL
  list(chunk_log = chunk_log, wall_elapsed_sec = wall_elapsed_sec, results = results)
}

collect_sim3_chunk_results <- function(out_dir, chunk_grid) {
  method_list <- vector("list", nrow(chunk_grid))
  diag_list <- vector("list", nrow(chunk_grid))
  for (i in seq_len(nrow(chunk_grid))) {
    ch <- chunk_grid[i, , drop = FALSE]
    path <- chunk_rds_path_s3(out_dir, ch$scenario_set[[1]], ch$scenario_id[[1]], ch$chunk_id[[1]])
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

write_sim3_summaries <- function(out_dir,
                                 config,
                                 scenario_grid,
                                 metric_library,
                                 method_results,
                                 pair_diagnostics = NULL) {
  method_summary <- summarize_method_performance_s3(method_results, config)
  design_keep <- unique(scenario_grid[, intersect(
    c(
      "scenario_id", "scenario_set", "pair_id", "pair_class", "J", "J0", "J1",
      "N_total", "K_total", "alignment_pattern", "alignment_strength",
      "nk_mode", "rho_true", "tau", "omega", "I_C_star", "J0_star", "J1_star",
      "nu_C_star", "predicted_variance", "predicted_se", "metric_status"
    ),
    names(scenario_grid)
  ), drop = FALSE])
  method_summary <- merge(method_summary, design_keep, by = "scenario_id", all.x = TRUE, suffixes = c("", "_design"))

  primary_metrics <- metric_library[metric_library$is_primary %in% TRUE, , drop = FALSE]
  var_pred <- summarize_variance_prediction_s3(method_summary, metric_library)
  df_pred <- summarize_df_prediction_s3(method_summary, metric_library)
  pwr_pred <- summarize_power_prediction_s3(method_summary, metric_library)
  df_che <- summarize_df_prediction_s3(method_summary, metric_library, method = "che_rve")
  pwr_che <- summarize_power_prediction_s3(method_summary, metric_library, method = "che_rve")

  rho_rows <- list()
  for (vid in c("sampling_rho_0.20", "sampling_rho_0.50", "sampling_rho_0.80", "sampling_rho_oracle")) {
    vp <- summarize_variance_prediction_s3(method_summary, metric_library, variant_id = vid)
    pp <- summarize_power_prediction_s3(method_summary, metric_library, variant_id = vid)
    rho_rows[[vid]] <- rbind(
      cbind(data.frame(quantity = "variance", stringsAsFactors = FALSE), vp$overall),
      cbind(data.frame(quantity = "power", stringsAsFactors = FALSE), pp$overall)
    )
  }
  rho_sens <- do.call(rbind, rho_rows)

  het_rows <- list()
  het_ids <- unique(metric_library$variant_id[metric_library$weight_type == "heterogeneity_adjusted"])
  for (vid in het_ids) {
    vp <- summarize_variance_prediction_s3(method_summary, metric_library, variant_id = vid)
    het_rows[[vid]] <- vp$overall
  }
  het_sens <- do.call(rbind, het_rows)

  pair_sum <- summarize_targeted_pairs_s3(method_summary, pair_diagnostics)

  mc <- method_summary[, c(
    "scenario_id", "beta1_true", "method", "n_success", "reject_rate",
    "mcse_reject_rate", "empirical_sd", "mean_df"
  )]

  utils::write.csv(design_keep, file.path(out_dir, "summaries", "scenario_design_summary.csv"), row.names = FALSE)
  utils::write.csv(primary_metrics, file.path(out_dir, "summaries", "scenario_metric_summary.csv"), row.names = FALSE)
  utils::write.csv(method_summary, file.path(out_dir, "summaries", "condition_method_summary.csv"), row.names = FALSE)
  utils::write.csv(var_pred$overall, file.path(out_dir, "summaries", "variance_prediction_summary.csv"), row.names = FALSE)
  utils::write.csv(df_pred$overall, file.path(out_dir, "summaries", "df_prediction_summary.csv"), row.names = FALSE)
  utils::write.csv(pwr_pred$overall, file.path(out_dir, "summaries", "power_prediction_summary.csv"), row.names = FALSE)
  utils::write.csv(rho_sens, file.path(out_dir, "summaries", "rho_sensitivity_summary.csv"), row.names = FALSE)
  utils::write.csv(het_sens, file.path(out_dir, "summaries", "heterogeneity_sensitivity_summary.csv"), row.names = FALSE)
  variant_power <- write_sim3_metric_variant_power_summaries(
    out_dir = out_dir,
    method_summary = method_summary,
    metric_library = metric_library
  )
  if (!is.null(pair_sum)) {
    utils::write.csv(pair_sum, file.path(out_dir, "summaries", "targeted_pair_summary.csv"), row.names = FALSE)
  }
  utils::write.csv(mc, file.path(out_dir, "summaries", "monte_carlo_uncertainty.csv"), row.names = FALSE)
  utils::write.csv(var_pred$by_scenario, file.path(out_dir, "summaries", "variance_prediction_by_scenario.csv"), row.names = FALSE)
  utils::write.csv(pwr_pred$by_scenario, file.path(out_dir, "summaries", "power_prediction_by_scenario.csv"), row.names = FALSE)

  list(
    method_summary = method_summary,
    variance_prediction = var_pred,
    df_prediction = df_pred,
    power_prediction = pwr_pred,
    df_prediction_che = df_che,
    power_prediction_che = pwr_che,
    metric_variant_power = variant_power,
    pair_summary = pair_sum
  )
}

write_sim3_metric_variant_power_summaries <- function(out_dir,
                                                      method_summary,
                                                      metric_library,
                                                      method = "mlma_cr2") {
  summaries_dir <- file.path(out_dir, "summaries")
  dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
  vp <- summarize_metric_variant_power_s3(
    method_summary = method_summary,
    metric_library = metric_library,
    method = method
  )
  slim_cols <- intersect(
    c(
      "scenario_id", "beta1_true", "metric_variant", "variant_id",
      "power", "predicted_power", "power_error", "predicted_se", "nu_C_star",
      "weight_type", "rho_plan", "tau_plan", "omega_plan",
      "is_primary", "is_default_het", "is_oracle_sampling", "is_oracle_het",
      "rho_true", "tau", "omega",
      "abs_rho_misspec", "abs_tau_misspec", "abs_omega_misspec",
      "signed_rho_misspec", "signed_tau_misspec", "signed_omega_misspec",
      "empirical_power_bin", "status"
    ),
    names(vp$by_scenario)
  )
  utils::write.csv(
    vp$overall,
    file.path(summaries_dir, "metric_variant_power_prediction_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    vp$by_empirical_range,
    file.path(summaries_dir, "metric_variant_power_prediction_by_empirical_range.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    vp$misspecification,
    file.path(summaries_dir, "metric_variant_power_prediction_misspecification.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    vp$misspecification_association,
    file.path(summaries_dir, "metric_variant_power_prediction_misspecification_association.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    vp$by_scenario[, slim_cols, drop = FALSE],
    file.path(summaries_dir, "metric_variant_power_prediction_by_scenario.csv"),
    row.names = FALSE
  )
  vp
}

# Simulation Study 3 — post-production heterogeneity-sensitivity analysis
#
# Exploratory / sensitivity layer only. Does not redefine the preregistered
# primary metric (sampling_rho_0.50), does not modify the frozen scenario
# library, and does not rerun replications. Assumed tau and omega vary
# independently on a fixed planning grid.

sim3_het_sens_analysis_label <- function() {
  "post_production_exploratory_heterogeneity_sensitivity"
}

sim3_het_sens_tau_grid <- function() c(0, 0.05, 0.10, 0.15, 0.20)
sim3_het_sens_omega_grid <- function() c(0, 0.05, 0.10, 0.15, 0.20)
sim3_het_sens_rho_grid <- function() c(0.20, 0.50, 0.80)
sim3_het_sens_primary_rho <- function() 0.50
sim3_het_sens_near_tol <- function() 0.025

sim3_het_sens_planning_spec <- function(tau_grid = sim3_het_sens_tau_grid(),
                                        omega_grid = sim3_het_sens_omega_grid(),
                                        rho_grid = sim3_het_sens_rho_grid(),
                                        primary_rho = sim3_het_sens_primary_rho()) {
  spec <- expand.grid(
    tau_plan = as.numeric(tau_grid),
    omega_plan = as.numeric(omega_grid),
    rho_plan = as.numeric(rho_grid),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  spec$tau_plan <- round(spec$tau_plan, 2)
  spec$omega_plan <- round(spec$omega_plan, 2)
  spec$rho_plan <- round(spec$rho_plan, 2)
  spec$rho_plan_role <- ifelse(
    abs(spec$rho_plan - primary_rho) < 1e-12,
    "primary",
    "secondary"
  )
  spec$tied_tau_omega <- abs(spec$tau_plan - spec$omega_plan) < 1e-12
  spec$is_sampling_equivalent <- spec$tau_plan < 1e-12 & spec$omega_plan < 1e-12
  spec$variant_id <- sprintf(
    "het_tau%.2f_omega%.2f_rho%.2f",
    spec$tau_plan, spec$omega_plan, spec$rho_plan
  )
  spec[order(spec$rho_plan, spec$tau_plan, spec$omega_plan), , drop = FALSE]
}

prediction_error_metrics_with_pearson_s3 <- function(observed, predicted) {
  err <- prediction_error_metrics_s3(observed, predicted)
  ok <- is.finite(observed) & is.finite(predicted)
  pearson <- if (sum(ok) >= 3L &&
                 stats::sd(observed[ok]) > 0 &&
                 stats::sd(predicted[ok]) > 0) {
    as.numeric(stats::cor(observed[ok], predicted[ok], method = "pearson"))
  } else {
    NA_real_
  }
  err$pearson <- pearson
  err
}

het_sens_direction_s3 <- function(plan, true, tol = sim3_het_sens_near_tol()) {
  d <- as.numeric(plan) - as.numeric(true)
  out <- rep(NA_character_, length(d))
  ok <- is.finite(d)
  out[ok & abs(d) <= tol] <- "near"
  out[ok & d < -tol] <- "under"
  out[ok & d > tol] <- "over"
  out
}

het_sens_joint_class_s3 <- function(tau_dir, omega_dir) {
  tau_dir <- as.character(tau_dir)
  omega_dir <- as.character(omega_dir)
  out <- rep(NA_character_, length(tau_dir))
  ok <- !is.na(tau_dir) & !is.na(omega_dir)
  both_near <- ok & tau_dir == "near" & omega_dir == "near"
  both_under <- ok & tau_dir == "under" & omega_dir == "under"
  both_over <- ok & tau_dir == "over" & omega_dir == "over"
  tau_near_omega_not <- ok & tau_dir == "near" & omega_dir != "near"
  omega_near_tau_not <- ok & omega_dir == "near" & tau_dir != "near"
  mixed <- ok & tau_dir != "near" & omega_dir != "near" & tau_dir != omega_dir
  out[both_near] <- "both_near"
  out[both_under] <- "both_under"
  out[both_over] <- "both_over"
  out[tau_near_omega_not] <- paste0("tau_near_omega_", omega_dir[tau_near_omega_not])
  out[omega_near_tau_not] <- paste0("omega_near_tau_", tau_dir[omega_near_tau_not])
  out[mixed] <- paste0("mixed_tau_", tau_dir[mixed], "_omega_", omega_dir[mixed])
  out
}

het_sens_one_correct_class_s3 <- function(tau_dir, omega_dir) {
  tau_dir <- as.character(tau_dir)
  omega_dir <- as.character(omega_dir)
  out <- rep(NA_character_, length(tau_dir))
  ok <- !is.na(tau_dir) & !is.na(omega_dir)
  out[ok & tau_dir == "near" & omega_dir == "near"] <- "both_near"
  out[ok & tau_dir == "near" & omega_dir != "near"] <- "tau_near_omega_misspecified"
  out[ok & omega_dir == "near" & tau_dir != "near"] <- "omega_near_tau_misspecified"
  out[ok & tau_dir != "near" & omega_dir != "near"] <- "neither_near"
  out
}

het_sens_abs_bin_s3 <- function(x) {
  abs_misspec_bin_s3(
    abs(as.numeric(x)),
    breaks = c(0, 0.025, 0.05, 0.10, 0.15, Inf),
    labels = c("[0,0.025)", "[0.025,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,Inf)")
  )
}

het_sens_signed_bin_s3 <- function(x) {
  abs_misspec_bin_s3(
    as.numeric(x),
    breaks = c(-Inf, -0.15, -0.10, -0.05, -0.025, 0.025, 0.05, 0.10, 0.15, Inf),
    labels = c(
      "< -0.15", "[-0.15,-0.10)", "[-0.10,-0.05)", "[-0.05,-0.025)",
      "near [-0.025,0.025)",
      "[0.025,0.05)", "[0.05,0.10)", "[0.10,0.15)", ">= 0.15"
    )
  )
}

rbind_fill_s3 <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1]]) && !is.data.frame(args[[1]])) {
    args <- args[[1]]
  }
  args <- args[!vapply(args, is.null, logical(1))]
  args <- args[vapply(args, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (!length(args)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  nms <- unique(unlist(lapply(args, names), use.names = FALSE))
  aligned <- lapply(args, function(d) {
    missing <- setdiff(nms, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, nms, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

het_sens_error_row_s3 <- function(d, meta) {
  ok <- is.finite(d$power) & is.finite(d$predicted_power)
  err <- prediction_error_metrics_with_pearson_s3(d$power, d$predicted_power)
  n_scen <- if ("scenario_id" %in% names(d)) {
    length(unique(d$scenario_id[ok]))
  } else {
    NA_integer_
  }
  cbind(
    meta,
    data.frame(
      n_scenarios = as.integer(n_scen),
      n_rows = err$n,
      rmse_power = err$rmse,
      mae_power = err$mae,
      mean_signed_error = err$mean_signed_error,
      spearman = err$spearman,
      pearson = err$pearson,
      calibration_intercept = err$calibration_intercept,
      calibration_slope = err$calibration_slope,
      r2 = err$r2,
      max_absolute_error = err$max_abs_error,
      median_ape = err$median_ape,
      mean_abs_tau_misspec = if ("abs_tau_misspec" %in% names(d)) {
        mean(d$abs_tau_misspec[ok], na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_abs_omega_misspec = if ("abs_omega_misspec" %in% names(d)) {
        mean(d$abs_omega_misspec[ok], na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_signed_tau_misspec = if ("signed_tau_misspec" %in% names(d)) {
        mean(d$signed_tau_misspec[ok], na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_signed_omega_misspec = if ("signed_omega_misspec" %in% names(d)) {
        mean(d$signed_omega_misspec[ok], na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  )
}

compute_heterogeneity_sensitivity_metrics_s3 <- function(study_structure,
                                                         scenario_ids = NULL,
                                                         spec = sim3_het_sens_planning_spec()) {
  ss <- study_structure
  if (!"scenario_id" %in% names(ss)) {
    stop("study_structure must include scenario_id.", call. = FALSE)
  }
  if (is.null(scenario_ids)) {
    scenario_ids <- sort(unique(as.integer(ss$scenario_id)))
  } else {
    scenario_ids <- as.integer(scenario_ids)
  }
  rows <- vector("list", length(scenario_ids) * nrow(spec))
  k <- 0L
  for (sid in scenario_ids) {
    ss_one <- ss[ss$scenario_id == sid, , drop = FALSE]
    for (i in seq_len(nrow(spec))) {
      v <- spec[i, , drop = FALSE]
      info <- compute_contrast_information_s3(
        study_structure = ss_one,
        rho_plan = v$rho_plan[[1]],
        tau_plan = v$tau_plan[[1]],
        omega_plan = v$omega_plan[[1]],
        weight_type = "heterogeneity_adjusted"
      )
      k <- k + 1L
      rows[[k]] <- cbind(
        data.frame(
          scenario_id = as.integer(sid),
          variant_id = as.character(v$variant_id[[1]]),
          rho_plan_role = as.character(v$rho_plan_role[[1]]),
          tied_tau_omega = as.logical(v$tied_tau_omega[[1]]),
          is_sampling_equivalent = as.logical(v$is_sampling_equivalent[[1]]),
          analysis_label = sim3_het_sens_analysis_label(),
          stringsAsFactors = FALSE
        ),
        info
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_heterogeneity_sensitivity_grid_s3 <- function(study_structure,
                                                    scenario_grid,
                                                    method_summary,
                                                    spec = sim3_het_sens_planning_spec(),
                                                    method = "mlma_cr2",
                                                    near_tol = sim3_het_sens_near_tol()) {
  metrics <- compute_heterogeneity_sensitivity_metrics_s3(
    study_structure = study_structure,
    scenario_ids = scenario_grid$scenario_id,
    spec = spec
  )
  pwr <- method_summary[
    method_summary$method == method & method_summary$beta1_true > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(pwr)) {
    stop("method_summary has no non-null rows for method ", method, ".", call. = FALSE)
  }
  design_cols <- intersect(
    c(
      "scenario_id", "scenario_set", "pair_id", "pair_class", "J", "J0", "J1",
      "N_total", "K_total", "alignment_pattern", "alignment_strength",
      "nk_mode", "rho_true", "tau", "omega", "n_cv", "k_cv",
      "minority_prop", "I_C_star", "nu_C_star"
    ),
    names(scenario_grid)
  )
  design <- unique(scenario_grid[, design_cols, drop = FALSE])
  names(design)[names(design) == "I_C_star"] <- "I_C_star_primary"
  names(design)[names(design) == "nu_C_star"] <- "nu_C_star_primary"
  if ("tau" %in% names(design)) names(design)[names(design) == "tau"] <- "tau_true"
  if ("omega" %in% names(design)) names(design)[names(design) == "omega"] <- "omega_true"

  keep_emp <- intersect(
    c(
      "scenario_id", "beta1_true", "method", "power", "n_attempted", "n_success",
      "error_rate", "empirical_variance", "mean_df"
    ),
    names(pwr)
  )
  pwr <- pwr[, keep_emp, drop = FALSE]
  pwr <- merge(pwr, design, by = "scenario_id", all.x = TRUE)

  pwr <- drop_overlapping_metric_cols_s3(pwr, metrics)
  merged <- merge(pwr, metrics, by = "scenario_id", suffixes = c("", "_metric"))
  merged$predicted_power <- predict_contrast_power_s3(
    beta1 = merged$beta1_true,
    predicted_se = merged$predicted_se,
    nu_C_star = merged$nu_C_star
  )
  merged$power_error <- merged$predicted_power - merged$power
  merged$signed_tau_misspec <- merged$tau_plan - merged$tau_true
  merged$signed_omega_misspec <- merged$omega_plan - merged$omega_true
  merged$signed_rho_misspec <- merged$rho_plan - merged$rho_true
  merged$abs_tau_misspec <- abs(merged$signed_tau_misspec)
  merged$abs_omega_misspec <- abs(merged$signed_omega_misspec)
  merged$abs_rho_misspec <- abs(merged$signed_rho_misspec)
  merged$tau_dir <- het_sens_direction_s3(merged$tau_plan, merged$tau_true, near_tol)
  merged$omega_dir <- het_sens_direction_s3(merged$omega_plan, merged$omega_true, near_tol)
  merged$joint_misspec_class <- het_sens_joint_class_s3(merged$tau_dir, merged$omega_dir)
  merged$one_correct_class <- het_sens_one_correct_class_s3(merged$tau_dir, merged$omega_dir)
  merged$abs_tau_bin <- het_sens_abs_bin_s3(merged$abs_tau_misspec)
  merged$abs_omega_bin <- het_sens_abs_bin_s3(merged$abs_omega_misspec)
  merged$signed_tau_bin <- het_sens_signed_bin_s3(merged$signed_tau_misspec)
  merged$signed_omega_bin <- het_sens_signed_bin_s3(merged$signed_omega_misspec)
  merged$empirical_power_bin <- as.character(empirical_power_bin_s3(merged$power))
  merged$analysis_label <- sim3_het_sens_analysis_label()
  merged[order(
    merged$rho_plan, merged$tau_plan, merged$omega_plan,
    merged$scenario_id, merged$beta1_true
  ), , drop = FALSE]
}

summarize_heterogeneity_sensitivity_s3 <- function(grid) {
  beta_groups <- c(0.05, 0.10, 0.20)
  strata <- c("all", "broad", "targeted")
  keys <- unique(grid[, c("rho_plan", "tau_plan", "omega_plan", "variant_id",
                          "rho_plan_role", "tied_tau_omega", "is_sampling_equivalent"),
                      drop = FALSE])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, , drop = FALSE]
    d0 <- grid[
      abs(grid$rho_plan - k$rho_plan) < 1e-12 &
        abs(grid$tau_plan - k$tau_plan) < 1e-12 &
        abs(grid$omega_plan - k$omega_plan) < 1e-12,
      ,
      drop = FALSE
    ]
    for (st in strata) {
      d_st <- if (identical(st, "all")) d0 else d0[d0$scenario_set == st, , drop = FALSE]
      if (!nrow(d_st)) next
      meta_base <- data.frame(
        analysis_family = "planning_cell",
        analysis_label = sim3_het_sens_analysis_label(),
        variant_id = k$variant_id,
        rho_plan = k$rho_plan,
        tau_plan = k$tau_plan,
        omega_plan = k$omega_plan,
        rho_plan_role = k$rho_plan_role,
        tied_tau_omega = k$tied_tau_omega,
        is_sampling_equivalent = k$is_sampling_equivalent,
        scenario_set = st,
        method = "mlma_cr2",
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- het_sens_error_row_s3(
        d_st,
        cbind(meta_base, data.frame(beta1_group = "all_non_null", stringsAsFactors = FALSE))
      )
      for (b in beta_groups) {
        d_b <- d_st[abs(d_st$beta1_true - b) < 1e-12, , drop = FALSE]
        if (!nrow(d_b)) next
        rows[[length(rows) + 1L]] <- het_sens_error_row_s3(
          d_b,
          cbind(meta_base, data.frame(beta1_group = sprintf("%.2f", b), stringsAsFactors = FALSE))
        )
      }
    }
  }
  range_sum <- summarize_heterogeneity_power_ranges_s3(grid)
  range_rows <- range_sum$summary
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  list(planning_cells = out, power_ranges = range_rows, power_range_by_scenario = range_sum$by_scenario)
}

summarize_heterogeneity_power_ranges_s3 <- function(grid,
                                                    primary_rho = sim3_het_sens_primary_rho()) {
  d <- grid[abs(grid$rho_plan - primary_rho) < 1e-12, , drop = FALSE]
  split_key <- paste(d$scenario_id, d$beta1_true, sep = ":")
  pieces <- split(d, split_key)
  by_scen <- do.call(rbind, lapply(pieces, function(x) {
    emp <- x$power[[1]]
    full <- x$predicted_power
    nozero <- x$predicted_power[!(x$is_sampling_equivalent %in% TRUE)]
    neighborhood <- x$predicted_power[x$abs_tau_misspec <= 0.075 & x$abs_omega_misspec <= 0.075]
    conservative <- x$predicted_power[x$tau_plan >= 0.10 - 1e-12 & x$omega_plan >= 0.10 - 1e-12]
    max_het <- x$predicted_power[abs(x$tau_plan - 0.20) < 1e-12 & abs(x$omega_plan - 0.20) < 1e-12]
    design_only <- x$predicted_power[x$is_sampling_equivalent %in% TRUE]
    default_het <- x$predicted_power[abs(x$tau_plan - 0.10) < 1e-12 & abs(x$omega_plan - 0.10) < 1e-12]
    range_stats <- function(pred, label) {
      pred <- pred[is.finite(pred)]
      if (!length(pred) || !is.finite(emp)) {
        return(data.frame(
          range_definition = label,
          n_assumptions = length(pred),
          min_predicted_power = NA_real_,
          max_predicted_power = NA_real_,
          predicted_power_range = NA_real_,
          empirical_inside_range = NA,
          stringsAsFactors = FALSE
        ))
      }
      mn <- min(pred)
      mx <- max(pred)
      data.frame(
        range_definition = label,
        n_assumptions = length(pred),
        min_predicted_power = mn,
        max_predicted_power = mx,
        predicted_power_range = mx - mn,
        empirical_inside_range = emp >= mn - 1e-12 & emp <= mx + 1e-12,
        stringsAsFactors = FALSE
      )
    }
    defs <- rbind(
      range_stats(full, "full_5x5_grid"),
      range_stats(nozero, "exclude_tau_omega_zero"),
      range_stats(neighborhood, "neighborhood_abs_misspec_le_0.075"),
      range_stats(conservative, "conservative_tau_omega_ge_0.10"),
      range_stats(max_het, "point_tau_omega_0.20")
    )
    cbind(
      data.frame(
        scenario_id = x$scenario_id[[1]],
        scenario_set = x$scenario_set[[1]],
        pair_id = if ("pair_id" %in% names(x)) x$pair_id[[1]] else NA_integer_,
        pair_class = if ("pair_class" %in% names(x)) x$pair_class[[1]] else NA_character_,
        beta1_true = x$beta1_true[[1]],
        tau_true = x$tau_true[[1]],
        omega_true = x$omega_true[[1]],
        J = if ("J" %in% names(x)) x$J[[1]] else NA_real_,
        I_C_star_primary = if ("I_C_star_primary" %in% names(x)) x$I_C_star_primary[[1]] else NA_real_,
        empirical_power = emp,
        predicted_power_sampling_equivalent = if (length(design_only)) design_only[[1]] else NA_real_,
        predicted_power_default_het = if (length(default_het)) default_het[[1]] else NA_real_,
        predicted_power_tau_omega_0.20 = if (length(max_het)) max_het[[1]] else NA_real_,
        stringsAsFactors = FALSE
      ),
      defs
    )
  }))
  rownames(by_scen) <- NULL

  sum_one <- function(x, st) {
    ok <- is.finite(x$predicted_power_range) & is.finite(x$empirical_power)
    inside <- x$empirical_inside_range[ok]
    data.frame(
      analysis_family = "power_range",
      analysis_label = sim3_het_sens_analysis_label(),
      variant_id = NA_character_,
      rho_plan = primary_rho,
      tau_plan = NA_real_,
      omega_plan = NA_real_,
      rho_plan_role = "primary",
      tied_tau_omega = NA,
      is_sampling_equivalent = NA,
      scenario_set = st,
      method = "mlma_cr2",
      beta1_group = "all_non_null",
      range_definition = x$range_definition[[1]],
      n_scenarios = length(unique(x$scenario_id[ok])),
      n_rows = as.integer(sum(ok)),
      mean_predicted_range = mean(x$predicted_power_range[ok]),
      median_predicted_range = stats::median(x$predicted_power_range[ok]),
      pct_empirical_inside_range = if (length(inside)) mean(inside) else NA_real_,
      mean_min_predicted = mean(x$min_predicted_power[ok]),
      mean_max_predicted = mean(x$max_predicted_power[ok]),
      mean_empirical_power = mean(x$empirical_power[ok]),
      stringsAsFactors = FALSE
    )
  }
  summary_rows <- list()
  for (def in unique(by_scen$range_definition)) {
    xd <- by_scen[by_scen$range_definition == def, , drop = FALSE]
    summary_rows[[length(summary_rows) + 1L]] <- sum_one(xd, "all")
    for (st in c("broad", "targeted")) {
      xs <- xd[xd$scenario_set == st, , drop = FALSE]
      if (nrow(xs)) summary_rows[[length(summary_rows) + 1L]] <- sum_one(xs, st)
    }
  }
  summary_out <- do.call(rbind, summary_rows)
  rownames(summary_out) <- NULL
  list(by_scenario = by_scen, summary = summary_out)
}

summarize_heterogeneity_misspecification_s3 <- function(grid,
                                                        primary_rho = sim3_het_sens_primary_rho()) {
  d <- grid[abs(grid$rho_plan - primary_rho) < 1e-12, , drop = FALSE]
  add_rows <- function(df, family, quantity, bin, extra = NULL) {
    if (!nrow(df)) return(NULL)
    meta <- data.frame(
      analysis_family = family,
      analysis_label = sim3_het_sens_analysis_label(),
      rho_plan = primary_rho,
      rho_plan_role = "primary",
      method = "mlma_cr2",
      beta1_group = "all_non_null",
      scenario_set = if ("scenario_set_scope" %in% names(extra)) extra$scenario_set_scope else "all",
      misspec_quantity = quantity,
      misspec_bin = bin,
      stringsAsFactors = FALSE
    )
    if (!is.null(extra)) {
      extra$scenario_set_scope <- NULL
      if (ncol(extra)) meta <- cbind(meta, extra)
    }
    het_sens_error_row_s3(df, meta)
  }

  rows <- list()
  for (st in c("all", "broad", "targeted")) {
    d_st <- if (identical(st, "all")) d else d[d$scenario_set == st, , drop = FALSE]
    extra_st <- data.frame(scenario_set_scope = st, stringsAsFactors = FALSE)

    for (bn in unique(d_st$signed_tau_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$signed_tau_bin == bn, , drop = FALSE],
        "signed_tau_misspec", "tau_plan_minus_true", bn, extra_st
      )
    }
    for (bn in unique(d_st$signed_omega_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$signed_omega_bin == bn, , drop = FALSE],
        "signed_omega_misspec", "omega_plan_minus_true", bn, extra_st
      )
    }
    for (bn in unique(d_st$abs_tau_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$abs_tau_bin == bn, , drop = FALSE],
        "abs_tau_misspec", "abs_tau_plan_minus_true", bn, extra_st
      )
    }
    for (bn in unique(d_st$abs_omega_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$abs_omega_bin == bn, , drop = FALSE],
        "abs_omega_misspec", "abs_omega_plan_minus_true", bn, extra_st
      )
    }
    for (bn in unique(d_st$joint_misspec_class)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$joint_misspec_class == bn, , drop = FALSE],
        "joint_direction", "joint_tau_omega_direction", bn, extra_st
      )
    }
    for (bn in unique(d_st$one_correct_class)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_st[d_st$one_correct_class == bn, , drop = FALSE],
        "one_component_correct", "one_component_near", bn, extra_st
      )
    }

    tau_bins <- unique(d_st$abs_tau_bin)
    omega_bins <- unique(d_st$abs_omega_bin)
    for (tb in tau_bins) {
      if (is.na(tb)) next
      for (ob in omega_bins) {
        if (is.na(ob)) next
        dd <- d_st[d_st$abs_tau_bin == tb & d_st$abs_omega_bin == ob, , drop = FALSE]
        if (!nrow(dd)) next
        rows[[length(rows) + 1L]] <- add_rows(
          dd, "abs_tau_x_abs_omega", "abs_tau_and_abs_omega",
          paste(tb, ob, sep = " x "), extra_st
        )
      }
    }

    d_omega_near <- d_st[d_st$omega_dir == "near", , drop = FALSE]
    for (bn in unique(d_omega_near$abs_tau_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_omega_near[d_omega_near$abs_tau_bin == bn, , drop = FALSE],
        "tau_given_omega_near", "abs_tau_plan_minus_true | omega_near",
        bn, extra_st
      )
    }
    d_tau_near <- d_st[d_st$tau_dir == "near", , drop = FALSE]
    for (bn in unique(d_tau_near$abs_omega_bin)) {
      if (is.na(bn)) next
      rows[[length(rows) + 1L]] <- add_rows(
        d_tau_near[d_tau_near$abs_omega_bin == bn, , drop = FALSE],
        "omega_given_tau_near", "abs_omega_plan_minus_true | tau_near",
        bn, extra_st
      )
    }
  }

  rbind_fill_s3(rows)
}

summarize_heterogeneity_sensitivity_by_design_slice_s3 <- function(
  grid,
  primary_rho = sim3_het_sens_primary_rho()
) {
  d <- grid[abs(grid$rho_plan - primary_rho) < 1e-12, , drop = FALSE]
  slice_one <- function(df, dim_name, dim_value, family, extra = NULL) {
    if (!nrow(df)) return(NULL)
    meta <- data.frame(
      analysis_family = family,
      analysis_label = sim3_het_sens_analysis_label(),
      rho_plan = primary_rho,
      rho_plan_role = "primary",
      method = "mlma_cr2",
      beta1_group = "all_non_null",
      design_dimension = dim_name,
      dimension_value = as.character(dim_value[[1]]),
      stringsAsFactors = FALSE
    )
    if (!is.null(extra)) meta <- cbind(meta, extra)
    het_sens_error_row_s3(df, meta)
  }
  add_cut <- function(df, dim_name, breaks, labels, family, extra = NULL) {
    df$bin <- cut(df[[dim_name]], breaks = breaks, labels = labels,
                  include.lowest = TRUE, right = TRUE)
    out <- lapply(split(df, df$bin, drop = TRUE), function(x) {
      slice_one(x, dim_name, as.character(x$bin[[1]]), family, extra)
    })
    out
  }

  rows <- list()
  planning_keys <- unique(d[, c("tau_plan", "omega_plan", "variant_id",
                                "tied_tau_omega", "is_sampling_equivalent"),
                            drop = FALSE])
  slice_dims <- function(df, family, extra) {
    out <- list(
      slice_one(df, "all", "all", family, extra)
    )
    if ("scenario_set" %in% names(df)) {
      out <- c(out, lapply(split(df, df$scenario_set), function(x) {
        slice_one(x, "scenario_set", x$scenario_set[[1]], family, extra)
      }))
    }
    if ("beta1_true" %in% names(df)) {
      out <- c(out, lapply(split(df, df$beta1_true), function(x) {
        slice_one(x, "beta1_true", x$beta1_true[[1]], family, extra)
      }))
    }
    if ("J" %in% names(df)) {
      out <- c(out, add_cut(df, "J", c(0, 36, 60, 96, Inf),
                            c("J<=36", "37-60", "61-96", "J>96"), family, extra))
    }
    if ("tau_true" %in% names(df)) {
      out <- c(out, add_cut(df, "tau_true", c(0, 0.05, 0.10, 0.15, Inf),
                            c("tau<0.05", "[0.05,0.10)", "[0.10,0.15)", "tau>=0.15"),
                            family, extra))
    }
    if ("omega_true" %in% names(df)) {
      out <- c(out, add_cut(df, "omega_true", c(0, 0.05, 0.10, 0.15, Inf),
                            c("omega<0.05", "[0.05,0.10)", "[0.10,0.15)", "omega>=0.15"),
                            family, extra))
    }
    if ("I_C_star_primary" %in% names(df)) {
      out <- c(out, add_cut(df, "I_C_star_primary", c(0, 2000, 8000, 20000, Inf),
                            c("I_C*<2000", "[2000,8000)", "[8000,20000)", "I_C*>=20000"),
                            family, extra))
    }
    if ("nu_C_star_primary" %in% names(df)) {
      out <- c(out, add_cut(df, "nu_C_star_primary", c(-Inf, 5, 10, 20, Inf),
                            c("nu_C*<5", "[5,10)", "[10,20)", "nu_C*>=20"),
                            family, extra))
    }
    if ("n_cv" %in% names(df)) {
      out <- c(out, add_cut(df, "n_cv", c(-0.01, 0.3, 0.6, Inf),
                            c("n_cv<=0.3", "0.3<n_cv<=0.6", "n_cv>0.6"), family, extra))
    }
    if ("alignment_pattern" %in% names(df)) {
      out <- c(out, lapply(split(df, df$alignment_pattern), function(x) {
        slice_one(x, "alignment_pattern", x$alignment_pattern[[1]], family, extra)
      }))
    }
    if ("nk_mode" %in% names(df)) {
      out <- c(out, lapply(split(df, df$nk_mode), function(x) {
        slice_one(x, "nk_mode", x$nk_mode[[1]], family, extra)
      }))
    }
    out
  }

  for (i in seq_len(nrow(planning_keys))) {
    k <- planning_keys[i, , drop = FALSE]
    dk <- d[
      abs(d$tau_plan - k$tau_plan) < 1e-12 &
        abs(d$omega_plan - k$omega_plan) < 1e-12,
      ,
      drop = FALSE
    ]
    extra <- data.frame(
      variant_id = k$variant_id,
      tau_plan = k$tau_plan,
      omega_plan = k$omega_plan,
      tied_tau_omega = k$tied_tau_omega,
      is_sampling_equivalent = k$is_sampling_equivalent,
      misspec_quantity = NA_character_,
      misspec_bin = NA_character_,
      stringsAsFactors = FALSE
    )
    rows <- c(rows, slice_dims(dk, "planning_cell_by_design", extra))
  }

  # Misspecification magnitude within design regions, pooling the 5x5 grid.
  extra_pool <- data.frame(
    variant_id = NA_character_,
    tau_plan = NA_real_,
    omega_plan = NA_real_,
    tied_tau_omega = NA,
    is_sampling_equivalent = NA,
    misspec_quantity = NA_character_,
    misspec_bin = NA_character_,
    stringsAsFactors = FALSE
  )
  for (bn in unique(d$abs_tau_bin)) {
    if (is.na(bn)) next
    extra_pool$misspec_quantity <- "abs_tau_plan_minus_true"
    extra_pool$misspec_bin <- bn
    rows <- c(rows, slice_dims(
      d[d$abs_tau_bin == bn, , drop = FALSE],
      "abs_tau_misspec_by_design",
      extra_pool
    ))
  }
  extra_omega <- extra_pool
  for (bn in unique(d$abs_omega_bin)) {
    if (is.na(bn)) next
    extra_omega$misspec_quantity <- "abs_omega_plan_minus_true"
    extra_omega$misspec_bin <- bn
    rows <- c(rows, slice_dims(
      d[d$abs_omega_bin == bn, , drop = FALSE],
      "abs_omega_misspec_by_design",
      extra_omega
    ))
  }

  rbind_fill_s3(rows)
}

write_heterogeneity_sensitivity_outputs_s3 <- function(grid,
                                                       out_dir,
                                                       include_power_range_by_scenario = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary_obj <- summarize_heterogeneity_sensitivity_s3(grid)
  misspec <- summarize_heterogeneity_misspecification_s3(grid)
  by_slice <- summarize_heterogeneity_sensitivity_by_design_slice_s3(grid)

  grid_cols <- intersect(
    c(
      "analysis_label", "scenario_id", "scenario_set", "pair_id", "pair_class",
      "beta1_true", "method", "variant_id", "rho_plan", "rho_plan_role",
      "tau_plan", "omega_plan", "tied_tau_omega", "is_sampling_equivalent",
      "tau_true", "omega_true", "rho_true",
      "signed_tau_misspec", "signed_omega_misspec", "signed_rho_misspec",
      "abs_tau_misspec", "abs_omega_misspec", "abs_rho_misspec",
      "tau_dir", "omega_dir", "joint_misspec_class", "one_correct_class",
      "abs_tau_bin", "abs_omega_bin", "signed_tau_bin", "signed_omega_bin",
      "J", "J0", "J1", "n_cv", "alignment_pattern", "nk_mode",
      "I_C_star_primary", "nu_C_star_primary",
      "I_C_star", "predicted_variance", "predicted_se",
      "J0_star", "J1_star", "nu_C_star", "status",
      "power", "predicted_power", "power_error",
      "n_attempted", "n_success", "error_rate", "empirical_power_bin"
    ),
    names(grid)
  )
  grid_out <- grid[, grid_cols, drop = FALSE]

  summary_out <- rbind_fill_s3(summary_obj$planning_cells, summary_obj$power_ranges)

  utils::write.csv(
    grid_out,
    file.path(out_dir, "heterogeneity_sensitivity_grid.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary_out,
    file.path(out_dir, "heterogeneity_sensitivity_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    misspec,
    file.path(out_dir, "heterogeneity_misspecification_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    by_slice,
    file.path(out_dir, "heterogeneity_sensitivity_by_design_slice.csv"),
    row.names = FALSE
  )
  if (isTRUE(include_power_range_by_scenario)) {
    utils::write.csv(
      summary_obj$power_range_by_scenario,
      file.path(out_dir, "heterogeneity_power_range_by_scenario.csv"),
      row.names = FALSE
    )
  }

  list(
    grid = grid_out,
    summary = summary_out,
    misspecification = misspec,
    by_design_slice = by_slice,
    power_range_by_scenario = summary_obj$power_range_by_scenario
  )
}

run_heterogeneity_sensitivity_analysis_s3 <- function(root,
                                                      method_summary,
                                                      out_dir,
                                                      spec = sim3_het_sens_planning_spec()) {
  lib <- read_frozen_scenario_library_s3(root)
  grid <- build_heterogeneity_sensitivity_grid_s3(
    study_structure = lib$study_structure,
    scenario_grid = lib$grid,
    method_summary = method_summary,
    spec = spec
  )
  write_heterogeneity_sensitivity_outputs_s3(grid, out_dir)
}
