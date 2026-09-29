# Simulation Study 4 runner.
#
# Published run (5,000 replications, chunk size 500; requires an output directory):
#   Rscript scripts/04_run_study4.R --production --output-dir reruns/study4_production
#   Rscript scripts/04_run_study4.R --production --workers 15 --output-dir reruns/study4_production
#
# Execution check (does not write results/study4/ or the frozen design):
#   Rscript scripts/04_run_study4.R --n-rep 1 --workers 1 --output-dir results/study4_check
#
# The run reads data/study4_design/. It does not construct or freeze scenarios.
# Outcome seeds come from make_sim4_config("production"):
# seed_root_outcomes, scenario_stride, and replication_stride.
# Paths are relative to reproducibility/ unless absolute and still inside it.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]))
} else {
  normalizePath(file.path("scripts", "04_run_study4.R"), mustWork = TRUE)
}
repro_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(repro_root, "R", "study1_functions.R"))
source(file.path(repro_root, "R", "study2_design.R"))
source(file.path(repro_root, "R", "study2_functions.R"))
source(file.path(repro_root, "R", "estimators.R"))
source(file.path(repro_root, "R", "contrast_information.R"))
source(file.path(repro_root, "R", "study3_functions.R"))
source(file.path(repro_root, "R", "study4_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
flag_value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) {
    return(NULL)
  }
  if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
    stop(flag, " requires a value.", call. = FALSE)
  }
  args[[i + 1L]]
}

production <- "--production" %in% args
n_rep_arg <- flag_value("--n-rep")
workers_arg <- flag_value("--workers")
output_arg <- flag_value("--output-dir")
chunk_arg <- flag_value("--chunk-size")

consumed <- c(
  if (production) "--production",
  if (!is.null(n_rep_arg)) c("--n-rep", n_rep_arg),
  if (!is.null(workers_arg)) c("--workers", workers_arg),
  if (!is.null(output_arg)) c("--output-dir", output_arg),
  if (!is.null(chunk_arg)) c("--chunk-size", chunk_arg)
)
unknown <- setdiff(args, consumed)
if (length(unknown)) {
  stop("Unrecognized argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

config <- make_sim4_config("production")

if (production) {
  if (is.null(output_arg)) {
    stop(
      "--production requires --output-dir. Example: --output-dir reruns/study4_production",
      call. = FALSE
    )
  }
  if (!is.null(n_rep_arg)) {
    stop(
      "--production does not take --n-rep. The published run uses ",
      config$n_rep_total, " replications.",
      call. = FALSE
    )
  }
  if (!is.null(chunk_arg) && as.integer(chunk_arg) != as.integer(config$chunk_size)) {
    stop("--production uses chunk size ", config$chunk_size, ".", call. = FALSE)
  }
  n_rep <- as.integer(config$n_rep_total)
  chunk_size <- as.integer(config$chunk_size)
  output_dir <- output_arg
} else {
  if (is.null(n_rep_arg) || is.null(output_arg)) {
    stop(
      "Pass --production for the published 5,000-replication run, or ",
      "--n-rep and --output-dir for a check run.",
      call. = FALSE
    )
  }
  n_rep <- as.integer(n_rep_arg)
  chunk_size <- if (is.null(chunk_arg)) n_rep else as.integer(chunk_arg)
  output_dir <- output_arg
}

workers <- if (is.null(workers_arg)) 1L else as.integer(workers_arg)

if (!is.finite(n_rep) || n_rep < 1L) {
  stop("Replication count must be a positive integer.", call. = FALSE)
}
if (!is.finite(chunk_size) || chunk_size < 1L) {
  stop("Chunk size must be a positive integer.", call. = FALSE)
}
if (n_rep %% chunk_size != 0L) {
  stop("n_rep must be divisible by chunk size.", call. = FALSE)
}
if (!is.finite(workers) || workers < 1L) {
  stop("Worker count must be a positive integer.", call. = FALSE)
}

if (!grepl("^(/|[A-Za-z]:)", output_dir)) {
  output_dir <- file.path(repro_root, output_dir)
}
output_dir <- normalizePath(output_dir, mustWork = FALSE)
repro_prefix <- paste0(repro_root, "/")
if (!startsWith(output_dir, repro_prefix)) {
  stop("Output directory must be a subdirectory of reproducibility/.", call. = FALSE)
}
protected <- c(
  "results/study1", "results/study2", "results/study3", "results/study4",
  "data", "manuscript"
)
for (rel in protected) {
  root <- normalizePath(file.path(repro_root, rel), mustWork = FALSE)
  if (identical(output_dir, root) || startsWith(output_dir, paste0(root, "/"))) {
    stop(
      "Refusing to write into ", rel, ". Choose an output directory outside the archived results, data, and manuscript directories.",
      call. = FALSE
    )
  }
}

needed <- c(
  "metafor", "clubSandwich", "MASS", "withr", "future", "future.apply", "digest"
)
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!exists("fit_mlma_base_s2", mode = "function") ||
    !exists("extract_cr2_s2", mode = "function") ||
    !exists("compute_contrast_information_s3", mode = "function") ||
    !exists("predict_contrast_power_s3", mode = "function") ||
    !exists("comparator_quantities_s3", mode = "function") ||
    !exists("make_sampling_Sigma", mode = "function") ||
    !exists("make_sim2_config", mode = "function") ||
    !exists("run_with_sim3_seed", mode = "function")) {
  stop("Study 4 dependencies were not loaded from reproducibility/R.", call. = FALSE)
}

lib_check <- validate_frozen_scenario_library_s4(repro_root, config)
if (!isTRUE(lib_check$ok)) {
  stop(
    "Frozen Study 4 library failed validation:\n",
    paste(lib_check$issues, collapse = "\n"),
    call. = FALSE
  )
}
lib <- lib_check$library
grid <- lib$grid
if (nrow(grid) != as.integer(config$n_scenarios) || anyDuplicated(grid$scenario_id)) {
  stop("Frozen library must contain ", config$n_scenarios, " unique scenarios.", call. = FALSE)
}
required_grid <- c(
  "scenario_id", "scenario_set", "targeted_class", "pair_id", "pair_arm",
  "family_id", "tau_family_level", "J", "J0", "J1", "N_total", "K_total",
  "rho_true", "tau", "omega", "outcome_scenario_id", "I_C_S_star", "I_C_H_star",
  "R_H", "predicted_power_H_0.10", "delta_P_H_0.10"
)
required_structure <- c("scenario_id", "study", "X", "n_j", "k_j")
required_metrics <- c(
  "scenario_id", "regime", "variant_id", "tau_plan", "omega_plan", "rho_plan",
  "I_C_star", "predicted_power_0.10"
)
required_pairs <- c(
  "pair_class", "pair_id", "scenario_id_a", "scenario_id_b", "pass_all_tolerances"
)
required_families <- c(
  "family_id", "profile", "F_I", "delta_P_H_0.10",
  "scenario_id_tau_0.05", "scenario_id_tau_0.12",
  "scenario_id_tau_0.20", "scenario_id_tau_0.30"
)
missing_cols <- list(
  grid = setdiff(required_grid, names(grid)),
  structure = setdiff(required_structure, names(lib$study_structure)),
  metrics = setdiff(required_metrics, names(lib$metrics)),
  pairs = setdiff(required_pairs, names(lib$pair_diagnostics)),
  families = setdiff(required_families, names(lib$family_diagnostics))
)
missing_cols <- missing_cols[lengths(missing_cols) > 0L]
if (length(missing_cols)) {
  detail <- vapply(names(missing_cols), function(nm) {
    paste0(nm, ": ", paste(missing_cols[[nm]], collapse = ", "))
  }, character(1))
  stop("Frozen library is missing columns. ", paste(detail, collapse = "; "), call. = FALSE)
}
if (!setequal(grid$scenario_id, unique(lib$study_structure$scenario_id))) {
  stop("Study-structure scenario IDs do not match the scenario library.", call. = FALSE)
}
if (!setequal(grid$scenario_id, unique(lib$metrics$scenario_id))) {
  stop("Diagnostic-library scenario IDs do not match the scenario library.", call. = FALSE)
}
pair_ids <- c(lib$pair_diagnostics$scenario_id_a, lib$pair_diagnostics$scenario_id_b)
targeted_ids <- grid$scenario_id[grid$scenario_set == "targeted"]
if (!all(pair_ids %in% targeted_ids)) {
  stop("Targeted-pair diagnostics do not match the targeted scenarios.", call. = FALSE)
}
family_ids <- unlist(lib$family_diagnostics[, c(
  "scenario_id_tau_0.05", "scenario_id_tau_0.12",
  "scenario_id_tau_0.20", "scenario_id_tau_0.30"
), drop = FALSE])
class_e_ids <- grid$scenario_id[!is.na(grid$targeted_class) & grid$targeted_class == "E"]
if (!setequal(as.integer(family_ids), as.integer(class_e_ids))) {
  stop("Class E family diagnostics do not match the Class E scenarios.", call. = FALSE)
}
hash_n <- lib$hashes$sha256[lib$hashes$file == "n_scenarios"]
if (length(hash_n) != 1L || !identical(as.character(hash_n), as.character(config$n_scenarios))) {
  stop("design_hashes.csv does not record ", config$n_scenarios, " scenarios.", call. = FALSE)
}

exec_config <- config
exec_config$n_rep_total <- n_rep
exec_config$chunk_size <- chunk_size
exec_config$workers <- workers
if (!production) {
  exec_config$run_mode <- "check"
}

chunk_grid <- make_chunk_grid_s4(grid, n_rep, chunk_size)
coverage <- assert_chunk_coverage_s4(chunk_grid, n_rep, n_scenarios = nrow(grid))
if (!isTRUE(coverage$ok)) {
  stop("Chunk grid is incomplete. ", paste(coverage$issues, collapse = " "), call. = FALSE)
}

dir.create(file.path(output_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "summaries"), recursive = TRUE, showWarnings = FALSE)

run_start <- Sys.time()
message(
  if (production) "Study 4 production run. " else "Study 4 check run. ",
  "n_rep=", n_rep, " chunk_size=", chunk_size, " workers=", workers,
  " scenarios=", nrow(grid)
)
message("Output directory: ", output_dir)
message("Frozen design: ", lib$dir)

run_res <- run_chunks_s4(
  chunk_grid = chunk_grid,
  scenario_grid = grid,
  study_structure = lib$study_structure,
  root = repro_root,
  out_dir = output_dir,
  config = exec_config
)
n_failed <- sum(run_res$chunk_log$status == "failed")
if (n_failed > 0L) {
  msgs <- unique(na.omit(run_res$chunk_log$error_message[run_res$chunk_log$status == "failed"]))
  stop(n_failed, " chunk(s) failed. ", paste(msgs, collapse = " | "), call. = FALSE)
}

collected <- collect_sim4_chunk_results(output_dir, chunk_grid)
summaries <- write_sim4_summaries(
  out_dir = output_dir,
  config = exec_config,
  scenario_grid = grid,
  metric_library = lib$metrics,
  method_results = collected$method_results,
  pair_diagnostics = lib$pair_diagnostics,
  family_diagnostics = lib$family_diagnostics
)
unlink(file.path(output_dir, "raw_chunks"), recursive = TRUE)

invisible(write_study4_derived_tables(
  out_dir = file.path(output_dir, "summaries"),
  scenario_library = grid,
  condition_method_summary = summaries$method_summary,
  family_diagnostics = lib$family_diagnostics,
  diagnostic_library = lib$metrics
))

message(
  "Study 4 finished in ",
  round(as.numeric(difftime(Sys.time(), run_start, units = "secs")), 1),
  " seconds."
)
