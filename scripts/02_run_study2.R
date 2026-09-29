# Simulation Study 2 runner.
#
# Published run (5,000 replications, chunk size 500; requires an output directory):
#   Rscript scripts/02_run_study2.R --production --output-dir reruns/study2_production
#   Rscript scripts/02_run_study2.R --production --workers 14 --output-dir reruns/study2_production
#
# Execution check (does not write results/study2/):
#   Rscript scripts/02_run_study2.R --n-rep 1 --workers 1 --output-dir results/study2_check
#
# Paths are relative to reproducibility/ unless absolute and still inside it.
# Replication count and chunk size for --production come from
# make_sim2_config("production"), which matches the archived production run.

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
  normalizePath(file.path("scripts", "02_run_study2.R"), mustWork = TRUE)
}
repro_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(repro_root, "R", "study1_functions.R"))
source(file.path(repro_root, "R", "study2_design.R"))
source(file.path(repro_root, "R", "study2_functions.R"))
source(file.path(repro_root, "R", "estimators.R"))

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

config <- make_sim2_config("production")

if (production) {
  if (is.null(output_arg)) {
    stop(
      "--production requires --output-dir. Example: --output-dir reruns/study2_production",
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

needed <- c("metafor", "clubSandwich", "MASS", "withr", "future", "future.apply")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

condition_grid <- make_unique_run_grid_s2(config)
required_cols <- c(
  "condition_id", "condition_set", "J", "k_imbalance", "moderator_split",
  "n_imbalance", "beta1_true", "allocation_pattern"
)
missing_cols <- setdiff(required_cols, names(condition_grid))
if (length(missing_cols)) {
  stop("Design grid is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
n_conditions <- as.integer(config$unique_run_conditions)
if (nrow(condition_grid) != n_conditions) {
  stop(
    "Expected ", n_conditions, " design conditions; found ",
    nrow(condition_grid), ".",
    call. = FALSE
  )
}
ids <- as.integer(condition_grid$condition_id)
if (anyDuplicated(ids) || !identical(sort(ids), seq_len(n_conditions))) {
  stop("Condition IDs must be the unique values 1:", n_conditions, ".", call. = FALSE)
}
if (sum(condition_grid$condition_set == "core") != as.integer(config$core_conditions) ||
    sum(condition_grid$condition_set == "alignment") != as.integer(config$alignment_new_conditions)) {
  stop("Core and alignment condition counts do not match the production configuration.", call. = FALSE)
}
if (!identical(sort(unique(as.integer(condition_grid$J))), sort(as.integer(config$J_levels)))) {
  stop("Design J values do not match the production configuration.", call. = FALSE)
}
if (!identical(sort(unique(condition_grid$k_imbalance)), sort(config$k_imbalance_levels))) {
  stop("Design k_imbalance values do not match the production configuration.", call. = FALSE)
}
if (!identical(sort(unique(condition_grid$moderator_split)), sort(config$split_levels))) {
  stop("Design moderator splits do not match the production configuration.", call. = FALSE)
}
if (!identical(sort(unique(condition_grid$n_imbalance)), sort(config$n_imbalance_levels))) {
  stop("Design n_imbalance values do not match the production configuration.", call. = FALSE)
}
beta1_found <- sort(unique(as.numeric(condition_grid$beta1_true)))
beta1_expected <- sort(as.numeric(config$beta1_levels))
if (!isTRUE(all.equal(beta1_found, beta1_expected, tolerance = 1e-12))) {
  stop("Design beta1 values do not match the production configuration.", call. = FALSE)
}

align_file <- "condition_grid_alignment_analysis.csv"
align_path <- file.path(sim2_fixture_dir(repro_root), align_file)
if (!file.exists(align_path)) {
  stop("Missing alignment fixture: ", align_path, call. = FALSE)
}
alignment <- read_sim2_fixture(align_file, root = repro_root)
align_cols <- c(
  "stress_cell_id", "condition_id", "reuses_core_condition", "J", "k_imbalance",
  "moderator_split", "n_imbalance", "beta1_true", "allocation_pattern"
)
if (!identical(names(alignment), align_cols)) {
  stop("Alignment fixture columns do not match the expected structure.", call. = FALSE)
}
n_stress <- as.integer(config$stress_analysis_cells)
if (nrow(alignment) != n_stress ||
    !identical(sort(as.integer(alignment$stress_cell_id)), seq_len(n_stress))) {
  stop("Alignment fixture must contain stress cells 1:", n_stress, ".", call. = FALSE)
}
reuse_text <- tolower(as.character(alignment$reuses_core_condition))
if (!all(reuse_text %in% c("true", "false"))) {
  stop("Alignment fixture reuses_core_condition values must be true or false.", call. = FALSE)
}
align_keys <- c(
  "condition_id", "J", "k_imbalance", "moderator_split",
  "n_imbalance", "beta1_true", "allocation_pattern"
)
joined_ids <- merge(
  alignment,
  condition_grid,
  by = align_keys,
  all.x = TRUE,
  sort = FALSE
)
if (nrow(joined_ids) != n_stress || anyNA(joined_ids$condition_set)) {
  stop("Alignment fixture did not match ", n_stress, " production conditions.", call. = FALSE)
}

exec_config <- config
exec_config$n_rep_total <- n_rep
exec_config$chunk_size <- chunk_size
exec_config$workers <- workers
if (!production) {
  exec_config$run_mode <- "check"
}

chunk_grid <- make_chunk_grid_s2(
  condition_grid = condition_grid,
  n_rep_total = n_rep,
  chunk_size = chunk_size
)
n_chunks_expected <- as.integer(n_conditions * (n_rep / chunk_size))
if (nrow(chunk_grid) != n_chunks_expected) {
  stop("Chunk grid has ", nrow(chunk_grid), " rows; expected ", n_chunks_expected, ".", call. = FALSE)
}

dir.create(file.path(output_dir, "raw_chunks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

run_start <- Sys.time()
message(
  if (production) "Study 2 production run. " else "Study 2 check run. ",
  "n_rep=", n_rep, " chunk_size=", chunk_size, " workers=", workers,
  " conditions=", n_conditions
)
message("Output directory: ", output_dir)

run_res <- run_chunks_s2(
  chunk_grid = chunk_grid,
  root = repro_root,
  out_dir = output_dir,
  config = exec_config
)
n_failed <- sum(run_res$chunk_log$status == "failed")
if (n_failed > 0L) {
  msgs <- unique(na.omit(run_res$chunk_log$error_message[run_res$chunk_log$status == "failed"]))
  stop(
    n_failed, " chunk(s) failed. ",
    paste(msgs, collapse = " | "),
    call. = FALSE
  )
}

collected <- collect_sim2_chunk_results(output_dir, chunk_grid)
method_summary <- summarize_method_performance_s2(
  results = collected$method_results,
  config = exec_config,
  condition_grid = condition_grid
)
design_summary <- summarize_design_diagnostics_s2_by_condition(
  diagnostics = collected$dataset_diagnostics,
  condition_grid = condition_grid
)
alignment_summary <- summarize_alignment_stress_s2(
  method_summary = method_summary,
  root = repro_root
)

summaries_dir <- file.path(output_dir, "summaries")
utils::write.csv(
  method_summary,
  file.path(summaries_dir, "condition_method_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  design_summary,
  file.path(summaries_dir, "condition_design_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  alignment_summary,
  file.path(summaries_dir, "alignment_stress_summary.csv"),
  row.names = FALSE
)
mc_cols <- c(
  "condition_id", "method", "bias", "mcse_bias", "empirical_sd",
  "mcse_empirical_sd", "coverage", "mcse_coverage", "reject_rate",
  "mcse_reject_rate", "n_success"
)
mcse <- method_summary[, intersect(mc_cols, names(method_summary)), drop = FALSE]
utils::write.csv(
  mcse,
  file.path(summaries_dir, "monte_carlo_uncertainty.csv"),
  row.names = FALSE
)
invisible(write_sim2_tables(
  method_summary = method_summary,
  design_summary = design_summary,
  alignment_summary = alignment_summary,
  condition_grid = condition_grid,
  tables_dir = file.path(output_dir, "tables"),
  config = exec_config,
  figure_data_dir = NULL,
  condition_runtime = NULL
))

unlink(file.path(output_dir, "raw_chunks"), recursive = TRUE)

message(
  "Study 2 finished in ",
  round(as.numeric(difftime(Sys.time(), run_start, units = "secs")), 1),
  " seconds."
)
