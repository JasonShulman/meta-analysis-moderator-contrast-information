# Simulation Study 1 runner.
#
# Published run (5,000 replications; requires an output directory):
#   Rscript scripts/01_run_study1.R --production --output-dir reruns/study1_production
#
# Execution check (does not write results/study1/):
#   Rscript scripts/01_run_study1.R --n-rep 1 --workers 1 --output-dir results/study1_check
#
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
  normalizePath(file.path("scripts", "01_run_study1.R"), mustWork = TRUE)
}
repro_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(repro_root, "R", "study1_functions.R"))

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

J_LEVELS <- c(20L, 40L, 80L)
IMBALANCE_LEVELS <- c("balanced", "moderate", "severe")
SPLIT_LEVELS <- c("50/50", "70/30", "85/15")
BETA1_LEVELS <- c(0, 0.10, 0.30)
BETA0 <- 0.20
RHO <- 0.50
TAU <- 0.10
OMEGA <- 0.10
N_J <- 100
N_CONDITIONS <- 81L

if (production) {
  if (is.null(output_arg)) {
    stop(
      "--production requires --output-dir. Example: --output-dir reruns/study1_production",
      call. = FALSE
    )
  }
  if (!is.null(n_rep_arg)) {
    stop("--production does not take --n-rep. The published run uses 5,000 replications.", call. = FALSE)
  }
  if (!is.null(chunk_arg) && as.integer(chunk_arg) != 500L) {
    stop("--production uses chunk size 500.", call. = FALSE)
  }
  n_rep <- 5000L
  chunk_size <- 500L
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

needed <- c("metafor", "clubSandwich", "MASS", "future", "future.apply")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

condition_grid <- make_study1_condition_grid(
  J = J_LEVELS,
  imbalance = IMBALANCE_LEVELS,
  split = SPLIT_LEVELS,
  beta1 = BETA1_LEVELS,
  beta0 = BETA0,
  rho = RHO,
  tau = TAU,
  omega = OMEGA,
  n_j = N_J,
  n_rep = n_rep
)

required_cols <- c(
  "condition_id", "J", "imbalance", "split", "beta1",
  "beta0", "rho", "tau", "omega", "n_j", "n_rep"
)
missing_cols <- setdiff(required_cols, names(condition_grid))
if (length(missing_cols)) {
  stop("Design grid is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (nrow(condition_grid) != N_CONDITIONS) {
  stop("Expected ", N_CONDITIONS, " design conditions; found ", nrow(condition_grid), ".", call. = FALSE)
}
if (anyDuplicated(condition_grid$condition_id) ||
    !identical(sort(as.integer(condition_grid$condition_id)), seq_len(N_CONDITIONS))) {
  stop("Condition IDs must be the unique values 1:", N_CONDITIONS, ".", call. = FALSE)
}
if (!identical(sort(unique(as.integer(condition_grid$J))), J_LEVELS)) {
  stop("Design J values do not match the published Study 1 levels.", call. = FALSE)
}
if (!identical(sort(unique(as.character(condition_grid$imbalance))), sort(IMBALANCE_LEVELS))) {
  stop("Design imbalance values do not match the published Study 1 levels.", call. = FALSE)
}
if (!identical(sort(unique(as.character(condition_grid$split))), sort(SPLIT_LEVELS))) {
  stop("Design splits do not match the published Study 1 levels.", call. = FALSE)
}
beta1_ok <- length(unique(condition_grid$beta1)) == length(BETA1_LEVELS) &&
  all(vapply(BETA1_LEVELS, function(b) any(abs(condition_grid$beta1 - b) < 1e-12), logical(1)))
if (!beta1_ok) {
  stop("Design beta1 values must be 0, 0.10, and 0.30.", call. = FALSE)
}
if (!isTRUE(all.equal(unique(condition_grid$beta0), BETA0)) ||
    !isTRUE(all.equal(unique(condition_grid$rho), RHO)) ||
    !isTRUE(all.equal(unique(condition_grid$tau), TAU)) ||
    !isTRUE(all.equal(unique(condition_grid$omega), OMEGA)) ||
    !isTRUE(all.equal(unique(as.numeric(condition_grid$n_j)), as.numeric(N_J)))) {
  stop("Design beta0, rho, tau, omega, or n_j does not match the published Study 1 values.", call. = FALSE)
}
if (any(as.integer(condition_grid$n_rep) != n_rep)) {
  stop("Every condition must use the requested replication count.", call. = FALSE)
}

chunk_grid <- make_chunk_grid(
  condition_grid = condition_grid,
  n_rep_total = n_rep,
  chunk_size = chunk_size
)
n_chunks_expected <- as.integer(N_CONDITIONS * (n_rep / chunk_size))
if (nrow(chunk_grid) != n_chunks_expected) {
  stop("Chunk grid has ", nrow(chunk_grid), " rows; expected ", n_chunks_expected, ".", call. = FALSE)
}

ensure_sim_study_dirs(output_dir)
extended_dir <- file.path(output_dir, "extended_summaries")
dir.create(extended_dir, recursive = TRUE, showWarnings = FALSE)

run_start <- Sys.time()
config <- data.frame(
  parameter = c(
    "run_mode", "n_rep_total", "chunk_size", "workers", "n_conditions",
    "j_levels", "imbalance_levels", "split_levels", "beta1_levels",
    "beta0", "rho", "tau", "omega", "n_j", "shuffle_kj", "seed_scheme"
  ),
  value = c(
    if (production) "production" else "check",
    as.character(n_rep),
    as.character(chunk_size),
    as.character(workers),
    as.character(nrow(condition_grid)),
    paste(J_LEVELS, collapse = ","),
    paste(IMBALANCE_LEVELS, collapse = ","),
    paste(SPLIT_LEVELS, collapse = ","),
    paste(BETA1_LEVELS, collapse = ","),
    as.character(BETA0),
    as.character(RHO),
    as.character(TAU),
    as.character(OMEGA),
    as.character(N_J),
    "TRUE",
    "condition_seed_base(condition_id) + rep_id (chunk-invariant)"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(config, file.path(output_dir, "run_config.csv"), row.names = FALSE)
utils::write.csv(chunk_grid, file.path(output_dir, "logs", "chunk_grid.csv"), row.names = FALSE)

message(
  if (production) "Study 1 production run. " else "Study 1 check run. ",
  "n_rep=", n_rep, " chunk_size=", chunk_size, " workers=", workers
)
message("Output directory: ", output_dir)

run_res <- run_chunks_parallel(
  chunk_grid = chunk_grid,
  root = repro_root,
  out_dir = output_dir,
  workers = workers,
  overwrite = TRUE
)
utils::write.csv(
  run_res$chunk_log,
  file.path(output_dir, "sim_study_1_chunk_log.csv"),
  row.names = FALSE
)
n_failed <- sum(run_res$chunk_log$status == "failed")
if (n_failed > 0L) {
  stop(n_failed, " chunk(s) failed. See sim_study_1_chunk_log.csv.", call. = FALSE)
}

raw_all <- combine_chunk_rds(output_dir, chunk_grid = chunk_grid)
summary_tab <- summarize_study1_performance(raw_all)
utils::write.csv(
  summary_tab,
  file.path(output_dir, "sim_study_1_summary.csv"),
  row.names = FALSE
)
condition_log <- aggregate_condition_runtime(run_res$chunk_log)
invisible(write_sim_study_tables(
  summary_tab = summary_tab,
  condition_grid = condition_grid,
  condition_log = condition_log,
  tables_dir = file.path(output_dir, "tables"),
  figure_data_dir = NULL
))

seed_lookup <- stats::setNames(
  vapply(seq_len(N_CONDITIONS), condition_seed_base, integer(1)),
  as.character(seq_len(N_CONDITIONS))
)
write_study1_extended_summaries(raw_all, extended_dir, seed_lookup)

message(
  "Study 1 finished in ",
  round(as.numeric(difftime(Sys.time(), run_start, units = "secs")), 1),
  " seconds."
)
