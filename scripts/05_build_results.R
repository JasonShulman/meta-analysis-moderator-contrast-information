# Reconstruct derived manuscript quantities from archived summaries.
#
# This step does not simulate data or fit models. It reads the archived
# Study 1–4 results and the frozen Study 3–4 designs, then writes only
# reproducibility/results/derived/.
#
# Study 1 extended summaries and Study 2 summaries are checked and left in
# place. Study 3 heterogeneity sensitivity and Study 4 Class E / planning-range
# tables are recomputed.
#
#   Rscript scripts/05_build_results.R

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]))
} else {
  normalizePath(file.path("scripts", "05_build_results.R"), mustWork = TRUE)
}
repro_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repro_root, "R", "contrast_information.R"))
source(file.path(repro_root, "R", "study3_functions.R"))
source(file.path(repro_root, "R", "study4_functions.R"))

results_root <- file.path(repro_root, "results")
derived_root <- file.path(results_root, "derived")
study3_out <- file.path(derived_root, "study3_heterogeneity")
study4_out <- file.path(derived_root, "study4")

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing archived input: ", path, call. = FALSE)
  }
  path
}

require_columns <- function(df, columns, label) {
  missing <- setdiff(columns, names(df))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(df)
}

require_unique <- function(df, columns, label) {
  if (anyDuplicated(df[, columns, drop = FALSE])) {
    stop(label, " has duplicate rows for ", paste(columns, collapse = ", "), ".", call. = FALSE)
  }
  invisible(df)
}

archived_inputs <- c(
  file.path(results_root, "study1", "sim_study_1_summary.csv"),
  file.path(results_root, "study1", "extended_summaries", "sim_study_1_summary_extended.csv"),
  file.path(results_root, "study1", "extended_summaries", "sim_study_1_failure_diagnostics.csv"),
  file.path(results_root, "study1", "extended_summaries", "table_se_calibration.csv"),
  file.path(results_root, "study1", "extended_summaries", "table_monte_carlo_uncertainty.csv"),
  file.path(results_root, "study2", "summaries", "condition_method_summary.csv"),
  file.path(results_root, "study2", "summaries", "alignment_stress_summary.csv"),
  file.path(results_root, "study3", "summaries", "condition_method_summary.csv"),
  file.path(results_root, "study4", "summaries", "condition_method_summary.csv")
)
invisible(lapply(archived_inputs, require_file))

cms3_path <- file.path(results_root, "study3", "summaries", "condition_method_summary.csv")
cms4_path <- file.path(results_root, "study4", "summaries", "condition_method_summary.csv")
cms3 <- utils::read.csv(cms3_path, stringsAsFactors = FALSE, check.names = FALSE)
cms4 <- utils::read.csv(cms4_path, stringsAsFactors = FALSE, check.names = FALSE)

lib3 <- read_frozen_scenario_library_s3(repro_root)
lib4 <- read_frozen_scenario_library_s4(repro_root)
n3 <- as.integer(make_sim3_config("production")$n_broad + make_sim3_config("production")$n_targeted)
n4 <- as.integer(make_sim4_config("production")$n_scenarios)

require_columns(
  lib3$grid,
  c("scenario_id", "scenario_set"),
  "Frozen Study 3 scenario grid"
)
require_columns(
  lib3$study_structure,
  c("scenario_id", "study", "X", "n_j", "k_j"),
  "Frozen Study 3 study structure"
)
require_columns(
  cms3,
  c(
    "scenario_id", "method", "beta1_true", "power", "n_attempted", "n_success",
    "error_rate", "empirical_variance", "mean_df"
  ),
  "Archived Study 3 condition-method summary"
)
require_unique(lib3$grid, "scenario_id", "Frozen Study 3 scenario grid")
require_unique(cms3, c("scenario_id", "beta1_true", "method"), "Archived Study 3 condition-method summary")
if (nrow(lib3$grid) != n3) {
  stop("Frozen Study 3 library has ", nrow(lib3$grid), " scenarios; expected ", n3, ".", call. = FALSE)
}
if (!identical(sort(unique(as.integer(cms3$scenario_id))), sort(as.integer(lib3$grid$scenario_id)))) {
  stop("Study 3 scenario IDs differ between the archived summary and the frozen library.", call. = FALSE)
}
if (!setequal(as.integer(lib3$grid$scenario_id), as.integer(lib3$study_structure$scenario_id))) {
  stop("Study 3 study-structure IDs do not match the frozen scenario grid.", call. = FALSE)
}

require_columns(
  lib4$grid,
  c(
    "scenario_id", "targeted_class", "family_id", "tau", "omega",
    "delta_P_H_0.10", "predicted_power_tau_L_0.10", "predicted_power_tau_U_0.10",
    "predicted_power_H_0.10"
  ),
  "Frozen Study 4 scenario library"
)
require_columns(
  lib4$family_diagnostics,
  c(
    "family_id", "profile", "profile_label", "base_candidate_id", "F_I",
    "delta_P_H_0.10", "omega_family",
    "scenario_id_tau_0.05", "scenario_id_tau_0.12",
    "scenario_id_tau_0.20", "scenario_id_tau_0.30"
  ),
  "Frozen Study 4 family diagnostics"
)
require_columns(
  lib4$metrics,
  c("scenario_id", "regime", "tau_plan", "omega_plan", "rho_plan", "I_C_star", "predicted_power_0.10"),
  "Frozen Study 4 diagnostic library"
)
require_columns(
  cms4,
  c("scenario_id", "method", "beta1_true", "power", "family_id", "tau_family_level", "empirical_variance", "n_success", "n_attempted"),
  "Archived Study 4 condition-method summary"
)
require_unique(lib4$grid, "scenario_id", "Frozen Study 4 scenario library")
require_unique(lib4$family_diagnostics, "family_id", "Frozen Study 4 family diagnostics")
require_unique(cms4, c("scenario_id", "beta1_true", "method"), "Archived Study 4 condition-method summary")
if (nrow(lib4$grid) != n4) {
  stop("Frozen Study 4 library has ", nrow(lib4$grid), " scenarios; expected ", n4, ".", call. = FALSE)
}
if (!identical(sort(unique(as.integer(cms4$scenario_id))), sort(as.integer(lib4$grid$scenario_id)))) {
  stop("Study 4 scenario IDs differ between the archived summary and the frozen library.", call. = FALSE)
}

# Replace the derived directories so a rerun cannot mix stale files.
unlink(study3_out, recursive = TRUE)
unlink(study4_out, recursive = TRUE)
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)

message("Recomputing Study 3 heterogeneity sensitivity.")
invisible(run_heterogeneity_sensitivity_analysis_s3(
  root = repro_root,
  method_summary = cms3,
  out_dir = study3_out
))

message("Recomputing Study 4 Class E and planning-range tables.")
invisible(write_study4_derived_tables(
  out_dir = study4_out,
  scenario_library = lib4$grid,
  condition_method_summary = cms4,
  family_diagnostics = lib4$family_diagnostics,
  diagnostic_library = lib4$metrics
))

message("Derived results written to ", derived_root)
