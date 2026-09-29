#!/usr/bin/env Rscript
# Manuscript and supplement displays from the public numerical record.
#
# Presentation only. Reads archived Study 1-4 summaries, results/derived/,
# and frozen design files. Does not simulate, refit, or rewrite those inputs.
#
# Display map (include -> numerical source -> output under manuscript/):
#
#   table1_simulation_program
#     Study 1 run_config + Study 2-4 production configs
#     + condition/scenario counts in the archived summaries
#     -> tables/table1_simulation_program.tex
#   table2_information_profile_landscape
#     Study 4 planning constants (sim4_tau_L, sim4_tau_U, sim4_fragility_omega,
#     sim4_rho_plan_primary, sim4_alpha_default, confirmatory beta1)
#     -> tables/table2_information_profile_landscape.tex
#   table3_study1_estimator_summary
#     Study 1 extended summary, unweighted condition means
#     -> tables/table3_study1_estimator_summary.tex
#   table4_study4_confirmation_revised
#     Study 4 power_calibration_summary (broad primary)
#     + retention, support, comparator, variance, targeted-pair, Type I summaries
#     + derived Class E fragility
#     -> tables/table4_study4_confirmation_revised.tex
#   figure1b_study1_df_and_power.pdf
#     Study 1 extended summary, MLMA+CR2 means by J and moderator allocation
#   figure2c_study2_alignment_j40.pdf
#     Study 2 alignment_stress_summary, J = 40
#   figure3_power_calibration.pdf
#     Study 3 metric-variant power by scenario (sampling-only and default heterogeneity)
#     + Study 4 scenario library predictions joined to MLMA+CR2 power
#   tableS1  Study 3/4 scenario libraries and study-structure rows
#   tableS2  Study 1 make_kj_template()
#   tableS3  Study 2 table_n_templates.csv
#   tableS4  Study 1 extended summary and Monte Carlo intervals
#   tableS5  Study 2 condition_method_summary, core factorial and one stress cell
#   tableS6  Study 2 alignment_stress_summary, null rows
#   tableS7  Study 3 variance_prediction_summary and df_prediction_summary
#   tableS8  Study 3 metric_variant_power_prediction_summary, all non-null
#   tableS9  Study 3 targeted_pair_summary, power joined from condition_method_summary
#   tableS10 derived Study 3 heterogeneity_misspecification_summary, absolute tau, rho = .50
#   tableS11 Study 4 power_calibration_summary, excluding the repeated robustness-band rows
#   tableS12 Study 4 variance_calibration_summary
#   tableS13 Study 4 targeted_pair_summary
#   tableS14 Study 4 comparator_performance_summary and support_validation_summary
#   tableS15 derived Study 4 heterogeneity_sensitivity_range_documentation
#   tableS16 derived Study 4 fragility_family_summary_derived; profile labels from the official family file
#   tableS17 Monte Carlo uncertainty and failure counts in the Study 1-4 summaries
#   sfig1  Study 1 table_type1_error.csv
#   sfig3  derived Study 3 heterogeneity_sensitivity_summary, planning cells at rho = .50
#   sfig5  derived Class E family powers at tau = .05, .12, .20, .30
#
# Output names follow the current main-text includes, including
# table2_information_profile_landscape and table4_study4_confirmation_revised.

options(warn = 1, scipen = 999)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run via Rscript reproducibility/scripts/06_build_manuscript_outputs.R", call. = FALSE)
}
repro <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
man_out <- file.path(repro, "manuscript")
tab_dir <- file.path(man_out, "tables")
fig_dir <- file.path(man_out, "figures")
stab_dir <- file.path(man_out, "supplement_tables")
sfig_dir <- file.path(man_out, "supplement_figures")
for (d in c(tab_dir, fig_dir, stab_dir, sfig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

source(file.path(repro, "R", "study1_functions.R"))
source(file.path(repro, "R", "study2_design.R"))
source(file.path(repro, "R", "study3_functions.R"))
source(file.path(repro, "R", "study4_functions.R"))

read_pub <- function(rel) {
  path <- file.path(repro, rel)
  if (!file.exists(path)) stop("Missing public input: ", rel, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "NaN", ""))
}

write_tex <- function(dir, name, lines) {
  path <- file.path(dir, name)
  if (!startsWith(normalizePath(dirname(path), mustWork = FALSE), normalizePath(man_out))) {
    stop("Refusing to write outside manuscript/: ", path, call. = FALSE)
  }
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

fmt_dec <- function(x, digits) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop("fmt_dec expected one finite number.", call. = FALSE)
  }
  sign <- if (x < 0) "-" else ""
  body <- formatC(abs(x), format = "f", digits = digits)
  if (abs(x) < 1) body <- sub("^0\\.", ".", body)
  paste0(sign, body)
}

fmt_fixed <- function(x, digits) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

fmt_math <- function(x, digits) {
  s <- fmt_dec(x, digits)
  if (startsWith(s, "-")) paste0("$-", substring(s, 2L), "$") else s
}

fmt_signed <- function(x, digits) {
  s <- fmt_dec(x, digits)
  if (startsWith(s, "-") || s %in% c("0", paste0(".", paste(rep("0", digits), collapse = "")))) {
    s
  } else {
    paste0("+", s)
  }
}

fmt_int <- function(x) {
  format(as.integer(round(as.numeric(x))), big.mark = "{,}", scientific = FALSE, trim = TRUE)
}

fmt_trim <- function(x, digits) {
  x <- as.numeric(x)
  s <- formatC(abs(x), format = "f", digits = digits)
  s <- sub("(\\.\\d*?)0+$", "\\1", s)
  s <- sub("\\.$", "", s)
  if (x < 0) s <- paste0("-", s)
  if (abs(x) < 1 && x != 0) s <- sub("^(-?)0\\.", "\\1.", s)
  add_thousands(s)
}

add_thousands <- function(s) {
  sign <- ""
  if (startsWith(s, "-")) {
    sign <- "-"
    s <- substring(s, 2L)
  }
  parts <- strsplit(s, ".", fixed = TRUE)[[1L]]
  int <- parts[[1L]]
  if (nchar(int) > 3L) {
    int <- gsub("(\\d)(?=(\\d{3})+$)", "\\1{,}", int, perl = TRUE)
  }
  if (length(parts) == 2L) paste0(sign, int, ".", parts[[2L]]) else paste0(sign, int)
}

fmt_beta <- function(bits) {
  paste(vapply(bits, function(x) {
    if (abs(x) < 1e-12) "$0$" else paste0("$", fmt_dec(x, 2L), "$")
  }, character(1)), collapse = ", ")
}

fmt_ic <- function(x) {
  x <- as.numeric(x)
  if (abs(x - round(x)) < 1e-6) {
    add_thousands(format(round(x), scientific = FALSE, trim = TRUE))
  } else {
    add_thousands(formatC(x, format = "f", digits = 1))
  }
}

mean_of <- function(x) mean(as.numeric(x), na.rm = TRUE)
cfg_vec <- function(x) as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])

method_label <- c(
  naive_re = "Independent RE",
  aggregate_hksj = "Aggregation + KH",
  mlma = "MLMA",
  mlma_cr2 = "MLMA+CR2",
  che_rve = "CHE-RVE"
)

# ---- Public inputs ----------------------------------------------------------

s1_cfg <- read_pub("results/study1/run_config.csv")
s1_cfg <- stats::setNames(s1_cfg$value, s1_cfg$parameter)
s1_ext <- read_pub("results/study1/extended_summaries/sim_study_1_summary_extended.csv")
s1_sum <- read_pub("results/study1/sim_study_1_summary.csv")
s1_mc <- read_pub("results/study1/extended_summaries/table_monte_carlo_uncertainty.csv")
s1_type1 <- read_pub("results/study1/tables/table_type1_error.csv")

s2_cfg <- make_sim2_config("production")
s2_cm <- read_pub("results/study2/summaries/condition_method_summary.csv")
s2_align <- read_pub("results/study2/summaries/alignment_stress_summary.csv")
s2_n <- read_pub("results/study2/tables/table_n_templates.csv")
s2_mc <- read_pub("results/study2/summaries/monte_carlo_uncertainty.csv")

s3_cfg <- make_sim3_config("production")
s3_scen <- read_pub("results/study3/summaries/scenario_design_summary.csv")
s3_ss <- read_pub("data/study3_design/scenario_study_structure.csv")
s3_cm <- read_pub("results/study3/summaries/condition_method_summary.csv")
s3_var <- read_pub("results/study3/summaries/variance_prediction_summary.csv")
s3_df <- read_pub("results/study3/summaries/df_prediction_summary.csv")
s3_pow <- read_pub("results/study3/summaries/metric_variant_power_prediction_summary.csv")
s3_by <- read_pub("results/study3/summaries/metric_variant_power_prediction_by_scenario.csv")
s3_pairs <- read_pub("results/study3/summaries/targeted_pair_summary.csv")
s3_mc <- read_pub("results/study3/summaries/monte_carlo_uncertainty.csv")
s3_het <- read_pub("results/derived/study3_heterogeneity/heterogeneity_sensitivity_summary.csv")
s3_miss <- read_pub("results/derived/study3_heterogeneity/heterogeneity_misspecification_summary.csv")

s4_cfg <- make_sim4_config("production")
s4_lib <- read_pub("data/study4_design/scenario_library.csv")
s4_ss <- read_pub("data/study4_design/study_structure.csv")
s4_cm <- read_pub("results/study4/summaries/condition_method_summary.csv")
s4_pow <- read_pub("results/study4/summaries/power_calibration_summary.csv")
s4_var <- read_pub("results/study4/summaries/variance_calibration_summary.csv")
s4_pairs <- read_pub("results/study4/summaries/targeted_pair_summary.csv")
s4_comp <- read_pub("results/study4/summaries/comparator_performance_summary.csv")
s4_sup <- read_pub("results/study4/summaries/support_validation_summary.csv")
s4_t1 <- read_pub("results/study4/summaries/type1_robustness_summary.csv")
s4_ret <- read_pub("results/study4/summaries/retention_validation_summary.csv")
s4_frag_off <- read_pub("results/study4/summaries/fragility_family_summary.csv")
s4_frag <- read_pub("results/derived/study4/fragility_family_summary_derived.csv")
s4_range <- read_pub("results/derived/study4/heterogeneity_sensitivity_range_documentation.csv")
s4_mc <- read_pub("results/study4/summaries/monte_carlo_uncertainty.csv")

stopifnot(
  as.integer(s1_cfg[["n_conditions"]]) == 81L,
  as.integer(s1_cfg[["n_rep_total"]]) == 5000L,
  s2_cfg$unique_run_conditions == 232L,
  s2_cfg$core_conditions == 216L,
  s2_cfg$alignment_new_conditions == 16L,
  s2_cfg$n_rep_total == 5000L,
  s3_cfg$n_broad + s3_cfg$n_targeted == nrow(s3_scen),
  s4_cfg$n_scenarios == nrow(s4_lib),
  s3_cfg$n_rep_total == 5000L,
  s4_cfg$n_rep_total == 5000L
)

# ---- Table 1 ----------------------------------------------------------------

j_tex <- function(x, sep = ", ") paste(as.integer(x), collapse = sep)
range_tex <- function(x) paste0(fmt_dec(x[[1L]], 2L), "$--$", fmt_dec(x[[2L]], 2L))

s1_beta <- cfg_vec(s1_cfg[["beta1_levels"]])
s1_J <- cfg_vec(s1_cfg[["j_levels"]])
n1 <- length(unique(s1_sum$condition_id))
if (n1 != as.integer(s1_cfg[["n_conditions"]])) {
  stop("Study 1 summary condition count does not match run_config.", call. = FALSE)
}
if (length(unique(s2_cm$condition_id)) != s2_cfg$unique_run_conditions) {
  stop("Study 2 summary condition count does not match the production config.", call. = FALSE)
}

t1_note <- paste0(
  "  \\textit{Note.} For Studies 3 and 4, reported $\\rho$, $\\tau$, and $\\omega$ values are rounded scenario-generation ranges; exact realized ranges are provided in the supplement. ",
  "Study 4 confirms the diagnostics at $\\beta_1=", fmt_dec(s4_cfg$confirmatory_beta1, 2L),
  "$ and uses the prespecified Type I error band ",
  fmt_dec(s4_cfg$type1_band[[1L]], 3L), "--", fmt_dec(s4_cfg$type1_band[[2L]], 3L), ". ",
  "RE, random effects; KH, Knapp--Hartung; MLMA, three-level multilevel meta-analysis; CR2, cluster-robust variance estimation with Satterthwaite degrees of freedom; CHE-RVE, correlated-and-hierarchical-effects robust variance estimation."
)

t1 <- c(
  "\\begin{sidewaystable}",
  "\\centering",
  "  \\caption{Simulation program}",
  "  \\label{tab:simulation-program}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{0.45em}",
  "  \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}cp{0.16\\linewidth}p{0.11\\linewidth}p{0.18\\linewidth}p{0.42\\linewidth}@{}}",
  "    \\toprule",
  "    Study & Design size & $\\beta_1$ & Estimators & Design \\\\",
  "    \\midrule",
  paste0(
    "    1 & ", n1, " conditions $\\times$ ", fmt_int(s1_cfg[["n_rep_total"]]),
    " & ", fmt_beta(s1_beta),
    " & Independent RE, MLMA, MLMA+CR2 & Dependence and moderator imbalance at a fixed study size. \\newline \\textit{Varied:} number of studies $J\\in\\{",
    j_tex(s1_J),
    "\\}$; effect-count balance (balanced, moderate, severe); moderator allocation (50/50, 70/30, 85/15). \\newline \\textit{Held fixed:} $\\beta_0=",
    fmt_dec(s1_cfg[["beta0"]], 2L), "$; $\\rho=", fmt_dec(s1_cfg[["rho"]], 2L),
    "$; $\\tau=", fmt_dec(s1_cfg[["tau"]], 2L), "$; $\\omega=", fmt_dec(s1_cfg[["omega"]], 2L),
    "$; $n_j=", as.integer(s1_cfg[["n_j"]]), "$. \\\\"
  ),
  "    \\addlinespace",
  paste0(
    "    2 & ", s2_cfg$unique_run_conditions, " conditions (", s2_cfg$core_conditions,
    " core $+$ ", s2_cfg$alignment_new_conditions, " alignment) $\\times$ ",
    fmt_int(s2_cfg$n_rep_total), " & ", fmt_beta(s2_cfg$beta1_levels),
    " & Aggregation + KH, MLMA, MLMA+CR2, CHE-RVE & Sample-size imbalance and alignment of informative studies. \\newline",
    "    \\textit{Varied:} $J\\in\\{", j_tex(s2_cfg$J_levels, ","),
    "\\}$; effect-count balance; moderator allocation (50/50, 70/30, 85/15); sample-size balance; alignment of large-$n$/high-$k$ studies with moderator group. \\newline",
    "    \\textit{Held fixed:} $\\beta_0=", fmt_dec(s2_cfg$beta0, 2L),
    "$; $\\rho=", fmt_dec(s2_cfg$rho_true, 2L),
    "$; $\\tau=\\omega=", fmt_dec(s2_cfg$tau, 2L),
    "$; mean study size ", as.integer(s2_cfg$mean_n), ". \\\\"
  ),
  "    \\addlinespace",
  paste0(
    "    3 & ", nrow(s3_scen), " scenarios (", s3_cfg$n_broad, " broad $+$ ", s3_cfg$n_targeted,
    " targeted) $\\times$ ", fmt_int(s3_cfg$n_rep_total), " & ", fmt_beta(s3_cfg$beta1_levels),
    " & Aggregation + KH, MLMA, MLMA+CR2, CHE-RVE & Development of the contrast-specific diagnostics. \\newline",
    "    \\textit{Varied:} $J\\in\\{", j_tex(s3_cfg$J_levels, ","),
    "\\}$; $\\rho\\approx", range_tex(s3_cfg$rho_true_range),
    "$; $\\tau,\\omega\\approx", range_tex(s3_cfg$tau_range),
    "$; $n$--$k$ relationship; alignment with moderator group. \\newline",
    "    \\textit{Held fixed:} $\\beta_0=", fmt_dec(s3_cfg$beta0, 2L),
    "$; planning $\\rho=", fmt_dec(s3_cfg$rho_plan_primary, 2L),
    "$; $\\alpha=", fmt_dec(s3_cfg$alpha, 2L), "$. \\\\"
  ),
  "    \\addlinespace",
  paste0(
    "    4 & ", nrow(s4_lib), " scenarios (", s4_cfg$n_broad, " broad $+$ ", s4_cfg$n_targeted,
    " targeted) $\\times$ ", fmt_int(s4_cfg$n_rep_total), " & ", fmt_beta(s4_cfg$beta1_levels),
    " & MLMA, MLMA+CR2 & Prospective confirmation on a new scenario library. \\newline",
    "    \\textit{Varied:} $J\\in\\{", j_tex(s4_cfg$J_levels, ","),
    "\\}$; $\\rho\\approx", range_tex(s4_cfg$rho_true_range),
    "$; $\\tau,\\omega\\approx", range_tex(s4_cfg$tau_range),
    "$; $n$--$k$ relationship; alignment with moderator group. \\newline",
    "    \\textit{Held fixed:} $\\beta_0=", fmt_dec(s4_cfg$beta0, 2L),
    "$; planning $\\rho=", fmt_dec(s4_cfg$rho_plan_primary, 2L),
    "$; $\\alpha=", fmt_dec(s4_cfg$alpha, 2L), "$. \\\\"
  ),
  "    \\bottomrule",
  "  \\end{tabular*}%",
  "  \\vspace{0.4em}",
  "  {\\footnotesize\\raggedright",
  t1_note,
  "  \\par}",
  "\\end{sidewaystable}"
)
if (!isTRUE(all.equal(s3_cfg$tau_range, s3_cfg$omega_range)) || !isTRUE(all.equal(s4_cfg$tau_range, s4_cfg$omega_range))) {
  stop("Table 1 combines tau and omega ranges; those ranges are not equal.", call. = FALSE)
}
if (!isTRUE(all.equal(s2_cfg$tau, s2_cfg$omega))) {
  stop("Table 1 prints one Study 2 heterogeneity value for tau and omega.", call. = FALSE)
}
write_tex(tab_dir, "table1_simulation_program.tex", t1)

# ---- Table 2 ----------------------------------------------------------------

het_step <- unique(round(diff(sort(sim4_het_grid)), 6L))
if (length(het_step) != 1L) stop("Study 4 heterogeneity grid is not evenly spaced.", call. = FALSE)
if (!isTRUE(all.equal(sim4_fragility_rho, sim4_rho_plan_primary))) {
  stop("Fragility correlation and primary planning correlation differ.", call. = FALSE)
}
d2 <- function(x) fmt_dec(x, 2L)
t2_def <- paste0(
  "The confirmatory power width is $\\Delta P_H(", d2(s4_cfg$confirmatory_beta1), ")$, ",
  "the change in predicted power for $\\beta_1=", d2(s4_cfg$confirmatory_beta1),
  "$ as $\\tau_p$ moves from ", d2(sim4_tau_L), " to ", d2(sim4_tau_U),
  " with $\\omega_p=", d2(sim4_fragility_omega), "$ and $\\rho_p=", d2(sim4_fragility_rho), "$."
)
t2_note <- paste0(
  "  \\textit{Note.} Every quantity is a function of the study design and the planning values $(\\tau_p,\\omega_p,\\rho_p)$. None is computed from simulated outcomes, fitted models, or realized CR2 degrees of freedom. ",
  "The heterogeneity-adjusted study weight is $w_j^{(H)}=[\\tau_p^2+\\omega_p^2/k_j+\\bar{v}_j]^{-1}$, where $\\bar{v}_j=[1+(k_j-1)\\rho_p]/[k_j(n_j-3)]$ is the sampling variance of the study-mean Fisher $z$. ",
  "The sampling-only weight is $w_j^{(S)}=1/\\bar{v}_j$ and does not use $\\tau_p$ or $\\omega_p$. ",
  "Group weight $W_g$ sums the chosen study weight over studies with moderator value $g$. ",
  "$I_{C,H}^*$ uses heterogeneity-adjusted weights and $I_{C,S}^*$ uses sampling-only weights, both at planning correlation $",
  d2(sim4_rho_plan_primary), "$. ",
  "$B_C=1$ when $W_0=W_1$, and $0<B_C\\le 1$ when both groups have positive weight. ",
  "$J_g^*$ equals the number of studies in the group when their weights are equal, and it approaches 1 when one study carries the group's weight. ",
  "$J_{\\min}^*=\\min(J_0^*,J_1^*)$. ",
  "$\\nu_C^*$ is an anticipated Satterthwaite-type support quantity, not the CR2 degrees of freedom returned by a fitted model. ",
  "The predicted standard error is $\\mathrm{SE}_C^*=2/\\sqrt{I_C^*}$. ",
  "Predicted power is a two-sided noncentral $t$ probability at level $\\alpha=", d2(sim4_alpha_default),
  "$, with noncentrality $\\delta_C^*=\\beta_1/\\mathrm{SE}_C^*$. ",
  "$\\Delta P_H(", d2(s4_cfg$confirmatory_beta1), ")$ is predicted power at $\\tau_p=", d2(sim4_tau_L),
  "$ minus predicted power at $\\tau_p=", d2(sim4_tau_U), "$, with $\\omega_p=", d2(sim4_fragility_omega),
  "$ and $\\rho_p=", d2(sim4_fragility_rho), "$ held fixed. ",
  "The companion information ratio on that same contrast is $F_I=1-I_{C,H}^*(\\tau_p=", d2(sim4_tau_U),
  ")/I_{C,H}^*(\\tau_p=", d2(sim4_tau_L), ")$."
)
t2 <- c(
  "\\begin{table*}",
  "  \\caption{Contrast-specific information profile.}",
  "  \\label{tab:information-profile}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{0.5em}",
  "  \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}p{0.20\\linewidth}p{0.26\\linewidth}p{0.46\\linewidth}@{}}",
  "    \\toprule",
  "    Feature & Question & Definition \\\\",
  "    \\midrule",
  "    Information magnitude & How much weight supports the contrast in total? & $T_C = W_0 + W_1$. Magnitude grows when either moderator group gains weight. \\\\",
  "    \\addlinespace",
  "    Balance & How evenly is that weight divided between the two moderator groups? & $B_C = 4 W_0 W_1 / (W_0 + W_1)^2$. Contrast information factors as $I_{C,H}^* = T_C B_C = 4 W_0 W_1 / (W_0 + W_1)$. Two designs can share $T_C$ and still differ in $I_{C,H}^*$ when the split of weight differs. \\\\",
  "    \\addlinespace",
  "    Independent-study\\newline support & How widely is the weight spread across studies, rather than concentrated in a few? & $J_g^* = W_g^2 / \\sum_{j:X_j=g} w_j^2$. Anticipated degrees of freedom use $V_g = 1/W_g$ and $\\nu_C^* = (V_0 + V_1)^2 / [V_0^2/(J_0^* - 1) + V_1^2/(J_1^* - 1)]$. Support is not a second copy of contrast information: designs matched on $I_{C,H}^*$ can still differ in $\\nu_C^*$. \\\\",
  "    \\addlinespace",
  paste0("    Heterogeneity\\newline sensitivity & How much sampling information remains after heterogeneity, and how far does predicted power move? & Retention is $R_H = I_{C,H}^* / I_{C,S}^*$. ", t2_def, " \\\\"),
  "    \\bottomrule",
  "  \\end{tabular*}%",
  "  \\vspace{0.4em}",
  "  {\\footnotesize\\raggedright",
  t2_note,
  "  \\par}",
  "\\end{table*}"
)
write_tex(tab_dir, "table2_information_profile_landscape.tex", t2)

# ---- Table 3 ----------------------------------------------------------------

s1_order <- c("naive_re", "mlma", "mlma_cr2")
if (!setequal(s1_ext$method, s1_order)) stop("Unexpected Study 1 estimators.", call. = FALSE)
s1_rows <- lapply(s1_order, function(m) {
  d <- s1_ext[s1_ext$method == m, , drop = FALSE]
  null <- d[abs(d$beta1_true) < 1e-12, , drop = FALSE]
  p10 <- d[abs(d$beta1_true - 0.10) < 1e-8, , drop = FALSE]
  p30 <- d[abs(d$beta1_true - 0.30) < 1e-8, , drop = FALSE]
  data.frame(
    method = method_label[[m]],
    n = nrow(d),
    n_null = nrow(null),
    max_abs_bias = max(abs(d$bias)),
    mean_se_ratio = mean_of(d$se_calibration_ratio),
    mean_type1 = mean_of(null$type_I_error),
    mean_coverage = mean_of(null$coverage),
    power10 = mean_of(p10$power),
    power30 = mean_of(p30$power),
    ci = mean_of(d$mean_ci_width),
    stringsAsFactors = FALSE
  )
})
s1_tab <- do.call(rbind, s1_rows)
if (length(unique(s1_tab$n)) != 1L || length(unique(s1_tab$n_null)) != 1L) {
  stop("Study 1 estimator rows do not share one condition count.", call. = FALSE)
}
if (abs(s1_tab$max_abs_bias[2L] - s1_tab$max_abs_bias[3L]) > 1e-12) {
  stop("MLMA and MLMA+CR2 no longer share the maximum absolute bias.", call. = FALSE)
}
t3_body <- vapply(seq_len(nrow(s1_tab)), function(i) {
  paste0(
    "    ", s1_tab$method[[i]], " & ",
    fmt_dec(s1_tab$max_abs_bias[[i]], 4L), " & ",
    fmt_dec(s1_tab$mean_se_ratio[[i]], 3L), " & ",
    fmt_dec(s1_tab$mean_type1[[i]], 4L), " & ",
    fmt_dec(s1_tab$mean_coverage[[i]], 3L), " & ",
    fmt_dec(s1_tab$power10[[i]], 3L), " & ",
    fmt_dec(s1_tab$power30[[i]], 3L), " & ",
    fmt_dec(s1_tab$ci[[i]], 3L), " \\\\"
  )
}, character(1))
t3 <- c(
  "\\begin{table*}",
  "  \\caption{Study 1 estimator summary.}",
  "  \\label{tab:study1-estimators}",
  "  \\small",
  "  \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrrrrr@{}}",
  "    \\toprule",
  "    Estimator & Max.\\ $|$bias$|$ & SE ratio & Type I & Coverage & Power, $\\beta_1=.10$ & Power, $\\beta_1=.30$ & CI width \\\\",
  "    \\midrule",
  t3_body,
  "    \\bottomrule",
  "  \\end{tabular*}%",
  "  \\vspace{0.4em}",
  "  {\\footnotesize\\raggedright",
  paste0(
    "  \\textit{Note.} Entries are unweighted means of the condition-level Study 1 summary. ",
    "Maximum absolute bias, the SE ratio, and mean confidence-interval width use all ", s1_tab$n[[1L]], " conditions. ",
    "Type I error and coverage use the ", s1_tab$n_null[[1L]], " null conditions. ",
    "Power uses the conditions at the stated $\\beta_1$. ",
    "MLMA and MLMA+CR2 share point estimates, so they share the maximum absolute bias. ",
    "The SE ratio is the mean ratio of the model standard error to the empirical standard deviation of the moderator estimate. ",
    "Type I error is shown to four decimals. ",
    "Independent RE treats effect sizes as independent. ",
    "MLMA is three-level REML with model-based Wald intervals. ",
    "MLMA+CR2 uses CR2 standard errors and Satterthwaite degrees of freedom. ",
    "SE, standard error; CI, confidence interval."
  ),
  "  \\par}",
  "\\end{table*}"
)
write_tex(tab_dir, "table3_study1_estimator_summary.tex", t3)

# ---- Table 4 ----------------------------------------------------------------

cal4 <- s4_pow[s4_pow$stratum == "broad_primary" & s4_pow$method == "mlma_cr2", , drop = FALSE]
regime_order <- c("sampling_only", "near_correct", "oracle")
regime_label <- c(sampling_only = "Sampling-only", near_correct = "Near-correct", oracle = "Oracle")
cal4 <- cal4[match(regime_order, cal4$regime), , drop = FALSE]
if (anyNA(cal4$regime)) stop("Broad-primary calibration regimes were not the expected three.", call. = FALSE)
t4_cal <- vapply(seq_len(nrow(cal4)), function(i) {
  row <- cal4[i, , drop = FALSE]
  cells <- c(
    fmt_math(row$mae, 3L), fmt_math(row$rmse, 3L), fmt_math(row$mean_signed_error, 3L),
    fmt_math(row$spearman, 3L), fmt_math(row$calibration_slope, 3L)
  )
  paste0("    ", regime_label[[row$regime]], " & ", fmt_int(row$n_cells), " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1))

comp_rho <- function(name) {
  hit <- s4_comp[s4_comp$comparator == name, , drop = FALSE]
  if (nrow(hit) != 1L) stop("Missing comparator row: ", name, call. = FALSE)
  hit
}
sup_rho <- function(name) {
  hit <- s4_sup[s4_sup$diagnostic == name & s4_sup$target == "mean_df", , drop = FALSE]
  if (nrow(hit) != 1L) stop("Missing support row: ", name, call. = FALSE)
  hit
}
var_row <- function(regime) {
  hit <- s4_var[s4_var$regime == regime, , drop = FALSE]
  if (nrow(hit) != 1L) stop("Missing variance regime: ", regime, call. = FALSE)
  hit
}
assoc <- function(diagnostic, comparison, n, estimate, digits) {
  paste0("    ", diagnostic, " & ", comparison, " & ", fmt_int(n), " & ", fmt_math(estimate, digits), " \\\\")
}
beta10 <- fmt_dec(s4_cfg$confirmatory_beta1, 2L)
pow_cmp <- paste0("Spearman with power at $\\beta_1=", beta10, "$")
class_letters <- c("A", "B", "C", "D")
class_ok <- vapply(class_letters, function(letter) {
  sub <- s4_pairs[s4_pairs$pair_class == letter, , drop = FALSE]
  nrow(sub) == 4L && all(sub$direction_supported %in% TRUE)
}, logical(1))
if (!all(class_ok)) {
  stop("Targeted classes A-D are not each 4/4 supported, so the collapsed confirmation row would be wrong.", call. = FALSE)
}
dir_phrase <- c(
  lower_R_H_has_larger_VIF_S = "lower retention has the larger sampling-only variance inflation factor",
  higher_B_C_has_lower_empirical_variance = "higher balance has the smaller empirical variance",
  greater_J_min_and_B_J_has_greater_mean_df = "greater independent-study support has the larger mean CR2 degrees of freedom",
  higher_I_C_H_has_lower_empirical_variance = "higher contrast information has the smaller empirical variance"
)
class_note <- paste(vapply(class_letters, function(letter) {
  code <- unique(s4_pairs$expected_direction[s4_pairs$pair_class == letter])
  if (length(code) != 1L || !code %in% names(dir_phrase)) {
    stop("Unexpected targeted-class direction for class ", letter, call. = FALSE)
  }
  paste0("Class ", letter, ", ", dir_phrase[[code]], " (4/4)")
}, character(1)), collapse = "; ")
if (any(!is.na(s4_frag_off$power_0.10_tau_0.05))) {
  stop("Official Class E empirical columns are no longer missing; the confirmation note would be wrong.", call. = FALSE)
}
if (nrow(s4_t1) != 1L || s4_t1$n_in_band != s4_t1$n_scenarios) {
  stop("Type I band count does not match the confirmation note.", call. = FALSE)
}
n_broad <- s4_cfg$n_broad
n_nonnull <- sum(abs(s4_cfg$beta1_levels) > 1e-12)
if (n_broad * n_nonnull != cal4$n_cells[[1L]]) {
  stop("Broad-primary cell count is not broad scenarios times non-null effects.", call. = FALSE)
}

t4 <- c(
  "\\begin{table*}",
  "  \\caption{Study 4 confirmation.}",
  "  \\label{tab:study4-confirmation}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{0.4em}",
  "  \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrrrr@{}}",
  "    \\toprule",
  "    \\multicolumn{7}{@{}l}{\\textit{A. Power calibration for MLMA+CR2 on the broad primary cells}} \\\\",
  "    \\midrule",
  "    Planning regime & Cells & MAE & RMSE & Signed error & Spearman & Slope \\\\",
  "    \\midrule",
  t4_cal,
  "    \\bottomrule",
  "  \\end{tabular*}%",
  "  \\vspace{0.8em}",
  "  \\setlength{\\tabcolsep}{0.45em}",
  "  \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}p{0.30\\linewidth}p{0.48\\linewidth}rr@{}}",
  "    \\toprule",
  "    \\multicolumn{4}{@{}l}{\\textit{B. Diagnostic confirmation}} \\\\",
  "    \\midrule",
  "    Diagnostic & Comparison & $n$ & Estimate \\\\",
  "    \\midrule",
  assoc("Retention, $R_H$", "Spearman with the log sampling-only variance inflation factor", s4_ret$n_broad_null, s4_ret$spearman_R_H_logVIF_S, 3L),
  "    \\addlinespace",
  assoc("Anticipated degrees of freedom, $\\nu_C^*$", "Spearman with mean CR2 degrees of freedom", sup_rho("nu_C_star")$n, sup_rho("nu_C_star")$spearman, 3L),
  assoc("Minimum group support, $J_{\\min}^*$", "Spearman with mean CR2 degrees of freedom", sup_rho("J_min_star")$n, sup_rho("J_min_star")$spearman, 3L),
  "    \\addlinespace",
  assoc("Near-correct contrast information, $I_{C,H}^*$", pow_cmp, comp_rho("I_C_H_star")$n, comp_rho("I_C_H_star")$spearman_vs_power_0.10, 3L),
  assoc("Sampling-only contrast information, $I_{C,S}^*$", pow_cmp, comp_rho("I_C_S_star")$n, comp_rho("I_C_S_star")$spearman_vs_power_0.10, 3L),
  assoc("Total sample size, $N$", pow_cmp, comp_rho("N_total")$n, comp_rho("N_total")$spearman_vs_power_0.10, 3L),
  assoc("Number of studies, $J$", pow_cmp, comp_rho("J")$n, comp_rho("J")$spearman_vs_power_0.10, 3L),
  assoc("Total number of effects, $K$", pow_cmp, comp_rho("K_total")$n, comp_rho("K_total")$spearman_vs_power_0.10, 3L),
  "    \\addlinespace",
  assoc("Sampling-only variance ratio", "Mean predicted variance divided by empirical variance at $\\beta_1=0$", var_row("sampling_only")$n_cells, var_row("sampling_only")$mean_predicted_empirical_ratio, 3L),
  assoc("Near-correct variance ratio", "Mean predicted variance divided by empirical variance at $\\beta_1=0$", var_row("near_correct")$n_cells, var_row("near_correct")$mean_predicted_empirical_ratio, 3L),
  "    \\addlinespace",
  "    Targeted classes A--D & Prespecified direction, 4 pairs in each class & 4 & 4/4 \\\\",
  "    \\addlinespace",
  assoc(paste0("Class E fragility, $\\Delta P_H(", beta10, ")$"), paste0("MAE versus empirical power loss at $\\beta_1=", beta10, "$"), s4_frag$n_families[[1L]], s4_frag$family_mae[[1L]], 3L),
  assoc(paste0("Class E fragility, $\\Delta P_H(", beta10, ")$"), paste0("RMSE versus empirical power loss at $\\beta_1=", beta10, "$"), s4_frag$n_families[[1L]], s4_frag$family_rmse[[1L]], 3L),
  assoc(paste0("Class E fragility, $\\Delta P_H(", beta10, ")$"), paste0("Spearman with empirical power loss at $\\beta_1=", beta10, "$"), s4_frag$n_families[[1L]], s4_frag$spearman_delta_P_H_vs_emp[[1L]], 2L),
  "    \\addlinespace",
  assoc("MLMA+CR2 Type I error", "Mean rejection rate at $\\beta_1=0$", s4_t1$n_scenarios, s4_t1$mean_type1, 4L),
  "    \\bottomrule",
  "  \\end{tabular*}%",
  "  \\vspace{0.4em}",
  "  {\\footnotesize\\raggedright",
  paste0(
    "  \\textit{Note.} Panel A reports the stored calibration of predicted MLMA+CR2 power. ",
    "Signed error is predicted minus empirical. ",
    "Those ", fmt_int(cal4$n_cells[[1L]]), " cells are the ", n_broad, " broad scenarios at each of the ", n_nonnull, " non-null moderator effects. ",
    "Near-correct planning rounds the true between-study and within-study heterogeneity to the nearest $", fmt_dec(het_step, 2L),
    "$ and uses planning correlation ", d2(sim4_rho_plan_primary), ". ",
    "Oracle planning uses the true heterogeneity and the true within-study correlation; it is a model-consistency benchmark. ",
    "Panel B reports official Study 4 summary rows, except the Class E fragility rows. ",
    "Those three rows are derived from the frozen scenario library and the MLMA+CR2 rows of condition\\_method\\_summary.csv. ",
    "Empirical power loss is power at $\\tau=.05$ minus power at $\\tau=.30$, the same direction as prespecified $\\Delta P_H(", beta10, ")$. ",
    "The official family summary still leaves the empirical columns missing. ",
    "The retention correlation is the stored Spearman association with the log variance inflation factor, on the broad null scenarios. ",
    "A secondary Spearman association between $R_H$ and sampling-only power optimism at $\\beta_1=", beta10, "$ was ",
    fmt_math(s4_ret$spearman_R_H_sampling_power_optimism, 3L), " on the same ", fmt_int(s4_ret$n_broad_beta10), " scenarios. ",
    "$\\nu_C^*$ and $J_{\\min}^*$ track realized CR2 degrees of freedom; they are not those degrees of freedom. ",
    "The support correlations use all ", fmt_int(sup_rho("nu_C_star")$n), " null scenarios. ",
    "The power-ranking correlations use the ", fmt_int(comp_rho("I_C_H_star")$n), " broad scenarios. ",
    class_note, ". ",
    "The mean Type I error is shown to four decimals. ",
    "All ", fmt_int(s4_t1$n_scenarios), " scenarios fell inside the prespecified band ",
    fmt_dec(s4_t1$band_low, 3L), " to ", fmt_dec(s4_t1$band_high, 3L),
    " (observed range ", fmt_dec(s4_t1$min_type1, 4L), " to ", fmt_dec(s4_t1$max_type1, 4L), "). ",
    "MAE, mean absolute error; RMSE, root mean squared error."
  ),
  "  \\par}",
  "\\end{table*}"
)
write_tex(tab_dir, "table4_study4_confirmation_revised.tex", t4)

# ---- Supplement S1 ----------------------------------------------------------

range_qs <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  stats::quantile(x, c(0, 0.5, 1), names = FALSE, type = 7)
}
fmt_count_q <- function(x) {
  if (abs(x - round(x)) > 1e-6) stop("Expected an integer design-range endpoint, got ", x, call. = FALSE)
  fmt_int(x)
}
fmt_cont_q <- function(x) fmt_dec(x, 3L)

s3_ss2 <- merge(s3_ss, s3_scen[, c("scenario_id", "scenario_set")], by = "scenario_id")
s4_ss2 <- merge(s4_ss, s4_lib[, c("scenario_id", "scenario_set")], by = "scenario_id")
if (nrow(s3_ss2) != nrow(s3_ss) || nrow(s4_ss2) != nrow(s4_ss)) {
  stop("Study structure did not join one-to-one onto scenarios.", call. = FALSE)
}
subset_block <- function(study_label, scen, ss) {
  sets <- c(Broad = "broad", Targeted = "targeted", Combined = "all")
  qty_scen <- c(
    "Number of studies, $J$" = "J",
    "Sampling correlation, $\\rho$" = "rho_true",
    "Between-study heterogeneity, $\\tau$" = "tau",
    "Within-study heterogeneity, $\\omega$" = "omega"
  )
  lines <- character()
  for (i in seq_along(sets)) {
    key <- sets[[i]]
    lab <- names(sets)[[i]]
    sub_s <- if (key == "all") scen else scen[scen$scenario_set == key, , drop = FALSE]
    sub_u <- if (key == "all") ss else ss[ss$scenario_set == key, , drop = FALSE]
    specs <- list(
      list("Number of studies, $J$", range_qs(sub_s$J), "count"),
      list("Study sample size, $n_j$", range_qs(sub_u$n_j), "count"),
      list("Effects per study, $k_j$", range_qs(sub_u$k_j), "count"),
      list("Sampling correlation, $\\rho$", range_qs(sub_s$rho_true), "cont"),
      list("Between-study heterogeneity, $\\tau$", range_qs(sub_s$tau), "cont"),
      list("Within-study heterogeneity, $\\omega$", range_qs(sub_s$omega), "cont")
    )
    for (r in seq_along(specs)) {
      qs <- specs[[r]][[2L]]
      fmt <- if (specs[[r]][[3L]] == "count") fmt_count_q else fmt_cont_q
      end <- if (r == length(specs)) " \\\\" else " \\\\*"
      lead <- if (r == 1L) lab else ""
      lines <- c(lines, paste0(
        if (r == 1L) "    " else "     ",
        lead, " & ", specs[[r]][[1L]], " & ",
        fmt(qs[[1L]]), " & ", fmt(qs[[2L]]), " & ", fmt(qs[[3L]]), end
      ))
    }
    if (i < length(sets)) lines <- c(lines, "    \\addlinespace")
  }
  lines
}
s1_lines <- c(
  "\\begin{longtable}{@{}llrrr@{}}",
  "  \\caption{Realized Study 3 and Study 4 scenario-library ranges.}",
  "  \\label{tab:s1}\\\\",
  "  \\toprule",
  "  Library & Quantity & Minimum & Median & Maximum \\\\",
  "  \\midrule",
  "  \\endfirsthead",
  "  \\caption[]{Realized Study 3 and Study 4 scenario-library ranges (continued).}\\\\",
  "  \\toprule",
  "  Library & Quantity & Minimum & Median & Maximum \\\\",
  "  \\midrule",
  "  \\endhead",
  "  \\midrule",
  "  \\multicolumn{5}{r@{}}{\\footnotesize Continued on the next page}\\\\",
  "  \\endfoot",
  "  \\bottomrule",
  "  \\endlastfoot",
  "  \\multicolumn{5}{@{}l}{\\textit{A. Study 3}}\\\\",
  "  \\midrule",
  subset_block("3", s3_scen, s3_ss2),
  "  \\midrule",
  "  \\multicolumn{5}{@{}l}{\\textit{B. Study 4}}\\\\",
  "  \\midrule",
  subset_block("4", s4_lib, s4_ss2),
  "\\end{longtable}",
  "\\vspace{-0.6em}",
  "{\\footnotesize\\raggedright",
  paste0(
    "\\textit{Note.} Entries are the exact observed minimum, median, and maximum in the realized scenario libraries. ",
    "$J$, $\\rho$, $\\tau$, and $\\omega$ are scenario-level. Study sample size $n_j$ and effects per study $k_j$ are the corresponding values across study rows, not scenario means. ",
    "The combined library is the broad library plus the targeted library. These observed ranges are distinct from the scenario-generation windows. ",
    "Study 3 scenarios were drawn with $J\\in\\{", j_tex(s3_cfg$J_levels), "\\}$, $\\rho$ approximately ",
    fmt_dec(s3_cfg$rho_true_range[[1L]], 2L), " to ", fmt_dec(s3_cfg$rho_true_range[[2L]], 2L),
    ", and $\\tau$ and $\\omega$ approximately ", fmt_dec(s3_cfg$tau_range[[1L]], 2L), " to ", fmt_dec(s3_cfg$tau_range[[2L]], 2L), ". ",
    "Study 4 scenarios were drawn with $J\\in\\{", j_tex(s4_cfg$J_levels), "\\}$, $\\rho$ approximately ",
    fmt_dec(s4_cfg$rho_true_range[[1L]], 2L), " to ", fmt_dec(s4_cfg$rho_true_range[[2L]], 2L),
    ", and $\\tau$ and $\\omega$ approximately ", fmt_dec(s4_cfg$tau_range[[1L]], 2L), " to ", fmt_dec(s4_cfg$tau_range[[2L]], 2L), ". ",
    "A printed endpoint such as .300 can fall on a generation bound. In Study 3 the targeted library lies inside the observed broad endpoints of $\\rho$, $\\tau$, and $\\omega$. In Study 4 the targeted library reaches a higher observed $\\tau$ and a lower observed $\\omega$ than the broad library. Continuous entries are shown to three decimals."
  ),
  "\\par}"
)
write_tex(stab_dir, "tableS1_design_ranges.tex", s1_lines)

# ---- Supplement S2 ----------------------------------------------------------

kj_phrase <- function(J, imbalance) {
  kj <- make_kj_template(J, imbalance)
  if (imbalance == "balanced") return("Every study has $k_j=4$")
  tab <- as.data.frame(table(kj), stringsAsFactors = FALSE)
  tab$kj <- as.integer(as.character(tab$kj))
  tab <- tab[order(tab$kj), , drop = FALSE]
  if (imbalance == "moderate") {
    if (length(unique(tab$Freq)) != 1L) stop("Moderate template is not balanced across k.", call. = FALSE)
    return(sprintf("%d studies at each of $k_j=%s$", tab$Freq[[1L]], paste(tab$kj, collapse = ", ")))
  }
  sprintf("%s studies with $k_j=%s$", paste(tab$Freq, collapse = ", "), paste(tab$kj, collapse = ", "))
}
s2_rows <- character()
imb_lab <- c(balanced = "Balanced", moderate = "Moderate", severe = "Severe")
for (J in c(20L, 40L, 80L)) {
  imbs <- c("balanced", "moderate", "severe")
  for (i in seq_along(imbs)) {
    lead <- if (i == 1L) as.character(J) else ""
    s2_rows <- c(s2_rows, paste0("    ", lead, " & ", imb_lab[[imbs[[i]]]], " & ", kj_phrase(J, imbs[[i]]), " \\\\"))
  }
  if (J != 80L) s2_rows <- c(s2_rows, "    \\addlinespace")
}
write_tex(stab_dir, "tableS2_study1_effect_counts.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Study 1 effect-count templates.}",
  "  \\label{tab:s2}",
  "  \\small",
  "  ",
  "\\begin{tabular}{@{}llp{0.62\\linewidth}@{}}",
  "    \\toprule",
  "    $J$ & Imbalance & Effect-count pattern \\\\",
  "    \\midrule",
  s2_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} In every condition the total number of effect sizes is $K=4J$, so the mean number of effects per study is four. The severe templates place many effects in a few studies. Templates were randomly permuted across study labels within each replication. $J$, number of studies; $k_j$, number of effects contributed by a study.",
  "  \\par}",
  "\\end{table}"
))

# ---- Supplement S3 ----------------------------------------------------------

s2_n$J <- as.integer(s2_n$J)
n_order <- c("balanced", "moderate", "severe", "extreme")
n_lab <- c(balanced = "Balanced", moderate = "Moderate", severe = "Severe", extreme = "Extreme")
s3n_rows <- character()
for (J in c(20L, 40L, 80L)) {
  for (i in seq_along(n_order)) {
    hit <- s2_n[s2_n$J == J & s2_n$n_imbalance == n_order[[i]], , drop = FALSE]
    if (nrow(hit) != 1L) stop("Missing Study 2 sample-size template.", call. = FALSE)
    lead <- if (i == 1L) as.character(J) else ""
    cells <- c(
      fmt_trim(hit$sigma_log_n, 1L), fmt_trim(hit$n_min, 2L), fmt_trim(hit$n_q1, 2L),
      fmt_trim(hit$n_median, 2L), fmt_trim(hit$n_mean, 2L), fmt_trim(hit$n_q3, 2L),
      fmt_trim(hit$n_max, 2L), fmt_trim(hit$n_total, 0L), fmt_trim(hit$cv_sample, 3L)
    )
    s3n_rows <- c(s3n_rows, paste0("    ", lead, " & ", n_lab[[n_order[[i]]]], " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  if (J != 80L) s3n_rows <- c(s3n_rows, "    \\addlinespace")
}
write_tex(stab_dir, "tableS3_study2_sample_sizes.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Study 2 sample-size template characteristics.}",
  "  \\label{tab:s3}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.35em}",
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}llrrrrrrrrr@{}}",
  "    \\toprule",
  "    $J$ & Imbalance & $\\sigma_{\\log n}$ & Min. & Q1 & Median & Mean & Q3 & Max. & Total $N$ & CV \\\\",
  "    \\midrule",
  s3n_rows,
  "    \\bottomrule",
  "  \\end{tabular*}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  paste0(
    "  \\textit{Note.} Templates are deterministic and bounded between ", as.integer(s2_cfg$n_min),
    " and ", fmt_int(s2_cfg$n_max), " participants. The mean study size is ", as.integer(s2_cfg$mean_n),
    " at every $J$, and the total number of participants is $N=", as.integer(s2_cfg$mean_n), "J$. ",
    "Increasing imbalance redistributes that fixed total. $\\sigma_{\\log n}$ is the log-scale dispersion used to build the template. ",
    "CV is the coefficient of variation of the template study sizes. At $J=20$ the extreme template reaches a maximum of ",
    fmt_trim(s2_n$n_max[s2_n$J == 20 & s2_n$n_imbalance == "extreme"], 0L),
    " rather than the upper bound of ", fmt_int(s2_cfg$n_max), ". $J$, number of studies."
  ),
  "  \\par}",
  "\\end{table}"
))

# ---- Supplement S4 ----------------------------------------------------------

se_rows <- lapply(c("balanced", "moderate", "severe"), function(imb) {
  d <- s1_ext[s1_ext$method == "naive_re" & s1_ext$imbalance == imb, , drop = FALSE]
  data.frame(imbalance = c(balanced = "Balanced", moderate = "Moderate", severe = "Severe")[[imb]], n = nrow(d), se = mean_of(d$se_calibration_ratio))
})
se_tab <- do.call(rbind, se_rows)
mc0 <- s1_mc[abs(s1_mc$beta1_true) < 1e-12, , drop = FALSE]
mc_rows <- lapply(s1_order, function(m) {
  d <- mc0[mc0$method == m, , drop = FALSE]
  inside <- d$rejection_mc_lower <= 0.05 & d$rejection_mc_upper >= 0.05
  data.frame(method = method_label[[m]], n = nrow(d), n_in = sum(inside))
})
mc_tab <- do.call(rbind, mc_rows)
write_tex(stab_dir, "tableS4_study1_operating.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Additional Study 1 operating characteristics.}",
  "  \\label{tab:s4}",
  "  \\small",
  "  ",
  "\\begin{tabular}{@{}lrr@{}}",
  "    \\toprule",
  "    \\multicolumn{3}{@{}l}{\\textit{A. Independent-RE standard-error calibration by effect-count imbalance}}\\\\",
  "    \\midrule",
  "    Imbalance & Conditions & Mean SE ratio \\\\",
  "    \\midrule",
  vapply(seq_len(nrow(se_tab)), function(i) paste0("    ", se_tab$imbalance[[i]], " & ", se_tab$n[[i]], " & ", fmt_dec(se_tab$se[[i]], 3L), " \\\\"), character(1)),
  "    \\midrule",
  "    \\multicolumn{3}{@{}l}{\\textit{B. Number of null conditions whose Monte Carlo interval contained .05}}\\\\",
  "    \\midrule",
  "    Estimator & Null conditions & Interval contained .05 \\\\",
  "    \\midrule",
  vapply(seq_len(nrow(mc_tab)), function(i) paste0("    ", mc_tab$method[[i]], " & ", mc_tab$n[[i]], " & ", mc_tab$n_in[[i]], " \\\\"), character(1)),
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  paste0(
    "  \\textit{Note.} In Panel A the SE ratio is the mean ratio of the model standard error to the empirical standard deviation of the moderator estimate. ",
    "Each imbalance level averages ", se_tab$n[[1L]], " conditions, crossing $J\\in\\{", j_tex(s1_J, ","),
    "\\}$, the three moderator allocations, and $\\beta_1=0,.10,.30$. ",
    "Independent RE is random-effects meta-regression that treats effect sizes as independent. ",
    "In Panel B each count is out of the ", mc_tab$n[[1L]], " null conditions. ",
    "The interval is a 95\\% normal-approximation Monte Carlo interval, calculated as the rejection rate $\\pm 1.96$ times its Monte Carlo standard error. ",
    "MLMA is three-level REML with model-based Wald intervals. MLMA+CR2 uses CR2 standard errors and Satterthwaite degrees of freedom. SE, standard error."
  ),
  "  \\par}",
  "\\end{table}"
))

# ---- Supplement S5 ----------------------------------------------------------

core <- s2_cm[s2_cm$condition_set == "core", , drop = FALSE]
split_lab <- c(`50_50` = "50/50", `70_30` = "70/30", `85_15` = "85/15")
k_labs <- c(balanced = "Balanced", moderate = "Moderate", severe = "Severe")
one_cell <- function(k, sp) {
  d <- core[core$method == "mlma_cr2" & abs(core$beta1_true - 0.10) < 1e-8 & core$J == 40 &
              core$n_imbalance == "balanced" & core$allocation_pattern == "independent" &
              core$k_imbalance == k & core$moderator_split == sp, , drop = FALSE]
  if (nrow(d) != 1L) stop("Study 2 Panel A cell is not unique.", call. = FALSE)
  d
}
panel_a_power <- character()
panel_a_df <- character()
for (k in names(k_labs)) {
  vals_p <- vapply(names(split_lab), function(sp) fmt_dec(one_cell(k, sp)$power, 3L), character(1))
  vals_d <- vapply(names(split_lab), function(sp) fmt_fixed(one_cell(k, sp)$mean_df, 2L), character(1))
  panel_a_power <- c(panel_a_power, paste0("    ", k_labs[[k]], " & ", paste(vals_p, collapse = " & "), " \\\\"))
  panel_a_df <- c(panel_a_df, paste0("    ", k_labs[[k]], " & ", paste(vals_d, collapse = " & "), " \\\\"))
}
panel_b <- lapply(n_order, function(ni) {
  out <- lapply(c("mlma_cr2", "che_rve"), function(m) {
    d <- core[core$method == m & core$n_imbalance == ni, , drop = FALSE]
    null <- d[abs(d$beta1_true) < 1e-12, , drop = FALSE]
    p10 <- d[abs(d$beta1_true - 0.10) < 1e-8, , drop = FALSE]
    if (nrow(null) != 27L || nrow(p10) != 27L) stop("Panel B cell count is not 27.", call. = FALSE)
    c(fmt_dec(mean_of(null$type1_error), 4L), fmt_dec(mean_of(p10$power), 3L))
  })
  paste0("    ", n_lab[[ni]], " & ", paste(unlist(out), collapse = " & "), " \\\\")
})
panel_c <- lapply(c("mlma_cr2", "che_rve"), function(m) {
  d <- core[core$method == m, , drop = FALSE]
  vals <- vapply(names(split_lab), function(sp) fmt_dec(mean_of(d$mean_df[d$moderator_split == sp]), 2L), character(1))
  paste0("    ", method_label[[m]], " & ", paste(vals, collapse = " & "), " \\\\")
})
dem <- core[core$J == 20 & core$moderator_split == "85_15" & core$k_imbalance == "severe" &
              core$n_imbalance == "extreme" & core$method %in% c("mlma_cr2", "che_rve"), , drop = FALSE]
panel_d <- character()
for (m in c("mlma_cr2", "che_rve")) {
  for (b in c(0, 0.10)) {
    d <- dem[dem$method == m & abs(dem$beta1_true - b) < 1e-8, , drop = FALSE]
    if (nrow(d) != 1L) stop("Panel D cell is not unique.", call. = FALSE)
    rate <- if (abs(b) < 1e-12) d$type1_error else d$power
    panel_d <- c(panel_d, paste0(
      "    ", method_label[[m]], " & ", if (abs(b) < 1e-12) "$0$" else paste0("$", fmt_dec(b, 2L), "$"),
      " & ", fmt_dec(d$mean_df, 2L), " & ", fmt_dec(d$prop_df_lt_4, 3L), " & ", fmt_dec(rate, 3L), " \\\\"
    ))
  }
}
write_tex(stab_dir, "tableS5_study2_operating.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Additional Study 2 operating characteristics.}",
  "  \\label{tab:s5}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{0.45em}",
  "\\begin{tabular}{@{}lccc@{}}",
  "    \\toprule",
  "    \\multicolumn{4}{@{}l}{\\textit{A. Power across effect-count imbalance}}\\\\",
  "    \\midrule",
  "    \\multicolumn{4}{@{}l}{MLMA+CR2 power at $\\beta_1=.10$}\\\\",
  "    Effect-count imbalance & 50/50 & 70/30 & 85/15 \\\\",
  "    \\midrule",
  panel_a_power,
  "    \\addlinespace",
  "    \\multicolumn{4}{@{}l}{Mean degrees of freedom}\\\\",
  "    Effect-count imbalance & 50/50 & 70/30 & 85/15 \\\\",
  "    \\midrule",
  panel_a_df,
  "    \\bottomrule",
  "  \\end{tabular}\\par\\bigskip",
  "  \\begin{tabular}{@{}lrrrr@{}}",
  "    \\toprule",
  "    \\multicolumn{5}{@{}l}{\\textit{B. Type I error and power across sample-size imbalance}}\\\\",
  "    \\midrule",
  "    Sample-size imbalance & \\multicolumn{2}{c}{MLMA+CR2} & \\multicolumn{2}{c}{CHE-RVE} \\\\",
  "    \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "    & Type I & Power & Type I & Power \\\\",
  "    \\midrule",
  unlist(panel_b),
  "    \\bottomrule",
  "  \\end{tabular}\\par\\bigskip",
  "  {\\textit{C. Mean CR2/Satterthwaite degrees of freedom across moderator allocations}\\par}",
  "  \\vspace{0.35em}",
  "  \\begin{tabular}{@{}lccc@{}}",
  "    \\toprule",
  "    Estimator & 50/50 & 70/30 & 85/15 \\\\",
  "    \\midrule",
  unlist(panel_c),
  "    \\bottomrule",
  "  \\end{tabular}\\par\\bigskip",
  "  \\begin{tabular}{@{}lrrrr@{}}",
  "    \\toprule",
  "    \\multicolumn{5}{@{}l}{\\textit{D. Most demanding $J=20$ condition}}\\\\",
  "    \\midrule",
  "    Estimator & $\\beta_1$ & Mean df & Prop.\\ df $<4$ & Rejection rate \\\\",
  "    \\midrule",
  panel_d,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Panel A holds $J=40$, balanced study sizes, and independent assignment of study size, effect multiplicity, and moderator membership. Every cell has $K=160$ and mean total $N=12{,}000$. Panel B averages the core factorial within each sample-size imbalance: Type I error over the 27 null cells, and power over the 27 cells at $\\beta_1=.10$. Panel C averages degrees of freedom over both values of $\\beta_1$ in the core factorial; the null-only means differ only in the third decimal. Panel D is the single core cell with $J=20$, an 85/15 allocation, severe effect-count imbalance, and extreme sample-size imbalance. Its rejection rate is the Type I error at $\\beta_1=0$ and power at $\\beta_1=.10$. MLMA+CR2 uses CR2 standard errors and Satterthwaite degrees of freedom. CHE-RVE is the correlated-and-hierarchical-effects robust-variance estimator. df, degrees of freedom.",
  "  \\par}",
  "\\end{table}"
))

# ---- Supplement S6 ----------------------------------------------------------

if (!all(s2_align$k_imbalance == "severe" & s2_align$moderator_split == "85_15" & s2_align$n_imbalance == "extreme")) {
  stop("Alignment summary is not the severe, 85/15, extreme stress design.", call. = FALSE)
}
alloc_lab <- c(
  independent = "Independent",
  n_majority = "Large $n$ in majority",
  n_minority = "Large $n$ in minority",
  nk_majority = "Large $n$ and $k$ in majority",
  nk_minority = "Large $n$ and $k$ in minority"
)
s2_methods <- c("aggregate_hksj", "mlma", "mlma_cr2", "che_rve")
fmt_df_cell <- function(method, x) {
  if (method == "mlma") {
    if (!is.na(x)) stop("MLMA alignment row has a degrees-of-freedom value.", call. = FALSE)
    return("---")
  }
  if (method == "aggregate_hksj") {
    if (abs(x - round(x)) > 1e-8) stop("Aggregation df is not an integer.", call. = FALSE)
    return(as.character(as.integer(round(x))))
  }
  fmt_dec(x, 2L)
}
s6_rows <- character()
for (J in c(20, 40)) {
  for (a in names(alloc_lab)) {
    first_a <- TRUE
    for (m in s2_methods) {
      d <- s2_align[s2_align$J == J & s2_align$allocation_pattern == a & s2_align$method == m & abs(s2_align$beta1_true) < 1e-12, , drop = FALSE]
      if (nrow(d) != 1L) stop("Alignment null row is not unique.", call. = FALSE)
      jcell <- if (a == names(alloc_lab)[[1L]] && first_a) as.character(J) else ""
      acell <- if (first_a) alloc_lab[[a]] else ""
      first_a <- FALSE
      s6_rows <- c(s6_rows, paste0(
        "    ", jcell, " & ", acell, " & ", method_label[[m]], " & ",
        fmt_dec(d$type1_error, 3L), " & ", fmt_dec(d$coverage, 3L), " & ",
        fmt_dec(d$se_ratio, 3L), " & ", fmt_df_cell(m, d$mean_df), " \\\\"
      ))
    }
  }
  if (J == 20) s6_rows <- c(s6_rows, "    \\addlinespace")
}
write_tex(stab_dir, "tableS6_alignment.tex", c(
  "{\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{longtable}{@{}lllcrrr@{}}",
  "  \\caption{Alignment-stress inferential performance under the null.}",
  "  \\label{tab:s6}\\\\",
  "  \\toprule",
  "  $J$ & Allocation & Estimator & Type I & Coverage & SE ratio & Mean df \\\\",
  "  \\midrule",
  "  \\endfirsthead",
  "  \\caption[]{Alignment-stress inferential performance under the null (continued).}\\\\",
  "  \\toprule",
  "  $J$ & Allocation & Estimator & Type I & Coverage & SE ratio & Mean df \\\\",
  "  \\midrule",
  "  \\endhead",
  "  \\bottomrule",
  "  \\endfoot",
  s6_rows,
  "\\end{longtable}",
  "}",
  "\\vspace{-0.2em}",
  "{\\footnotesize\\raggedright",
  "\\textit{Note.} All rows are null alignment-stress cells with severe effect-count imbalance, an 85/15 allocation, and extreme sample-size imbalance. Allocation changes where the same study sizes and effect counts fall across moderator groups. Aggregation + KH aggregates dependent effects and applies the Knapp--Hartung adjustment. MLMA has no Satterthwaite df. MLMA+CR2 uses CR2 standard errors and Satterthwaite df. CHE-RVE is correlated-and-hierarchical-effects robust variance estimation. The SE ratio is the model standard error divided by the empirical standard deviation of the moderator estimate. df, degrees of freedom.",
  "\\par}"
))

# ---- Supplement S7-S12 ------------------------------------------------------

if (nrow(s3_var) != 1L || nrow(s3_df) != 1L) stop("Study 3 variance or df summary is not one row.", call. = FALSE)
write_tex(stab_dir, "tableS7_study3_variance_df.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Prediction of variance and finite-sample degrees of freedom in Study 3.}",
  "  \\label{tab:s7}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.4em}",
  "\\begin{tabular}{@{}lrrrrr@{}}",
  "    \\toprule",
  "    Target & $n$ & RMSE & MAE & Signed error & Spearman \\\\",
  "    \\midrule",
  paste0("    Moderator-estimate variance & ", s3_var$n, " & ", fmt_dec(s3_var$rmse, 5L), " & ", fmt_dec(s3_var$mae, 5L), " & ", fmt_dec(s3_var$mean_signed_error, 5L), " & ", fmt_dec(s3_var$spearman, 3L), " \\\\"),
  paste0("    Mean CR2/Satterthwaite df & ", s3_df$n, " & ", fmt_dec(s3_df$rmse, 2L), " & ", fmt_dec(s3_df$mae, 2L), " & ", fmt_dec(s3_df$mean_signed_error, 2L), " & ", fmt_dec(s3_df$spearman, 3L), " \\\\"),
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Both rows use the sampling-only contrast information at the planning correlation $\\rho_p=.50$, evaluated on the 100 Study 3 scenarios under the null. The variance target is the empirical variance of the MLMA moderator estimate. The degrees-of-freedom target is the mean CR2/Satterthwaite degrees of freedom. Signed error is predicted minus empirical, so a negative value means the diagnostic was too small. RMSE, root mean squared error; MAE, mean absolute error; df, degrees of freedom.",
  "  \\par}",
  "\\end{table}"
))

variant_order <- c(
  "sampling_rho_0.50",
  "het_tau0.10_omega0.10_rho0.50",
  "sampling_rho_oracle",
  "het_oracle",
  "het_tau0.05_omega0.05_rho0.20",
  "het_tau0.10_omega0.10_rho0.20",
  "het_tau0.20_omega0.20_rho0.20",
  "het_tau0.05_omega0.05_rho0.50",
  "het_tau0.20_omega0.20_rho0.50",
  "het_tau0.05_omega0.05_rho0.80",
  "het_tau0.10_omega0.10_rho0.80",
  "het_tau0.20_omega0.20_rho0.80"
)
variant_label <- c(
  sampling_rho_0.50 = "Sampling only ($\\rho_p=.50$)",
  het_tau0.10_omega0.10_rho0.50 = "Default ($\\tau_p=\\omega_p=.10$, $\\rho_p=.50$)",
  sampling_rho_oracle = "Oracle sampling ($\\rho$ known)",
  het_oracle = "Oracle heterogeneity ($\\tau$, $\\omega$, and $\\rho$ known)",
  het_tau0.05_omega0.05_rho0.20 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.05$, $\\rho_p=.20$)",
  het_tau0.10_omega0.10_rho0.20 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.10$, $\\rho_p=.20$)",
  het_tau0.20_omega0.20_rho0.20 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.20$, $\\rho_p=.20$)",
  het_tau0.05_omega0.05_rho0.50 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.05$, $\\rho_p=.50$)",
  het_tau0.20_omega0.20_rho0.50 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.20$, $\\rho_p=.50$)",
  het_tau0.05_omega0.05_rho0.80 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.05$, $\\rho_p=.80$)",
  het_tau0.10_omega0.10_rho0.80 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.10$, $\\rho_p=.80$)",
  het_tau0.20_omega0.20_rho0.80 = "Heterogeneity-adjusted ($\\tau_p=\\omega_p=.20$, $\\rho_p=.80$)"
)
pow_nn <- s3_pow[s3_pow$beta1_group == "all_non_null", , drop = FALSE]
if (!setequal(pow_nn$variant_id, variant_order)) stop("Study 3 power variants do not match the supplement table.", call. = FALSE)
s8_rows <- vapply(variant_order, function(id) {
  row <- pow_nn[pow_nn$variant_id == id, , drop = FALSE]
  paste0(
    "    ", variant_label[[id]], " & ", row$n_rows, " & ",
    fmt_dec(row$rmse_power, 3L), " & ", fmt_dec(row$mae_power, 3L), " & ",
    fmt_signed(row$mean_signed_error, 3L), " & ", fmt_dec(row$spearman, 3L), " & ",
    fmt_dec(row$calibration_slope, 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS8_study3_power_variants.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Study 3 power-prediction variants.}",
  "  \\label{tab:s8}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.35em}",
  "\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{0.44\\linewidth}rrrrrr@{}}",
  "    \\toprule",
  "    Specification & $n$ & RMSE & MAE & Signed error & Spearman & Slope \\\\",
  "    \\midrule",
  s8_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Each row predicts MLMA+CR2 power for the 300 non-null Study 3 conditions (100 scenarios at $\\beta_1=.05,.10,.20$). Sampling-only specifications use sampling precision and do not include heterogeneity. Heterogeneity-adjusted specifications use the stated planning values, with $\\tau_p$ tied to $\\omega_p$. The default row is the common planning specification $\\tau_p=\\omega_p=.10$ and $\\rho_p=.50$. Oracle sampling uses the scenario's true $\\rho$ and still omits heterogeneity. Oracle heterogeneity uses the true $\\tau$, $\\omega$, and $\\rho$. Signed error is predicted power minus empirical power. The slope is the calibration slope from regressing empirical power on predicted power. RMSE, root mean squared error; MAE, mean absolute error.",
  "  \\par}",
  "\\end{table}"
))

pow10 <- s3_cm[s3_cm$method == "mlma_cr2" & abs(s3_cm$beta1_true - 0.10) < 1e-8, c("scenario_id", "power")]
s9_rows <- vapply(seq_len(nrow(s3_pairs)), function(i) {
  pa <- pow10$power[match(s3_pairs$scenario_id_a[[i]], pow10$scenario_id)]
  pb <- pow10$power[match(s3_pairs$scenario_id_b[[i]], pow10$scenario_id)]
  if (!is.finite(pa) || !is.finite(pb)) stop("Study 3 pair power join failed.", call. = FALSE)
  paste0(
    "    ", s3_pairs$pair_class[[i]], " & ", s3_pairs$pair_id[[i]], " & ",
    s3_pairs$scenario_id_a[[i]], "/", s3_pairs$scenario_id_b[[i]], " & ",
    fmt_ic(s3_pairs$I_C_star_a[[i]]), " & ", fmt_ic(s3_pairs$I_C_star_b[[i]]), " & ",
    fmt_dec(s3_pairs$mean_df_cr2_a[[i]], 2L), " & ", fmt_dec(s3_pairs$mean_df_cr2_b[[i]], 2L), " & ",
    fmt_dec(pa, 3L), " & ", fmt_dec(pb, 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS9_study3_targeted_pairs.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Study 3 targeted matched-design comparisons.}",
  "  \\label{tab:s9}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.35em}",
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}ccr rr rr rr@{}}",
  "    \\toprule",
  "    & & & \\multicolumn{2}{c}{$I_{C,S}^{*}$} & \\multicolumn{2}{c}{Mean CR2 df} & \\multicolumn{2}{c}{Power} \\\\",
  "    \\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "    Class & Pair & Scenarios & Arm A & Arm B & Arm A & Arm B & Arm A & Arm B \\\\",
  "    \\midrule",
  s9_rows,
  "    \\bottomrule",
  "  \\end{tabular*}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Classes A and B match $J$, $N$, $K$, and the minority-group count while $I_{C,S}^{*}$ differs; Class B also separates the arms on the minority weight share. Class C matches $I_{C,S}^{*}$ and varies independent-study support. Class D matches anticipated degrees of freedom and varies $I_{C,S}^{*}$. Class E matches the study counts and their concentration and varies $I_{C,S}^{*}$. Mean CR2 df is the realized CR2/Satterthwaite degrees of freedom. Power is the MLMA+CR2 rejection rate at $\\beta_1=.10$. Those power columns were taken from the MLMA+CR2 condition summaries because the targeted-pair summary does not store power. $I_{C,S}^{*}$, sampling-only contrast information; df, degrees of freedom.",
  "  \\par}",
  "\\end{table}"
))

bin_order <- c("[0,0.025)", "[0.025,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,Inf)")
bin_label <- c(
  "[0,0.025)" = "$|\\tau_p-\\tau|<.025$",
  "[0.025,0.05)" = "$.025\\le|\\tau_p-\\tau|<.05$",
  "[0.05,0.10)" = "$.05\\le|\\tau_p-\\tau|<.10$",
  "[0.10,0.15)" = "$.10\\le|\\tau_p-\\tau|<.15$",
  "[0.15,Inf)" = "$|\\tau_p-\\tau|\\ge.15$"
)
miss <- s3_miss[s3_miss$analysis_family == "abs_tau_misspec" & abs(s3_miss$rho_plan - 0.50) < 1e-8 & s3_miss$scenario_set == "all", , drop = FALSE]
if (!setequal(miss$misspec_bin, bin_order)) stop("Tau-misspecification bins do not match the supplement table.", call. = FALSE)
s10_rows <- vapply(bin_order, function(bin) {
  row <- miss[miss$misspec_bin == bin, , drop = FALSE]
  paste0(
    "    ", bin_label[[bin]], " & ", fmt_int(row$n_rows), " & ",
    fmt_dec(row$mae_power, 3L), " & ", fmt_dec(row$rmse_power, 3L), " & ",
    fmt_signed(row$mean_signed_error, 3L), " & ", fmt_dec(row$spearman, 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS10_tau_misspecification.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Error by between-study heterogeneity misspecification.}",
  "  \\label{tab:s10}",
  "  \\small",
  "  ",
  "\\begin{tabular}{@{}lrrrrr@{}}",
  "    \\toprule",
  "    $|\\tau_p-\\tau|$ & $n$ & MAE & RMSE & Signed error & Spearman \\\\",
  "    \\midrule",
  s10_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Rows are the absolute discrepancy between the planning value $\\tau_p$ and the true between-study heterogeneity, at the planning correlation $\\rho_p=.50$. The outcome is MLMA+CR2 power. $n$ is the number of scenario-by-effect-size-by-planning-value evaluation cells in each misspecification bin. Signed error is predicted power minus empirical power. MAE, mean absolute error; RMSE, root mean squared error.",
  "  \\par}",
  "\\end{table}"
))

stratum_order <- c("overall_nonnull", "beta1_0.05", "beta1_0.10", "beta1_0.20", "broad_primary", "targeted_secondary")
stratum_label <- c(
  overall_nonnull = "All non-null conditions",
  beta1_0.05 = "$\\beta_1=.05$",
  beta1_0.10 = "$\\beta_1=.10$",
  beta1_0.20 = "$\\beta_1=.20$",
  broad_primary = "Broad scenarios",
  targeted_secondary = "Targeted scenarios"
)
regime_long <- c(
  sampling_only = "Sampling only",
  near_correct = "Near-correct heterogeneity",
  oracle = "Oracle heterogeneity"
)
s11_rows <- character()
for (reg in regime_order) {
  first <- TRUE
  for (st in stratum_order) {
    row <- s4_pow[s4_pow$regime == reg & s4_pow$stratum == st, , drop = FALSE]
    if (nrow(row) != 1L) stop("Missing Study 4 calibration stratum.", call. = FALSE)
    lead <- if (first) regime_long[[reg]] else ""
    first <- FALSE
    s11_rows <- c(s11_rows, paste0(
      "    ", lead, " & ", stratum_label[[st]], " & ", row$n_cells, " & ",
      fmt_dec(row$rmse, 3L), " & ", fmt_dec(row$mae, 3L), " & ",
      fmt_signed(row$mean_signed_error, 3L), " & ", fmt_dec(row$spearman, 3L), " & ",
      fmt_dec(row$calibration_slope, 3L), " \\\\"
    ))
  }
  if (reg != "oracle") s11_rows <- c(s11_rows, "    \\addlinespace")
}
write_tex(stab_dir, "tableS11_study4_power_calibration.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Power calibration across Study 4 strata.}",
  "  \\label{tab:s11}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.3em}",
  "\\begin{tabular}{@{}p{0.24\\linewidth}p{0.24\\linewidth}rrrrrr@{}}",
  "    \\toprule",
  "    Planning regime & Stratum & $n$ & RMSE & MAE & Signed error & Spearman & Slope \\\\",
  "    \\midrule",
  s11_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} The outcome is MLMA+CR2 power. Sampling only ignores heterogeneity. Near-correct heterogeneity uses each scenario's $\\tau$ and $\\omega$ rounded to the nearest .05, with $\\rho_p=.50$. Oracle heterogeneity uses the true $\\tau$, $\\omega$, and $\\rho$. Broad and targeted rows use the non-null conditions in that part of the library. Signed error is predicted power minus empirical power. The slope is the calibration slope. A robustness-band summary that repeated the all-non-null metrics is not reprinted. RMSE, root mean squared error; MAE, mean absolute error.",
  "  \\par}",
  "\\end{table}"
))
s12_rows <- vapply(regime_order, function(reg) {
  row <- var_row(reg)
  paste0(
    "    ", regime_long[[reg]], " & ", fmt_dec(row$mean_predicted_empirical_ratio, 3L), " & ",
    row$n_cells, " & ", fmt_dec(row$rmse, 6L), " & ", fmt_dec(row$mae, 6L), " & ",
    fmt_dec(row$mean_signed_error, 6L), " & ", fmt_dec(row$spearman, 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS12_study4_variance.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Variance calibration in Study 4.}",
  "  \\label{tab:s12}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.4em}",
  "\\begin{tabular}{@{}lrrrrrr@{}}",
  "    \\toprule",
  "    Planning regime & Ratio & $n$ & RMSE & MAE & Signed error & Spearman \\\\",
  "    \\midrule",
  s12_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} The ratio is the mean of predicted variance divided by the empirical variance of the MLMA+CR2 moderator estimate, over the 120 null scenarios. RMSE, MAE, and signed error are on that variance scale. Signed error is predicted minus empirical. Planning regimes match Table~\\ref{tab:s11}. RMSE, root mean squared error; MAE, mean absolute error.",
  "  \\par}",
  "\\end{table}"
))

# ---- Supplement S13-S17 -----------------------------------------------------

outcome_label <- c(
  A = "Sampling-only variance inflation factor",
  B = "Empirical variance $\\times 10^{3}$",
  C = "Mean CR2/Satterthwaite df",
  D = "Empirical variance $\\times 10^{3}$"
)
expect_label <- c(
  A = "Lower information retention, greater variance inflation",
  B = "Greater contrast balance, lower variance",
  C = "Greater independent-study support $\\rightarrow$ larger df",
  D = "Greater heterogeneity-adjusted contrast information, lower variance"
)
s13_rows <- vapply(seq_len(nrow(s4_pairs)), function(i) {
  letter <- s4_pairs$pair_class[[i]]
  if (letter == "A") {
    va <- fmt_fixed(s4_pairs$VIF_S_a[[i]], 3L)
    vb <- fmt_fixed(s4_pairs$VIF_S_b[[i]], 3L)
  } else if (letter == "C") {
    va <- fmt_fixed(s4_pairs$mean_df_a[[i]], 2L)
    vb <- fmt_fixed(s4_pairs$mean_df_b[[i]], 2L)
  } else {
    va <- fmt_fixed(s4_pairs$empirical_variance_a[[i]] * 1000, 3L)
    vb <- fmt_fixed(s4_pairs$empirical_variance_b[[i]] * 1000, 3L)
  }
  supported <- if (isTRUE(s4_pairs$direction_supported[[i]])) "Yes" else "No"
  paste0(
    "    ", letter, " & ", s4_pairs$pair_id[[i]], " & ", s4_pairs$scenario_id_a[[i]], "/", s4_pairs$scenario_id_b[[i]],
    " & ", outcome_label[[letter]], " & ", va, " & ", vb, " & ", expect_label[[letter]], " & ", supported, " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS13_study4_targeted_pairs.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Targeted matched-pair confirmation in Study 4.}",
  "  \\label{tab:s13}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.32em}",
  "\\begin{tabularx}{\\linewidth}{@{}cc l >{\\raggedright\\arraybackslash}X rr >{\\raggedright\\arraybackslash}X c@{}}",
  "    \\toprule",
  "    Class & Pair & Scenarios & Prespecified outcome & Arm A & Arm B & Expected direction & Supported \\\\",
  "    \\midrule",
  s13_rows,
  "    \\bottomrule",
  "  \\end{tabularx}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Classes A--D each contain four matched pairs. The prespecified outcome is the sampling-only variance inflation factor for Class A, empirical variance of the moderator estimate for Classes B and D, and mean CR2/Satterthwaite degrees of freedom for Class C. Variance entries are multiplied by $10^{3}$. The expected direction is the ordering the pair was built to produce, and Supported records whether that ordering was observed. $R_H$ is information retention, $B_C$ is contrast balance, $J_{\\min}^{*}$ is independent-study support on the weaker side of the contrast, and $I_{C,H}^{*}$ is heterogeneity-adjusted contrast information. Class C pairs were also constructed so that the arm with greater independent-study support has the more even division of that support across the two moderator groups. df, degrees of freedom.",
  "  \\par}",
  "\\end{table}"
))

comp_label <- c(
  N_total = "Total sample size, $N$",
  J = "Number of studies, $J$",
  K_total = "Total effect count, $K$",
  mean_n_target = "Mean study size",
  minority_study_count = "Minority study count",
  J_prec = "Precision-effective number of studies",
  J_k = "Effect-count-effective number of studies",
  I_C_S_star = "Sampling-only contrast information, $I_{C,S}^{*}$",
  I_C_H_star = "Heterogeneity-adjusted contrast information, $I_{C,H}^{*}$",
  R_H = "Information retention, $R_H$",
  B_C = "Contrast balance, $B_C$",
  J_min_star = "Minimum independent-study support, $J_{\\min}^{*}$",
  B_J = "Support balance, $B_J$"
)
if (!identical(s4_comp$comparator, names(comp_label))) {
  stop("Comparator row order does not match the supplement table.", call. = FALSE)
}
s14a <- vapply(seq_len(nrow(s4_comp)), function(i) {
  paste0("    ", comp_label[[s4_comp$comparator[[i]]]], " & ", fmt_dec(s4_comp$spearman_vs_power_0.10[[i]], 3L), " \\\\")
}, character(1))
sup_label <- c(
  J_min_star = "Minimum independent-study support, $J_{\\min}^{*}$",
  B_J = "Support balance, $B_J$",
  nu_C_star = "Anticipated degrees of freedom, $\\nu_C^{*}$"
)
sup_df <- s4_sup[s4_sup$target == "mean_df", , drop = FALSE]
if (!identical(sup_df$diagnostic, names(sup_label))) stop("Support-row order does not match the supplement table.", call. = FALSE)
s14b <- vapply(seq_len(nrow(sup_df)), function(i) {
  paste0("    ", sup_label[[sup_df$diagnostic[[i]]]], " & ", fmt_dec(sup_df$spearman[[i]], 3L), " & ", fmt_dec(sup_df$pearson[[i]], 3L), " \\\\")
}, character(1))
write_tex(stab_dir, "tableS14_comparators_support.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Comparator and independent-study-support associations in Study 4.}",
  "  \\label{tab:s14}",
  "  \\small",
  "  ",
  "\\begin{tabular}{@{}p{0.78\\linewidth}r@{}}",
  "    \\toprule",
  "    \\multicolumn{2}{@{}l}{\\textit{A. Comparator rank correlations with empirical power at $\\beta_1=.10$}}\\\\",
  "    \\midrule",
  "    Comparator & Spearman \\\\",
  "    \\midrule",
  s14a,
  "    \\bottomrule",
  "  \\end{tabular}\\par\\bigskip",
  "  \\begin{tabular}{@{}p{0.62\\linewidth}rr@{}}",
  "    \\toprule",
  "    \\multicolumn{3}{@{}l}{\\textit{B. Independent-study-support diagnostics versus realized CR2/Satterthwaite df}}\\\\",
  "    \\midrule",
  "    Diagnostic & Spearman & Pearson \\\\",
  "    \\midrule",
  s14b,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Panel A uses the 64 broad scenarios and the MLMA+CR2 rejection rate at $\\beta_1=.10$. The precision-effective and effect-count-effective numbers of studies compress unequal study sizes and unequal effect counts, respectively, toward an equivalent count of equal studies. Panel B uses the 120 null scenarios. $B_J$ is the balance of independent-study support across the two moderator groups, constructed in the same way as contrast balance $B_C$. df, degrees of freedom.",
  "  \\par}",
  "\\end{table}"
))

range_lab <- c(
  tau_plan = "Between-study heterogeneity, $\\tau_p$",
  omega_plan = "Within-study heterogeneity, $\\omega_p$",
  rho_plan = "Sampling correlation, $\\rho_p$"
)
held_lab <- c(
  tau_plan = "$\\omega_p$ and $\\rho_p$",
  omega_plan = "$\\tau_p$ and $\\rho_p$",
  rho_plan = "$\\tau_p$ and $\\omega_p$"
)
if (!identical(s4_range$varied_component, names(range_lab))) stop("Planning-range rows are not tau, omega, rho.", call. = FALSE)
s15_rows <- vapply(seq_len(nrow(s4_range)), function(i) {
  key <- s4_range$varied_component[[i]]
  paste0(
    "    ", range_lab[[key]], " & ", held_lab[[key]], " & ",
    s4_range$n_groups[[i]], " & ", fmt_dec(s4_range$mean_range_predicted_power_0.10[[i]], 3L), " & ",
    fmt_dec(s4_range$median_range_predicted_power_0.10[[i]], 3L), " & ",
    fmt_dec(s4_range$max_range_predicted_power_0.10[[i]], 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS15_planning_sensitivity.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Predicted-power sensitivity to planning assumptions.}",
  "  \\label{tab:s15}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.55em}",
  "\\begin{tabular}{@{}llrrrr@{}}",
  "    \\toprule",
  "    Component varied & Held fixed & Evaluation groups & Mean range & Median range & Maximum range \\\\",
  "    \\midrule",
  s15_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Each entry summarizes the range of predicted power at $\\beta_1=.10$ obtained by varying one planning component while holding the other two fixed. An evaluation group is one scenario at one fixed setting of the held components. There are 480 groups when $\\tau_p$ varies, because each of the 120 scenarios is crossed with four values of $\\omega_p$ at $\\rho_p=.50$, and 840 when $\\omega_p$ or $\\rho_p$ varies, because each scenario is crossed with the seven values of $\\tau_p$. The range is a design and planning sensitivity range. It is not a confidence interval and it is not a prediction interval.",
  "  \\par}",
  "\\end{table}"
))

profile_prose <- c(
  low_I_C_S_balanced_contrast_support = "Low $I_{C,S}^{*}$, balanced support",
  moderate_I_C_S_balanced_contrast_support = "Moderate $I_{C,S}^{*}$, balanced support",
  high_I_C_S_balanced_contrast_support = "High $I_{C,S}^{*}$, balanced support",
  high_I_C_S_poor_contrast_balance = "High $I_{C,S}^{*}$, poor balance",
  high_I_C_S_bottleneck_support = "High $I_{C,S}^{*}$, bottleneck support",
  high_I_C_S_strong_alignment = "High $I_{C,S}^{*}$, strong alignment"
)
s4_frag <- s4_frag[order(s4_frag$family_id), , drop = FALSE]
prof <- s4_frag_off$profile_label[match(s4_frag$family_id, s4_frag_off$family_id)]
if (any(!prof %in% names(profile_prose))) stop("Unknown Class E profile label.", call. = FALSE)
s16_rows <- vapply(seq_len(nrow(s4_frag)), function(i) {
  row <- s4_frag[i, , drop = FALSE]
  paste0(
    "    ", row$family_id, " & ", profile_prose[[prof[[i]]]], " & ",
    fmt_dec(row$delta_P_H_0.10, 3L), " & ", fmt_dec(row$delta_P_emp_0.10, 3L), " & ",
    fmt_signed(row$signed_error, 3L), " \\\\"
  )
}, character(1))
write_tex(stab_dir, "tableS16_class_e.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Class E heterogeneity-sensitivity validation.}",
  "  \\label{tab:s16}",
  "  \\small",
  "  ",
  "\\begin{tabular}{@{}clrrr@{}}",
  "    \\toprule",
  "    Family & Profile & Predicted $\\Delta P_H(.10)$ & Empirical power loss & Signed error \\\\",
  "    \\midrule",
  s16_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  paste0(
    "  \\textit{Note.} Each family holds the study structure fixed and varies $\\tau$ from .05 to .30, with $\\omega=.15$. ",
    "Predicted $\\Delta P_H(.10)$ is the change in predicted MLMA+CR2 power for $\\beta_1=.10$ over that range. ",
    "Empirical power loss is the corresponding drop in the simulated rejection rate. Signed error is predicted change minus empirical loss. ",
    "Across the six families, MAE = ", fmt_dec(s4_frag$family_mae[[1L]], 3L),
    ", RMSE = ", fmt_dec(s4_frag$family_rmse[[1L]], 3L),
    ", and Spearman correlation = ", fmt_dec(s4_frag$spearman_delta_P_H_vs_emp[[1L]], 2L),
    ". $I_{C,S}^{*}$, sampling-only contrast information; MAE, mean absolute error; RMSE, root mean squared error."
  ),
  "  \\par}",
  "\\end{table}"
))

mcse_summ <- function(x) {
  x <- as.numeric(x)
  c(n = length(x), mean = mean(x), median = stats::median(x), min = min(x), max = max(x))
}
mc_line <- function(study, x) {
  sm <- mcse_summ(x)
  paste0(
    "    ", study, " & MLMA+CR2 & ", sm[["n"]], " & ",
    fmt_dec(sm[["mean"]], 5L), " & ", fmt_dec(sm[["median"]], 5L), " & ",
    fmt_dec(sm[["min"]], 5L), " & ", fmt_dec(sm[["max"]], 5L), " \\\\"
  )
}
fmt_fail <- function(rate) {
  if (!is.finite(rate) || abs(rate) < 1e-15) "0" else paste0(formatC(100 * rate, format = "f", digits = 4), "\\%")
}
fail_line <- function(study, method, attempted, fail) {
  paste0(
    "    ", study, " & ", method_label[[method]], " & ", fmt_int(attempted), " & ",
    fmt_int(fail), " & ", fmt_fail(fail / attempted), " \\\\"
  )
}
fail_rows <- character()
for (m in sort(unique(s1_sum$method))) {
  sub <- s1_sum[s1_sum$method == m, , drop = FALSE]
  fail_rows <- c(fail_rows, fail_line("1", m, sum(sub$n_attempted), sum(sub$n_fail)))
}
for (stud in list(
  list(id = "2", df = s2_cm, col = "n_error"),
  list(id = "3", df = s3_cm, col = "n_error"),
  list(id = "4", df = s4_cm, col = "n_error")
)) {
  for (m in sort(unique(stud$df$method))) {
    sub <- stud$df[stud$df$method == m, , drop = FALSE]
    fail_rows <- c(fail_rows, fail_line(stud$id, m, sum(sub$n_attempted), sum(sub[[stud$col]])))
  }
}
write_tex(stab_dir, "tableS17_monte_carlo.tex", c(
  "\\begin{table}[H]",
  "  \\centering",
  "  \\caption{Monte Carlo uncertainty and convergence.}",
  "  \\label{tab:s17}",
  "  \\footnotesize",
  "  \\setlength{\\tabcolsep}{0.4em}",
  "\\begin{tabular}{@{}llrrrrr@{}}",
  "    \\toprule",
  "    \\multicolumn{7}{@{}l}{\\textit{A. Rejection-rate Monte Carlo standard errors by study}}\\\\",
  "    \\midrule",
  "    Study & Estimator & $n$ & Mean & Median & Minimum & Maximum \\\\",
  "    \\midrule",
  mc_line("1", s1_mc$mcse_rejection_rate[s1_mc$method == "mlma_cr2"]),
  mc_line("2", s2_mc$mcse_reject_rate[s2_mc$method == "mlma_cr2"]),
  mc_line("3", s3_mc$mcse_reject_rate[s3_mc$method == "mlma_cr2"]),
  mc_line("4", s4_mc$mcse_reject_rate[s4_mc$method == "mlma_cr2"]),
  "    \\bottomrule",
  "  \\end{tabular}\\par\\bigskip",
  "  \\begin{tabular}{@{}llrrr@{}}",
  "    \\toprule",
  "    \\multicolumn{5}{@{}l}{\\textit{B. Attempted fits and failures by study and estimator}}\\\\",
  "    \\midrule",
  "    Study & Estimator & Attempted & Failures & Failure rate \\\\",
  "    \\midrule",
  fail_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\par\\vspace{0.35em}",
  "  {\\footnotesize\\raggedright",
  "  \\textit{Note.} Panel A reports Monte Carlo standard errors of MLMA+CR2 rejection rates. Each production condition used 5{,}000 replications, so a rejection probability near .05 has a Monte Carlo standard error of about .003 and a probability near .50 has one of about .007. Across all reported simulation summaries, the maximum rejection-rate Monte Carlo standard error did not exceed approximately .0071. Panel B counts attempted model fits and fits that did not converge. The failure rate is failures divided by attempted fits, shown as a percentage. Independent RE treats effect sizes as independent. Aggregation + KH aggregates dependent effects and applies the Knapp--Hartung adjustment. MLMA and MLMA+CR2 use the same fitted multilevel coefficients; the stored summaries report identical failure counts for the two approaches. CHE-RVE is the correlated-and-hierarchical-effects robust-variance estimator.",
  "  \\par}",
  "\\end{table}"
))

# ---- Figures ----------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.", call. = FALSE)
measure_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(measure_pdf, width = 8, height = 6, useDingbats = FALSE)
on.exit({
  while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  unlink(measure_pdf)
  unlink(file.path(repro, "Rplots.pdf"))
  unlink("Rplots.pdf")
}, add = TRUE)

theme_main <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text = ggplot2::element_text(size = 8.5, face = "bold", margin = ggplot2::margin(4, 3, 4, 3)),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9),
      legend.key.width = grid::unit(1.2, "lines"),
      legend.key.height = grid::unit(0.9, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 8, colour = "grey15"),
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 8, 2, 4)
    )
}
theme_supp <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position = "bottom",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}
bottom_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p + ggplot2::theme(legend.position = "bottom"))
  idx <- which(g$layout$name == "guide-box-bottom")
  if (length(idx) != 1L) stop("Could not find the bottom legend.", call. = FALSE)
  g$grobs[[idx]]
}
save_two_panel <- function(left, right, legend_plot, path, width, height, legend_height) {
  left <- left + ggplot2::theme(legend.position = "none")
  right <- right + ggplot2::theme(legend.position = "none")
  leg <- bottom_legend(legend_plot)
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = 2, ncol = 2,
      heights = grid::unit(c(1, legend_height), "null"),
      widths = grid::unit(c(1, 1), "null")
    )))
    print(left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid::grid.draw(leg)
    grid::popViewport()
    grid::popViewport()
  }
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  draw()
  grDevices::dev.off()
  invisible(path)
}

fig1_rows <- list()
for (J in c(20, 40, 80)) {
  for (sp in c("50/50", "70/30", "85/15")) {
    d <- s1_ext[s1_ext$method == "mlma_cr2" & s1_ext$J == J & s1_ext$split == sp, , drop = FALSE]
    p10 <- d[abs(d$beta1_true - 0.10) < 1e-8, , drop = FALSE]
    fig1_rows[[length(fig1_rows) + 1L]] <- data.frame(
      J = J, moderator_split = sp, mean_df = mean_of(d$mean_df),
      mean_power_beta10 = mean_of(p10$power), stringsAsFactors = FALSE
    )
  }
}
fig1 <- do.call(rbind, fig1_rows)
fig1$J <- factor(fig1$J, levels = c(20, 40, 80))
fig1$moderator_split <- factor(fig1$moderator_split, levels = c("50/50", "70/30", "85/15"))
j_colors <- c("20" = "#0072B2", "40" = "#009E73", "80" = "#D55E00")
j_shapes <- c("20" = 16, "40" = 17, "80" = 15)
j_lines <- c("20" = "solid", "40" = "longdash", "80" = "dotted")
p1_base <- function(y) {
  ggplot2::ggplot(fig1, ggplot2::aes(x = moderator_split, y = .data[[y]], color = J, shape = J, linetype = J, group = J)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = j_colors, name = "Studies (J)") +
    ggplot2::scale_shape_manual(values = j_shapes, name = "Studies (J)") +
    ggplot2::scale_linetype_manual(values = j_lines, name = "Studies (J)") +
    ggplot2::labs(x = "Moderator allocation") +
    theme_main()
}
p1_df <- p1_base("mean_df") + ggplot2::labs(y = "Mean CR2 / Satterthwaite df", title = "A. Finite-sample support")
p1_pw <- p1_base("mean_power_beta10") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(y = expression("Mean power, " * beta[1] * " = .10"), title = "B. Power")
save_two_panel(p1_df, p1_pw, p1_df, file.path(fig_dir, "figure1b_study1_df_and_power.pdf"), 7.2, 3.95, 0.14)

alloc_raw <- c(
  independent = "Independent",
  n_majority = "Large n\nin majority",
  n_minority = "Large n\nin minority",
  nk_majority = "Large n,k\nin majority",
  nk_minority = "Large n,k\nin minority"
)
alloc_plot <- c(
  "Independent" = "Independent",
  "Large n\nin majority" = "Large n\nin majority",
  "Large n\nin minority" = "Large n\nin minority",
  "Large n,k\nin majority" = "Large n and k\nin majority",
  "Large n,k\nin minority" = "Large n and k\nin minority"
)
fig2_rows <- list()
for (p in names(alloc_raw)) {
  for (m in s2_methods) {
    d0 <- s2_align[s2_align$J == 40 & s2_align$allocation_pattern == p & s2_align$method == m & abs(s2_align$beta1_true) < 1e-12, , drop = FALSE]
    d1 <- s2_align[s2_align$J == 40 & s2_align$allocation_pattern == p & s2_align$method == m & abs(s2_align$beta1_true - 0.10) < 1e-8, , drop = FALSE]
    if (nrow(d0) != 1L || nrow(d1) != 1L) stop("Figure 2 alignment cell is not unique.", call. = FALSE)
    fig2_rows[[length(fig2_rows) + 1L]] <- data.frame(
      allocation = alloc_raw[[p]], method = method_label[[m]],
      type1 = d0$type1_error, type1_mcse = d0$mcse_reject_rate,
      power = d1$power, power_mcse = d1$mcse_reject_rate,
      stringsAsFactors = FALSE
    )
  }
}
fig2 <- do.call(rbind, fig2_rows)
fig2$allocation <- factor(unname(alloc_plot[fig2$allocation]), levels = unname(alloc_plot))
fig2$method <- factor(fig2$method, levels = unname(method_label[s2_methods]))
method_colors <- c("Aggregation + KH" = "#56B4E9", "MLMA" = "#009E73", "MLMA+CR2" = "#0072B2", "CHE-RVE" = "#D55E00")
method_shapes <- c("Aggregation + KH" = 16, "MLMA" = 17, "MLMA+CR2" = 15, "CHE-RVE" = 18)
dodge <- ggplot2::position_dodge(width = 0.72)
p2_base <- ggplot2::ggplot(fig2, ggplot2::aes(x = allocation, color = method, shape = method)) +
  ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
  ggplot2::scale_shape_manual(values = method_shapes, drop = FALSE) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 1.7, linewidth = 0.4)),
    shape = ggplot2::guide_legend(nrow = 1)
  ) +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
  ggplot2::labs(x = NULL, color = NULL, shape = NULL) +
  theme_main() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 7.5, lineheight = 0.92, colour = "grey15"),
    plot.margin = ggplot2::margin(4, 6, 8, 4)
  )
p2_t1 <- p2_base +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.4, colour = "grey25") +
  ggplot2::geom_pointrange(
    ggplot2::aes(y = type1, ymin = type1 - 1.96 * type1_mcse, ymax = type1 + 1.96 * type1_mcse),
    position = dodge, linewidth = 0.35, size = 0.65
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 0.14)) +
  ggplot2::labs(y = "Type I error", title = "A. Null rejection rate")
p2_pw <- p2_base +
  ggplot2::geom_pointrange(
    ggplot2::aes(y = power, ymin = power - 1.96 * power_mcse, ymax = power + 1.96 * power_mcse),
    position = dodge, linewidth = 0.35, size = 0.65
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 0.6)) +
  ggplot2::labs(y = expression("Power, " * beta[1] * " = .10"), title = "B. Power")
save_two_panel(p2_t1, p2_pw, p2_pw, file.path(fig_dir, "figure2c_study2_alignment_j40.pdf"), 7.2, 4.7, 0.13)

s3_plot <- s3_by[s3_by$metric_variant %in% c("primary_design_only", "default_heterogeneity_adjusted"), , drop = FALSE]
s3_plot$panel <- ifelse(
  s3_plot$metric_variant == "primary_design_only",
  "Study 3: Sampling-only",
  "Study 3: Default heterogeneity (tau = omega = 0.10)"
)
pred_map <- data.frame(
  regime = rep(c("sampling_only", "near_correct", "oracle"), each = 3),
  beta1 = rep(c(0.05, 0.10, 0.20), times = 3),
  column = c(
    "predicted_power_S_0.05", "predicted_power_S_0.10", "predicted_power_S_0.20",
    "predicted_power_H_0.05", "predicted_power_H_0.10", "predicted_power_H_0.20",
    "predicted_power_oracle_0.05", "predicted_power_oracle_0.10", "predicted_power_oracle_0.20"
  ),
  stringsAsFactors = FALSE
)
emp4 <- s4_cm[s4_cm$method == "mlma_cr2", , drop = FALSE]
join_pieces <- vector("list", nrow(pred_map))
for (i in seq_len(nrow(pred_map))) {
  b <- pred_map$beta1[[i]]
  e <- emp4[abs(emp4$beta1_true - b) < 1e-8, , drop = FALSE]
  e$predicted <- s4_lib[[pred_map$column[[i]]]][match(as.character(e$scenario_id), as.character(s4_lib$scenario_id))]
  e$regime <- pred_map$regime[[i]]
  e$beta1 <- b
  join_pieces[[i]] <- e
}
s4_join <- do.call(rbind, join_pieces)
s4_join <- s4_join[s4_join$regime %in% c("sampling_only", "near_correct"), , drop = FALSE]
s4_panel <- c(sampling_only = "Study 4: Sampling-only", near_correct = "Study 4: Near-correct heterogeneity")
fig3 <- rbind(
  data.frame(panel = s3_plot$panel, beta1 = s3_plot$beta1_true, predicted = s3_plot$predicted_power, empirical = s3_plot$power, stringsAsFactors = FALSE),
  data.frame(panel = unname(s4_panel[s4_join$regime]), beta1 = s4_join$beta1, predicted = s4_join$predicted, empirical = s4_join$power, stringsAsFactors = FALSE)
)
panel_expr <- c(
  "Study 3: Sampling-only" = "\"Study 3: Sampling-only\"",
  "Study 3: Default heterogeneity (tau = omega = 0.10)" = "\"Study 3: Default heterogeneity (\" * tau == omega * \" = .10)\"",
  "Study 4: Sampling-only" = "\"Study 4: Sampling-only\"",
  "Study 4: Near-correct heterogeneity" = "\"Study 4: Near-correct heterogeneity\""
)
fig3$panel_expr <- factor(unname(panel_expr[fig3$panel]), levels = unname(panel_expr))
fig3$beta_key <- factor(formatC(as.numeric(fig3$beta1), format = "f", digits = 2), levels = c("0.05", "0.10", "0.20"))
beta_labels <- c(
  expression(beta[1] * " = .05"),
  expression(beta[1] * " = .10"),
  expression(beta[1] * " = .20")
)
beta_colors <- c("0.05" = "#0072B2", "0.10" = "#009E73", "0.20" = "#D55E00")
beta_shapes <- c("0.05" = 16, "0.10" = 17, "0.20" = 15)
p3 <- ggplot2::ggplot(fig3, ggplot2::aes(x = predicted, y = empirical, color = beta_key, shape = beta_key)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey35") +
  ggplot2::geom_point(alpha = 0.45, size = 1.45) +
  ggplot2::facet_wrap(~ panel_expr, nrow = 2, labeller = ggplot2::label_parsed) +
  ggplot2::scale_color_manual(values = beta_colors, labels = beta_labels, name = NULL) +
  ggplot2::scale_shape_manual(values = beta_shapes, labels = beta_labels, name = NULL) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2.4)),
    shape = ggplot2::guide_legend(nrow = 1)
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "Predicted power", y = "Empirical power") +
  theme_main() +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 8, 10, 4))
ggplot2::ggsave(file.path(fig_dir, "figure3_power_calibration.pdf"), p3, width = 7.2, height = 6.5, device = grDevices::pdf, useDingbats = FALSE)

type1 <- s1_type1
type1$J_n <- as.numeric(type1$J)
type1$type1_n <- as.numeric(type1$type_I_error)
type1$method <- factor(type1$method, levels = c("naive_re", "mlma", "mlma_cr2"), labels = c("Independent RE", "MLMA", "MLMA+CR2"))
type1$imbalance <- factor(type1$imbalance, levels = c("balanced", "moderate", "severe"), labels = c("Balanced", "Moderate", "Severe"))
type1$split <- factor(type1$split, levels = c("50/50", "70/30", "85/15"))
p_s1 <- ggplot2::ggplot(
  type1,
  ggplot2::aes(x = J_n, y = type1_n, color = method, linetype = method, shape = imbalance, group = interaction(method, imbalance))
) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_wrap(~ split, nrow = 1) +
  ggplot2::scale_color_manual(values = c("Independent RE" = "#D55E00", "MLMA" = "#E69F00", "MLMA+CR2" = "#0072B2"), name = "Estimator") +
  ggplot2::scale_linetype_manual(values = c("Independent RE" = "solid", "MLMA" = "longdash", "MLMA+CR2" = "dotted"), name = "Estimator") +
  ggplot2::scale_shape_manual(values = c("Balanced" = 16, "Moderate" = 17, "Severe" = 15), name = "Effect-count imbalance") +
  ggplot2::scale_x_continuous(breaks = c(20, 40, 80)) +
  ggplot2::guides(
    color = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 0.9)),
    linetype = ggplot2::guide_legend(order = 1),
    shape = ggplot2::guide_legend(order = 2, override.aes = list(color = "black", linetype = 0, size = 2.2))
  ) +
  ggplot2::labs(x = "J", y = "Type I error") +
  theme_supp() +
  ggplot2::theme(
    legend.box = "vertical",
    legend.spacing.y = grid::unit(1, "pt"),
    legend.key.width = grid::unit(1.35, "cm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.text = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_text(size = 9)
  )
ggplot2::ggsave(file.path(sfig_dir, "sfig1_study1_type1.pdf"), p_s1, width = 8.8, height = 4.9, device = grDevices::pdf)

het <- s3_het[
  s3_het$analysis_family == "planning_cell" & s3_het$beta1_group == "all_non_null" &
    s3_het$scenario_set == "all" & abs(s3_het$rho_plan - 0.5) < 1e-8,
  , drop = FALSE
]
het$tau_n <- as.numeric(het$tau_plan)
het$omega_n <- as.numeric(het$omega_plan)
het$mae_n <- as.numeric(het$mae_power)
p_het <- ggplot2::ggplot(het, ggplot2::aes(x = factor(omega_n), y = factor(tau_n), fill = mae_n)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = formatC(mae_n, format = "f", digits = 3)), size = 3.1) +
  ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", name = "Power-prediction MAE") +
  ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = grid::unit(5.2, "cm"), barheight = grid::unit(0.35, "cm"), title.position = "top")) +
  ggplot2::labs(
    x = expression("Within-study heterogeneity planning value, " * omega[p]),
    y = expression("Between-study heterogeneity planning value, " * tau[p])
  ) +
  theme_supp() +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 9), plot.margin = ggplot2::margin(8, 14, 4, 10))
ggplot2::ggsave(file.path(sfig_dir, "sfig3_study3_heterogeneity_mae.pdf"), p_het, width = 7.8, height = 5.8, device = grDevices::pdf)

class_e_labels <- c(
  "1" = "1. Low sampling information, balanced support",
  "2" = "2. Moderate sampling information, balanced support",
  "3" = "3. High sampling information, balanced support",
  "4" = "4. High sampling information, poor balance",
  "5" = "5. High sampling information, bottleneck support",
  "6" = "6. High sampling information, strong alignment"
)
frag_long <- do.call(rbind, lapply(seq_len(nrow(s4_frag)), function(i) {
  data.frame(
    family_id = s4_frag$family_id[[i]],
    tau = c(0.05, 0.12, 0.20, 0.30),
    empirical_power = as.numeric(c(s4_frag$power_0.10_tau_0.05[[i]], s4_frag$power_0.10_tau_0.12[[i]], s4_frag$power_0.10_tau_0.20[[i]], s4_frag$power_0.10_tau_0.30[[i]])),
    predicted_power = as.numeric(c(s4_frag$predicted_power_H_0.10_tau_0.05[[i]], s4_frag$predicted_power_H_0.10_tau_0.12[[i]], s4_frag$predicted_power_H_0.10_tau_0.20[[i]], s4_frag$predicted_power_H_0.10_tau_0.30[[i]])),
    stringsAsFactors = FALSE
  )
}))
frag_long$family_lab <- unname(class_e_labels[as.character(frag_long$family_id)])
frag_series <- rbind(
  data.frame(family_lab = frag_long$family_lab, tau = frag_long$tau, power = frag_long$empirical_power, series = "Empirical MLMA+CR2", stringsAsFactors = FALSE),
  data.frame(family_lab = frag_long$family_lab, tau = frag_long$tau, power = frag_long$predicted_power, series = "Predicted", stringsAsFactors = FALSE)
)
frag_series$family_lab <- factor(frag_series$family_lab, levels = unname(class_e_labels))
frag_series$series <- factor(frag_series$series, levels = c("Empirical MLMA+CR2", "Predicted"))
p_frag <- ggplot2::ggplot(frag_series, ggplot2::aes(x = tau, y = power, color = series, linetype = series, shape = series, group = series)) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~ family_lab, ncol = 2) +
  ggplot2::scale_color_manual(values = c("Empirical MLMA+CR2" = "#0072B2", "Predicted" = "grey25"), name = NULL) +
  ggplot2::scale_linetype_manual(values = c("Empirical MLMA+CR2" = "solid", "Predicted" = "dashed"), name = NULL) +
  ggplot2::scale_shape_manual(values = c("Empirical MLMA+CR2" = 16, "Predicted" = 1), name = NULL) +
  ggplot2::scale_x_continuous(breaks = c(0.05, 0.12, 0.20, 0.30), labels = c(".05", ".12", ".20", ".30")) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  ggplot2::labs(x = expression("Between-study heterogeneity, " * tau), y = expression("Power, " * beta[1] * " = .10")) +
  theme_supp() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"))
ggplot2::ggsave(file.path(sfig_dir, "sfig5_study4_fragility.pdf"), p_frag, width = 8.6, height = 8.2, device = grDevices::pdf)

# ---- Confirm the generated displays -----------------------------------------

snames <- c(
  "tableS1_design_ranges.tex", "tableS2_study1_effect_counts.tex", "tableS3_study2_sample_sizes.tex",
  "tableS4_study1_operating.tex", "tableS5_study2_operating.tex", "tableS6_alignment.tex",
  "tableS7_study3_variance_df.tex", "tableS8_study3_power_variants.tex", "tableS9_study3_targeted_pairs.tex",
  "tableS10_tau_misspecification.tex", "tableS11_study4_power_calibration.tex", "tableS12_study4_variance.tex",
  "tableS13_study4_targeted_pairs.tex", "tableS14_comparators_support.tex", "tableS15_planning_sensitivity.tex",
  "tableS16_class_e.tex", "tableS17_monte_carlo.tex"
)
expected <- c(
  "tables/table1_simulation_program.tex",
  "tables/table3_study1_estimator_summary.tex",
  "figures/figure1b_study1_df_and_power.pdf",
  "figures/figure2c_study2_alignment_j40.pdf",
  "tables/table2_information_profile_landscape.tex",
  "figures/figure3_power_calibration.pdf",
  "tables/table4_study4_confirmation_revised.tex",
  file.path("supplement_tables", snames),
  "supplement_figures/sfig1_study1_type1.pdf",
  "supplement_figures/sfig3_study3_heterogeneity_mae.pdf",
  "supplement_figures/sfig5_study4_fragility.pdf"
)
missing <- expected[!file.exists(file.path(man_out, expected))]
if (length(missing)) stop("Missing generated displays: ", paste(missing, collapse = ", "), call. = FALSE)
message("Regenerated ", length(expected), " manuscript displays.")
