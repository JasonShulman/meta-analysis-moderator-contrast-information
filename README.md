# Reproducibility materials

These materials accompany the manuscript *How Much Information Supports a Moderator Contrast? A Framework for Meta-Analysis With Dependent Effects*. They contain the simulation code for four studies, the frozen design files, and the archived production summaries.

## Quick reproduction of manuscript results

The manuscript displays can be regenerated from the archived summaries and frozen designs. The Monte Carlo simulations do not need to be rerun.

From this directory:

```
Rscript scripts/05_build_results.R
Rscript scripts/06_build_manuscript_outputs.R
```

Script 05 recomputes the Study 3 heterogeneity-sensitivity results and the Study 4 Class E and planning-range quantities. It reads the archived production summaries and the frozen designs, and it writes `results/derived/`. Script 06 reads those archives, the derived tables, and the frozen designs. The second command creates `manuscript/` and writes the regenerated main-text and supplement tables and figures there. Neither script simulates data or refits the Monte Carlo models.

Tables are written to `manuscript/tables/` and `manuscript/supplement_tables/`. Figures are written to `manuscript/figures/` and `manuscript/supplement_figures/`.

## Repository contents

`R/` contains the functions used by the simulations and by the post-processing that rebuilds derived quantities.

`scripts/` contains six entry points. Scripts 01-04 run the four simulation studies. Script 05 rebuilds derived tables from the archived summaries. Script 06 regenerates the manuscript and supplement displays.

`data/` contains the frozen scenario libraries and Study 2 design input.

`results/study1/` through `results/study4/` contain the archived outputs of the original production simulations. `results/derived/` contains the Study 3 heterogeneity-sensitivity tables and the Study 4 Class E and planning-range tables, which are inexpensive calculations from those archives and the frozen designs.

Script 06 creates `manuscript/` and writes the regenerated main-text and supplement tables and figures there.

`environment/` records production software versions and checksums.

## Full simulation reproduction

Each study can be rerun on its own. From this directory, the published 5,000-replication designs are:

```
Rscript scripts/01_run_study1.R --production --output-dir reruns/study1_production
Rscript scripts/02_run_study2.R --production --output-dir reruns/study2_production
Rscript scripts/03_run_study3.R --production --output-dir reruns/study3_production
Rscript scripts/04_run_study4.R --production --output-dir reruns/study4_production
```

`--production` is required for the published replication count, and it also requires `--output-dir`. It uses 5,000 replications and a chunk size of 500. It does not accept `--n-rep`. If `--workers` is omitted, the run uses one worker. The original production runs used `--workers 10` for Study 1, `--workers 14` for Studies 2 and 3, and `--workers 15` for Study 4.

A smaller execution check requires `--n-rep` and `--output-dir`. `--workers` is optional and defaults to 1. `--chunk-size` is optional and defaults to `--n-rep`; when supplied, it must divide `--n-rep` evenly. For example:

```
Rscript scripts/01_run_study1.R --n-rep 1 --workers 1 --output-dir results/study1_check
Rscript scripts/02_run_study2.R --n-rep 1 --workers 1 --output-dir results/study2_check
Rscript scripts/03_run_study3.R --n-rep 1 --workers 1 --output-dir results/study3_check
Rscript scripts/04_run_study4.R --n-rep 1 --workers 1 --output-dir results/study4_check
```

`results/study1/` through `results/study4/` contain the archived outputs from the original production runs. A rerun is written to the directory named by `--output-dir`, for example `reruns/study1_production/`. The runners refuse to write into those archived directories, into `data/`, or into `manuscript/`. A check uses the same restriction and still requires its own `--output-dir`. Relative output paths are resolved inside this directory.

## Study designs

The manuscript and supplement give the full design. The production settings below are the ones used by the runners.

Study 1 has 81 conditions: number of studies J = 20, 40, or 80; effect-count imbalance (balanced, moderate, severe); moderator allocation (50/50, 70/30, 85/15); and moderator effect beta1 = 0, 0.10, or 0.30. The intercept is beta0 = 0.20, the sampling correlation is rho = 0.50, between- and within-study heterogeneity are tau = omega = 0.10, and each study has n_j = 100 participants. Each condition uses 5,000 replications. The estimators are independent random effects, multilevel meta-analysis (MLMA), and MLMA with CR2 standard errors.

Study 2 has 232 conditions: 216 in the core factorial and 16 additional alignment conditions. The moderator effect is beta1 = 0 or 0.10. Relative to Study 1, the core design adds sample-size imbalance, and the alignment conditions place large studies, or large and effect-rich studies, with the majority or minority moderator group. Each condition uses 5,000 replications. The estimators are aggregation with the Knapp-Hartung adjustment, MLMA, MLMA with CR2 standard errors, and CHE-RVE.

Study 3 reads a frozen library of 100 scenarios, 80 broad and 20 targeted. The moderator effect is beta1 = 0, 0.05, 0.10, or 0.20. Each scenario uses 5,000 replications. The public runner reads `data/study3_design/` and does not rebuild the library. The estimators match Study 2.

Study 4 reads a frozen library of 120 scenarios, 64 broad and 56 targeted. The targeted set contains paired comparisons in Classes A-D (eight scenarios in each class) and six Class E fragility families (24 scenarios). Classes A-D hold selected contrast-information or support quantities close within a pair and separate others. Class E varies between-study heterogeneity within a family. The moderator effect is beta1 = 0, 0.05, 0.10, or 0.20, with the confirmatory value beta1 = 0.10. Each scenario uses 5,000 replications. The public runner reads `data/study4_design/` and does not rebuild the library. The estimators are MLMA and MLMA with CR2 standard errors.

## Frozen scenario libraries

Studies 3 and 4 use the scenario libraries stored in `data/study3_design/` and `data/study4_design/`. The public runners read those files. They do not redraw the libraries.

The Study 3 and Study 4 configuration objects still contain the seeds used when the libraries were constructed (`scenario_seed_root` and `candidate_seed_root` for Study 3; `scenario_seed_root`, `pair_candidate_seed`, and `pair_candidate_seed_expanded` for Study 4). Those seeds are separate from the seeds that generate simulated outcomes. Redrawing a Latin-hypercube or candidate pool in another R session is a different operation from the outcome simulation, and these materials do not treat that redraw as a substitute for the deposited libraries.

## Random-number generation

Each replication receives a seed from its condition or scenario identifier and its replication number. Parallel workers do not supply that seed. The runners call `future.apply::future_lapply()` with `future.seed = NULL`, and the seed is applied inside the replication.

Study 1 uses `seed_base = 202700000 + condition_id * 100000` and calls `set.seed(seed_base + rep_id)` inside the replication. Study 2 uses seed root 600000000, a condition stride of 1000000, a replication stride of 100, and a small component code. Study 3 uses outcome seed root 900000000, a scenario stride of 1000000, a replication stride of 100, and a stream code. Study 4 uses outcome seed root 1100000000 with the same strides and stream codes. Studies 2-4 apply those seeds through `withr::with_seed()` with the Mersenne-Twister, inversion normal, and rejection sampling kinds.

## Archived results

`results/study1/` through `results/study4/` are the outputs of the original production simulations. They contain summary and table files. They do not contain the replication-level raw chunk files. The manuscript tables and figures can therefore be rebuilt without rerunning the simulations. A new replication-level Monte Carlo file requires rerunning the relevant study.

Study 1's archived results include the extended summaries used for the manuscript, so its raw chunks are not required for the quick workflow above.

## Computational requirements

The times below are the wall times of the original production runs.

- Study 1: 25,926 seconds (7.2 hours) on 10 workers.
- Study 2: 88,425 seconds (24.6 hours) on 14 workers.
- Study 3: 205,118 seconds (57 hours) on 14 workers.
- Study 4: 163,902 seconds (45.5 hours) on 15 workers.

## Software environment

The production runs used R 4.6.1 on 64-bit Windows. Core fitting packages included metafor 5.0-1, clubSandwich 0.7.0, MASS 7.3-66, and Matrix 1.7-5. 

The complete study-specific environment record is in `environment/production_versions.csv`.

## Main dependencies

Simulation scripts 01-04 require:

- metafor
- clubSandwich
- MASS
- future
- future.apply

Scripts 02-04 also require withr. Scripts 03 and 04 also require digest. Matrix is part of the recorded fitting session because metafor uses it; the scripts do not call Matrix directly.

Display generation in script 06 requires ggplot2 and grid. Script 05 uses base R together with the functions in `R/`.

## Output locations

Simulation output is written only to the directory named by `--output-dir`. The archived directories `results/study1/` through `results/study4/` are not output locations for the public runners.

Derived results are written to `results/derived/study3_heterogeneity/` and `results/derived/study4/`.

Script 06 creates the generated directories `manuscript/tables/`, `manuscript/supplement_tables/`, `manuscript/figures/`, and `manuscript/supplement_figures/`.

## Notes on numerical reproduction

The runners keep the original seed formulas and the original model-fitting calls. The archived summaries were produced under the Windows, R, and package versions recorded in `environment/production_versions.csv`. A rerun on another platform can differ slightly where optimization stops near a convergence boundary. The archived summaries are included so the reported tables and figures can be reconstructed from those results directly.

## License

These materials are available under the MIT License. See `LICENSE` for details.