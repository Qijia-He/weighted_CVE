# Dealing with positivity violations in mediation analysis via weighted controlled effects, 
with application to assessing immune correlates of protection in antigen-experienced participants

This repository contains the simulation code for paper 'Dealing with positivity violations in mediation analysis via weighted
  controlled effects, with application to assessing immune correlates of
  protection in antigen-experienced participants'.

## Structure

- `simulation/simulation_cor`: STWCR (CoR) simulations
- `simulation/simulation_cve`: STWCRVE (CVE) simulations
- `application`: Stage 1 application analysis

## Run CoR simulations

Run these commands from `simulation/simulation_cor`.

### 1. Compute ground truth

```bash
Rscript ground_truth_calculation.R
```

This step computes the ground-truth target values and saves them as an `.rds` file. The output is used later when summarizing simulation bias and coverage.

### 2. Run simulation jobs

```bash
sbatch simu_once_cluster.sh
```

This step submits the SLURM array job for the Monte Carlo simulations. It generates replicate-level result files for each simulation setting.

### 3. Generate final results

```bash
Rscript analysis_results_CoR.R
```

This step combines the simulation outputs from Step 2 and compares them with the ground truth from Step 1. It produces the final summary files for the CoR analysis.

## Run CVE simulations

Run these commands from `simulation/simulation_cve`.

### 1. Compute ground truth

```bash
Rscript ground_truth_calculation.R
```

This step computes the ground-truth target values and saves them as an `.rds` file. The output is used later in the final simulation summary.

### 2. Run simulation jobs

```bash
sbatch simu_once_cve_cluster.sh
```

This step submits the SLURM array job for the CVE Monte Carlo simulations. It generates replicate-level result files for each simulation setting.

### 3. Generate final results

```bash
Rscript analysis_results_CVE_log.R
```

This step combines the simulation outputs from Step 2 and compares them with the ground truth from Step 1. It produces the final summary files for the CVE analysis.

## Run application analysis

Run this command from `application`.

### 1. Run the main analysis

```bash
Rscript rerun_stage1_cve_log_ci.R
```

This script reruns the Stage 1 CVE application analysis used for the main manuscript figures and tables. It produces the application outputs, including the Figure 3 and Figure 4 results and the Appendix Table S1 results.

### Supporting file

```bash
functions_application_superlearner_stage1_cde11.R
```

This file contains the estimation functions used by the main application script.

## Notes

- The current scripts are configured by commenting and uncommenting lines for different simulation scenarios.
- If you switch scenarios, update the settings consistently in the ground-truth, simulation, and analysis scripts within the same folder.
- Some output paths are currently hard-coded for the cluster environment and may need to be adjusted before running.
