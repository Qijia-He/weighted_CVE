# ================================================================
# Simulation Script for est.cve.eif.smooth.function  (EIF-based CVE)
# ================================================================

touse <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(dplyr)
library(DescTools)

source('./data_generation.R')
source('./simulation_based_estimation.R')

# --------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------
t0       = 0.1        # Trimming threshold
epsilon  = 0.1        # Smoothing epsilon
h        = 0.1           # Bandwidth for kernel smoothing

# Evaluation points (for s0 and s1 combinations)
s0.eval  = seq(7, 9, by = 1)
s1.eval  = seq(8, 10, by = 1)

# Kernel support sequences
#s0.seq   = seq(4, 11, by = 0.1)
#s1.seq   = seq(4, 11, by = 0.1)
s0.seq   = seq(4, 13, by = 0.1)
s1.seq   = seq(4, 13, by = 0.1)
# Sample sizes
n_list <- c(1000, 2000, 5000)

# Create all (s0, s1) combinations
combos <- expand.grid(s0 = s0.eval, s1 = s1.eval)
n_combo <- nrow(combos)

# --------------------------------------------------------------------------------
# Result storage
# --------------------------------------------------------------------------------
total_rows <- length(n_list) * n_combo
results <- data.frame(
  rep              = integer(total_rows),
  n                = integer(total_rows),
  s0               = numeric(total_rows),
  s1               = numeric(total_rows),
  eif_est          = numeric(total_rows),
  ci_lower         = numeric(total_rows),
  ci_upper         = numeric(total_rows),
  ci_ratio_log_to_cve_lower = numeric(total_rows),
  ci_ratio_log_to_cve_upper = numeric(total_rows),
  true_value       = numeric(total_rows),
  se               = numeric(total_rows)
)

row_idx <- 1L

# --------------------------------------------------------------------------------
# Main simulation loop
# --------------------------------------------------------------------------------
for (n in n_list) {
  cat("Simulating for n =", n, "\n")
  # data_i <- generate_data(n) # I
  data_i <- generate_data_B_continuous(n) # II
  # data_i <- generate_data_B_discrete_zero(n) # III

  for (i in seq_len(n_combo)) {
    s0_i <- combos$s0[i]
    s1_i <- combos$s1[i]
    
    cat("  → s0 =", s0_i, ", s1 =", s1_i, "\n")
    
    res <- est.cve.eif.smooth.function(
      data     = data_i,
      s0       = s0_i,
      s1       = s1_i,
      a0       = 0,
      a1       = 1,
      t        = t0,
      h        = h,
      s0.seq   = s0.seq,
      s1.seq   = s1.seq,
      epsilon  = epsilon
    )
    
    results$rep[row_idx]          <- touse
    results$n[row_idx]            <- n
    results$s0[row_idx]           <- s0_i
    results$s1[row_idx]           <- s1_i
    results$eif_est[row_idx]      <- res$cve_est
    results$ci_lower[row_idx]     <- res$ci[1]
    results$ci_upper[row_idx]     <- res$ci[2]
    results$ci_ratio_log_to_cve_lower[row_idx] <- res$ci_ratio_log_to_cve[1]
    results$ci_ratio_log_to_cve_upper[row_idx] <- res$ci_ratio_log_to_cve[2]
    results$true_value[row_idx]   <- res$cve_true
    results$se[row_idx]           <- res$se
    
    row_idx <- row_idx + 1L
  }
}

# --------------------------------------------------------------------------------
# Save results
# --------------------------------------------------------------------------------
dir_name <- 'results_II'
if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)

FileSave <- file.path(dir_name, paste0('CVE_est_II_', touse, '.csv'))
write.table(results, file = FileSave, row.names = FALSE, sep = ',')

cat("Simulation finished. Results saved to:", FileSave, "\n")
