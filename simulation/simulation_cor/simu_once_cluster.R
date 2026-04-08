# --------------------------------------------------------------------------------
# Simulation Script (Smoothed Estimator Only)
# --------------------------------------------------------------------------------
# Purpose: Run simulation using est.tau.smoothed.function only (ignore est.tau.function)
# --------------------------------------------------------------------------------

touse <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(DescTools)  

# --------------------------------------------------------------------------------
# Source functions
# --------------------------------------------------------------------------------
source('./data_generation.R')              # contains generate_data()
source('./simulation_based_estimation.R')  # contains est.tau.smoothed.function()

# --------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------
t0       = 0.1        # trimming threshold (still needed as smoothing param input)
epsilon  = 0.1        # smoothing epsilon
h        = 0.1          # bandwidth for kernel smoothing

# Evaluation points and smoothing grid
# When data is generated from generate_data
s.eval   = 7:10
s0.seq   = seq(4, 11, by = 0.05)

# when data is generated from generate_data_B_continuous
s.eval   = 7:10
s0.seq   = seq(4, 13, by = 0.05)

# Sample sizes
n_list <- c(1000, 2000, 5000)

# --------------------------------------------------------------------------------
# Prepare results data.frame (smoothed estimator only)
# --------------------------------------------------------------------------------
total_rows <- length(n_list) * length(s.eval)
results    <- data.frame(
  rep             = integer(total_rows),
  n               = integer(total_rows),
  s               = numeric(total_rows),
  tau_h_est       = numeric(total_rows),
  tau_h_ci_l      = numeric(total_rows),
  tau_h_ci_u      = numeric(total_rows)
)

row_idx <- 1L

# --------------------------------------------------------------------------------
# Main simulation loop
# --------------------------------------------------------------------------------
for (n in n_list) {
  #data_i <- generate_data(n)
  #data_i <- generate_data_B_continuous(n)
  data_i <- generate_data_B_discrete_zero(n)
  for (s_j in s.eval) {
    # ---- Smoothed estimator + analytic CI ----
    smooth_res <- est.tau.smoothed.function(
      data      = data_i,
      s         = s_j,
      t         = t0,
      epsilon   = epsilon,
      h         = h,
      s0.seq    = s0.seq,
      variance  = TRUE
    )
    
    # Store results
    results$rep[row_idx]       <- touse
    results$n[row_idx]         <- n
    results$s[row_idx]         <- s_j
    results$tau_h_est[row_idx] <- smooth_res$tau.h.est
    results$tau_h_ci_l[row_idx]<- smooth_res$CI_lower
    results$tau_h_ci_u[row_idx]<- smooth_res$CI_upper
    
    row_idx <- row_idx + 1L
  }
}

# --------------------------------------------------------------------------------
# Save results
# --------------------------------------------------------------------------------
dir_name = '/home/qhe2/weighted_cve/simulation/simu_cor_new/results_B_III'
FileSave = paste0(dir_name, '/CoR_B_III_', touse, ".csv")
write.table(results, file = FileSave, row.names = FALSE)
