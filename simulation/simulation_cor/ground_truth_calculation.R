# -------------------------------------------------------------------------
# Compute Ground Truth for Smoothed τ (τ_smooth_true)
# -------------------------------------------------------------------------
source('./data_generation.R')
source('./super_population_estimation.R')  # needed for getSt(), gaussianKernel(), etc.

compute_truth <- function(Npop, s.eval, t0, epsilon, h, s0.seq) {
  #data_pop <- generate_data(Npop)
  #data_pop <- generate_data_B_continuous(Npop)
  data_pop <- generate_data_B_discrete_zero(Npop)
  
  n_eval <- length(s.eval)
  tau_smooth_true <- numeric(n_eval)
  
  for (j in seq_along(s.eval)) {
    s_j <- s.eval[j]
    message("Computing τ_smooth_true for s = ", s_j)
    
    smooth_res <- tau.smoothed.function(
      data      = data_pop,
      s         = s_j,
      t         = t0,
      h         = h,
      s0.seq    = s0.seq,
      epsilon   = epsilon,
      additional_info = FALSE
    )
    
    tau_smooth_true[j] <- smooth_res$tau.h
    print(smooth_res$tau.h)
  }
  
  return(tau_smooth_true)
}

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------
Npop     <- 2e5          # population size for numerical ground truth
t0       <- 0.1          # smoothing threshold t
epsilon  <- 0.1          # smoothing epsilon
h        <- 0.1            # kernel bandwidth
s.eval   = 7:10
#s0.seq   = seq(4, 11, by = 0.05)
s0.seq   = seq(4, 13, by = 0.05)
s0.len   <- s0.seq[2] - s0.seq[1]

# -------------------------------------------------------------------------
# Run computation
# -------------------------------------------------------------------------
message("Running τ_smooth_true computation...")

tau_smooth_true <- compute_truth(
  Npop     = Npop,
  s.eval   = s.eval,
  t0       = t0,
  epsilon  = epsilon,
  h        = h,
  s0.seq   = s0.seq
)

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------
#saveRDS(tau_smooth_true, file = "ground_truth_tau_h01.rds")
saveRDS(tau_smooth_true, file = "ground_truth_tau_h01_B_III.rds")

message("Done! Saved τ_smooth_true to ground_truth_tau_h01_B_III.rds")
