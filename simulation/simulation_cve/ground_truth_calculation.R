#####################################################
# Superpopulation Truth Estimation for CVE
#   - Computes true CVE, numerator, and denominator
#   - Based on cve.smooth.function()
#####################################################

source('./data_generation.R')
source('./super_population_estimation.R')
# Requires: cve.smooth.function()

compute_truth_cve = function(Npop,
                             s0.eval,
                             s1.eval,
                             a0 = 0, a1 = 1,
                             t0, epsilon, h,
                             s0.seq, s1.seq) {
  
  # Generate superpopulation
  #data_pop = generate_data(Npop)
  data_pop = generate_data_B_discrete_zero(Npop)

  # Create all (s0, s1) combinations
  combo = expand.grid(s0 = s0.eval, s1 = s1.eval)
  n_combo = nrow(combo)
  
  # Storage vectors
  cve_true      = numeric(n_combo)
  cve_num_true  = numeric(n_combo)
  cve_den_true  = numeric(n_combo)
  
  for (i in seq_len(n_combo)) {
    s0_i = combo$s0[i]
    s1_i = combo$s1[i]
    cat("Computing for s0 =", s0_i, ", s1 =", s1_i, "\n")
    
    res = cve.smooth.function(
      data     = data_pop,
      s0       = s0_i,
      s1       = s1_i,
      a0       = a0,
      a1       = a1,
      t        = t0,
      h        = h,
      s0.seq   = s0.seq,
      s1.seq   = s1.seq,
      epsilon  = epsilon
    )

    print(names(res))   
    
    cve_true[i]     = res$cve_true
    cve_num_true[i] = res$cve_num_true
    cve_den_true[i] = res$cve_den_true
  }
  
  # Return tidy data frame for easy plotting / comparison
  truth_df = data.frame(
    s0 = combo$s0,
    s1 = combo$s1,
    cve_true      = cve_true,
    cve_num_true  = cve_num_true,
    cve_den_true  = cve_den_true
  )
  
  return(truth_df)
}

# --------------------------------------------------------------------
# Run superpopulation computation
# --------------------------------------------------------------------
Npop     = 2e5
t0       = 0.1        # Trimming threshold
epsilon  = 0.1       # Smoothing epsilon
h        = 0.1          # Bandwidth for kernel smoothing

# Evaluation points (for s0 and s1 combinations)
s0.eval  = seq(7, 9, by = 1)
s1.eval  = seq(8, 10, by = 1)

# Kernel support sequences
#s0.seq   = seq(4, 11, by = 0.1)
#s1.seq   = seq(4, 11, by = 0.1)
s0.seq   = seq(4, 13, by = 0.1)
s1.seq   = seq(4, 13, by = 0.1)

cat("Running superpopulation truth estimation for CVE...\n")

truth_df = compute_truth_cve(
  Npop     = Npop,
  s0.eval  = s0.eval,
  s1.eval  = s1.eval,
  a0       = 0,
  a1       = 1,
  t0       = t0,
  epsilon  = epsilon,
  h        = h,
  s0.seq   = s0.seq,
  s1.seq   = s1.seq
)

cat("Superpopulation truth estimation completed.\n")
print(head(truth_df))

# Example output:
#   s0 s1 cve_true cve_num_true cve_den_true
# 1  5  6  0.83452    0.03122     0.04709
# 2  6  6  0.84103    0.02855     0.04492
# 3  7  6  0.85487    0.02104     0.03692
# 4  5  7  0.82241    0.03416     0.04815

# Save results
saveRDS(truth_df, "ground_truth_cve_all_pairs_B_discrete_zero.rds")
cat("Saved ground truth CVE estimates to 'ground_truth_cve_all_pairs_B_discrete_zero.rds'\n")
