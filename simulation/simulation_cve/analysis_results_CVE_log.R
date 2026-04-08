#############################################################
# Performance Evaluation for CVE Estimator
# -----------------------------------------------------------
# 1. Combine 1000 “CVE_est_<i>.csv” files
# 2. Load superpopulation truth (“ground_truth_cve_all_pairs.rds”)
# 3. Compute bias, coverage, empirical SD, and mean SE
#############################################################

library(dplyr)
library(DescTools)

# -----------------------------------------------------------
# Step 1. Combine all simulation results
# -----------------------------------------------------------

combined_df <- data.frame()
for (i in 1:1000) {
  # file_path = paste0("results_III/CVE_est_B_discrete_zero_", i, ".csv")
  # file_path = paste0("results_I/CVE_est_I_", i, ".csv")
  file_path = paste0("results_II/CVE_est_II_", i, ".csv")
  if (file.exists(file_path)) {
    res = read.csv(file_path, header = TRUE)
    combined_df = rbind(combined_df, res)
  }
}

# write.csv(combined_df, "combined_CVE_results_III.csv", row.names = FALSE)
# write.csv(combined_df, "combined_CVE_results_I.csv", row.names = FALSE)
write.csv(combined_df, "combined_CVE_results_II.csv", row.names = FALSE)
# cat("✅ Combined 1000 CVE result files into 'combined_CVE_results_III.csv'\n")
# cat("✅ Combined 1000 CVE result files into 'combined_CVE_results_I.csv'\n")
cat("✅ Combined 1000 CVE result files into 'combined_CVE_results_II.csv'\n")

# -----------------------------------------------------------
# Step 2. Load ground truth for CVE
# -----------------------------------------------------------
# truth_df <- readRDS("ground_truth_cve_all_pairs_B_discrete_zero.rds") # III
# truth_df <- readRDS("ground_truth_cve_all_pairs.rds") # I
truth_df <- readRDS("ground_truth_cve_all_pairs_B_conti.rds") # II

# Expected columns: s0, s1, cve_true, cve_num_true, cve_den_true
str(truth_df)

# -----------------------------------------------------------
# Step 3. Basic parameters and data split
# -----------------------------------------------------------
h        <- 0.1
epsilon  <- 0.1
t0       <- 0.1
n_list   <- c(1000, 2000, 5000)

# combined_df <- read.csv("combined_CVE_results_III.csv") # III
# combined_df <- read.csv("combined_CVE_results_I.csv") # I
combined_df <- read.csv("combined_CVE_results_II.csv") # II

# Split by n
df_list <- lapply(n_list, function(nv) combined_df %>% filter(n == nv))
names(df_list) <- as.character(n_list)

# -----------------------------------------------------------
# Step 4. Define summarizing function
# -----------------------------------------------------------
summarize_for_n_cve <- function(df_sub, truth_df) {
  nsims <- length(unique(df_sub$rep))
  combo <- unique(df_sub[, c("s0", "s1")])
  n_combo <- nrow(combo)
  
  bias_vec        <- numeric(n_combo)
  cover_vec       <- numeric(n_combo)
  cover_log_vec   <- numeric(n_combo)
  sd_vec          <- numeric(n_combo)
  mean_se_vec     <- numeric(n_combo)

  for (i in seq_len(n_combo)) {
    s0_i <- combo$s0[i]
    s1_i <- combo$s1[i]

    df_pair <- df_sub %>%
      filter(s0 == s0_i, s1 == s1_i) %>%
      arrange(rep)

    # True value lookup
    cve_true_i <- truth_df %>%
      filter(s0 == s0_i, s1 == s1_i) %>%
      pull(cve_true)

    est_vals    <- df_pair$eif_est
    ci_l        <- df_pair$ci_lower
    ci_u        <- df_pair$ci_upper
    ci_log_l    <- df_pair$ci_ratio_log_to_cve_lower
    ci_log_u    <- df_pair$ci_ratio_log_to_cve_upper

    # Bias (%)
    winsorized_est_vals <- Winsorize(est_vals, val = c(0.01, 0.99))
    bias_vec[i] <- mean(winsorized_est_vals - cve_true_i) / cve_true_i * 100

    # Coverage (original CI)
    cover_vec[i] <- mean(ci_l <= cve_true_i & ci_u >= cve_true_i, na.rm = TRUE)

    # Coverage (log-ratio CI mapped back to CVE)
    cover_log_vec[i] <- mean(ci_log_l <= cve_true_i & ci_log_u >= cve_true_i, na.rm = TRUE)

    # Empirical SD
    sd_vec[i] <- sd(est_vals, na.rm = TRUE)

    # Mean estimated SE (from half-width of original CI)
    se_est <- (ci_u - ci_l) / (2 * qnorm(0.975))
    mean_se_vec[i] <- mean(se_est, na.rm = TRUE)
  }

  n_val <- unique(df_sub$n)

  data.frame(
    n            = rep(n_val, n_combo),
    h            = rep(h, n_combo),
    epsilon      = rep(epsilon, n_combo),
    s0           = combo$s0,
    s1           = combo$s1,
    bias_pct     = round(bias_vec, 2),
    cover95      = round(cover_vec, 3),
    cover95_log  = round(cover_log_vec, 3),
    sd_est       = round(sd_vec, 3),
    mean_se      = round(mean_se_vec, 3),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------
# Step 5. Apply summary function for each n
# -----------------------------------------------------------
summary_list <- lapply(df_list, summarize_for_n_cve, truth_df = truth_df)
final_summary <- bind_rows(summary_list, .id = "n_label")

# -----------------------------------------------------------
# Step 6. Save and print
# -----------------------------------------------------------
# write.csv(final_summary, "summary_CVE_performance_III.csv", row.names = FALSE) # III
# write.csv(final_summary, "summary_CVE_performance_I.csv", row.names = FALSE) # I
write.csv(final_summary, "summary_CVE_performance_II.csv", row.names = FALSE) # II
# cat("✅ Saved final summary to 'summary_CVE_performance_III.csv'\n") # III
# cat("✅ Saved final summary to 'summary_CVE_performance_I.csv'\n") # I
cat("✅ Saved final summary to 'summary_CVE_performance_II.csv'\n") # II

print(final_summary)

