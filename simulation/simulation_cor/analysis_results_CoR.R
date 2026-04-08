# ------------------------------------------------------------------------
# Combine results from all CoR_<i>.csv files (Smoothed τ_h only)
# ------------------------------------------------------------------------
library(DescTools)
library(dplyr)

combined_df <- data.frame()
for (i in 1:1000) {
  #file_path <- paste0('results/CoR_smoothed_only_', i, '.csv')
  #file_path <- paste0('results_B_conti/CoR_smoothed_only_', i, '.csv')
  file_path <- paste0('results_B_III/CoR_B_III_', i, '.csv')
  if (file.exists(file_path)) {
    res <- read.table(file_path, header = TRUE)
    combined_df <- rbind(combined_df, res)
  }
}
#write.csv(combined_df, "result_data.csv", row.names = FALSE)
#write.csv(combined_df, "result_data_B_conti.csv", row.names = FALSE)
write.csv(combined_df, "result_data_B_III.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Simulation settings (must match estimation script)
# ------------------------------------------------------------------------
#s.eval   <- 7:10
s.eval   <- 7:10
n_eval   <- length(s.eval)
h        <- 0.1
epsilon  <- 0.1
t0       <- 0.1

# ------------------------------------------------------------------------
# Read combined results
# ------------------------------------------------------------------------
#combined_df <- read.csv("result_data.csv", sep = ',')
#combined_df <- read.csv("result_data_B_conti.csv", sep = ',')
combined_df <- read.csv("result_data_B_III.csv", sep = ',')
str(combined_df)

# Keep only n ∈ {1000, 2000, 5000}
df_n1000 <- combined_df %>% filter(n == 1000)
df_n2000 <- combined_df %>% filter(n == 2000)
df_n5000 <- combined_df %>% filter(n == 5000)

df_list <- list(
  "1000" = df_n1000,
  "2000" = df_n2000,
  "5000" = df_n5000
)

# ------------------------------------------------------------------------
# Load ground-truth τ_smooth_true (RDS contains only that vector)
# ------------------------------------------------------------------------
#tau_smooth_true <- readRDS("ground_truth_tau_h01.rds")
#tau_smooth_true <- readRDS("ground_truth_tau_h01_B_continuous.rds")
tau_smooth_true <- readRDS("ground_truth_tau_h01_B_III.rds")

# ------------------------------------------------------------------------
# Summarize τ_h results for each n
# ------------------------------------------------------------------------
summarize_tau_h <- function(df_sub, tau_true_vec, s.eval, h, epsilon) {
  nsims  <- length(unique(df_sub$rep))
  n_eval <- length(s.eval)
  
  tau_est_mat  <- matrix(NA_real_, nrow = nsims, ncol = n_eval)
  tau_ci_l_mat <- matrix(NA_real_, nrow = nsims, ncol = n_eval)
  tau_ci_u_mat <- matrix(NA_real_, nrow = nsims, ncol = n_eval)
  
  for (j in seq_along(s.eval)) {
    s_j <- s.eval[j]
    df_sj <- df_sub %>%
      filter(s == s_j) %>%
      arrange(rep)
    
    tau_est_mat[ , j]  <- df_sj$tau_h_est
    tau_ci_l_mat[ , j] <- df_sj$tau_h_ci_l
    tau_ci_u_mat[ , j] <- df_sj$tau_h_ci_u
  }
  
  # Winsorize each column at 1%/99%
  tau_est_wins <- apply(tau_est_mat, 2, function(col) {
    cutoffs <- quantile(col, probs = c(0.01, 0.98), na.rm = TRUE)
    Winsorize(col, val = cutoffs)
  })
  
  # Bias (%)
  bias_vec <- (colMeans(tau_est_wins) - tau_true_vec) / tau_true_vec * 100
  #bias_vec <- (apply(tau_est_wins, 2, median) - tau_true_vec) / tau_true_vec * 100

  # Coverage probability (95%)
  cover_vec <- colMeans(
    (tau_ci_l_mat <= matrix(tau_true_vec, nrow = nsims, ncol = n_eval, byrow = TRUE)) &
      (tau_ci_u_mat >= matrix(tau_true_vec, nrow = nsims, ncol = n_eval, byrow = TRUE)),
    na.rm = TRUE
  )
  
  # Empirical SD of raw τ_h estimates
  sd_vec <- apply(tau_est_mat, 2, sd, na.rm = TRUE)
  
  # Mean estimated SE from CI half-width
  se_mat <- (tau_ci_u_mat - tau_ci_l_mat) / (2 * qnorm(0.975))
  mean_se_vec <- apply(se_mat, 2, mean, na.rm = TRUE)
  
  n_val <- unique(df_sub$n)
  df_summary <- data.frame(
    n        = rep(n_val,      n_eval),
    h        = rep(h,          n_eval),
    epsilon  = rep(epsilon,    n_eval),
    s        = s.eval,
    bias_pct = round(bias_vec,   2),
    cover95  = round(cover_vec,  3),
    sd_est   = round(sd_vec,     3),
    mean_se  = round(mean_se_vec, 3),
    stringsAsFactors = FALSE
  )
  
  return(df_summary)
}

# Compute summaries for n = 1000, 2000, 5000
summary_n1000 <- summarize_tau_h(df_list[["1000"]], tau_smooth_true, s.eval, h, epsilon)
summary_n2000 <- summarize_tau_h(df_list[["2000"]], tau_smooth_true, s.eval, h, epsilon)
summary_n5000 <- summarize_tau_h(df_list[["5000"]], tau_smooth_true, s.eval, h, epsilon)

final_summary <- bind_rows(summary_n1000, summary_n2000, summary_n5000)

# ------------------------------------------------------------------------
# Print or save final summary
# ------------------------------------------------------------------------
print(final_summary)
#write.csv(final_summary, "summary_tau_h_results.csv", row.names = FALSE)
#write.csv(final_summary, "summary_tau_h_results_B_conti.csv", row.names = FALSE)
write.csv(final_summary, "summary_tau_h_results_B_III.csv", row.names = FALSE)

