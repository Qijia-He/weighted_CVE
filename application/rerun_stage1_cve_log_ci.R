library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)

source("./functions_application_superlearner_stage1_cde11.R")

prepare_stage1_data <- function() {
  data <- read.csv("covail_data_processed_20240313.csv")

  dt <- data %>%
    filter(TrtonedosemRNA == 1, stage == 1) %>%
    select(
      COVIDIndD22toD181,
      TrtC,
      Bpseudoneutid50_BA.1,
      Day15pseudoneutid50_BA.1,
      naive,
      Age
    )

  names(dt) <- c("Y", "A", "B", "S", "X1", "X2")
  dt <- na.omit(dt)

  dt$A <- as.numeric(dt$A)
  dt$Y <- as.numeric(dt$Y)
  dt$X1 <- as.numeric(dt$X1)

  dt <- dt[dt$S - dt$B > 0, ]
  rownames(dt) <- NULL

  dt
}

run_cve_points <- function(dt, s0_vals, s1_vals, estimator_fun, a0, a1,
                           t, epsilon, h, s0.seq, s1.seq, label) {
  stopifnot(length(s0_vals) == length(s1_vals))

  results <- Map(
    function(s0, s1) {
      out <- estimator_fun(
        data = dt,
        s0 = s0,
        s1 = s1,
        a0 = a0,
        a1 = a1,
        t = t,
        h = h,
        s0.seq = s0.seq,
        s1.seq = s1.seq,
        epsilon = epsilon
      )

      data.frame(
        analysis = label,
        s0 = s0,
        s1 = s1,
        cve = out$cve_est,
        cve_num = out$cve_num,
        cve_den = out$cve_den,
        se = out$se,
        ci_low_normal = out$ci[1],
        ci_high_normal = out$ci[2],
        ci_low = out$ci_ratio_log_to_cve[1],
        ci_high = out$ci_ratio_log_to_cve[2]
      )
    },
    s0_vals,
    s1_vals
  )

  bind_rows(results)
}

plot_cve_curve <- function(df, x_var, ylab_text, x_breaks,
                           y_limits = c(-1, 1),
                           y_breaks = seq(-1, 1, by = 0.2)) {
  ggplot(df, aes(x = .data[[x_var]], y = cve)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.8
    ) +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.02,
      linewidth = 0.9,
      color = "#2166AC"
    ) +
    geom_point(
      size = 4,
      color = "#2166AC"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      x = if (x_var == "s1") expression(s[1]) else expression(s),
      y = ylab_text,
      title = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      panel.grid.major = element_line(linewidth = 0.4, color = "gray90"),
      panel.grid.minor = element_blank()
    )
}

fmt_pct1 <- function(x) sprintf("%.1f%%", 100 * x)
fmt_num3 <- function(x) sprintf("%.3f", x)

dt <- prepare_stage1_data()

s.eval <- seq(3.5, 4, by = 0.05)
s0.seq <- seq(3.3, 4.3, by = 0.01)
s1.seq <- seq(3.3, 4.3, by = 0.01)
t <- 0.1
epsilon <- 0.1
h <- 0.1

cat("Stage 1 analysis sample size:", nrow(dt), "\n")
cat("Arm counts after filtering:\n")
print(table(dt$A))

figure3_left <- run_cve_points(
  dt = dt,
  s0_vals = rep(3.5, length(seq(3.55, 4.00, by = 0.05))),
  s1_vals = seq(3.55, 4.00, by = 0.05),
  estimator_fun = est.cve.eif.smooth.function.new,
  a0 = 1,
  a1 = 1,
  t = t,
  epsilon = epsilon,
  h = h,
  s0.seq = s0.seq,
  s1.seq = s1.seq,
  label = "Figure 3 left"
)

figure3_right <- run_cve_points(
  dt = dt,
  s0_vals = rep(3.7, length(seq(3.75, 4.00, by = 0.05))),
  s1_vals = seq(3.75, 4.00, by = 0.05),
  estimator_fun = est.cve.eif.smooth.function.new,
  a0 = 1,
  a1 = 1,
  t = t,
  epsilon = epsilon,
  h = h,
  s0.seq = s0.seq,
  s1.seq = s1.seq,
  label = "Figure 3 right"
)

figure4_df <- run_cve_points(
  dt = dt,
  s0_vals = s.eval,
  s1_vals = s.eval,
  estimator_fun = est.cve.eif.smooth.function,
  a0 = 0,
  a1 = 1,
  t = t,
  epsilon = epsilon,
  h = h,
  s0.seq = s0.seq,
  s1.seq = s1.seq,
  label = "Figure 4"
) %>%
  mutate(s = s0)

table_s1_full <- run_cve_points(
  dt = dt,
  s0_vals = c(3.50, 3.50, 3.50, 3.75, 3.75, 3.75, 4.00, 4.00, 4.00),
  s1_vals = c(3.50, 3.75, 4.00, 3.50, 3.75, 4.00, 3.50, 3.75, 4.00),
  estimator_fun = est.cve.eif.smooth.function,
  a0 = 0,
  a1 = 1,
  t = t,
  epsilon = epsilon,
  h = h,
  s0.seq = s0.seq,
  s1.seq = s1.seq,
  label = "Table S1"
)

table_s1_offdiag <- table_s1_full %>%
  filter(s0 != s1) %>%
  mutate(
    cve_pct = fmt_pct1(cve),
    cve_num_pct = fmt_pct1(cve_num),
    cve_den_pct = fmt_pct1(cve_den),
    se_fmt = fmt_num3(se),
    ci_normal_fmt = sprintf(
      "(%s, %s)",
      fmt_pct1(ci_low_normal),
      fmt_pct1(ci_high_normal)
    ),
    ci_fmt = sprintf("(%s, %s)", fmt_pct1(ci_low), fmt_pct1(ci_high)),
    latex_row = sprintf(
      "%.2f & %.2f & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\",
      s0, s1, cve_pct, cve_num_pct, cve_den_pct, se_fmt, ci_fmt
    )
  )

figure3_left_plot <- plot_cve_curve(
  figure3_left,
  x_var = "s1",
  ylab_text = "STWCRVE(1, 1, s1, s0=3.5)",
  x_breaks = figure3_left$s1
)

figure3_right_plot <- plot_cve_curve(
  figure3_right,
  x_var = "s1",
  ylab_text = "STWCRVE(1, 1, s1, s0=3.7)",
  x_breaks = figure3_right$s1
)

figure4_plot <- plot_cve_curve(
  figure4_df,
  x_var = "s",
  ylab_text = "Controlled Direct Effect",
  x_breaks = figure4_df$s,
  y_limits = c(-1.2, 1.2),
  y_breaks = seq(-1.2, 1.2, by = 0.2)
)

ggsave(
  filename = "STWCRVE_s0_3.5.pdf",
  plot = figure3_left_plot,
  width = 10,
  height = 8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "STWCRVE_s0_3.7.pdf",
  plot = figure3_right_plot,
  width = 10,
  height = 8,
  units = "in",
  device = cairo_pdf
)

pdf("STWCRVE_two_panels.pdf", width = 20, height = 8)
grid.arrange(figure3_left_plot, figure3_right_plot, ncol = 2)
dev.off()

ggsave(
  filename = "controlled_direct_effect.pdf",
  plot = figure4_plot,
  width = 10,
  height = 8,
  units = "in",
  device = cairo_pdf
)

write.csv(figure3_left, "figure3_left_log_ci.csv", row.names = FALSE)
write.csv(figure3_right, "figure3_right_log_ci.csv", row.names = FALSE)
write.csv(figure4_df, "figure4_controlled_direct_effect_log_ci.csv", row.names = FALSE)
write.csv(table_s1_full, "table_S1_full_log_ci.csv", row.names = FALSE)
write.csv(table_s1_offdiag, "table_S1_offdiag_log_ci.csv", row.names = FALSE)
writeLines(table_s1_offdiag$latex_row, "table_S1_offdiag_log_ci_rows.tex")

saveRDS(
  list(
    sample_n = nrow(dt),
    arm_counts = table(dt$A),
    figure3_left = figure3_left,
    figure3_right = figure3_right,
    figure4 = figure4_df,
    table_s1_full = table_s1_full,
    table_s1_offdiag = table_s1_offdiag
  ),
  file = "stage1_cve_log_ci_outputs.rds"
)

cat("\nFigure 3 left:\n")
print(figure3_left)

cat("\nFigure 3 right:\n")
print(figure3_right)

cat("\nFigure 4:\n")
print(figure4_df)

cat("\nTable S1 off-diagonal rows:\n")
print(table_s1_offdiag)
