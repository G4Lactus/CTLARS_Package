rm(list = ls()) # clear the workspace

library(ctlars)
library(ggplot2)

set.seed(824)
t_vec <- c(1, 2, 5, 10, 20, 50, 100)
sim_res <- simulate_model_TPR_FDR(
  t_vec,
  MC = 500,
  n = 100,
  p = 300,
  num_actives = 10,
  num_dummies = 300,
  snr_is_linear = TRUE,
  snr_val_linear = 1
)

# data frame containing data to be plotted (FDR and TPR in %)
# -------------------------------------------------------------
plot_data <- data.frame(T_vec = t_vec,
                        FDR = 100 * sim_res$FDR,
                        TPR = 100 * sim_res$TPR)

# FDR vs. T
# ------------
fdr_vs_t <- ggplot(plot_data, aes(x = t_vec, y = FDR)) +
  labs(x = "T", y = "FDR") +
  scale_x_continuous(breaks = t_vec[-2], minor_breaks = c(2),
                     limits = c(t_vec[1], t_vec[length(t_vec)])) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     minor_breaks = c(), limits = c(0, 100)) +
  geom_line(linewidth = 1.5, colour = "#336C68") +
  geom_point(size = 2.5, colour = "#336C68") +
  theme_bw(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white",
                                    color = "black",
                                    linewidth = 1)
  ) +
  coord_fixed(ratio =  0.85 * t_vec[length(t_vec)] / (100 - 0))
plot(fdr_vs_t)


# TPR vs. T
# ------------
tpr_vs_t <- ggplot(plot_data, aes(x = t_vec, y = TPR)) +
  labs(x = "T", y = "TPR") +
  scale_x_continuous(breaks = t_vec[-2], minor_breaks = c(2),
                     limits = c(t_vec[1], t_vec[length(t_vec)])) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     minor_breaks = c(), limits = c(0, 100)) +
  geom_line(linewidth = 1.5, colour = "#336C68") +
  geom_point(size = 2.5, colour = "#336C68") +
  theme_bw(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white",
                                    color = "black",
                                    linewidth = 1)
  ) +
  coord_fixed(ratio =  0.85 * t_vec[length(t_vec)] / (100 - 0))
plot(tpr_vs_t)


# TPR vs. FDR
# ------------
tpr_vs_fdr <- ggplot(plot_data, aes(x = FDR, y = TPR)) +
  labs(x = "FDR",
       y = "TPR") +
  geom_line(linewidth = 1.5, colour = "#336C68") +
  geom_point(size = 2.5, colour = "#336C68") +
  theme_bw(base_size = 16) +
  theme(panel.background = element_rect(fill = "white",
                                        color = "black",
                                        linewidth = 1))
plot(tpr_vs_fdr)
