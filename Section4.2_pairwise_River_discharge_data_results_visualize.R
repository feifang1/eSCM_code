rm(list=ls())

library(tidyr)
library(ggplot2)

load("eSCM_code/datasets/river_data_pairwise_semi_alpha_2_lambda_2_fit_marginal")

# Function to compute the error rate across all k
error_rate <- function(x) {
  est_cause_cols <- grep("^est_cause_k_", names(x), value = TRUE)
  sapply(est_cause_cols, function(col) {
    mean(x[[col]] != x$truth_cause)
  })
}

# Compute error rates for different methods.
AAC_NMT_error_rate <- error_rate(fin_res$AAC_NMT_order)
ease_CTC_error_rate <- error_rate(fin_res$ease_CTC_order)
AAC_MT_error_rate <- error_rate(fin_res$AAC_MT_order)

# Construct a data frame and reshape to long format for plotting.
error_df <- data.frame(
  k = fin_res[[9]],
  "AAC (NMT)" = AAC_NMT_error_rate,
  "CTC" = ease_CTC_error_rate,
  "AAC (MT)" = AAC_MT_error_rate,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = -k,
    names_to = "Method",
    values_to = "ErrorRate",
    names_transform = list(Method = as.character)
  )

# Specify the ordering of methods in the legend.
error_df$Method <- factor(error_df$Method,
                          levels = c("CTC",
                                     "AAC (MT)",
                                     "AAC (NMT)"))

# Plot error rates as a function of k for different methods.
ggplot(error_df, aes(x = k, y = ErrorRate, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "k", y = "Causal Direction Error Rate") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c(
    "CTC" = "#1f77b4",
    "AAC (MT)" = "#ff7f0e",
    "AAC (NMT)" = "#2ca02c"
  ))

# Save the plot to file.
filename <- paste0("eSCM_code/plots/pairwise_river_dataset_lambda_2_fit_tail.pdf")
ggsave(filename, width = 6.5, height = 4.5, units = "in")
