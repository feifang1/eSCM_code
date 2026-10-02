rm(list=ls())

library(tidyr)
library(ggplot2)

load("datasets/Section_4.2/river_data_pairwise_alpha_2_lambda_4_a_0.7_fit_marginal_m_sqrt_n_alpha_2")

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

error_df

# Specify the ordering of methods in the legend.
error_df$Method <- factor(error_df$Method,
                          levels = c("CTC",
                                     "AAC (MT)",
                                     "AAC (NMT)"))

ggplot(
  error_df,
  aes(
    x = k,
    y = ErrorRate,
    color = Method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # Emphasize CTC at k = 29 and AAC at k = 103
  geom_point(
    data = error_df %>%
      filter(
        (Method == "CTC" & k == 29) |
          (Method %in% c("AAC (MT)", "AAC (NMT)") & k == 366)
      ),
    size = 4,
    show.legend = FALSE
  ) +
  
  coord_cartesian(
    ylim = c(0, 0.8)
  ) +
  labs(
    x = expression(k),
    y = "Causal direction error rate"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(
    values = c(
      "CTC" = "#1f77b4",
      "AAC (MT)" = "#ff7f0e",
      "AAC (NMT)" = "#2ca02c"
    )
  )

# # Save the plot to file.
filename <- paste0("plots/Section_4.2_pairwise_violation_rate_river_dataset_lambda_4_m_sqrt_n_a_0.7_alpha_2.pdf")
ggsave(filename, width = 6.5, height = 4.5, units = "in")
