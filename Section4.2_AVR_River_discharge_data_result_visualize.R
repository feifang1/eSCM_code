rm(list=ls())

library(tidyr)
library(ggplot2)
library(tidyverse)

load("eSCM_code/datasets/river_data_lambda_2_parametric_fit_upper_50")

# Map indices to the corresponding river station names.
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)
station_names[res_rec$ease_order_rec[2,]]
station_names[res_rec$AAC_order_NO_transform_rec[8,]]
station_names[res_rec$AAC_order_NO_transform_rec[2,]]

n<-4600
q_vec<-seq(0.001,0.2,by=0.005)
k_vec<-round(n*q_vec)

# Create a summary table of error rates across different methods.
res_tab<-data.frame(k_vec=k_vec,ease_error_rate=res_rec$ease_error_rec,ACC_error_rate_with_transform=res_rec$AAC_error_transform_rec,ACC_error_rate_without_transform=res_rec$AAC_error_NO_transform_rec)

# Convert the table to long format for plotting.
res_tab_long <- res_tab %>%
  pivot_longer(
    cols = c(ease_error_rate, ACC_error_rate_with_transform, ACC_error_rate_without_transform),
    names_to = "metric",
    values_to = "error_rate"
  )

# Set the order of methods in the legend.
res_tab_long$metric <- factor(res_tab_long$metric,
                              levels = c("ease_error_rate", "ACC_error_rate_with_transform", "ACC_error_rate_without_transform")
)

# Plot error rates as a function of k for different methods.
ggplot(res_tab_long, aes(x = k_vec, y = error_rate, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "k", y = "Ancestral Violation Rate") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("ease_error_rate" = "#1f77b4",
               "ACC_error_rate_with_transform" = "#ff7f0e",
               "ACC_error_rate_without_transform" = "#2ca02c"),
    labels = c(
      expression("EASE with CTC"),
      expression("EASE with AAC (MT)"),
      expression("EASE with AAC (NMT)")
    )
  )

# Save the plot to file.
filename <- paste0("eSCM_code/plots/river_dataset_lambda_2_para_k_fit_tail_adaptive.pdf")
ggsave(filename, width = 6.5, height = 4.5, units = "in")