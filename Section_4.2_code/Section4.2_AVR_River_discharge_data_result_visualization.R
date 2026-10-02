rm(list=ls())
library(tidyr)
library(ggplot2)
library(tidyverse)

load("datasets/Section_4.2/river_data_whole_causal_order_lambda_4_a_0.7_parametric_fit_upper_fixed_cut_m_sqrt_n_alpha_2")

# Map indices to the corresponding river station names.
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)
station_names[res_rec$ease_order_rec[1,]]
station_names[res_rec$AAC_order_NO_transform_rec[1,]]
station_names[res_rec$AAC_order_NO_transform_rec[1,]]

n<-4600
q_vec<-seq(0.001,0.2,by=0.005)
k_vec<-round(n*q_vec)

k_vec<-c(k_vec[1:2],29,k_vec[3:16],366,k_vec[17:length(k_vec)])

#Create a summary table of error rates across different methods.
res_tab<-data.frame(k_vec=k_vec,ease_error_rate=res_rec$ease_error_rec,ACC_error_rate_with_transform=res_rec$AAC_error_transform_rec,ACC_error_rate_without_transform=res_rec$AAC_error_NO_transform_rec)

# Convert the table to long format for plotting.
res_tab_long <- res_tab %>%
  pivot_longer(
    cols = c(ease_error_rate, ACC_error_rate_with_transform, ACC_error_rate_without_transform),
    names_to = "metric",
    values_to = "error_rate"
  )

res_tab_long


# Set the order of methods in the legend.
res_tab_long$metric <- factor(res_tab_long$metric,
                              levels = c("ease_error_rate", "ACC_error_rate_with_transform", "ACC_error_rate_without_transform")
)

# Plot error rates as a function of k for different methods.
ggplot(
  res_tab_long,
  aes(
    x = k_vec,
    y = error_rate,
    color = metric
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # Emphasize selected k values:
  # CTC at k = 29, AAC at k = 366
  geom_point(
    data = res_tab_long %>%
      filter(
        (metric == "ease_error_rate" & k_vec == 29) |
          (metric %in% c(
            "ACC_error_rate_with_transform",
            "ACC_error_rate_without_transform"
          ) & k_vec == 366)
      ),
    size = 4,
    show.legend = FALSE
  )+
  
  coord_cartesian(
    ylim = c(0, 0.8)
  ) +
  labs(
    x = expression(k),
    y = "Ancestral violation rate"
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
      "ease_error_rate" = "#1f77b4",
      "ACC_error_rate_with_transform" = "#ff7f0e",
      "ACC_error_rate_without_transform" = "#2ca02c"
    ),
    labels = c(
      expression("EASE with CTC"),
      expression("EASE with AAC (MT)"),
      expression("EASE with AAC (NMT)")
    )
  )


# Save the plot to file.
filename <- paste0("plots/Section_4.2_river_dataset_AVR_lambda_4_a_0.7_m_sqrt_n_alpha_2.pdf")
ggsave(filename, width = 6.5, height = 4.5, units = "in")
