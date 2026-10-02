rm(list=ls())

library(tidyr)
library(ggplot2)
load("datasets/Section_4.3/pairwise_datasets_lambda_4_a_0.7_m_sqrt_n_alpha_2_with_four_quadrant_fit_with_tail_strong_stop_adjust_m")

lambda<-4

# Add the original and normalized weights for each dataset.
weight_tab<-read.table("datasets/pairs/pairmeta.txt", header = F)
colnames(weight_tab)[6]<-"weight"

# Dataset pairs to exclude:
# Exclude pairs with categorical variables; exclude high-dimensional datasets.
remov_ind<-c(33:37,47,52:55,70:71,105,107)

weight<-weight_tab$weight[-remov_ind]
nor_weight<-(weight_tab$weight[-remov_ind]/sum(weight_tab$weight[-remov_ind]))

names_AAC_order_list<-names(fin_res$AAC_order_list)

fin_res$AAC_order_list <- lapply(seq_along(fin_res$AAC_order_list), function(i) {
  df <- fin_res$AAC_order_list[[i]]
  df$weight <- weight
  df$nor_weight <- nor_weight
  df
})

names(fin_res$AAC_order_list) <- names_AAC_order_list

length(fin_res$AAC_order_list)

weighted_correct_rate_first_quad<-function(x){
  
  return( sum((x$truth_cause==x$first_quad_cause)*x$nor_weight) )
  
}

weighted_correct_rate_second_quad<-function(x){
  
  return( sum((x$truth_cause==x$second_quad_cause)*x$nor_weight) )
  
}

weighted_correct_rate_third_quad<-function(x){
  
  return( sum((x$truth_cause==x$third_quad_cause)*x$nor_weight) )
  
}

weighted_correct_rate_fourth_quad<-function(x){
  
  return( sum((x$truth_cause==x$fourth_quad_cause)*x$nor_weight) )
  
}

correct_rate_first_quad<-unlist(lapply( fin_res$AAC_order_list, weighted_correct_rate_first_quad ))
correct_rate_second_quad<-unlist(lapply( fin_res$AAC_order_list, weighted_correct_rate_second_quad ))
correct_rate_third_quad<-unlist(lapply( fin_res$AAC_order_list, weighted_correct_rate_third_quad ))
correct_rate_fourth_quad<-unlist(lapply( fin_res$AAC_order_list, weighted_correct_rate_fourth_quad ))

# Estimate the weighted standard error for the normal approximation.
weighted_SE_first_quad<-function(x,p_hat){
  
  trial_vec<-(x$truth_cause==x$first_quad_cause)
  se_p_hat<-sqrt(sum((x$nor_weight)^2 * (trial_vec - p_hat)^2) / (sum(x$nor_weight)^2))
  return( se_p_hat )
  
}

# Estimate the weighted standard error for the normal approximation.
weighted_SE_second_quad<-function(x,p_hat){
  
  trial_vec<-(x$truth_cause==x$second_quad_cause)
  se_p_hat<-sqrt(sum((x$nor_weight)^2 * (trial_vec - p_hat)^2) / (sum(x$nor_weight)^2))
  return( se_p_hat )
  
}

# Estimate the weighted standard error for the normal approximation.
weighted_SE_third_quad<-function(x,p_hat){
  
  trial_vec<-(x$truth_cause==x$third_quad_cause)
  se_p_hat<-sqrt(sum((x$nor_weight)^2 * (trial_vec - p_hat)^2) / (sum(x$nor_weight)^2))
  return( se_p_hat )
  
}

# Estimate the weighted standard error for the normal approximation.
weighted_SE_fourth_quad<-function(x,p_hat){
  
  trial_vec<-(x$truth_cause==x$fourth_quad_cause)
  se_p_hat<-sqrt(sum((x$nor_weight)^2 * (trial_vec - p_hat)^2) / (sum(x$nor_weight)^2))
  return( se_p_hat )
  
}


compute_p_val<-function(p_hat,p_0,se_p_hat){
  
  # Calculate the Z statistic.
  z_stat <- (p_hat - p_0) / se_p_hat
  
  # Calculate the p-value for a one-sided test.
  # Upper-tail alternative: p > p0.
  p_value_right <- 1 - pnorm(z_stat)
  
  return(p_value_right)
  
}


# Apply weighted_SE_first_quad to each result table and its corresponding accuracy estimate.
weighted_SE_first_quad <- unlist(Map(weighted_SE_first_quad, 
                                     fin_res$AAC_order_list, 
                                     correct_rate_first_quad))

# Apply weighted_SE_second_quad to each result table and its corresponding accuracy estimate.
weighted_SE_second_quad <- unlist(Map(weighted_SE_second_quad, 
                                      fin_res$AAC_order_list, 
                                      correct_rate_second_quad))

# Apply weighted_SE_third_quad to each result table and its corresponding accuracy estimate.
weighted_SE_third_quad <- unlist(Map(weighted_SE_third_quad, 
                                     fin_res$AAC_order_list, 
                                     correct_rate_third_quad))

# Apply weighted_SE_fourth_quad to each result table and its corresponding accuracy estimate.
weighted_SE_fourth_quad <- unlist(Map(weighted_SE_fourth_quad, 
                                      fin_res$AAC_order_list, 
                                      correct_rate_fourth_quad))


sig_level<-0.05

upper_bound_first_quad<-correct_rate_first_quad+qnorm(1-sig_level/2)*weighted_SE_first_quad
lower_bound_first_quad<-correct_rate_first_quad-qnorm(1-sig_level/2)*weighted_SE_first_quad

upper_bound_second_quad<-correct_rate_second_quad+qnorm(1-sig_level/2)*weighted_SE_second_quad
lower_bound_second_quad<-correct_rate_second_quad-qnorm(1-sig_level/2)*weighted_SE_second_quad

upper_bound_third_quad<-correct_rate_third_quad+qnorm(1-sig_level/2)*weighted_SE_third_quad
lower_bound_third_quad<-correct_rate_third_quad-qnorm(1-sig_level/2)*weighted_SE_third_quad

upper_bound_fourth_quad<-correct_rate_fourth_quad+qnorm(1-sig_level/2)*weighted_SE_fourth_quad
lower_bound_fourth_quad<-correct_rate_fourth_quad-qnorm(1-sig_level/2)*weighted_SE_fourth_quad

res_rec<-data.frame(correct_rate_first_quad=correct_rate_first_quad,lower_bound_first_quad=lower_bound_first_quad,upper_bound_first_quad=upper_bound_first_quad,
                    correct_rate_second_quad=correct_rate_second_quad,lower_bound_second_quad=lower_bound_second_quad,upper_bound_second_quad=upper_bound_second_quad,
                    correct_rate_third_quad=correct_rate_third_quad,lower_bound_third_quad=lower_bound_third_quad,upper_bound_third_quad=upper_bound_third_quad,
                    correct_rate_fourth_quad=correct_rate_fourth_quad,lower_bound_fourth_quad=lower_bound_fourth_quad,upper_bound_fourth_quad=upper_bound_fourth_quad)

library(tidyverse)

# Define the horizontal offset for the quadrant-specific points and error bars.
pdodge <- position_dodge(width = 0.6)

# Loop over the rows of res_rec.
for (i in 1:nrow(res_rec)) {
  
  # Extract the i-th row.
  row_data <- res_rec[i, ]
  
  # Use the row name as the method label.
  method_label <- rownames(res_rec)[i]
  
  # Create a data frame in long format for the current row.
  long_df <- tibble(
    method = rep(method_label, 4),
    quadrant = factor(
      c("(x_1,x_2)", "(-x_1,x_2)", "(-x_1,-x_2)", "(x_1,-x_2)"),
      levels = c("(x_1,x_2)", "(-x_1,x_2)", "(-x_1,-x_2)", "(x_1,-x_2)")
    ),
    accuracy = c(
      row_data$correct_rate_first_quad,
      row_data$correct_rate_second_quad,
      row_data$correct_rate_third_quad,
      row_data$correct_rate_fourth_quad
    ),
    lower = c(
      row_data$lower_bound_first_quad,
      row_data$lower_bound_second_quad,
      row_data$lower_bound_third_quad,
      row_data$lower_bound_fourth_quad
    ),
    upper = c(
      row_data$upper_bound_first_quad,
      row_data$upper_bound_second_quad,
      row_data$upper_bound_third_quad,
      row_data$upper_bound_fourth_quad
    )
  )
  
  # Create the plot.
  p <- ggplot(long_df, aes(x = method, y = accuracy, color = quadrant, group = quadrant)) +
    geom_point(position = pdodge, size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = pdodge, width = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
    coord_cartesian(ylim = c(0.35, 0.85)) +
    labs(
      x = "",  # Remove the x-axis label.
      y = "Accuracy with 95% CI"
    ) +
    scale_color_discrete(
      labels = c(
        expression("(" * x[1] * "," ~ x[2] * ")"),
        expression("(-" * x[1] * "," ~ x[2] * ")"),
        expression("(-" * x[1] * "," ~ -x[2] * ")"),
        expression("(" * x[1] * "," ~ -x[2] * ")")
      )
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
  
  # Save the plot.
  filename <- paste0("plots/Section_4.3_pairwise_semi_", method_label,"_lambda_",lambda,"_a_0.7_fit_tail_strong_stop_adjust_m",".pdf")
  ggsave(filename, p, width = 6, height = 4.5, units = "in")
}

