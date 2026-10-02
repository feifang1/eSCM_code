rm(list=ls())

library(tidyr)
library(ggplot2)

cal_mean<-function(mat_list,round_digit=4){
  
  mat_list_num <- lapply(mat_list, function(df) {
    m <- as.matrix(df)              # Convert each data frame to a matrix with a single data type.
    storage.mode(m) <- "numeric"    # Convert numeric strings to numbers; nonnumeric strings become NA.
    m
  })
  
  # Combine the matrices into a 3D array: [row, column, simulation].
  array_data <- simplify2array(mat_list_num)
  
  # Calculate the mean of each entry across simulations.
  mean_matrix <- apply(array_data, c(1, 2), mean, na.rm = TRUE)
  
  # Calculate the standard deviation of each entry across simulations.
  sd_matrix <- apply(array_data, c(1, 2), sd, na.rm = TRUE)
  
  # Calculate the standard error using the total number of simulations.
  se_matrix <- sd_matrix / sqrt(dim(array_data)[3])
  
  # Convert the summary matrices to data frames.
  mean_df <- as.data.frame(mean_matrix)
  sd_df   <- as.data.frame(sd_matrix)
  se_df   <- as.data.frame(se_matrix)
  
  # Preserve the row and column names from the first input data frame.
  rownames(mean_df) <- rownames(mat_list[[1]])
  colnames(mean_df) <- colnames(mat_list[[1]])
  
  rownames(sd_df) <- rownames(mat_list[[1]])
  colnames(sd_df) <- colnames(mat_list[[1]])
  
  rownames(se_df) <- rownames(mat_list[[1]])
  colnames(se_df) <- colnames(mat_list[[1]])
  
  # Round the means and standard deviations to round_digit decimal places.
  mean_df <- round(mean_df, round_digit)
  sd_df   <- round(sd_df, round_digit)
  
  return(list(mean_df=mean_df,sd_df=sd_df) )
  
}

# To load results for a different shape parameter, uncomment the corresponding load() line and comment out the currently active line.

#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_5_n_1000_shape_1_M_1_500_alpha_2_m_sqrt_n")
load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_5_n_1000_shape_3_M_1_500_alpha_2_m_sqrt_n")
#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_5_n_1000_shape_5_M_1_500_alpha_2_m_sqrt_n")
length(fin_res$error_rate_list)
res_rec_p_5_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_5_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_5_CTC)


#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_10_n_1000_shape_1_M_1_500_alpha_2_m_sqrt_n")
load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_10_n_1000_shape_3_M_1_500_alpha_2_m_sqrt_n")
#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_10_n_1000_shape_5_M_1_500_alpha_2_m_sqrt_n")
res_rec_p_10_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_10_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_10_CTC)

#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_15_n_1000_shape_1_M_1_500_alpha_2_m_sqrt_n")
load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_15_n_1000_shape_3_M_1_500_alpha_2_m_sqrt_n")
#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_15_n_1000_shape_5_M_1_500_alpha_2_m_sqrt_n")
res_rec_p_15_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_15_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_15_CTC)

#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_30_n_1000_shape_1_M_1_500_alpha_2_m_sqrt_n")
load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_30_n_1000_shape_3_M_1_500_alpha_2_m_sqrt_n")
#load("datasets/Section_4.1/CTC_AAC_DAG_lambda_4_p_30_n_1000_shape_5_M_1_500_alpha_2_m_sqrt_n")
res_rec_p_30_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_30_AAC<-fin_res$error_rate_list$temp_rec_AAC

mean_CI_p_30_CTC<-t( (cal_mean(res_rec_p_30_CTC))$mean_df )
mean_CI_p_15_CTC<-t((cal_mean(res_rec_p_15_CTC))$mean_df)
mean_CI_p_10_CTC<-t((cal_mean(res_rec_p_10_CTC))$mean_df) 
mean_CI_p_5_CTC<-t((cal_mean(res_rec_p_5_CTC))$mean_df )

mean_CI_p_30_AAC<-t((cal_mean(res_rec_p_30_AAC))$mean_df)
mean_CI_p_15_AAC<-t((cal_mean(res_rec_p_15_AAC))$mean_df)
mean_CI_p_10_AAC<-t((cal_mean(res_rec_p_10_AAC))$mean_df) 
mean_CI_p_5_AAC<-t((cal_mean(res_rec_p_5_AAC))$mean_df )

mean_CI_tab_AAC<-rbind(mean_CI_p_5_AAC[1,],mean_CI_p_10_AAC[1,],mean_CI_p_15_AAC[1,], mean_CI_p_30_AAC[1,])
mean_CI_tab_CTC<-rbind(mean_CI_p_5_CTC[1,],mean_CI_p_10_CTC[1,],mean_CI_p_15_CTC[1,], mean_CI_p_30_CTC[1,])

mean_CI_tab_AAC<mean_CI_tab_CTC

data.frame(d=rep(c(5,10,15,30),each=1),
           c(mean_CI_p_5_AAC[,1],mean_CI_p_10_AAC[,1], mean_CI_p_15_AAC[,1],mean_CI_p_30_AAC[,1]),
           c(mean_CI_p_5_CTC[,1],mean_CI_p_10_CTC[,1], mean_CI_p_15_CTC[,1],mean_CI_p_30_CTC[,1]),
           c(mean_CI_p_5_AAC[,2],mean_CI_p_10_AAC[,2], mean_CI_p_15_AAC[,2],mean_CI_p_30_AAC[,2]),
           c(mean_CI_p_5_CTC[,2],mean_CI_p_10_CTC[,2], mean_CI_p_15_CTC[,2],mean_CI_p_30_CTC[,2]),
           c(mean_CI_p_5_AAC[,3],mean_CI_p_10_AAC[,3], mean_CI_p_15_AAC[,3],mean_CI_p_30_AAC[,3]),
           c(mean_CI_p_5_CTC[,3],mean_CI_p_10_CTC[,3], mean_CI_p_15_CTC[,3],mean_CI_p_30_CTC[,3]),
           c(mean_CI_p_5_AAC[,4],mean_CI_p_10_AAC[,4], mean_CI_p_15_AAC[,4],mean_CI_p_30_AAC[,4]),
           c(mean_CI_p_5_CTC[,4],mean_CI_p_10_CTC[,4], mean_CI_p_15_CTC[,4],mean_CI_p_30_CTC[,4])
)


