###### Visualize the table in Section 4.1 and Supplementary Section S.1.1

rm(list=ls())

library(tidyr)
library(ggplot2)

cal_mean<-function(mat_list,round_digit=4){
  
  mat_list_num <- lapply(mat_list, function(df) {
    m <- as.matrix(df)              # Force a single data type.
    storage.mode(m) <- "numeric"    # Convert "1.2" -> 1.2; non-numeric values -> NA.
    m
  })
  
  # Convert the list of matrices to a 3D array: [row, col, simulation].
  array_data <- simplify2array(mat_list_num)
  
  # Elementwise mean.
  mean_matrix <- apply(array_data, c(1, 2), mean, na.rm = TRUE)
  
  # Elementwise standard deviation.
  sd_matrix <- apply(array_data, c(1, 2), sd, na.rm = TRUE)
  
  # Elementwise standard error.
  se_matrix <- sd_matrix / sqrt(dim(array_data)[3])
  
  # Convert to data frames.
  mean_df <- as.data.frame(mean_matrix)
  sd_df   <- as.data.frame(sd_matrix)
  se_df   <- as.data.frame(se_matrix)
  
  # Preserve row and column names.
  rownames(mean_df) <- rownames(mat_list[[1]])
  colnames(mean_df) <- colnames(mat_list[[1]])
  
  rownames(sd_df) <- rownames(mat_list[[1]])
  colnames(sd_df) <- colnames(mat_list[[1]])
  
  rownames(se_df) <- rownames(mat_list[[1]])
  colnames(se_df) <- colnames(mat_list[[1]])
  
  # Round if needed.
  mean_df <- round(mean_df, round_digit)
  sd_df   <- round(sd_df, round_digit)
  
  return(list(mean_df=mean_df,sd_df=sd_df) )
  
}

load("eSCM_code/datasets/CTC_AAC_DAG_lambda_2_p_5_n_1000_w_min_0.04_w_max_0.4_shape_1_M_1_500_fixed_k_tail_para")
#load("eSCM_code/datasets/CTC_AAC_DAG_lambda_2_p_5_n_1000_w_min_0.04_w_max_0.4_shape_3_M_1_500_fixed_k_tail_para")
#load("eSCM_code/datasets/CTC_AAC_DAG_lambda_2_p_5_n_1000_w_min_0.04_w_max_0.4_shape_5_M_1_500_fixed_k_tail_para")
length(fin_res$error_rate_list)
res_rec_p_5_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_5_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_5_CTC)


load("result_record/CTC_AAC_DAG_lambda_2_p_10_n_1000_w_min_0.04_w_max_0.4_shape_1_M_1_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_10_n_1000_w_min_0.04_w_max_0.4_shape_3_M_1_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_10_n_1000_w_min_0.04_w_max_0.4_shape_5_M_1_500_fixed_k_tail_para")
res_rec_p_10_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_10_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_10_CTC)

load("result_record/CTC_AAC_DAG_lambda_2_p_15_n_1000_w_min_0.04_w_max_0.4_shape_1_M_1_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_15_n_1000_w_min_0.04_w_max_0.4_shape_3_M_1_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_15_n_1000_w_min_0.04_w_max_0.4_shape_5_M_1_500_fixed_k_tail_para")
res_rec_p_15_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_15_AAC<-fin_res$error_rate_list$temp_rec_AAC
length(res_rec_p_15_CTC)

load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_1_M_1_250_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_3_M_1_250_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_5_M_1_250_fixed_k_tail_para")
res_rec_p_30_CTC<-fin_res$error_rate_list$temp_rec_CTC
res_rec_p_30_AAC<-fin_res$error_rate_list$temp_rec_AAC

load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_1_M_251_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_3_M_251_500_fixed_k_tail_para")
#load("result_record/CTC_AAC_DAG_lambda_2_p_30_n_1000_w_min_0.04_w_max_0.4_shape_5_M_251_500_fixed_k_tail_para")
temp_CTC<-fin_res$error_rate_list$temp_rec_CTC
temp_AAC<-fin_res$error_rate_list$temp_rec_AAC
res_rec_p_30_CTC<<-c(res_rec_p_30_CTC,temp_CTC)
res_rec_p_30_AAC<<-c(res_rec_p_30_AAC,temp_AAC)

# Example: res_p_15_n_1000_w_min_0.004_w_max_0.008[498:500] is a list of matrices.
mean_CI_p_30_CTC<-t( (cal_mean(res_rec_p_30_CTC))$mean_df )
mean_CI_p_15_CTC<-t((cal_mean(res_rec_p_15_CTC))$mean_df)
mean_CI_p_10_CTC<-t((cal_mean(res_rec_p_10_CTC))$mean_df) 
mean_CI_p_5_CTC<-t((cal_mean(res_rec_p_5_CTC))$mean_df )

mean_CI_p_30_AAC<-t((cal_mean(res_rec_p_30_AAC))$mean_df)
mean_CI_p_15_AAC<-t((cal_mean(res_rec_p_15_AAC))$mean_df)
mean_CI_p_10_AAC<-t((cal_mean(res_rec_p_10_AAC))$mean_df) 
mean_CI_p_5_AAC<-t((cal_mean(res_rec_p_5_AAC))$mean_df )

mean_CI_tab_AAC<-rbind(mean_CI_p_5_AAC[1:3,],mean_CI_p_10_AAC[1:3,],mean_CI_p_15_AAC[1:3,], mean_CI_p_30_AAC[1:3,])
mean_CI_tab_CTC<-rbind(mean_CI_p_5_CTC[1:3,],mean_CI_p_10_CTC[1:3,],mean_CI_p_15_CTC[1:3,], mean_CI_p_30_CTC[1:3,])

mean_CI_tab_AAC<mean_CI_tab_CTC

data.frame(d=rep(c(5,10,15,30),each=3),k=rep(c(50,100,150),4),
           c(mean_CI_p_5_AAC[,1],mean_CI_p_10_AAC[,1], mean_CI_p_15_AAC[,1],mean_CI_p_30_AAC[,1]),
           c(mean_CI_p_5_CTC[,1],mean_CI_p_10_CTC[,1], mean_CI_p_15_CTC[,1],mean_CI_p_30_CTC[,1]),
           c(mean_CI_p_5_AAC[,2],mean_CI_p_10_AAC[,2], mean_CI_p_15_AAC[,2],mean_CI_p_30_AAC[,2]),
           c(mean_CI_p_5_CTC[,2],mean_CI_p_10_CTC[,2], mean_CI_p_15_CTC[,2],mean_CI_p_30_CTC[,2]),
           c(mean_CI_p_5_AAC[,3],mean_CI_p_10_AAC[,3], mean_CI_p_15_AAC[,3],mean_CI_p_30_AAC[,3]),
           c(mean_CI_p_5_CTC[,3],mean_CI_p_10_CTC[,3], mean_CI_p_15_CTC[,3],mean_CI_p_30_CTC[,3]),
           c(mean_CI_p_5_AAC[,4],mean_CI_p_10_AAC[,4], mean_CI_p_15_AAC[,4],mean_CI_p_30_AAC[,4]),
           c(mean_CI_p_5_CTC[,4],mean_CI_p_10_CTC[,4], mean_CI_p_15_CTC[,4],mean_CI_p_30_CTC[,4])
)
