# Section 4.1 simulation study
# Each job reads one (p, shape) combination from param_simulation.txt.
# p is the DAG size; shape is the Pareto shape parameter alpha_0.
# This script uses lambda = 4, 500 replications, and AAC exponent 0.7.
# Before running, adjust the source and output paths below.

rm(list=ls())
# Optional: uncomment and adjust if you use a custom R package library.
# .libPaths("/path/to/your/R/library")
library(graph)
library(pcalg)
library(EnvStats)
# Replace this path with the location of function_class.R.
source("/path/to/your/project/function_class.R")


library(dplyr)
library(Matrix)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(expm)
library(mev)

args <- commandArgs(trailingOnly = TRUE)

p <- as.numeric(args[1])
shape <- as.numeric(args[2])

cat("p =", p, "\n")
cat("shape =", shape, "\n")


# Auxiliary parameters used in the AAC calculation.
L<-c(3)
d_vec<-c(1)
r<-c(0.5)

# Model names in the code and in the paper: 
# LCM = SL0; MLC = ML0.
# ALCM = SL1; AMLC = ML1.
MOD_vec<-c("LCM","MLC","ALCM","AMLC")

n<-1000

# Extremal subsample sizes used by AAC and CTC.
k_vec_AAC<-round(n^0.7)
k_vec_CTC<-floor(n^0.4)
# Upper-tail subsample size used to fit the generalized Pareto distribution.
k_tail_vec<-round(n^0.5)

# p: DAG size d in the paper, d \in {5,10,15,30}
# Average node degree for random DAG generation.
deg<-3

# Replication indices; the paper uses replications 1 through 500.
M_lower<-1
M_upper<-500
M<-500

# Uniform edge-weight bounds for SL0 and ML0.
w_min<-0.04
w_max<-0.4

# Lognormal meanlog and sdlog for SL1 and ML1.
log_normal_par_1<--1.528578
log_normal_par_2<-0.4681941

# parameter lambda in g_n(s,t)
lambda<-4

error_rate_vec<-matrix(0,nrow=(M_upper-M_lower+1),ncol=length(k_vec_AAC))

# Create row labels for each model and auxiliary distance parameter.
name_vec<-rep(0,(length(MOD_vec)*length(d_vec)))

index <- 1
for (MOD_ind in 1:length(MOD_vec)) {
  for (d_ind in 1:length(d_vec)) {
    name <- paste0("MOD_", MOD_vec[MOD_ind], "_d_", d_vec[d_ind])
    name_vec[index] <- name
    index <- index + 1
  }
}

# Noise distribution parameters.
# shape is alpha_0: 3 in the main study and 1 or 5 in the Appendix.
#shape<-5
location<-1

# Tail index used for the marginal transformation.
alpha<-2

start_time <- Sys.time()  # Start time

# Use the CPUs allocated by Slurm; default to two cores outside Slurm.
no_cores<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))
cl<-makeCluster(no_cores,outfile="")
registerDoParallel(cl)

res_list <- foreach(m = M_lower:M_upper,
                    .combine = combine_custom, 
                    .packages = c("igraph", "Matrix", "pcalg", "expm", "EnvStats","mev")) %dopar% { 
                      
                      print( paste0("iteration m: ", m) )
                      
                      set.seed(m)
                      dag<-randDAG(n=p,d=deg,method="er", DAG=T, weighted=F)
                      adj_mat <- as(dag, "matrix")
                      
                      # Initialize result tables for this replication.
                      
                      temp_rec<-matrix(0,nrow=(length(MOD_vec)*length(d_vec)),ncol=length(k_vec_AAC) )
                      temp_rec_CTC<-temp_rec_AAC<-data.frame(temp_rec)
                      rownames(temp_rec_CTC)<-rownames(temp_rec_AAC)<-name_vec
                      colnames(temp_rec_AAC)<-k_vec_AAC
                      colnames(temp_rec_CTC)<-k_vec_CTC
                      
                      for (MOD_ind in 1:length(MOD_vec) ) {
                        
                        if(MOD_vec[MOD_ind] %in% c("LCM","MLC")){
                          
                          # Assign weights only to existing edges
                          adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1), min = w_min, max = w_max)
                          adj_mat<-Matrix(adj_mat)
                          
                          A<-matrix(0,p,p)
                          # Sum powers of the weighted adjacency matrix to obtain path weights.
                          for( i in 1:p ){
                            temp_mat<-as.matrix(adj_mat) %^% i
                            A <- A+temp_mat
                            if(sum(temp_mat)==0){
                              break 
                            }
                            
                          }  
                          
                          
                          A<-diag(p)+A
                          A<-t(A)
                          A_weighted <- A
                          
                          
                        }else if(MOD_vec[MOD_ind] %in% c("ALCM","AMLC")){
                          
                          # Assign weights only to existing edges
                          adj_mat[adj_mat == 1] <- rlnorm(sum(adj_mat == 1), meanlog = log_normal_par_1, sdlog = log_normal_par_2)
                          adj_mat<-Matrix(adj_mat)
                          
                          A<-matrix(0,p,p)
                          for( i in 1:p ){
                            temp_mat<-as.matrix(adj_mat) %^% i
                            A <- A+temp_mat
                            if(sum(temp_mat)==0){
                              break 
                            }
                            
                          }  
                          
                          A<-diag(p)+A
                          A<-t(A)
                          A_weighted <- A
                          
                        }
                        
                        # Generate data from the random DAG for the selected model.
                        if(MOD_vec[MOD_ind]%in% c("LCM","ALCM") ){
                          
                          eta_mat<-matrix(0,nrow=n, ncol=p)
                          
                          for(i in 1:p){
                            
                            eta_mat[,i]<-rpareto(n,location, shape)
                            
                          }
                          
                          X_mat<- t( A_weighted %*% t(eta_mat) )
                          R<-apply(X_mat,2,rank_top_down)
                          Z<-(n/R)^(1/alpha)
                          
                          
                        }else if(MOD_vec[MOD_ind]%in% c("MLC","AMLC")){
                          
                          eta_mat<-matrix(0,nrow=n, ncol=p)
                          
                          for(i in 1:p){
                            
                            eta_mat[,i]<-rpareto(n,location, shape)
                            
                          }
                          
                          X_mat<-plus_max_linear_mat_mat(A,eta_mat)
                          R<-apply(X_mat,2,rank_top_down)
                          Z<-(n/R)^(1/alpha)
                          
                        }
                        
                        
                        for(d_ind in 1:length(d_vec)){
                          
                          name <- paste0("MOD_", MOD_vec[MOD_ind], "_d_", d_vec[d_ind])
                          temp_ind_CTC<-which(rownames(temp_rec_CTC)==name)
                          temp_ind_AAC<-which(rownames(temp_rec_AAC)==name)
                          
                          for(k_ind in 1:length(k_vec_AAC) ){
                            
                            # Fit a generalized Pareto distribution to each marginal upper tail.
                            tau0 <- 1 - k_tail_vec[k_ind]/dim(X_mat)[1]
                            
                            u <- sapply(data.frame(X_mat), function(col) {
                              as.numeric(quantile(col, probs = tau0, type = 8))
                            })
                            
                            # If a threshold equals the maximum, lower it to obtain exceedances.
                            # Lower it again if the exceedances have zero variance.
                            ind_equal_max<-which(u==apply(data.frame(X_mat),2,max))
                            if(length(ind_equal_max)>0){
                              
                              for(ind_equal_max_index in ind_equal_max){
                                
                                data_exc_second_large<-which(X_mat[,ind_equal_max_index]<u[ind_equal_max_index])
                                u[ind_equal_max_index]<-max(X_mat[data_exc_second_large,ind_equal_max_index])
                                if(var(X_mat[,ind_equal_max_index][which(X_mat[,ind_equal_max_index]>u[ind_equal_max_index])])==0 ){
                                  data_exc_third_large<-which(X_mat[,ind_equal_max_index]<u[ind_equal_max_index])
                                  u[ind_equal_max_index]<-max(X_mat[data_exc_third_large,ind_equal_max_index])
                                }
                                
                                
                              }
                              
                              
                            }
                            
                            
                            fits <- mapply(
                              function(col, thr) {
                                fit.gpd(col, threshold = thr, method = "Grimshaw")
                              },
                              # Fit each column using its corresponding threshold.
                              as.data.frame(X_mat),
                              u,
                              SIMPLIFY = FALSE
                            )
                            
                            params <- lapply(fits, coef)
                            sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
                            gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
                            
                            exceed_idx <- mapply(
                              function(col, thr) which(col > thr),
                              as.data.frame(X_mat),
                              u,
                              SIMPLIFY = FALSE
                            )
                            
                            
                            F_x_u <- vector("list", p)
                            MT_upper_u <- vector("list", p)
                            
                            for (j in seq_len(p)) {
                              
                              idx <- exceed_idx[[j]]   # indices for column j
                              xexc <- X_mat[idx, j]     # exceedances
                              
                              z <- (xexc - u[j]) / sigma_u_hat[j]
                              
                              Fval <- 1 - (k_tail_vec[k_ind] / n) * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
                              
                              # Avoid CDF values of exactly one due to floating-point rounding.
                              if (any(Fval >= 1)) {
                                
                                ind <- which(Fval >= 1)
                                
                                if(sum(Fval >= 1)<length(Fval)){
                                  
                                  Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
                                }else{
                                  # If every CDF value is at least one, apply the same cap.
                                  Fval[ind]<-0.99999
                                }
                                
                              }
                              
                              F_x_u[[j]] <- Fval
                              MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
                            }
                            
                            for(j in 1:dim(Z)[2]){
                              
                              Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
                              
                            }
                            
                            temp_CTC<-ease(dat=Z,k=k_vec_CTC[k_ind],both_tails = F)
                            temp_rec_CTC[temp_ind_CTC,k_ind]<-order_error_measure(dag=as(dag, "matrix"),ab_dir=temp_CTC)
                            
                            temp_AAC<-AAC_ease_order(dat=Z, a_b_val=F, k = k_vec_AAC[k_ind], lambda=lambda, L_r=L, dis=d_vec[d_ind], r=r )
                            temp_rec_AAC[temp_ind_AAC,k_ind]<-order_error_measure(dag=as(dag, "matrix"),ab_dir=temp_AAC)
                            
                          }
                          
                        }
                        
                        
                      }
                      
                      return(list(temp_rec_CTC=temp_rec_CTC,temp_rec_AAC=temp_rec_AAC)  )
                      
                    }

stopCluster(cl)

end_time <- Sys.time()  # End time
elapsed_time <- end_time - start_time
elapsed_time

fin_res<-list(
  error_rate_list=res_list,
  k_vec_AAC=k_vec_AAC,
  k_vec_CTC=k_vec_CTC,
  k_tail_vec=k_tail_vec,
  n=n,
  alpha=alpha,
  alpha_0=shape,
  beta_satisfy=paste0("unif[",w_min,",",w_max,"]"),
  beta_NOT_satisfy=paste0("unif[",log_normal_par_1,",",log_normal_par_2,"]"),
  L=L,
  d_vec=d_vec,
  lambda=lambda,
  r=r, 
  M=M,
  p=p
)


# Replace this directory with your output directory; create it before running.
file_name<-paste0("/path/to/your/output/CTC_AAC_DAG_lambda_",lambda,"_p_",p,"_n_",n,"_shape_",shape, "_M_",M_lower,"_",M_upper,"_alpha_",alpha,"_m_sqrt_n")
save(fin_res, file=file_name)








