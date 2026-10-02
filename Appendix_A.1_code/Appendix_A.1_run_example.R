# Compare CTC and AAC across k = round(n^a), a = 0.20, 0.25, ..., 0.80.
# Full local version of Appendix_A.1_Simulation_choose_k.R.
# Edit the parameter block below, then source("run_example.R") in RStudio
# or run Rscript run_example.R. No command-line parameters are required.
# The original simulation calculations and result object are retained.
# As in the original script, rm(list=ls()) clears the current workspace.

rm(list=ls())

# ---- Choose one parameter combination and local computing settings ----
n <- 500               # Sample size: 500, 1000, or 2000 in params_v7.txt.
p <- 5                 # DAG size: 5, 10, 15, or 30.
shape <- 3             # Pareto noise shape: 1, 3, or 5.
lambda <- 4            # Penalty parameter, fixed at 4 in params_v7.txt.
quantile_val <- 1/2     # Exponent for m = round(n^quantile_val).
M_lower <- 1           # First Monte Carlo replication index.
M_upper <- 2           # Last index; use 50 for the paper's runs.
no_cores <- 1          # Number of local R workers.
output_folder <- file.path("datasets", "example")
# ---- Run the full simulation below with the values selected above ----

# Locate this file when sourced in RStudio or run with Rscript.
# For interactive execution of selected lines, set the working directory to
# the folder containing run_example.R before running the full file.
script_dir <- local({
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(x) is.character(x) && length(x) == 1L, source_files)
  if (length(source_files)) {
    dirname(normalizePath(source_files[[length(source_files)]], mustWork = TRUE))
  } else {
    script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(script_arg)) {
      dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
    } else {
      normalizePath(getwd(), mustWork = TRUE)
    }
  }
})
output_dir <- file.path(script_dir, output_folder)

function_file <- Sys.getenv("ESCM_FUNCTION_FILE", unset = file.path(script_dir, "function_class.R"))
if (!file.exists(function_file)) {
  stop("Missing helper file: ", function_file,
       ". Include the original eSCM_code/function_class.R, or set ESCM_FUNCTION_FILE.")
}


# Use standard R library locations; configure R_LIBS_USER externally if needed.
required_packages <- c("graph", "pcalg", "EnvStats", "dplyr", "Matrix", "igraph",
                       "foreach", "doParallel", "expm", "mev")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                            logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required R packages: ", paste(missing_packages, collapse = ", "))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#BiocManager::install("graph", update = FALSE, ask = FALSE)
#BiocManager::install("RBGL", update = FALSE, ask = FALSE)

library(RBGL)
library(graph)
library(pcalg)
library(EnvStats)
source(function_file)

library(dplyr)
library(Matrix)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(expm)
library(mev)


# Validate the selected combination before starting worker processes.
chosen_parameters <- list(n, p, shape, lambda, quantile_val, M_lower, M_upper, no_cores)
if (any(lengths(chosen_parameters) != 1L) ||
    !all(vapply(chosen_parameters, is.numeric, logical(1)))) {
  stop("Specify one numeric value for each parameter, replication index, and worker count.")
}
numeric_parameters <- unlist(chosen_parameters)
if (any(!is.finite(numeric_parameters)) || any(numeric_parameters <= 0) ||
    n != floor(n) || p != floor(p) || n < 2 || p < 2 || quantile_val >= 1) {
  stop("Use positive finite values, integer n and p >= 2, and 0 < quantile_val < 1.")
}
if (M_lower != floor(M_lower) || M_upper != floor(M_upper) || M_upper < M_lower ||
    no_cores != floor(no_cores) || no_cores > .Machine$integer.max) {
  stop("Use integer replication indices with M_upper >= M_lower and a positive integer worker count.")
}
M <- M_upper - M_lower + 1

cat("Running simulation with:\n")
cat("n =", n, "\n")
cat("p =", p, "\n")
cat("shape =", shape, "\n")
cat("lambda =", lambda, "\n")
cat("quantile_val =", quantile_val, "\n")


# nuisance parameters in D_k(s,t)
L<-c(3)
d_vec<-c(1)
r<-c(0.5)

# LCM: SL0 in the paper; MLC: ML0 in the paper 
# ALCM: SL1 in the paper; AMLC: ML1 in the paper
MOD_vec<-c("LCM","MLC","ALCM","AMLC")

#n<-1000
# Extremal subsample sizes used to estimate causal order.
nv_vec<-seq(0.2,0.8,by=0.05)
k_vec<-round(n^nv_vec)
# Subsample size m used to fit the marginal upper tails.
#k_tail_vec<-round(n^nv_vec)
k_tail_vec<-rep(round(n^quantile_val),length(nv_vec))

# p: DAG size d in the paper, d \in {5,10,15,30}
#p<-5
# graph degree 
deg<-3

# Replication indices and worker count are selected at the top of this file.

# parameters for SL0
w_min<-0.04
w_max<-0.4

log_normal_par_1<--1.528578
log_normal_par_2<-0.4681941

# lambda is selected at the top of this file.

error_rate_vec<-matrix(0,nrow=(M_upper-M_lower+1),ncol=length(k_vec))

# Generate name combinations
name_vec<-rep(0,(length(MOD_vec)*length(d_vec)))

index <- 1
for (MOD_ind in 1:length(MOD_vec)) {
  for (d_ind in 1:length(d_vec)) {
    name <- paste0("MOD_", MOD_vec[MOD_ind], "_d_", d_vec[d_ind])
    name_vec[index] <- name
    index <- index + 1
  }
}

# Data generation process parameters
# shape: alpha_0 in the paper, alpha_0 =3 in the paper, alpha_0 \in {1,5} in Appendix
#shape<-3
location<-1

# marginal transformation parameter
alpha<-2

start_time <- Sys.time()  # Start time

# Use the chosen local worker count, without a scheduler or shell launcher.
no_cores <- as.integer(no_cores)
cl<-makeCluster(no_cores,outfile="")
registerDoParallel(cl)

res_list <- foreach(m = M_lower:M_upper,
                    .combine = combine_custom, 
                    .packages = c("igraph", "Matrix", "pcalg", "expm", "EnvStats","mev")) %dopar% { 
                      
                      print( paste0("iteration m: ", m) )
                      
                      set.seed(m)
                      dag<-randDAG(n=p,d=deg,method="er", DAG=T, weighted=F)
                      adj_mat <- as(dag, "matrix")
                      
                      # Initialize weight matrix with zeros
                      
                      temp_rec<-matrix(0,nrow=(length(MOD_vec)*length(d_vec)),ncol=length(k_vec) )
                      temp_rec_CTC<-temp_rec_AAC<-data.frame(temp_rec)
                      rownames(temp_rec_CTC)<-rownames(temp_rec_AAC)<-name_vec
                      colnames(temp_rec_CTC)<-colnames(temp_rec_AAC)<-k_vec
                      
                      for (MOD_ind in 1:length(MOD_vec) ) {
                        
                        if(MOD_vec[MOD_ind] %in% c("LCM","MLC")){
                          
                          # Assign weights only to existing edges
                          adj_mat[adj_mat == 1] <- runif(sum(adj_mat == 1), min = w_min, max = w_max)
                          adj_mat<-Matrix(adj_mat)
                          
                          A<-matrix(0,p,p)
                          
                          ############ generate adjacency matrix with DAG ###################
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
                        
                        # Generate observations from the current model and weighted DAG.
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
                          
                          for(k_ind in 1:length(k_vec) ){
                            
                            # Fit a generalized Pareto distribution above each marginal threshold.
                            tau0 <- 1 - k_tail_vec[k_ind]/dim(X_mat)[1]
                            
                            u <- sapply(data.frame(X_mat), function(col) {
                              as.numeric(quantile(col, probs = tau0, type = 8))
                            })
                            
                            # Lower thresholds equal to the maximum to obtain exceedances.
                            # If their variance is zero, lower the threshold again.
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
                              # Fit each column separately at its corresponding threshold.
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
                              
                              # avoid exact 1 due to floating-point issues
                              if (any(Fval >= 1)) {
                                
                                ind <- which(Fval >= 1)
                                
                                if(sum(Fval >= 1)<length(Fval)){
                                  
                                  Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
                                }else{
                                  # All computed CDF values are at least 1.
                                  Fval[ind]<-0.99999
                                }
                                
                              }
                              
                              F_x_u[[j]] <- Fval
                              MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
                            }
                            
                            for(j in 1:dim(Z)[2]){
                              
                              Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
                              
                            }
                            
                            temp_CTC<-ease(dat=Z,k=k_vec[k_ind],both_tails = F)
                            temp_rec_CTC[temp_ind_CTC,k_ind]<-order_error_measure(dag=as(dag, "matrix"),ab_dir=temp_CTC)
                            
                            temp_AAC<-AAC_ease_order(dat=Z, a_b_val=F, k = k_vec[k_ind], lambda=lambda, L_r=L, dis=d_vec[d_ind], r=r )
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
  k_vec=k_vec,
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


# Save the original fin_res object and filename convention in a portable directory.
# output_folder is selected at the top and resolved relative to this file.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
file_name <- file.path(output_dir, paste0("CTC_AAC_DAG_lambda_",lambda,"_p_",p,"_n_",n,"_shape_",shape, "_M_",M_lower,"_",M_upper,"_m_sqrt_n"))
save(fin_res, file=file_name)
cat("Saved results to:", normalizePath(file_name), "\n")
