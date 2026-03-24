# Purpose: Apply the EASE algorithm with AAC and calculate causal direction accuracy across 94 datasets.

# Data source:
# J. M. Mooij, J. Peters, D. Janzing, J. Zscheischler, B. Schoelkopf:
# "Distinguishing cause from effect using observational data: methods and benchmarks"
# https://webdav.tuebingen.mpg.de/cause-effect/

rm(list=ls())
library(graph)
library(pcalg)
library(EnvStats)
library(igraph)
library(ismev)
library(dplyr)
library(Matrix)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(expm)
library(mev)
library(POT)
library(eva)
source("eSCM_code/function_class.R")

opt_k_est<-function(Z){
  
  B<-Z[,1]+Z[,2]
  cutoff<-fit_power_law(B)$xmin
  k<-sum(B>= cutoff)
  return(k)
  
}

########### Add functions specific to the pairwise datasets ##############
AAC_val_opt<-function(Z,a_b_val=F,k,lambda,L_r=3,dis,r=1/2){
  
  B<-Z[,1]+Z[,2]
  B<-data.frame(ind=1:length(B),B)
  sorted_B<-B[order(-B$B), ]
  first_k_Z<-Z[sorted_B$ind[1:k],]
  first_k_large_stat<-sorted_B$B[1:k]
  W<-(first_k_Z[,1]/first_k_large_stat)
  
  if(L_r==1){
    
    L_r_vec<-rep(1,k)
    
  }else if(L_r==2){
    
    L_r_vec<-first_k_large_stat
    
  }else if(L_r==3){
    
    L_r_vec<-first_k_large_stat/first_k_large_stat[length(first_k_large_stat)]*log(first_k_large_stat/first_k_large_stat[length(first_k_large_stat)])
    
  }else if(L_r==4){
    
    L_r_vec<-first_k_large_stat/first_k_large_stat[length(first_k_large_stat)]
    
  }
  
  if(dis==1){
    
    objective_function <- function(x, first_k_Z, first_k_large_stat, k, W,lambda,r) {
      a<-x[1]
      b<-x[2]
      
      # Input validation to avoid division by zero or negative values.
      if (a <= 0 || b <= 0) {
        stop("Both parameters (a and b) must be positive.")
      }
      
      # Compute component-wise differences.
      part_1 <- data.frame(
        comp_1 = L_r_vec*(W-b),
        comp_2 = L_r_vec*(a-W)
      )
      
      # Compute the element-wise maximum using pmax (vectorized and faster).
      d_x_y_C <- pmax(part_1$comp_1, part_1$comp_2, 0)
      
      D_minus_H<-sum(d_x_y_C)/k
      
      # Compute the final objective value.
      return((b - a) + lambda*(k^(r))* D_minus_H)
      
      
    }
    
    
    A <- matrix(c(1, 0,   # Ensure a > 0.
                  0, -1,  # Ensure b > 0.
                  -1, 1), # Optionally enforce b >= a.
                byrow = TRUE, ncol = 2)
    
    rhs_scale<-c(0,-1,0)
    
    # Choose an initial feasible point strictly inside the region.
    start_point <- c(0.5, 0.6)  # Ensure this satisfies the constraints.
    
    # Run constrained optimization.
    result <- constrOptim(
      theta = start_point,
      f = function(x) objective_function(x, first_k_Z, first_k_large_stat, k, W,lambda,r),
      grad = NULL,      # No gradient provided; let constrOptim() approximate it.
      ui = A,           # Constraint matrix.
      ci =rhs_scale,    # Constraint vector.
      method = "Nelder-Mead"
    )
    
    a_opt<-result$par[1]
    b_opt<-result$par[2]
    
    
  }else if(dis==2){
    
    
    objective_function <- function(x, first_k_Z, first_k_large_stat, k, W,lambda,r) {
      a<-x[1]
      b<-x[2]
      
      # Input validation to avoid division by zero or negative values.
      if (a <= 0 || b <= 0) {
        stop("Both parameters (a and b) must be positive.")
      }
      
      # Compute component-wise differences.
      part_1 <- data.frame(
        comp_1 = L_r_vec*(1-(1-W)/(1-b)),
        comp_2 = L_r_vec*(1-W/a)
      )
      
      # Compute the element-wise maximum using pmax (vectorized and faster).
      d_x_y_C <- pmax(part_1$comp_1, part_1$comp_2, 0)
      
      D_minus_H<-sum(d_x_y_C)/k
      
      #a_b_stat_mat[i,j]<-(b-a[i])+lambda*sqrt(k)*D_minus_H
      
      # Compute D_minus_H.
      # D_minus_H <- sum(d_x_y_C / first_k_large_stat[k] * log(first_k_large_stat / first_k_large_stat[k])) / k
      
      # Compute the final objective value.
      return((b - a) + lambda*(k^(r))* D_minus_H)
      
    }
    
    
    A <- matrix(c(1, 0,   # Ensure a > 0.
                  0, -1,  # Ensure b > 0.
                  -1, 1), # Optionally enforce b >= a.
                byrow = TRUE, ncol = 2)
    
    rhs_scale<-c(0,-1,0)
    
    # Choose an initial feasible point strictly inside the region.
    start_point <- c(0.5, 0.6)  # Ensure this satisfies the constraints.
    
    # Run constrained optimization.
    result <- constrOptim(
      theta = start_point,
      f = function(x) objective_function(x, first_k_Z, first_k_large_stat, k, W,lambda,r),
      grad = NULL,      # No gradient provided; let constrOptim() approximate it.
      ui = A,           # Constraint matrix.
      ci =rhs_scale,    # Constraint vector.
      method = "Nelder-Mead"
    )
    
    a_opt<-result$par[1]
    b_opt<-result$par[2]
    
  }
  
  if(a_b_val==F){
    
    return(1-(a_opt+b_opt))
    
  }else{
    return(c(a=a_opt,b=b_opt,stat=(1-(a_opt+b_opt))))
  }
  
}

find_k_for_fit_GPD <- function(x, n, power_order, C_vec_k_tail_fit) {
  
  k_fallback <- round((n)^(power_order) * C_vec_k_tail_fit)
  
  k <- tryCatch({
    tab <- gpdSeqTests(x, nextremes = seq(15, 90, by = 1))
    
    tab_ind <- max(which(tab$StrongStop >= 0.05))
    
    # Handle the case where which() returns integer(0).
    if (!is.finite(tab_ind) || length(tab_ind) == 0) stop("No StrongStop >= 0.05")
    
    tab$num.above[tab_ind]
  }, error = function(e) {
    k_fallback
  })
  
  return(k)
}

weight_tab<-read.table("eSCM_code/datasets/pairs/pairmeta.txt", header = F)
colnames(weight_tab)[6]<-"weight"

# Indices of the removed datasets.
remov_ind<-c(33:37,47,52:55,70,71,105,107)

# In the paper, we use lambda \in {0.5,1,2,3}.
lambda<-3
r<-1/2

# Constant for the extremal subsample size k.
C_vec<-1.2

# Nuisance parameter for the subsample size m.
C_vec_k_tail_fit<-0.6

# Nuisance parameters for the objective function D_k(s,t).
L_vec<-c(3)
d_vec<-c(1)


# Use both the ab method and the EASE method.
AAC_order<-matrix(0,nrow=108-length(remov_ind),ncol=10 )
AAC_order<-as.data.frame(AAC_order)

# Generate column names.
col_names <- c("truth_cause", "truth_effect", 
               "first_quad_cause","first_quad_effect",
               "second_quad_cause","second_quad_effect", 
               "third_quad_cause","third_quad_effect",
               "fourth_quad_cause","fourth_quad_effect")

# Assign the column names to AAC_order.
colnames(AAC_order) <- col_names

AAC_val<-matrix(0,nrow=108-length(remov_ind),ncol=8 )
AAC_val<-as.data.frame(AAC_val)

col_names <- c("first_quad_cause","first_quad_effect",
               "second_quad_cause","second_quad_effect", 
               "third_quad_cause","third_quad_effect",
               "fourth_quad_cause","fourth_quad_effect")

colnames(AAC_val)<-col_names

data_set_ind_kept<-setdiff(1:108,remov_ind)

# Encode the underlying truth.
trueth<-matrix(0,nrow=length(data_set_ind_kept),ncol=2)

trueth[,1]<-1
trueth[,2]<-2

data_set_ind_kept

rownames(trueth)<-data_set_ind_kept
ind_reverse_order<-which(rownames(trueth) %in% as.character(c(48,49:53,55,56:63,68:69,73,77,79,80,84,89,90,92,99,106,108)))
trueth[ind_reverse_order,1]<-2
trueth[ind_reverse_order,2]<-1

AAC_order[,1:2]<-trueth

num_obj_func<-length(C_vec)*length(L_vec)*length(d_vec)

AAC_order_list <- replicate(num_obj_func, AAC_order, simplify = FALSE)
AAC_val_list <- replicate(num_obj_func, AAC_val, simplify = FALSE)

# Create all parameter combinations.
param_grid <- expand.grid(C = C_vec, L = L_vec, d = d_vec)

# Generate list names.
list_names <- with(param_grid, paste0("C_", C, "_L_", L, "_d_", d))

# Assign names to the lists.
names(AAC_order_list) <- names(AAC_val_list)<-list_names

power_order<-1/2

size_rec<-rep(0,length(setdiff(1:108,remov_ind)))
k_pos_rec<-rep(0,length(setdiff(1:108,remov_ind)))

# Marginal transformation parameter.
alpha<-2

size_graph<-rep(0,length(setdiff(1:108,remov_ind)))

for(i in 1:length(setdiff(1:108,remov_ind)) ){
  
  print(paste0("dataset ",i))
  
  # Read and preprocess each dataset.
  original_string <- as.character(data_set_ind_kept[i])
  
  # Pad with zeros to make it 4 digits.
  padded_string <- sprintf("%04s", original_string)
  
  file_name<-paste0("eSCM_code/datasets/pairs/pair",padded_string,".txt")
  
  data <- read.table(file_name, header = F)
  
  if(dim(data)[2]==1){
    data <- read.table(file_name, header = F, sep = "\t")
  }else if(dim(data)[2]==3){
    
    if(data_set_ind_kept[i] %in% c(81,82,83)){
      data<-data[,c(1,2)]
    }else{
      data<-data[,c(1,3)]
    }
    
  }else if(dim(data)[2]==4){
    
    data<-data[,c(1,3)]
    
  }
  
  na_ind<-which(apply(data,1,function(x){sum(is.na(x))>0} ))
  if(length( na_ind==0)>0){
    data<-data[-na_ind,]
  }
  
  
  for(C_ind in 1:length(C_vec)){
    
    for(L_ind in 1:length(L_vec)){
      
      for(d_ind in 1:length(d_vec)){
        
        name <- paste0("C_", C_vec[C_ind], "_L_", L_vec[L_ind], "_d_", d_vec[d_ind])
        
        # Extremal subsample size k.
        k<-round((nrow(data))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        # Subsample size m used to fit the upper tail; selected adaptively
        # using StrongStop in gpdSeqTests.
        k_tail_fit<-apply(data,2, function(x){find_k_for_fit_GPD(x,n=dim(data)[1],power_order=power_order,C_vec_k_tail_fit=C_vec_k_tail_fit) })
        if(k_tail_fit[1]< 15){
          k_tail_fit[1]<-15
        }
        if(k_tail_fit[2]< 15){
          k_tail_fit[2]<-15
        }
        
        # Convert to a vector.
        tau0 <- 1 - k_tail_fit/dim(data)[1]
        
        u <- mapply(function(col, tau) {
          as.numeric(quantile(col, probs = tau, type = 8))
        }, data, tau0)
        
        # If u equals the maximum value, step down the threshold until the
        # exceedance sample has variation.
        ind_equal_max<-which(u==apply(data,2,max))
        if(length(ind_equal_max)>0){
          
          for(ind_equal_max_index in ind_equal_max){
            
            data_exc_second_large<-which(data[,ind_equal_max_index]<u[ind_equal_max_index])
            u[ind_equal_max_index]<-max(data[data_exc_second_large,ind_equal_max_index])
            # If there is still no variation in the upper tail, step down once more.
            if(var(data[,ind_equal_max_index][which(data[,ind_equal_max_index]>u[ind_equal_max_index])])==0 ){
              data_exc_third_large<-which(data[,ind_equal_max_index]<u[ind_equal_max_index])
              u[ind_equal_max_index]<-max(data[data_exc_third_large,ind_equal_max_index])
            }
            
            
          }
          
          
        }
        
        fits <- mapply(
          function(col, thr) {
            fit.gpd(col, threshold = thr, method = "Grimshaw")
          },
          # This makes each object correspond to one column.
          as.data.frame(data),
          u,
          SIMPLIFY = FALSE
        )
        
        params <- lapply(fits, coef)
        sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
        gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
        
        exceed_idx <- mapply(
          function(col, thr) which(col > thr),
          as.data.frame(data),
          u,
          SIMPLIFY = FALSE
        )
        
        p <- ncol(data)
        n <- nrow(data)
        
        F_x_u <- vector("list", p)
        MT_upper_u <- vector("list", p)
        
        for (j in seq_len(p)) {
          
          idx <- exceed_idx[[j]]   # Indices for column j.
          xexc <- data[idx, j]     # Exceedances.
          
          z <- (xexc - u[j]) / sigma_u_hat[j]
          
          Fval <- 1 - (k_tail_fit[j] / n) * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
          
          # Avoid exact 1 due to floating-point issues.
          if (any(Fval >= 1)) {
            
            ind <- which(Fval >= 1)
            
            if(sum(Fval >= 1)<length(Fval)){
              
              Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
            }else{
              # If all values are equal to 1.
              Fval[ind]<-0.99999
            }
            
          }
          
          F_x_u[[j]] <- Fval
          MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
        }
        
        # Marginal transformation.
        R<-apply(data,2,rank_top_down)
        Z<-(dim(data)[1]/R)^(1/alpha)
        
        for(j in 1:dim(Z)[2]){
          
          Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
          
        }
        
        k<-round((nrow(data))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        AAC_val_list[[name]][i,1:2]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,3:4]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
        data_trans<-cbind(-data[,1],data[,2])
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        AAC_val_list[[name]][i,3:4]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,5:6]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
        data_trans<-cbind(-data[,1],-data[,2])
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        AAC_val_list[[name]][i,5:6]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,7:8]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
        data_trans<-cbind(data[,1],-data[,2])
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        AAC_val_list[[name]][i,7:8]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,9:10]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
      }
      
      
    }
    
  }
  
  
}


fin_res<-list(AAC_order_list=AAC_order_list,AAC_val_list=AAC_val_list,lambda=lambda,r=r,C_vec=C_vec,C_vec_k_tail_fit=C_vec_k_tail_fit
              ,L_vec=L_vec,d_vec=d_vec, k_tail_fit= k_tail_fit)
file_name<-paste0("eSCM_code/datasets/pairwise_datasets_lambda_",lambda,"_fit_tail_adaptive")
save(fin_res,file=file_name)