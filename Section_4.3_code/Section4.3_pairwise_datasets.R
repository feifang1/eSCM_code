# Purpose: Apply EASE with AAC to estimate causal orderings across 94 datasets.

# Data source:
# J. M. Mooij, J. Peters, D. Janzing, J. Zscheischler, B. Schoelkopf:
# "Distinguishing cause from effect using observational data: methods and benchmarks"
# https://webdav.tuebingen.mpg.de/cause-effect/

rm(list=ls())
#install.packages("graph")
library(graph)
#install.packages("pcalg")
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
source("function_class.R")

opt_k_est<-function(Z){
  
  B<-Z[,1]+Z[,2]
  cutoff<-fit_power_law(B)$xmin
  k<-sum(B>= cutoff)
  return(k)
  
}

########### Functions specific to the pairwise datasets ##############
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
      
      # Reject nonpositive interval endpoints.
      if (a <= 0 || b <= 0) {
        stop("Both parameters (a and b) must be positive.")
      }
      
      # Compute the weighted penalties for falling outside [a, b].
      part_1 <- data.frame(
        comp_1 = L_r_vec*(W-b),
        comp_2 = L_r_vec*(a-W)
      )
      
      # Take the larger endpoint penalty, with zero penalty inside the interval.
      d_x_y_C <- pmax(part_1$comp_1, part_1$comp_2, 0)
      
      D_minus_H<-sum(d_x_y_C)/k
      
      # Minimize interval width plus the average penalty scaled by lambda * k^r.
      return((b - a) + lambda*(k^(r))* D_minus_H)
      
      
    }
    
    
    A <- matrix(c(1, 0,   # Constraint: a >= 0.
                  0, -1,  # Constraint: b <= 1.
                  -1, 1), # Constraint: b >= a.
                byrow = TRUE, ncol = 2)
    
    rhs_scale<-c(0,-1,0)
    
    # Start strictly inside the feasible region: 0 < a < b < 1.
    start_point <- c(0.5, 0.6)  # Strictly feasible starting values.
    
    # Minimize the constrained objective using Nelder-Mead.
    result <- constrOptim(
      theta = start_point,
      f = function(x) objective_function(x, first_k_Z, first_k_large_stat, k, W,lambda,r),
      grad = NULL,      # Nelder-Mead does not require a supplied gradient.
      ui = A,           # Rows define A %*% c(a, b) >= rhs_scale.
      ci =rhs_scale,    # Right-hand sides for the three constraints.
      method = "Nelder-Mead"
    )
    
    a_opt<-result$par[1]
    b_opt<-result$par[2]
    
    
  }else if(dis==2){
    
    
    objective_function <- function(x, first_k_Z, first_k_large_stat, k, W,lambda,r) {
      a<-x[1]
      b<-x[2]
      
      # Reject nonpositive interval endpoints.
      if (a <= 0 || b <= 0) {
        stop("Both parameters (a and b) must be positive.")
      }
      
      # Compute the weighted penalties for falling outside [a, b].
      part_1 <- data.frame(
        comp_1 = L_r_vec*(1-(1-W)/(1-b)),
        comp_2 = L_r_vec*(1-W/a)
      )
      
      # Take the larger endpoint penalty, with zero penalty inside the interval.
      d_x_y_C <- pmax(part_1$comp_1, part_1$comp_2, 0)
      
      D_minus_H<-sum(d_x_y_C)/k
      
      # Minimize interval width plus the average penalty scaled by lambda * k^r.
      return((b - a) + lambda*(k^(r))* D_minus_H)
      
    }
    
    
    A <- matrix(c(1, 0,   # Constraint: a >= 0.
                  0, -1,  # Constraint: b <= 1.
                  -1, 1), # Constraint: b >= a.
                byrow = TRUE, ncol = 2)
    
    rhs_scale<-c(0,-1,0)
    
    # Start strictly inside the feasible region: 0 < a < b < 1.
    start_point <- c(0.5, 0.6)  # Strictly feasible starting values.
    
    # Minimize the constrained objective using Nelder-Mead.
    result <- constrOptim(
      theta = start_point,
      f = function(x) objective_function(x, first_k_Z, first_k_large_stat, k, W,lambda,r),
      grad = NULL,      # Nelder-Mead does not require a supplied gradient.
      ui = A,           # Rows define A %*% c(a, b) >= rhs_scale.
      ci =rhs_scale,    # Right-hand sides for the three constraints.
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

# Select the upper-tail count m using sequential GPD goodness-of-fit tests with StrongStop.
find_k_for_fit_GPD <- function(
    x, n, power_order, C_vec_k_tail_fit,
    dataset_id, quadrant, column_id) {
  
  k_fallback <- round(n^power_order * C_vec_k_tail_fit)
  
  context <- sprintf(
    "Dataset %s | quadrant %s | column %s",
    dataset_id, quadrant, column_id
  )
  
  k <- tryCatch(
    withCallingHandlers({
      
      lower_k <- min(k_fallback, 15)
      by_raw <- (lower_k - k_fallback) / 50
      
      by_use <- if (by_raw > -1 && by_raw < 0) {
        -1L
      } else {
        as.integer(round(by_raw))
      }
      
      tab <- as.data.frame(
        gpdSeqTests(
          x,
          nextremes = seq(
            k_fallback,
            lower_k,
            by = by_use
          )
        )
      )
      
      if (any(!is.finite(tab$StrongStop))) {
        stop("Non-finite StrongStop values")
      }
      
      # Last candidate rejected by StrongStop at level 0.05; zero if none.
      tab_ind <- max(c(0L, which(tab$StrongStop <= 0.05)))
      
      if (tab_ind == nrow(tab)) {
        stop("All candidate thresholds rejected")
      }
      
      tab$num.above[tab_ind + 1L]
      
    }, warning = function(w) {
      message(sprintf(
        "[WARNING] %s: %s",
        context, conditionMessage(w)
      ))
      invokeRestart("muffleWarning")
    }),
    
    error = function(e) {
      message(sprintf(
        "[ERROR / FALLBACK] %s: %s; using tail count %d",
        context, conditionMessage(e), k_fallback
      ))
      k_fallback
    }
  )
  
  return(k)
}


weight_tab<-read.table("datasets/pairs/pairmeta.txt", header = F)
colnames(weight_tab)[6]<-"weight"

# Exclude these dataset IDs from the 108-pair benchmark.
remov_ind<-c(33:37,47,52:55,70,71,105,107)

# Penalty multiplier lambda and exponent r used in this run.
lambda<-4
r<-1/2

# Multiplier C in the AAC tail size k = max(round(n^power_order * C), 15).
C_vec<-1

# Multiplier for the initial/fallback marginal GPD tail count m.
C_vec_k_tail_fit<-1

# Settings for the radial weights and penalty in the objective function.
L_vec<-c(3)
d_vec<-c(1)


# Store the true ordering and the AAC-based ordering for each quadrant.
AAC_order<-matrix(0,nrow=108-length(remov_ind),ncol=10 )
AAC_order<-as.data.frame(AAC_order)

# Label the true ordering and the four estimated orderings.
col_names <- c("truth_cause", "truth_effect", 
               "first_quad_cause","first_quad_effect",
               "second_quad_cause","second_quad_effect", 
               "third_quad_cause","third_quad_effect",
               "fourth_quad_cause","fourth_quad_effect")

# Apply the ordering labels to the result table.
colnames(AAC_order) <- col_names

AAC_val<-matrix(0,nrow=108-length(remov_ind),ncol=8 )
AAC_val<-as.data.frame(AAC_val)

col_names <- c("first_quad_cause","first_quad_effect",
               "second_quad_cause","second_quad_effect", 
               "third_quad_cause","third_quad_effect",
               "fourth_quad_cause","fourth_quad_effect")

colnames(AAC_val)<-col_names

data_set_ind_kept<-setdiff(1:108,remov_ind)

# Encode the true causal ordering.
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

# Enumerate all combinations of the AAC tail-size multiplier and objective settings.
param_grid <- expand.grid(C = C_vec, L = L_vec, d = d_vec)

# Construct a descriptive name for each parameter combination.
list_names <- with(param_grid, paste0("C_", C, "_L_", L, "_d_", d))

# Use the same parameter names for ordering and interval-endpoint results.
names(AAC_order_list) <- names(AAC_val_list)<-list_names

power_order<-0.7

size_rec<-rep(0,length(setdiff(1:108,remov_ind)))
k_pos_rec<-rep(0,length(setdiff(1:108,remov_ind)))

# Tail index alpha for the Pareto-type marginal transformation.
alpha<-2

size_graph<-rep(0,length(setdiff(1:108,remov_ind)))

for(i in 1:length(setdiff(1:108,remov_ind)) ){
  
  print(paste0("dataset ",i))
  
  # Read the retained benchmark pair and select its two analysis variables.
  original_string <- as.character(data_set_ind_kept[i])
  
  # Format the dataset ID for the four-character benchmark filename.
  padded_string <- sprintf("%04s", original_string)
  
  file_name<-paste0("datasets/pairs/pair",padded_string,".txt")
  
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
        
        # Number of largest radial sums used for AAC estimation; enforce k >= 15.
        k<-round((nrow(data))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        # Select the marginal tail counts m used to set the initial quantile thresholds.
        k_tail_fit <- vapply(seq_len(ncol(data)), function(j) {
          find_k_for_fit_GPD(
            x = data[, j],
            n = nrow(data),
            power_order = power_order,
            C_vec_k_tail_fit = C_vec_k_tail_fit,
            dataset_id = data_set_ind_kept[i],
            quadrant = 1,
            column_id = j
          )
        }, numeric(1))
        
        if(k_tail_fit[1]< 15){
          k_tail_fit[1]<-15
        }
        if(k_tail_fit[2]< 15){
          k_tail_fit[2]<-15
        }
        
        # Convert the selected marginal tail counts to quantile probabilities.
        tau0 <- 1 - k_tail_fit/dim(data)[1]
        
        u <- mapply(function(col, tau) {
          as.numeric(quantile(col, probs = tau, type = 8))
        }, as.data.frame(data), tau0)
        
        # Lower the threshold until there are at least three strict
        # exceedances and at least two distinct values above it.
        for (j in seq_len(ncol(data))) {
          x <- data[, j]
          if (any(!is.finite(x)) || !is.finite(u[j])) {
            stop(sprintf("Dataset %s, quadrant 1, column %d: non-finite data or threshold.",
                         original_string, j))
          }
          repeat {
            exc <- x[x > u[j]]
            if (length(exc) >= 3L && length(unique(exc)) >= 2L) {
              break
            }
            lower_values <- x[x < u[j]]
            if (length(lower_values) == 0L) {
              stop(sprintf("Dataset %s, quadrant 1, column %d: insufficient upper-tail variation for GPD fitting.",
                           original_string, j))
            }
            u[j] <- max(lower_values)
          }
        }
        
        fits <- mapply(
          function(col, thr) {
            fit.gpd(col, threshold = thr, method = "Grimshaw")
          },
          # Fit each marginal separately at its final threshold; retain a list of fits.
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
        
        # Actual strict-exceedance proportion at the final threshold.
        # Use the actual exceedance counts rather than the nominal tail counts.
        tail_prob <- lengths(exceed_idx) / n
        
        F_x_u <- vector("list", p)
        MT_upper_u <- vector("list", p)
        
        for (j in seq_len(p)) {
          
          idx <- exceed_idx[[j]]   # Row indices strictly above the threshold in margin j.
          xexc <- data[idx, j]     # Observed values above the threshold (before subtracting it).
          
          z <- (xexc - u[j]) / sigma_u_hat[j]
          
          Fval <- 1 - tail_prob[j] * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
          
          # Keep fitted CDF values below one so the marginal transformation stays finite.
          if (any(Fval >= 1)) {
            
            ind <- which(Fval >= 1)
            
            if(sum(Fval >= 1)<length(Fval)){
              
              Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
            }else{
              # If every fitted CDF value is at least one, use the fixed cap.
              Fval[ind]<-0.99999
            }
            
          }
          
          F_x_u[[j]] <- Fval
          MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
        }
        
        # Apply the rank-based marginal transformation and replace the fitted upper-tail values.
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
        
        # ---- Second quadrant: (-X1, X2) ----
        data_trans<-cbind(-data[,1],data[,2])
        
        k_tail_fit <- vapply(seq_len(ncol(data_trans)), function(j) {
          find_k_for_fit_GPD(
            x = data_trans[, j],
            n = nrow(data_trans),
            power_order = power_order,
            C_vec_k_tail_fit = C_vec_k_tail_fit,
            dataset_id = data_set_ind_kept[i],
            quadrant = 2,  # Quadrant label used in warning and error messages.
            column_id = j
          )
        }, numeric(1))
  
        if(k_tail_fit[1]< 15){
          k_tail_fit[1]<-15
        }
        if(k_tail_fit[2]< 15){
          k_tail_fit[2]<-15
        }
        
        
        # Convert the selected marginal tail counts to quantile probabilities.
        tau0 <- 1 - k_tail_fit/dim(data_trans)[1]
        
        u <- mapply(function(col, tau) {
          as.numeric(quantile(col, probs = tau, type = 8))
        }, as.data.frame(data_trans), tau0)
        
        # Lower the threshold until there are at least three strict
        # exceedances and at least two distinct values above it.
        for (j in seq_len(ncol(data_trans))) {
          x <- data_trans[, j]
          if (any(!is.finite(x)) || !is.finite(u[j])) {
            stop(sprintf("Dataset %s, quadrant 2, column %d: non-finite data or threshold.",
                         original_string, j))
          }
          repeat {
            exc <- x[x > u[j]]
            if (length(exc) >= 3L && length(unique(exc)) >= 2L) {
              break
            }
            lower_values <- x[x < u[j]]
            if (length(lower_values) == 0L) {
              stop(sprintf("Dataset %s, quadrant 2, column %d: insufficient upper-tail variation for GPD fitting.",
                           original_string, j))
            }
            u[j] <- max(lower_values)
          }
        }
        
        fits <- mapply(
          function(col, thr) {
            fit.gpd(col, threshold = thr, method = "Grimshaw")
          },
          # Fit each marginal separately at its final threshold; retain a list of fits.
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        params <- lapply(fits, coef)
        sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
        gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
        
        exceed_idx <- mapply(
          function(col, thr) which(col > thr),
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        p <- ncol(data_trans)
        n <- nrow(data_trans)
        
        # Actual strict-exceedance proportion at the final threshold.
        # Use the actual exceedance counts rather than the nominal tail counts.
        tail_prob <- lengths(exceed_idx) / n
        
        F_x_u <- vector("list", p)
        MT_upper_u <- vector("list", p)
        
        for (j in seq_len(p)) {
          
          idx <- exceed_idx[[j]]   # Row indices strictly above the threshold in margin j.
          xexc <- data_trans[idx, j]     # Observed values above the threshold (before subtracting it).
          
          z <- (xexc - u[j]) / sigma_u_hat[j]
          
          Fval <- 1 - tail_prob[j] * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
          
          # Keep fitted CDF values below one so the marginal transformation stays finite.
          if (any(Fval >= 1)) {
            
            ind <- which(Fval >= 1)
            
            if(sum(Fval >= 1)<length(Fval)){
              
              Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
            }else{
              # If every fitted CDF value is at least one, use the fixed cap.
              Fval[ind]<-0.99999
            }
            
          }
          
          F_x_u[[j]] <- Fval
          MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
        }
        
        
        # Apply the rank-based marginal transformation and replace the fitted upper-tail values.
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        
        for(j in 1:dim(Z)[2]){
          
          Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
          
        }
        
        k<-round((nrow(data_trans))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        AAC_val_list[[name]][i,3:4]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,5:6]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
        
        # ---- Third quadrant: (-X1, -X2) ----
        
        data_trans<-cbind(-data[,1],-data[,2])
        
        
        k_tail_fit <- vapply(seq_len(ncol(data_trans)), function(j) {
          find_k_for_fit_GPD(
            x = data_trans[, j],
            n = nrow(data_trans),
            power_order = power_order,
            C_vec_k_tail_fit = C_vec_k_tail_fit,
            dataset_id = data_set_ind_kept[i],
            quadrant = 3,  # Quadrant label used in warning and error messages.
            column_id = j
          )
        }, numeric(1))
        
        if(k_tail_fit[1]< 15){
          k_tail_fit[1]<-15
        }
        if(k_tail_fit[2]< 15){
          k_tail_fit[2]<-15
        }
        
        # Convert the selected marginal tail counts to quantile probabilities.
        tau0 <- 1 - k_tail_fit/dim(data_trans)[1]
        
        u <- mapply(function(col, tau) {
          as.numeric(quantile(col, probs = tau, type = 8))
        }, as.data.frame(data_trans), tau0)
        
        # Lower the threshold until there are at least three strict
        # exceedances and at least two distinct values above it.
        for (j in seq_len(ncol(data_trans))) {
          x <- data_trans[, j]
          if (any(!is.finite(x)) || !is.finite(u[j])) {
            stop(sprintf("Dataset %s, quadrant 3, column %d: non-finite data or threshold.",
                         original_string, j))
          }
          repeat {
            exc <- x[x > u[j]]
            if (length(exc) >= 3L && length(unique(exc)) >= 2L) {
              break
            }
            lower_values <- x[x < u[j]]
            if (length(lower_values) == 0L) {
              stop(sprintf("Dataset %s, quadrant 3, column %d: insufficient upper-tail variation for GPD fitting.",
                           original_string, j))
            }
            u[j] <- max(lower_values)
          }
        }
        
        fits <- mapply(
          function(col, thr) {
            fit.gpd(col, threshold = thr, method = "Grimshaw")
          },
          # Fit each marginal separately at its final threshold; retain a list of fits.
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        params <- lapply(fits, coef)
        sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
        gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
        
        exceed_idx <- mapply(
          function(col, thr) which(col > thr),
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        p <- ncol(data_trans)
        n <- nrow(data_trans)
        
        # Actual strict-exceedance proportion at the final threshold.
        # Use the actual exceedance counts rather than the nominal tail counts.
        tail_prob <- lengths(exceed_idx) / n
        
        F_x_u <- vector("list", p)
        MT_upper_u <- vector("list", p)
        
        for (j in seq_len(p)) {
          
          idx <- exceed_idx[[j]]   # Row indices strictly above the threshold in margin j.
          xexc <- data_trans[idx, j]     # Observed values above the threshold (before subtracting it).
          
          z <- (xexc - u[j]) / sigma_u_hat[j]
          
          Fval <- 1 - tail_prob[j] * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
          
          # Keep fitted CDF values below one so the marginal transformation stays finite.
          if (any(Fval >= 1)) {
            
            ind <- which(Fval >= 1)
            
            if(sum(Fval >= 1)<length(Fval)){
              
              Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
            }else{
              # If every fitted CDF value is at least one, use the fixed cap.
              Fval[ind]<-0.99999
            }
            
          }
          
          F_x_u[[j]] <- Fval
          MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
        }
        
        
        # Apply the rank-based marginal transformation and replace the fitted upper-tail values.
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        
        for(j in 1:dim(Z)[2]){
          
          Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
          
        }
        
        k<-round((nrow(data_trans))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        AAC_val_list[[name]][i,5:6]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,7:8]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
        # ---- Fourth quadrant: (X1, -X2) ----
        
        data_trans<-cbind(data[,1],-data[,2])
        
        k_tail_fit <- vapply(seq_len(ncol(data_trans)), function(j) {
          find_k_for_fit_GPD(
            x = data_trans[, j],
            n = nrow(data_trans),
            power_order = power_order,
            C_vec_k_tail_fit = C_vec_k_tail_fit,
            dataset_id = data_set_ind_kept[i],
            quadrant = 4,  # Quadrant label used in warning and error messages.
            column_id = j
          )
        }, numeric(1))
      
        if(k_tail_fit[1]< 15){
          k_tail_fit[1]<-15
        }
        if(k_tail_fit[2]< 15){
          k_tail_fit[2]<-15
        }
        
        # Convert the selected marginal tail counts to quantile probabilities.
        tau0 <- 1 - k_tail_fit/dim(data_trans)[1]
        
        u <- mapply(function(col, tau) {
          as.numeric(quantile(col, probs = tau, type = 8))
        }, as.data.frame(data_trans), tau0)
        
        # Lower the threshold until there are at least three strict
        # exceedances and at least two distinct values above it.
        for (j in seq_len(ncol(data_trans))) {
          x <- data_trans[, j]
          if (any(!is.finite(x)) || !is.finite(u[j])) {
            stop(sprintf("Dataset %s, quadrant 4, column %d: non-finite data or threshold.",
                         original_string, j))
          }
          repeat {
            exc <- x[x > u[j]]
            if (length(exc) >= 3L && length(unique(exc)) >= 2L) {
              break
            }
            lower_values <- x[x < u[j]]
            if (length(lower_values) == 0L) {
              stop(sprintf("Dataset %s, quadrant 4, column %d: insufficient upper-tail variation for GPD fitting.",
                           original_string, j))
            }
            u[j] <- max(lower_values)
          }
        }
        
        fits <- mapply(
          function(col, thr) {
            fit.gpd(col, threshold = thr, method = "Grimshaw")
          },
          # Fit each marginal separately at its final threshold; retain a list of fits.
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        params <- lapply(fits, coef)
        sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
        gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
        
        exceed_idx <- mapply(
          function(col, thr) which(col > thr),
          as.data.frame(data_trans),
          u,
          SIMPLIFY = FALSE
        )
        
        p <- ncol(data_trans)
        n <- nrow(data_trans)
        
        # Actual strict-exceedance proportion at the final threshold.
        # Use the actual exceedance counts rather than the nominal tail counts.
        tail_prob <- lengths(exceed_idx) / n
        
        F_x_u <- vector("list", p)
        MT_upper_u <- vector("list", p)
        
        for (j in seq_len(p)) {
          
          idx <- exceed_idx[[j]]   # Row indices strictly above the threshold in margin j.
          xexc <- data_trans[idx, j]     # Observed values above the threshold (before subtracting it).
          
          z <- (xexc - u[j]) / sigma_u_hat[j]
          
          Fval <- 1 - tail_prob[j] * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
          
          # Keep fitted CDF values below one so the marginal transformation stays finite.
          if (any(Fval >= 1)) {
            
            ind <- which(Fval >= 1)
            
            if(sum(Fval >= 1)<length(Fval)){
              
              Fval[ind] <- min(0.99999, 1.001 * max(Fval[-ind]))
            }else{
              # If every fitted CDF value is at least one, use the fixed cap.
              Fval[ind]<-0.99999
            }
            
          }
          
          F_x_u[[j]] <- Fval
          MT_upper_u[[j]] <- 1 / (1 - Fval)^(1 / alpha)
        }
        
        
        # Apply the rank-based marginal transformation and replace the fitted upper-tail values.
        R<-apply(data_trans,2,rank_top_down)
        Z<-(dim(data_trans)[1]/R)^(1/alpha)
        
        for(j in 1:dim(Z)[2]){
          
          Z[exceed_idx[[j]],j]<-MT_upper_u[[j]]
          
        }
        
        k<-round((nrow(data_trans))^(power_order)*C_vec[C_ind])
        if(k < 15){
          k<-15
        }
        
        AAC_val_list[[name]][i,7:8]<-AAC_val_opt(Z=Z,a_b_val=T,k=k,lambda,L_r=L_vec[L_ind],dis=d_vec[d_ind],r=r)[1:2]
        AAC_order_list[[name]][i,9:10]<-AAC_ease_order(dat=Z, a_b_val=F, k = k, lambda=lambda, L_r=L_vec[L_ind], dis=d_vec[d_ind], r=r )
        
      }
      
      
    }
    
  }
  
  
}


fin_res<-list(AAC_order_list=AAC_order_list,AAC_val_list=AAC_val_list,lambda=lambda,r=r,C_vec=C_vec,C_vec_k_tail_fit=C_vec_k_tail_fit
              ,L_vec=L_vec,d_vec=d_vec, k_tail_fit= k_tail_fit)
file_name<-paste0("datasets/Section_4.3/pairwise_datasets_lambda_",lambda,"_a_0.7_m_sqrt_n_alpha_",alpha,"_with_four_quadrant_fit_with_tail_strong_stop_adjust_m")
save(fin_res,file=file_name)




