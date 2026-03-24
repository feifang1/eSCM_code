rank_top_down<-function(x){
  
  # Input: vector to be ranked from largest to smallest.
  rank_vec<-rank(-x,ties.method = "max")
  return(rank_vec)
  
}


order_error_measure<-function(dag,ab_dir){
  
  # Compute the ancestral violation rate.
  
  p<-dim(dag)[1]
  
  num_path_count<-0
  count_error_ab<-0
  
  for( i in 1:p ){
    
    path_mat <- as.matrix(dag) %^% i
    if(sum(path_mat)==0){
      break 
    }
    
    indices <- which(path_mat > 0, arr.ind = TRUE)
    num_path_count<-num_path_count+nrow(indices)
    # Check whether the estimated order is consistent. 
    
    order_match_ab<-apply(indices,1,function(x){ match(x, ab_dir) }  )
    # True order for first-distance neighbors:
    # 1->2, 2->3, 1->3
    
    # Check whether the order in b is consistent with the order in a.
    count_error_ab<-count_error_ab+sum( (order_match_ab[1,]<order_match_ab[2,])==F )
  }
  
  # num_path_count is the total number of oracle direct and indirect paths.
  error_rate_ab<-count_error_ab/num_path_count
  
  return(error_rate_ab)
  
}


ease <- function(dat, k = floor(n ^ 0.4),
                 both_tails = TRUE){
  # Set up variables.
  n <- NROW(dat)
  d <- NCOL(dat)
  
  # Compute the causal tail matrix.
  causal_mat <- causal_tail_matrix(dat, k, both_tails)
  
  # Run Extremal Ancestral Search.
  current_order <- add <- which.min(apply(causal_mat, 2, max, na.rm = TRUE))
  for (k in 2:d){
    causal_mat[add, ] <- NA
    avail <- (1:d)[-current_order]
    
    add <- if (k < d){
      avail[which.min(apply(causal_mat, 2, max, na.rm = TRUE)[avail])]
    } else {
      avail
    }
    current_order <- c(current_order, add)
  }
  
  order <- current_order
  
  # Return the causal order.
  return(order=order)
}

AAC_val_opt<-function(Z,a_b_val=F,k,lambda,L_r=3,dis,r=1/2){
  
  # Obtain the optimal estimates of a and b by optimizing g_n(s,t).
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
      ci =rhs_scale,           # Constraint vector.
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
      ci =rhs_scale,           # Constraint vector.
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


causal_tail_matrix <- function(dat, k = floor(n ^ 0.4),
                               both_tails = TRUE){
  
  # Get the numbers of observations and variables.
  n   <- NROW(dat)
  p   <- NCOL(dat)
  
  # Rank the variables.
  ranked_dat <- apply(dat, 2, rank, ties.method = "first")
  
  # Compute the causal tail coefficient for all variable pairs.
  causal_mat <- sapply(1:p,
                       function(j){
                         sapply(1:p,
                                function(i){
                                  if (i == j){
                                    NA
                                  } else {
                                    causal_tail_coeff(ranked_dat[, i],
                                                      ranked_dat[, j], k,
                                                      to_rank = FALSE,
                                                      both_tails = both_tails)
                                  }
                                })
                       })
  
  # Return the matrix.
  return(causal_mat)
}

causal_tail_coeff <- function(v1, v2, k = floor(n ^ 0.4), to_rank = TRUE,
                              both_tails = TRUE){
  # Number of observations.
  n <- NROW(v1)
  
  # Check k.
  if (k <= 1 | k >= n) {
    stop("k must be greater than 1 and smaller than n.")
  }
  
  # Rank variables if needed.
  if (to_rank){
    r1 <- rank(v1, ties.method = "first") # Ranks of v1.
    r2 <- rank(v2, ties.method = "first") # Ranks of v2.
  } else{
    r1 <- v1
    r2 <- v2
  }
  
  # Compute the causal tail coefficient.
  if (both_tails){
    k <- (k %/% 2) * 2
    1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / 2))
  } else{
    1 / (k * n) * sum(r2[r1 > n - k])
  }
}


AAC_ease_order<-function(dat, a_b_val=F, k = floor(n ^ 0.4), lambda, L_r, dis, r ){
  
  # EASE algorithm based on AAC coefficients.
  
  # Set up variables.
  n <- NROW(dat)
  d <- NCOL(dat)
  
  causal_mat<-matrix(0,d,d)
  
  # Compute the pairwise AAC orientation matrix.
  for(i in 1:d){
    
    for(j in 1:d){
      
      causal_mat[i,j]<-AAC_val_opt(Z=dat[,c(i,j)],a_b_val=F,k=k,lambda=lambda,L_r=L_r,dis=dis,r=r)
      
    }
    
  }
  
  # Run Extremal Ancestral Search.
  current_order <- add <- which.min(apply(causal_mat, 2, max, na.rm = TRUE))
  for (h in 2:d){
    causal_mat[add, ] <- NA
    avail <- (1:d)[-current_order]
    
    add <- if (h < d){
      avail[which.min(apply(causal_mat, 2, max, na.rm = TRUE)[avail])]
    } else {
      avail
    }
    current_order <- c(current_order, add)
  }
  
  order <- current_order
  
  # Return the causal order.
  return(order=order)
}


combine_custom<-function( list_1, list_2) {
  # Save the results in parallel mode.
  
  if ( length(list_1$temp_rec_CTC) > 1 & class(list_1$temp_rec_CTC)=="list"     ) {
    # When list_1$temp_rec_CTC has already become a list.
    temp_rec_CTC<- c(list_1$temp_rec_CTC, list(list_2$temp_rec_CTC))
    
  } else {
    # list_1$temp_rec_CTC and list_2$temp_rec_CTC are still data frames.
    temp_rec_CTC <- c(list(list_1$temp_rec_CTC), list(list_2$temp_rec_CTC))
    
  }
  
  if ( length(list_1$temp_rec_AAC) > 1 & class(list_1$temp_rec_AAC)=="list"     ) {
    # When list_1$temp_rec_AAC has already become a list.
    temp_rec_AAC<- c(list_1$temp_rec_AAC, list(list_2$temp_rec_AAC))
    
  } else {
    # list_1$temp_rec_AAC and list_2$temp_rec_AAC are still data frames.
    temp_rec_AAC <- c(list(list_1$temp_rec_AAC), list(list_2$temp_rec_AAC))
    
  }
  
  # Note: the names must match those used in the list;
  # otherwise, the second combine step cannot find them.
  ls <- c( temp_rec_CTC=list(temp_rec_CTC), temp_rec_AAC=list(temp_rec_AAC) )
  
  return(ls)
  
}

plus_max_linear_mat_mat<-function(A,eta_mat){
  
  Y_mat<-matrix(0, nrow=dim(eta_mat)[1], ncol=dim(A)[1] )
  
  for(i in 1:dim(A)[1]){
    
    Y_mat[,i]<-plus_max_linear_vec_mat(A[i,],eta_mat)
    
  }
  
  return(Y_mat)
  
}

plus_max_linear_vec_mat<-function(a,eta_mat){
  
  return(apply(eta_mat,1, function(x){ unique(max(a*x)) } ))
  
}