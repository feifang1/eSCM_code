# Goal: Obtain the pairwise causal direction identification error rate for river discharge data.

rm(list=ls())
library(graph)
library(pcalg)
library(EnvStats)
source("function_class.R")

library(dplyr)
library(Matrix)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(expm)

library(Hmisc)
library(latex2exp)
require(reshape2)
library(tidyquant)
library(timetk)
library(tsibble)
library(ismev)
library(pracma)
library(tidyverse)

library(mev)

# Function to obtain the optimal estimates of a and b by minimizing g_n(s,t).
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
      grad = NULL,      # No gradient is provided; let constrOptim() approximate it.
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
      
      # Compute the element-wise maximum using pmax.
      d_x_y_C <- pmax(part_1$comp_1, part_1$comp_2, 0)
      
      D_minus_H<-sum(d_x_y_C)/k
      
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
      grad = NULL,      # No gradient is provided; let constrOptim() approximate it.
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


###################### Preprocess the dataset based on the code and data from the paper "EXTREMES ON RIVER NETWORKS" ########################
theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank()))

## Define constants ####
OUTPUT_FILE <- "output/river_results.txt"
tolPalette <- c(tolBlack = "#000000",
                tolBlue = "#4477AA",
                tolRed = "#EE6677",
                tolGreen = "#228833",
                tolYellow = "#CCBB44",
                tolCyan = "#66CCEE",
                tolPurple = "#AA3377",
                tolGrey = "#BBBBBB") %>%
  unname()



## Define functions ####
dfr2tibble <- function(dfr){
  ## dataframe -> tibble
  ## Convert a dataframe to a tibble and rename column "Date" to "date".
  
  dfr %>%
    mutate(Date = as_date(Date)) %>%
    rename(date = Date) %>%
    as_tibble()
}

PP.lik.linear <- function(Params, Data, u, NoYears, CovarMu, CovarSc, CovarXi){
  ## Function from the paper "EXTREMES ON RIVER NETWORKS"
  ## By PEIMAN ASADI, ANTHONY C. DAVISON and SEBASTIAN ENGELKE
  
  NoSt <- length(Data)
  NoParMu <- ncol(CovarMu)
  NoParSc <- ncol(CovarSc)
  NoParXi <- ncol(CovarXi)
  NoPar <- NoParMu + NoParSc + NoParXi
  Out <- 0
  
  for (i in 1:NoSt) {
    mu <- sum(Params[1:NoParMu]*(CovarMu[i,]))
    sc <- sum(Params[(NoParMu+1):(NoParMu+NoParSc)]*CovarSc[i,])
    xi <- sum(Params[(NoParMu+NoParSc+1):NoPar]*CovarXi[i,])
    
    y1 <- 1 + xi * ((u[i] - mu)/sc)
    y2 <- 1 + xi * ((Data[[i]]- mu)/sc)
    InterceptSc <- Params[NoParMu+1]
    InterceptMu <- Params[1]
    
    if ((sc <= 0) | (min(y1) <= 0) | (min(y2) <= 0) |
        InterceptSc <= 0 | InterceptMu < 0 ) {
      l <- 10^6
    } else {
      l <- NoYears[i] * (y1 ^ (-1 / xi)) +
        length(Data[[i]]) * log(sc) + (1 / xi + 1) * sum(log(y2))
    }
    Out <- Out + l
  }
  return(Out)
}

PPFit <- function(Data, u, NoYears, Init, CovarMu, CovarSc, CovarXi,
                  method = method, control = control) {
  ## Function from the paper "EXTREMES ON RIVER NETWORKS"
  ## By PEIMAN ASADI, ANTHONY C. DAVISON and SEBASTIAN ENGELKE
  
  x <- optim(Init, PP.lik.linear,
             Data = Data, u = u, NoYears = NoYears, CovarMu = CovarMu,
             CovarSc = CovarSc, CovarXi = CovarXi,
             method = method, control = control, hessian = TRUE)
  output <- x
  return(output)
  
}

split_stations <- function(dat, year = NULL, stations = NULL){
  ## tibble integer character_vector -> list
  ## Produce a named list of discharges for the stations in the given year.
  if (!is.null(year) & !is.null(stations)){
    dat_tmp <- dat %>%
      dplyr::filter(year(date) %in% year, station %in% stations)
    
    if (nrow(dat_tmp) == 0){
      return(list())
    }
  } else {
    dat_tmp <- dat
  }
  
  unique_names <- unique(dat_tmp$station)
  
  dat_tmp %>%
    select(station, disc) %>%
    group_split() %>%
    set_names(unique_names) %>%
    lapply(function(tbl){tbl$disc})
  
}

## Import dataset ####
load("datasets/River_raw_data_ease_algorithm/StsTSs.RData")
load("datasets/River_raw_data_ease_algorithm/StsInfo.RData")

class(StsTSs[[1]])
StsTSs[[1]][1:10,]

# Select the 31 stations used in the AOAS paper plus station 19 (renamed as 32).
StsChos <- c(c(1:47)[-c(16,30,31,34,42,43,44,45,46,47,3,1,2,29,18,19)], 19)

NoSt <- length(StsChos)

StsTSsChos <- StsTSs[StsChos]
StsInfoChos <- StsInfo[StsChos,] %>%
  mutate(id_old = StsChos,
         id = 1:NoSt)

# Clean the data and find common dates.
river_dat <- map(.x = StsTSsChos, .f = dfr2tibble) %>%
  reduce(.f = inner_join, by = "date") %>%
  dplyr::filter(month(date) %in% c(6, 7, 8),
                year(date) < 2010)

# Rename columns according to the AOAS paper.
colnames(river_dat)[2:(NoSt + 1)] <- sprintf("station_%02d", 1:NoSt)

# Select stations.
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)
river_dat <- river_dat[, c(1, station_names + 1)]
NoSt <- length(station_names)

station_info <- StsInfoChos %>%
  dplyr::filter(id %in% station_names) %>%
  rename(name = RivNames,
         lat = Lat,
         lon = Long,
         ave_vol = AveVol)

# Clean workspace.
rm(StsInfo, StsInfoChos, StsTSs, StsTSsChos, StsChos)

## Check model assumptions ####
method <- "BFGS"
control <- list(maxit = 5000, reltol=10^(-30), abstol=0.0001, trace=0)

# Take the maximum for each station within each year.
SummerMaxima <- river_dat %>%
  mutate(year = year(date)) %>%
  gather(key = "station", value = "disc", -date, -year) %>%
  group_by(year, station) %>%
  summarise(maxDisc = max(disc)) %>%
  spread(key = station, value = maxDisc) %>%
  select(year, sprintf("station_%02d", station_names))

# GEV analysis.
GevPars <- tibble(station = character(), location = double(),
                  scale = double(), shape = double())

help(gev.fit)
for (i in 1:NoSt){
  
  vec <- SummerMaxima[[i + 1]]
  station_nm <- colnames(SummerMaxima)[i + 1]
  tmpres <- gev.fit(vec, show = FALSE)
  GevPars <- bind_rows(GevPars,
                       tibble(station = station_nm,
                              location = tmpres$mle[1],
                              scale = tmpres$mle[2],
                              shape = tmpres$mle[3]))
}


# Compute threshold data.
q <- .9
threshold_dat <- river_dat %>%
  gather(key = "station", value = "disc", -date) %>%
  mutate(station = factor(station,
                          levels = sprintf("station_%02d", station_names))) %>%
  group_by(station) %>%
  mutate(thres = quantile(disc, q)) %>%
  dplyr::filter(disc > thres) %>%
  mutate(noyears = length(unique(year(date))))

river_dat

AllEvents <- split_stations(threshold_dat)
Threshold <- threshold_dat %>% distinct(station, thres) %>% deframe()
NoOfYears <- threshold_dat %>% distinct(station, noyears) %>% deframe()


# Fit a GEVD to each station by maximizing the joint Poisson
# process likelihood; see formula (21) in the paper by
# Asadi, Davison, and Engelke (2015).
ParsPP <- tibble(station = character(), location = double(),
                 scale = double(), shape = double())

for (i in 1:NoSt) {
  Init <- GevPars %>% select(-station) %>% slice(i) %>% as.double()
  station_nm <- GevPars %>% select(station) %>% slice(i) %>% as.character()
  
  tmpres <- PPFit(Data=AllEvents[i],
                  u=Threshold[i],
                  NoYears=NoOfYears[i],
                  Init=Init,
                  CovarMu=matrix(1,1,1),CovarSc=matrix(1,1,1),
                  CovarXi=matrix(1,1,1),method=method,control=control)
  ParsPP <- bind_rows(ParsPP,
                      tibble(station = station_nm,
                             location = tmpres$par[1],
                             scale = tmpres$par[2],
                             shape = tmpres$par[3]))
}



# Regionalized model from the paper by Asadi, Davison, and Engelke (2015).
Grp1 <- sprintf("station_%02d", c(11, 19, 21))
Grp2 <- sprintf("station_%02d", c(13, 28, 32))
Grp3 <- sprintf("station_%02d", c(1, 7, 9, 14))
Grp4 <- sprintf("station_%02d", c(23, 26))
Grp <- list(Grp1, Grp2, Grp3, Grp4)

# Estimate parameters: shape, location, and scale.
ParsCovModel <- tibble(station = character(), region = double(),
                       location = double(), scale = double(),
                       shape = double(), shape_se = double())

hessian_list <- list()

for (i in 1:length(Grp)) {
  GrSts <- Grp[[i]]
  
  n_st <- length(GrSts)
  
  CovarMuReg <- diag(n_st)
  CovarScReg <- diag(n_st)
  CovarXi <- matrix(1, nrow = n_st)
  
  station_nms <- ParsPP %>%
    dplyr::filter(station %in% GrSts) %>%
    select(station) %>%
    deframe()
  
  InitMu <- ParsPP %>% dplyr::filter(station %in% GrSts) %>% select(location) %>%
    deframe()
  InitSc <- ParsPP %>% dplyr::filter(station %in% GrSts) %>% select(scale) %>%
    deframe()
  InitXi <- ParsPP %>% dplyr::filter(station %in% GrSts) %>% select(shape) %>%
    summarise(mean(shape)) %>% deframe()
  Init <- c(InitMu, InitSc, InitXi)
  
  tmpres <- PPFit(Data=AllEvents[names(AllEvents) %in% GrSts],
                  u=Threshold[names(Threshold) %in% GrSts],
                  NoYears=NoOfYears[names(NoOfYears) %in% GrSts],
                  Init=Init,
                  CovarMu=CovarMuReg,
                  CovarSc=CovarScReg,
                  CovarXi=CovarXi,
                  method=method, control=control)
  
  
  MuCov <- tmpres$par[1:n_st]
  ScCov <- tmpres$par[1:n_st + n_st]
  XiCov <- tmpres$par[2 * n_st + 1]
  
  nms <- c(paste("loc_", station_nms, sep = ""),
           paste("scale_", station_nms, sep = ""),
           "shape")
  H <- tmpres$hessian
  dimnames(H) <- list(nms, nms)
  hessian_list[[i]] <- H
  
  
  sderror <- sqrt(solve(H)[2 * n_st + 1, 2 * n_st + 1])
  
  ParsCovModel <- bind_rows(ParsCovModel,
                            tibble(station = station_nms,
                                   region = i, location = MuCov, scale = ScCov,
                                   shape = XiCov, shape_se = sderror))
}



# Compute adjusted standard errors accounting for time dependence
# (see Fawcett and Walshaw, 2007 and 2012).
se_adj <- tibble(region = double(), shape_se_adj = double())

for (i in 1:length(Grp)) {
  GrSts <- Grp[[i]]
  
  my_grad <- NULL
  yrs <- unique(year(threshold_dat$date))
  
  for (j in seq_along(yrs)){
    y <- unique(year(threshold_dat$date))[j]
    AllEvents_year <- split_stations(threshold_dat, year = y, stations = GrSts)
    
    GrSts_year <- names(AllEvents_year)
    
    if(length(GrSts) != length(GrSts_year)){
      next()
    }
    
    cat(GrSts_year, "--- year =", y, "\n")
    
    
    n_st <- length(GrSts_year)
    
    # Initialize parameters.
    CovarMuReg <- diag(n_st)
    CovarScReg <- diag(n_st)
    CovarXi <- matrix(1, nrow = n_st)
    
    station_nms <- ParsPP %>%
      dplyr::filter(station %in% GrSts_year) %>%
      select(station) %>%
      deframe()
    
    InitMu <- ParsCovModel %>%
      dplyr::filter(station %in% GrSts_year) %>%
      select(location) %>%
      deframe()
    InitSc <- ParsCovModel %>%
      dplyr::filter(station %in% GrSts_year) %>%
      select(scale) %>%
      deframe()
    InitXi <- ParsCovModel %>%
      dplyr::filter(station %in% GrSts_year) %>%
      select(shape) %>%
      summarise(mean(shape)) %>%
      deframe()
    
    Init <- c(InitMu, InitSc, InitXi)
    
    # Compute the gradient at the optimal parameter value.
    tmpres <- pracma::grad(f = PP.lik.linear,
                           x0 = Init, Data=AllEvents_year[GrSts_year],
                           u=Threshold[GrSts_year],
                           NoYears=NoOfYears[GrSts_year],
                           CovarMu=CovarMuReg,CovarSc=CovarScReg,
                           CovarXi=CovarXi)
    
    my_grad <- rbind(my_grad, tmpres)
    
  }
  
  H <- hessian_list[[i]]
  V <- max(NoOfYears[GrSts_year]) * cov(my_grad)
  
  d <- ncol(my_grad)
  cov_theta <- solve(H) %*% V %*% solve(H)
  
  se_adj <- bind_rows(se_adj,
                      tibble(region = i,
                             shape_se_adj = sqrt(cov_theta[d, d])))
}


## Spatial structure ####
# True DAG.
true_dag <- rbind(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 0, 1, 0),
                  rep(0, 12),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 0, 1, 0))
colnames(true_dag) <- rownames(true_dag) <- station_names

g_true_dat <- igraph::graph_from_adjacency_matrix(true_dag)
igraph::V(g_true_dat)$color <- "white"
# igraph::tkplot(g)

# Data.
mat <- river_dat %>%  select(-date) %>% as.matrix()

# Ground truth.
ll <- list(
  dataset = mat,
  dag = true_dag,
  pos_confounders = integer(0)
)

true_dag
################################## Data pre-processing finishes #######################################

dat<-mat

# n: sample size
n<-dim(dat)[1]

# Set exceedance levels.
q_vec_AAC<-seq(0.001,0.2,by=0.005)
k_vec_AAC<-round(n*q_vec_AAC)

# set exceedance
q_vec_ease<-seq(0.001,0.2,by=0.005)
k_vec_ease<-round(n*q_vec_ease)

round(n^0.7)
floor(n^0.4)

k_vec_AAC<-k_vec_ease<-c(k_vec_AAC[1:2],29,k_vec_AAC[3:16],366,k_vec_AAC[17:length(k_vec_AAC)])

# Marginal transformation parameter.
alpha<-2
# Penalty parameter in the objective function D_k(s,t).
lambda<-4

L_r<-3
dis<-1
r<-1/2

# True direction from site to site.
causal_pair_truth<-matrix(0,nrow=18,ncol=2)
causal_pair_truth[1,]<-c(1,2)
causal_pair_truth[2,]<-c(1,4)
causal_pair_truth[3,]<-c(1,10)
causal_pair_truth[4,]<-c(2,4)
causal_pair_truth[5,]<-c(2,10)
causal_pair_truth[6,]<-c(4,10)
causal_pair_truth[7,]<-c(3,4)
causal_pair_truth[8,]<-c(3,10)
causal_pair_truth[9,]<-c(5,6)
causal_pair_truth[10,]<-c(5,10)
causal_pair_truth[11,]<-c(6,10)
causal_pair_truth[12,]<-c(7,10)
causal_pair_truth[13,]<-c(8,10)
causal_pair_truth[14,]<-c(9,11)
causal_pair_truth[15,]<-c(9,10)
causal_pair_truth[16,]<-c(11,10)
causal_pair_truth[17,]<-c(12,11)
causal_pair_truth[18,]<-c(12,10)

# k_for_fit_GPD_vec<-k_vec_AAC
#k_for_fit_GPD_vec<-pmax(50,k_vec_AAC)
k_for_fit_GPD_vec<-rep(50,length(k_vec_AAC) )

R_NO_transform<-dat
Z_NO_transform<-dat

data<-dat
data<-data.frame(data)

# Marginal transformation with parameter alpha.
R_transform<-apply(dat[,1:dim(dat)[2]],2,rank_top_down)
Z_transform<-(n/R_transform)^(1/alpha)

# Store EASE based on AAC from our paper and CTC from the paper
# "Causal discovery in heavy-tailed models".
# MT: marginal transform; NMT: no marginal transform.
AAC_MT_order<-AAC_NMT_order<-ease_CTC_order<-matrix(0,nrow=dim(causal_pair_truth)[1],ncol=(2+2*length(k_vec_AAC)) )
AAC_MT_order<-as.data.frame(AAC_MT_order)
AAC_NMT_order<-as.data.frame(AAC_NMT_order)
ease_CTC_order<-as.data.frame(ease_CTC_order)

# Construct the column names.
col_names <- c("truth_cause", "truth_effect")

for (k in k_vec_AAC) {
  col_names <- c(col_names, paste0("est_cause_k_", k), paste0("est_effect_k_", k))
}

# Assign the column names to ab_order.
colnames(AAC_MT_order) <-colnames(AAC_NMT_order) <- colnames(ease_CTC_order)<-col_names

AAC_MT_val<-AAC_NMT_val<-ease_CTC_val<-matrix(0,nrow=dim(causal_pair_truth)[1],ncol=2*length(k_vec_AAC) )
AAC_MT_val<-as.data.frame(AAC_MT_val)
AAC_NMT_val<-as.data.frame(AAC_NMT_val)
ease_CTC_val<-as.data.frame(ease_CTC_val)

col_names <- c()

for (k in k_vec_AAC) {
  col_names <- c(col_names, paste0("est_a_k_", k), paste0("est_b_k_", k))
}

colnames(AAC_MT_val)<-colnames(AAC_NMT_val)<-col_names

col_names <- c()

for (k in k_vec_AAC) {
  col_names <- c(col_names, paste0("Gamma_12_k_", k), paste0("Gamma_21_k_", k))
}

colnames(ease_CTC_val)<-col_names

AAC_MT_order[,1:2]<-AAC_NMT_order[,1:2]<-ease_CTC_order[,1:2]<-causal_pair_truth

for(k_ind in 1:length(k_vec_AAC)){
  
  print(paste0("k: ", k_vec_AAC[k_ind]))
  
  for(i in 1:dim(causal_pair_truth)[1]){
    
    Z_NO_transform_sub<-Z_NO_transform[,causal_pair_truth[i,]]
    
    k<-k_vec_AAC[k_ind]
    
    # Fit a generalized Pareto distribution to the upper tail with
    # estimated shape parameter \gamma and scale parameter \sigma.
    tau0 <- 1 - k_for_fit_GPD_vec[k_ind]/dim(data)[1]
    
    u <- mapply(function(col, tau) {
      as.numeric(quantile(col, probs = tau, type = 8))
    }, data, tau0)
    
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
      
      Fval <- 1 - (k / n) * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
      
      
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
    
    for(j in 1:dim(Z_transform)[2]){
      
      Z_transform[exceed_idx[[j]],j]<-MT_upper_u[[j]]
      
    }
    
    Z_transform_sub<-Z_transform[,causal_pair_truth[i,]]
    
    ##############################################################
    
    AAC_NMT_val[i, (2*(k_ind-1)+1):(2*(k_ind-1)+2) ]<-AAC_val_opt(Z=Z_NO_transform_sub,a_b_val=T,k=k,lambda,L_r=4,dis=dis,r=r)[1:2]
    temp_order<-AAC_ease_order(dat=Z_NO_transform_sub, a_b_val=F, k = k, lambda=lambda, L_r=L_r,dis=dis,r=r )
    AAC_NMT_order[i, (2+2*(k_ind-1)+1):(2+2*(k_ind-1)+2) ]<-causal_pair_truth[i,][temp_order]
    
    ease_CTC_val[i, (2*(k_ind-1)+1)]<-causal_tail_coeff(v1=Z_NO_transform_sub[,1],v2=Z_NO_transform_sub[,2], k=k_vec_ease[k_ind], both_tails=F)
    ease_CTC_val[i, (2*(k_ind-1)+2)]<-causal_tail_coeff(v1=Z_NO_transform_sub[,2],v2=Z_NO_transform_sub[,1], k=k_vec_ease[k_ind], both_tails=F)
    temp_order<-ease(Z_NO_transform_sub,k=k_vec_ease[k_ind],both_tails = F)
    ease_CTC_order[i, (2+2*(k_ind-1)+1):(2+2*(k_ind-1)+2) ]<-causal_pair_truth[i,][temp_order]
    
    AAC_MT_val[i, (2*(k_ind-1)+1):(2*(k_ind-1)+2) ]<-AAC_val_opt(Z=Z_transform_sub,a_b_val=T,k=k,lambda,L_r=4,dis=dis,r=r)[1:2]
    temp_order<-AAC_ease_order(dat=Z_transform_sub, a_b_val=F, k = k, lambda=lambda, L_r=L_r,dis=dis,r=r )
    AAC_MT_order[i, (2+2*(k_ind-1)+1):(2+2*(k_ind-1)+2) ]<-causal_pair_truth[i,][temp_order]
    
  }
  
}

fin_res<-list(AAC_NMT_val=AAC_NMT_val,AAC_NMT_order=AAC_NMT_order,
              ease_CTC_val=ease_CTC_val,ease_CTC_order=ease_CTC_order,
              AAC_MT_val=AAC_MT_val, AAC_MT_order=AAC_MT_order,
              lambda=lambda,r=r,k_vec_AAC=k_vec_AAC,k_for_fit_GPD_vec=k_for_fit_GPD_vec, n=n
              ,L_r=3,dis=dis)

file_name<-paste0("datasets/Section_4.2/river_data_pairwise_alpha_2_lambda_",lambda,"_a_0.7_fit_marginal_m_sqrt_n_alpha_",alpha)
save(fin_res,file=file_name)
