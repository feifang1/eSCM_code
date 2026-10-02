
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
library(eva)


###################### Preprocess the dataset based on the code and datasets from the paper "EXTREMES ON RIVER NETWORKS" ########################
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
  ## convert a dataframe to a tibble and rename the column "Date" to "date"
  
  dfr %>%
    mutate(Date = as_date(Date)) %>%
    rename(date = Date) %>%
    as_tibble()
}

PP.lik.linear <- function(Params, Data, u, NoYears, CovarMu, CovarSc, CovarXi){
  ## Function from the paper "EXTREMES ON RIVER NETWORKS"
  ## by Peiman Asadi, Anthony C. Davison, and Sebastian Engelke
  
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
  ## by Peiman Asadi, Anthony C. Davison, and Sebastian Engelke
  
  x <- optim(Init, PP.lik.linear,
             Data = Data, u = u, NoYears = NoYears, CovarMu = CovarMu,
             CovarSc = CovarSc, CovarXi = CovarXi,
             method = method, control = control, hessian = TRUE)
  output <- x
  return(output)
  
}

split_stations <- function(dat, year = NULL, stations = NULL){
  ## tibble integer character_vector -> list
  ## produce a named list with discharges for the stations in a given year
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

## Import datasets ####
getwd()
load("datasets/River_raw_data_ease_algorithm/StsTSs.RData")
load("datasets/River_raw_data_ease_algorithm/StsInfo.RData")

class(StsTSs[[1]])
StsTSs[[1]][1:10,]

# Select the 31 stations that are in the AOAS paper plus station 19 (renamed 32)
StsChos <- c(c(1:47)[-c(16,30,31,34,42,43,44,45,46,47,3,1,2,29,18,19)], 19)

NoSt <- length(StsChos)

StsTSsChos <- StsTSs[StsChos]
StsInfoChos <- StsInfo[StsChos,] %>%
  mutate(id_old = StsChos,
         id = 1:NoSt)

# clean the data and find common dates
river_dat <- map(.x = StsTSsChos, .f = dfr2tibble) %>%
  reduce(.f = inner_join, by = "date") %>%
  dplyr::filter(month(date) %in% c(6, 7, 8),
                year(date) < 2010)

# rename columns (according to the AOAS paper)
colnames(river_dat)[2:(NoSt + 1)] <- sprintf("station_%02d", 1:NoSt)

# select stations
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)
river_dat <- river_dat[, c(1, station_names + 1)]
NoSt <- length(station_names)

station_info <- StsInfoChos %>%
  dplyr::filter(id %in% station_names) %>%
  rename(name = RivNames,
         lat = Lat,
         lon = Long,
         ave_vol = AveVol)

# clean the workspace
rm(StsInfo, StsInfoChos, StsTSs, StsTSsChos, StsChos)

## Check model assumptions ####
method <- "BFGS"
control <- list(maxit = 5000, reltol=10^(-30), abstol=0.0001, trace=0)

# Take the maximum for each station within a year.
SummerMaxima <- river_dat %>%
  mutate(year = year(date)) %>%
  gather(key = "station", value = "disc", -date, -year) %>%
  group_by(year, station) %>%
  summarise(maxDisc = max(disc)) %>%
  spread(key = station, value = maxDisc) %>%
  select(year, sprintf("station_%02d", station_names))

# GEV analysis
GevPars <- tibble(station = character(), location = double(),
                  scale = double(), shape = double())

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


# Compute threshold data
q <- .9
threshold_dat <- river_dat %>%
  gather(key = "station", value = "disc", -date) %>%
  mutate(station = factor(station,
                          levels = sprintf("station_%02d", station_names))) %>%
  group_by(station) %>%
  mutate(thres = quantile(disc, q)) %>%
  dplyr::filter(disc > thres) %>%
  mutate(noyears = length(unique(year(date))))

AllEvents <- split_stations(threshold_dat)
Threshold <- threshold_dat %>% distinct(station, thres) %>% deframe()
NoOfYears <- threshold_dat %>% distinct(station, noyears) %>% deframe()


# Fit a GEVD to each station (by maximizing the joint Poisson
# process likelihood; cf. formula (21) in the paper by
# Asadi, Davison, and Engelke (2015))
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



# Regionalized model from the paper by Asadi, Davison, and Engelke (2015)
Grp1 <- sprintf("station_%02d", c(11, 19, 21))
Grp2 <- sprintf("station_%02d", c(13, 28, 32))
Grp3 <- sprintf("station_%02d", c(1, 7, 9, 14))
Grp4 <- sprintf("station_%02d", c(23, 26))
Grp <- list(Grp1, Grp2, Grp3, Grp4)

# Estimate parameters: shape, location, and scale
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



# Compute adjusted standard errors taking into account time dependence
# (see Fawcett and Walshaw, 2007 and 2012)
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
    
    # Initialize parameters
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
    
    # Compute the gradient at the optimal parameter value
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

ParsCovModel  %>% left_join(se_adj) %>%
  mutate(lb = round(shape - 2 * shape_se_adj, 3),
         ub = round(shape + 2 * shape_se_adj, 3)) %>%
  group_by(region) %>%
  summarise(shape = mean(shape),
            lb = mean(lb),
            ub = mean(ub)) %>%
  print()
cat("\n\n\n")
#sink()

# Encode the true DAG for river data
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

mat <- river_dat %>%  select(-date) %>% as.matrix()

# ground truth
ll <- list(
  dataset = mat,
  dag = true_dag,
  pos_confounders = integer(0)
)
################################## data pre-processing finishes #######################################

dat<-mat

# n: sample size
n<-dim(dat)[1]

# set exceedance
q_vec_AAC<-seq(0.001,0.2,by=0.005)
k_vec_AAC<-round(n*q_vec_AAC)

# set exceedance
q_vec_ease<-seq(0.001,0.2,by=0.005)
k_vec_ease<-round(n*q_vec_ease)

round(n^0.7)
floor(n^0.4)
k_vec_AAC<-k_vec_ease<-c(k_vec_AAC[1:2],29,k_vec_AAC[3:16],366,k_vec_AAC[17:length(k_vec_AAC)])

alpha<-2
lambda<-4

L_r<-3
dis<-1
r<-1/2

R_NO_transform<-dat
Z_NO_transform<-dat

R_transform<-apply(dat[,1:dim(dat)[2]],2,rank_top_down)
Z_transform<-(n/R_transform)^(1/alpha)

ease_error_rec<-AAC_error_transform_rec<-AAC_error_NO_transform_rec<-rep(0,length(k_vec_AAC))
ease_order_rec<-AAC_order_transform_rec<-AAC_order_NO_transform_rec<-matrix(0,nrow=length(k_vec_AAC), ncol=dim(Z_NO_transform)[2])

k_fit_tail_AAC<-rep(round(n^(1/2)),length(k_vec_AAC) )
k_fit_tail_ease<-rep(round(n^(1/2)),length(k_vec_ease) )

for(i in 1:length(k_vec_AAC)){
  
  print(paste0("iteration i:",i))
  
  #### AAC parametric fit for the tail
  k_AAC<-k_vec_AAC[i]
  
  tau0_AAC <- 1 - k_fit_tail_AAC[i]/dim(R_NO_transform)[1]
  
  u <- sapply(data.frame(R_NO_transform), function(col) {
    as.numeric(quantile(col, probs = tau0_AAC, type = 8))
  })
  
  fits <- mapply(
    function(col, thr) {
      fit.gpd(col, threshold = thr, method = "Grimshaw")
    },
    as.data.frame(R_NO_transform),
    u,
    SIMPLIFY = FALSE
  )
  
  params <- lapply(fits, coef)
  sigma_u_hat <- sapply(params, function(p) unname(p["scale"]))
  gamma_u_hat <- sapply(params, function(p) unname(p["shape"]))
  
  exceed_idx <- mapply(
    function(col, thr) which(col > thr),
    as.data.frame(R_NO_transform),
    u,
    SIMPLIFY = FALSE
  )
  
  p <- ncol(R_NO_transform)
  n <- nrow(R_NO_transform)
  
  F_x_u <- vector("list", p)
  MT_upper_u <- vector("list", p)
  
  for (j in seq_len(p)) {
    
    idx <- exceed_idx[[j]]   # indices for column j
    xexc <- R_NO_transform[idx, j]     # exceedances
    
    z <- (xexc - u[j]) / sigma_u_hat[j]
    
    Fval <- 1 - (k_fit_tail_AAC[i]/ n) * (1 + gamma_u_hat[j] * z)^(-1 / gamma_u_hat[j])
    
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
  
  Z_transform_AAC<-Z_transform
  
  for(j in 1:dim(Z_transform_AAC)[2]){
    
    Z_transform_AAC[exceed_idx[[j]],j]<-MT_upper_u[[j]]
    
  }

  AAC_order_NO_transform_rec[i,]<-AAC_ease_order(dat=Z_NO_transform, a_b_val=F, k = k_vec_AAC[i] , lambda=lambda, L_r=L_r, dis=dis, r=r )
  AAC_error_NO_transform_rec[i]<-order_error_measure(dag=true_dag,ab_dir=AAC_order_NO_transform_rec[i,])
  
  AAC_order_transform_rec[i,]<-AAC_ease_order(dat=Z_transform_AAC, a_b_val=F, k = k_vec_AAC[i], lambda=lambda, L_r=L_r, dis=dis, r=r )
  
  AAC_error_transform_rec[i]<-order_error_measure(dag=true_dag,ab_dir=AAC_order_transform_rec[i,])
  
  # function "ease" is based on the paper "CAUSAL DISCOVERY IN HEAVY-TAILED MODELS"
  # and their corresponding codes for function "causal_tail_matrix" and "causal_tail_coeff"
  ease_order_rec[i,]<-ease(Z_NO_transform,k=k_vec_ease[i],both_tails = F)
  ease_error_rec[i]<-order_error_measure(dag=true_dag,ab_dir=ease_order_rec[i,])
  
}

res_rec<-list(AAC_order_NO_transform_rec=AAC_order_NO_transform_rec
              ,AAC_error_NO_transform_rec=AAC_error_NO_transform_rec
              ,AAC_order_transform_rec=AAC_order_transform_rec
              ,AAC_error_transform_rec=AAC_error_transform_rec
              ,ease_order_rec=ease_order_rec
              ,ease_error_rec=ease_error_rec
              ,lambda=lambda
              ,alpha=alpha
              ,L=L_r
              ,d=dis
              ,r=r)

file_name<-paste0("datasets/Section_4.2/river_data_whole_causal_order_lambda_",lambda,"_a_0.7_parametric_fit_upper_fixed_cut_m_sqrt_n_alpha_",alpha )
save(res_rec,file=file_name)


