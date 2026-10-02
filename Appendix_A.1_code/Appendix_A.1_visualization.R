################################################################################
# Clean workspace and load packages
################################################################################

rm(list = ls())

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

################################################################################
# Result directory
################################################################################

result_dir <- "datasets/Appendix_A.1"

################################################################################
# 1. Find all relevant files
#
# Expected form:
#
# CTC_AAC_DAG_lambda_4_p_5_n_500_shape_1_M_1_50_m_sqrt_n
################################################################################

all_files <- list.files(
  path = result_dir,
  pattern = paste0(
    "^CTC_AAC_DAG_lambda_[0-9.]+",
    "_p_[0-9]+",
    "_n_[0-9]+",
    "_shape_[0-9]+",
    "_M_[0-9]+_[0-9]+",
    "_m_sqrt_n",
    "(\\.RData|\\.rda)?$"
  ),
  full.names = TRUE
)

cat("Total number of files found:", length(all_files), "\n")

head(basename(all_files))


################################################################################
# 2. Parse parameter values from file names
################################################################################

parse_file_info <- function(file_path) {
  
  file_name <- basename(file_path)
  
  pattern <- paste0(
    "^CTC_AAC_DAG_lambda_([0-9.]+)",
    "_p_([0-9]+)",
    "_n_([0-9]+)",
    "_shape_([0-9]+)",
    "_M_([0-9]+)_([0-9]+)",
    "_m_sqrt_n",
    "(?:\\.(?:RData|rda))?$"
  )
  
  parsed <- str_match(file_name, pattern)
  
  if (is.na(parsed[1, 1])) {
    stop(
      paste(
        "File name does not match expected pattern:",
        file_name
      )
    )
  }
  
  data.frame(
    file_path = file_path,
    file_name = file_name,
    lambda = as.numeric(parsed[1, 2]),
    p = as.numeric(parsed[1, 3]),
    n = as.numeric(parsed[1, 4]),
    shape = as.numeric(parsed[1, 5]),
    M_start = as.numeric(parsed[1, 6]),
    M_end = as.numeric(parsed[1, 7]),
    stringsAsFactors = FALSE
  )
}


file_info <- bind_rows(
  lapply(all_files, parse_file_info)
) %>%
  mutate(
    M_num = M_end - M_start + 1
  ) %>%
  arrange(
    lambda,
    n,
    p,
    shape,
    M_start
  )

print(file_info)


################################################################################
# 3. Check the number of files for each lambda
################################################################################

group_check <- file_info %>%
  count(lambda, name = "n_files") %>%
  arrange(lambda)

print(group_check)

cat(
  "\nNumber of distinct lambda values:",
  nrow(group_check),
  "\n"
)

# If:
#
# n     = {500, 1000, 2000}
# p     = {5, 10, 15, 30}
# shape = {1, 3, 5}
#
# there should be:
#
# 3 * 4 * 3 = 36
#
# files for each fixed lambda.

expected_files_per_group <- 3 * 4 * 3

if (any(group_check$n_files != expected_files_per_group)) {
  warning(
    "Some lambda groups do not contain exactly 36 files."
  )
}


################################################################################
# 4. Safely convert data frame / matrix to numeric matrix
################################################################################

to_numeric_matrix <- function(x) {
  
  if (is.data.frame(x)) {
    
    out <- lapply(
      x,
      function(z) as.numeric(as.character(z))
    )
    
    out <- as.matrix(as.data.frame(out))
    
    rownames(out) <- rownames(x)
    colnames(out) <- colnames(x)
    
  } else {
    
    x <- as.matrix(x)
    
    out <- matrix(
      as.numeric(x),
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = dimnames(x)
    )
  }
  
  return(out)
}


################################################################################
# 5. Elementwise mean across Monte Carlo repetitions
################################################################################

cal_mean <- function(mat_list, round_digit = NULL) {
  
  mat_list <- lapply(
    mat_list,
    to_numeric_matrix
  )
  
  reference_dim <- dim(mat_list[[1]])
  
  same_dim <- all(
    vapply(
      mat_list,
      function(x) identical(dim(x), reference_dim),
      logical(1)
    )
  )
  
  if (!same_dim) {
    stop("Matrices within the list have different dimensions.")
  }
  
  if (length(mat_list) == 1) {
    
    mean_matrix <- mat_list[[1]]
    
  } else {
    
    array_data <- simplify2array(mat_list)
    
    mean_matrix <- apply(
      array_data,
      c(1, 2),
      mean,
      na.rm = TRUE
    )
  }
  
  if (!is.null(round_digit)) {
    mean_matrix <- round(mean_matrix, round_digit)
  }
  
  mean_df <- as.data.frame(mean_matrix)
  
  rownames(mean_df) <- rownames(mat_list[[1]])
  colnames(mean_df) <- colnames(mat_list[[1]])
  
  return(mean_df)
}


################################################################################
# 6. Process one file
################################################################################

process_one_file <- function(file_path) {
  
  env <- new.env()
  
  load(file_path, envir = env)
  
  if (!exists("fin_res", envir = env)) {
    stop(
      paste(
        "fin_res not found in:",
        file_path
      )
    )
  }
  
  fin_res <- env$fin_res
  
  AAC_mean <- cal_mean(
    fin_res$error_rate_list$temp_rec_AAC
  )
  
  CTC_mean <- cal_mean(
    fin_res$error_rate_list$temp_rec_CTC
  )
  
  return(
    list(
      AAC = AAC_mean,
      CTC = CTC_mean
    )
  )
}


################################################################################
# 7. Process all files
################################################################################

all_results <- vector(
  "list",
  nrow(file_info)
)

for (i in seq_len(nrow(file_info))) {
  
  if (i %% 20 == 0 || i == 1 || i == nrow(file_info)) {
    
    cat(
      "Processing file",
      i,
      "of",
      nrow(file_info),
      "\n"
    )
  }
  
  all_results[[i]] <- process_one_file(
    file_info$file_path[i]
  )
}

file_info$result_id <- seq_len(nrow(file_info))


################################################################################
# 8. Elementwise weighted mean across files
################################################################################

elementwise_weighted_mean <- function(
    mat_list,
    weights = NULL,
    round_digit = NULL) {
  
  mat_list <- lapply(
    mat_list,
    to_numeric_matrix
  )
  
  if (is.null(weights)) {
    weights <- rep(1, length(mat_list))
  }
  
  if (length(weights) != length(mat_list)) {
    stop("Length of weights must equal length of mat_list.")
  }
  
  reference_dim <- dim(mat_list[[1]])
  
  same_dim <- all(
    vapply(
      mat_list,
      function(x) identical(dim(x), reference_dim),
      logical(1)
    )
  )
  
  if (!same_dim) {
    stop("Matrices across files have different dimensions.")
  }
  
  if (length(mat_list) == 1) {
    
    mean_matrix <- mat_list[[1]]
    
  } else {
    
    array_data <- simplify2array(mat_list)
    
    mean_matrix <- apply(
      array_data,
      c(1, 2),
      function(x) {
        weighted.mean(
          x,
          w = weights,
          na.rm = TRUE
        )
      }
    )
  }
  
  if (!is.null(round_digit)) {
    mean_matrix <- round(mean_matrix, round_digit)
  }
  
  mean_df <- as.data.frame(mean_matrix)
  
  rownames(mean_df) <- rownames(mat_list[[1]])
  colnames(mean_df) <- colnames(mat_list[[1]])
  
  return(mean_df)
}


################################################################################
# 9. Average separately for each fixed lambda
################################################################################

parameter_groups <- file_info %>%
  distinct(lambda) %>%
  arrange(lambda)

final_results <- vector(
  "list",
  nrow(parameter_groups)
)

names(final_results) <- paste0(
  "lambda_",
  parameter_groups$lambda
)


for (g in seq_len(nrow(parameter_groups))) {
  
  lambda_g <- parameter_groups$lambda[g]
  
  ind <- which(
    file_info$lambda == lambda_g
  )
  
  cat(
    "Averaging lambda =",
    lambda_g,
    ":",
    length(ind),
    "files\n"
  )
  
  AAC_list_g <- lapply(
    ind,
    function(j) all_results[[j]]$AAC
  )
  
  CTC_list_g <- lapply(
    ind,
    function(j) all_results[[j]]$CTC
  )
  
  weights_g <- file_info$M_num[ind]
  
  AAC_mean_g <- elementwise_weighted_mean(
    AAC_list_g,
    weights = weights_g
  )
  
  CTC_mean_g <- elementwise_weighted_mean(
    CTC_list_g,
    weights = weights_g
  )
  
  group_name <- paste0(
    "lambda_",
    lambda_g
  )
  
  final_results[[group_name]] <- list(
    lambda = lambda_g,
    AAC_mean_df = AAC_mean_g,
    CTC_mean_df = CTC_mean_g,
    files_used = file_info$file_name[ind],
    parameter_table = file_info[ind, ]
  )
}


################################################################################
# 10. Convert averaged results to plotting data
#
# k_n = n^a
#
# a = 0.20, 0.25, ..., 0.80
################################################################################

exponent_grid <- seq(
  0.20,
  0.80,
  by = 0.05
)

plot_list <- list()

counter <- 1

for (g in seq_len(nrow(parameter_groups))) {
  
  lambda_g <- parameter_groups$lambda[g]
  
  group_name <- paste0(
    "lambda_",
    lambda_g
  )
  
  AAC_mat <- as.matrix(
    final_results[[group_name]]$AAC_mean_df
  )
  
  CTC_mat <- as.matrix(
    final_results[[group_name]]$CTC_mean_df
  )
  
  if (ncol(AAC_mat) != length(exponent_grid)) {
    stop(
      paste(
        "Number of AAC columns does not match exponent grid for",
        group_name
      )
    )
  }
  
  AAC_average <- colMeans(
    AAC_mat,
    na.rm = TRUE
  )
  
  CTC_average <- colMeans(
    CTC_mat,
    na.rm = TRUE
  )
  
  plot_list[[counter]] <- data.frame(
    lambda = lambda_g,
    exponent = exponent_grid,
    AAC_error_rate = AAC_average,
    CTC_error_rate = CTC_average
  )
  
  counter <- counter + 1
}


plot_df <- bind_rows(plot_list)

print(plot_df)


plot_df_lambda4 <- plot_df %>%
  filter(lambda == 4)

p_AAC <- ggplot(
  plot_df_lambda4,
  aes(
    x = exponent,
    y = AAC_error_rate
  )
) +
  geom_line(
    linewidth = 1,
    color = "red"
  ) +
  geom_point(
    size = 2,
    color = "red"
  ) +
  geom_vline(
    xintercept = 0.7,
    linetype = "dashed",
    linewidth = 1,
    color = "darkgray"
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 0.8, by = 0.1)
  ) +
  labs(
    x = expression("Exponent " * a),
    y = "Average ancestral violation rate"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    panel.grid.minor = element_blank()
  )

p_AAC

ggsave(
  filename = "plots/Appendix_A.1_AAC_error_rate_lambda_4_alpha_2.pdf",
  plot = p_AAC,
  width = 6.5,
  height = 4.5,
  units = "in"
)


