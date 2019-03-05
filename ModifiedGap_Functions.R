# Functions for calculation of modified gap statistic.

## Generate gap statistics for simulated uniform distributions.

# Function to apply clusGap to a dataframe and return output as a tibble.
clusFunc <- function(df, K_max) {
  mat <- as.matrix(df)
  out <- cluster::clusGap(mat, kmeans, K.max = K_max, B = 50, nstart = 20, verbose = FALSE)
  as_tibble(out$Tab)
}

# Function to generate gap values for data simulated from a uniform reference distribution
gap_uniform_sim <- function(n_cells = 20, K.max = 6, n_sim = 20, val_min = 0, val_max = 1){
  
  sim_data <- matrix(rep(runif(n_cells, value_min, value_max), n_sim), nrow = n_cells, ncol = n_sim)
  sim_gaps <- matrix(nrow = n_sim, ncol = K.max)
  
  sim_tib <- tibble(sim_num = c(1:n_sim), num_cells = n_cells, v_min = value_min, v_max = value_max)
  
  # Make random numbers for each simulation, then calculate gap statistics for each set of random numbers, then calculate dela gap values.
  sim_tib <- mutate(sim_tib,
                    data = pmap(list(num_cells, v_min, v_max), runif),
                    gap_stats = map(data, clusFunc),
                    gaps = map(gap_stats, ~.$gap),
                    gap_delta = map(gaps, diff))
}


# Function to calculate percentile cut offs.
# Extract delta gap values into a matrix with columns for each k. Then calculate cut offs.
# df should be a frame containing delta gap outputs, e.g. sim_tib$gap_delta
calc_cutoffs <- function(df, cut_off = 0.9){
  n_sim = length(df)
  n_col = length(df[[1]])
  
  delta_gap_mat <- matrix(nrow = n_sim, ncol = n_col)
  for (i in 1:n_col) {
    delta_gap_mat[, i] <- sapply(df, function(x)
      x[[i]])
  }
  
  cut_off_vals <- vector(length = n_col)
  for (i in 1:n_col) {
    cut_off_vals[i] <- quantile(delta_gap_mat[, i], probs = (cut_off))
  }
  
  cut_off_vals
}
