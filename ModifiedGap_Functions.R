# Functions for calculation of modified gap statistic.

clusGap_Extra <- function(df, K_max = 6){
  clu_df <- clusFunc(df, K_max = K.max)
  
  clu_df$clus_num <- c(1:K.max)
  clu_df <- mutate(clu_df,
                     gap_min = gap - SE.sim,
                     gap_max = gap + SE.sim)
}


## Generate gap statistics for simulated uniform distributions.

# Function to apply clusGap to a dataframe and return output as a tibble.
clusFunc <- function(df, K_max) {
  df <- scale(df) # Center data. This shouldn't be necessary for 1d data.
  mat <- as.matrix(df)
  out <- cluster::clusGap(mat, kmeans, K.max = K_max, B = 50, nstart = 20, verbose = FALSE)
  as_tibble(out$Tab)
}

# Function to generate gap values for data simulated from a uniform reference distribution
gap_uniform_sim <- function(n_cells = 20, K_max = 6, n_sim = 20, value_min = 0, value_max = 1){
  
    sim_tib <- tibble(sim_num = c(1:n_sim), num_cells = n_cells, v_min = value_min, v_max = value_max)
    
    # Make random numbers for each simulation, then calculate gap statistics for each set of random numbers, then calculate dela gap values.
    sim_tib <- mutate(sim_tib,
                      data = pmap(list(num_cells, v_min, v_max), runif),
                      gap_stats = map2(data, K_max, clusFunc),
                      gaps = map(gap_stats, ~.$gap),
                      gap_delta = map(gaps, diff))
  }

# Function to calculate percentile cut offs.
# Extract delta gap values into a matrix with columns for each k. Then calculate cut offs.
# df should be a frame containing delta gap outputs, e.g. sim_tib$gap_delta
calc_thresholds <- function(df, cut_off = 0.9){
  n_sim = length(df)
  n_col = length(df[[1]])
  
  delta_gap_mat <- matrix(nrow = n_sim, ncol = n_col)
  for (i in 1:n_col) {
    delta_gap_mat[, i] <- sapply(df, function(x)
      x[[i]])
  }
  
  threshold_vals <- vector(length = n_col)
  for (i in 1:n_col) {
    threshold_vals[i] <- quantile(delta_gap_mat[, i], probs = (cut_off))
  }
  
  cut_off_vals
}


# Function to plot LogW, ELogW versus cluster number.
logW_plot <- function(df) {
  ggplot(df) +
    geom_point(aes(clus_num, logW), colour = "blue") +
    geom_point(aes(clus_num, E.logW), colour = "red")
}


# Function to plot gap statistic versus cluster number.
gap_plot <- function(df) {
  ggplot(df, aes(clus_num, gap)) +
    geom_point() +
    geom_errorbar(aes(ymin = gap_min, ymax = gap_max))
}

# Function to plot delta gap as a function of cluster number.
diff_plot <- function(df){
  ggplot(df) +
    geom_point(aes(clus_num, gap_diff), colour = "red") +
    geom_point(aes(clus_num, gap_thresh), colour = "blue")
}