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
# gap_delta a frame containing delta gap outputs found within df, e.g. sim_tib$sims[[n]]$gap_delta
calc_thresholds <- function(df, cut_off = 0.9){
  df_gd <- df$gap_delta
  n_sim = length(df_gd)
  n_col = length(df_gd[[1]])
  
  delta_gap_mat <- matrix(nrow = n_sim, ncol = n_col)
  for (i in 1:n_col) {
    delta_gap_mat[, i] <- sapply(df_gd, function(x)
      x[[i]])
  }
  
  thresh_vals <- vector(length = n_col)
  for (i in 1:n_col) {
    thresh_vals[i] <- quantile(delta_gap_mat[, i], probs = (cut_off))
  }
  
  thresh_vals
}

# Function to plot data for a given feature and if K_est>1 to add cluster colouring
data_plot <- function(df, feature_name){
  data_to_plot <- filter(df, property == feature_name)
  ggplot(data_to_plot, aes(dvloc, value)) +
    geom_point()
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

# Helper function for accessing clusGap using purr::map.
clusGap_Extra_helper <- function(df, K.max){
  df <- na.omit(df) # remove Nas
  clusGap_Extra(df$value, K_max = K.max) # Cluster only the value of the property
}

# Function for calculating gap differences.
# df should be test_data_r$gap
diff_calc <- function(df, K.max, tvs){
  gd <- tibble(clus_num = c(2:(K.max)), gap_diff = diff(df$gap))
  gd$gap_threshold <- tvs
  gd
}

# Function for calculating K_est.
# df is test_data_r$gap_diff
calc_K_est <- function(df){
  comp <- which(df$gap_diff > df$gap_threshold)
  if (sum(comp)> 0) {
    K_est <- min(comp)
  } else {
    K_est <- 1
  }
  K_est
}
