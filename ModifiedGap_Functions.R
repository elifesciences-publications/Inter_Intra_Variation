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
clusFunc <- function(df, K_max = 8) {
  df <- scale(df) # Center data. This shouldn't be necessary for 1d data.
  mat <- as.matrix(df)
  out <- cluster::clusGap(mat,
                          kmeans,
                          K.max = K_max,
                          B = 50,
                          d.power = 2,
                          spaceH0 = "original",
                          nstart = 20,
                          verbose = FALSE)
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
calc_thresholds <- function(df, cut_off = 0.99){
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

# Function to apply clusGap to the complete dataset.
# df should be all_data_r$data[[n]]
clusGap_AllData <- function(df, sims, K_max = 8) {
  sims <- sims$threshold_vals
  df <- select(df, vm:fi) %>%
    gather("property", "value", vm:fi) %>%
    group_by(property) %>%
    nest() %>%
    mutate(gap = map2(data, K_max, clusGap_Extra_helper),
           gap_diff = pmap(list(gap, K_max, sims), diff_calc),
           K_est = map_dbl(gap_diff, calc_K_est),
           Gap_best = map_dbl(gap, calc_K_gap)) %>%
    select(-data)
}


# Function to apply clusGap to a simulated dataset.
# df should be sim_data$data[[n]]
clusGap_SimData <- function(df, sims, K_max = 8) {
  sims <- sims$threshold_vals
  df <- select(df, x2) %>%
    mutate(gap = map2(data, K_max, clusGap_Extra_helper),
           gap_diff = pmap(list(gap, K_max, sims), diff_calc),
           K_est = map_dbl(gap_diff, calc_K_est),
           Gap_best = map_dbl(gap, calc_K_gap)) %>%
    select(-data)
}

# Function to return gap thresholds for a given n
return_gap_thresholds <- function(n_cells, thresh_ref_df = sim_tib) {
  filter(thresh_ref_df, cell_count == n_cells)
}

return_data <- function(df, mouse){
  data_to_return <- filter(df, id == mouse) %>% select(data) %>% unlist(recursive = FALSE)
  data_to_return <- data_to_return$data
}

# Function to plot data for a given feature and if K_est>1 to add colours to indicate cluster identity.
data_plot <- function(df, mouse, feature){
  # Extract data to plot into a data frame.
  data_to_plot <- filter(df, id == mouse) %>% select(data) %>% unlist(recursive = FALSE)
  data_to_plot <- data_to_plot$data
  
  # Extract K_est.
  K_est <- filter(df, id == mouse) %>%
    select(clusGap)
  K_est <- filter(K_est[[1]][[1]], property == feature) %>%
    select(K_est) %>%
    unlist()

  
  # If K_est > 1, assign clusters and colour code
  if (K_est > 1) {
    data_to_cluster <- data_to_plot %>% select(!! rlang::sym(feature))
    # Need to repeat K means clustering with  settings used by clusGap::cluster.
    # Possible caveate here is that clusters may still be different.
    km <- kmeans(as.matrix(data_to_cluster), centers = K_est, nstart = 20)
    data_to_plot$cluster <- as.factor(km$cluster)
  } else {
    data_to_plot$cluster <- 1
  }
  
  ggplot(data_to_plot, aes(dvloc,!!rlang::sym(feature)), colour = cluster) +
    geom_point(aes(colour = cluster)) +
    geom_rug(size = 0.2)
}


# Function to plot LogW, ELogW versus cluster number.
logW_plot <- function(df, mouse, feature) {
  data_to_plot <- filter(df, id == mouse) %>% select(clusGap) %>% unlist(recursive = FALSE)
  data_to_plot <- filter(data_to_plot[[1]], property == feature) %>% select(gap) %>% unlist(recursive = FALSE)
  data_to_plot <- select(data_to_plot[[1]], clus_num, logW, E.logW) %>%
    gather("Group", "logW", 2:3)
  ggplot(data_to_plot) +
    geom_point(aes(clus_num, logW, colour = Group)) +
    labs(x = "K", y = "logW") +
    theme(legend.title = element_blank())
}

# Function to plot gap statistic versus cluster number.
gap_plot <- function(df, mouse, feature) {
  data_to_plot <- filter(df, id == mouse) %>% select(clusGap) %>% unlist(recursive = FALSE)
  data_to_plot <- filter(data_to_plot[[1]], property == feature) %>% select(gap) %>% unlist(recursive = FALSE)
  ggplot(data_to_plot[[1]], aes(clus_num, gap)) +
    geom_point() +
    geom_errorbar(aes(ymin = gap_min, ymax = gap_max)) +
    labs(x = "K", y = "Gap")
}

# Function to plot delta gap as a function of cluster number.
diff_plot <- function(df, mouse, feature){
  data_to_plot <- filter(df, id == mouse) %>% select(clusGap) %>% unlist(recursive = FALSE)
  data_to_plot <- filter(data_to_plot[[1]], property == feature) %>% select(gap_diff) %>% unlist(recursive = FALSE)
  data_to_plot <- select(data_to_plot[[1]], clus_num, gap_diff, gap_threshold) %>%
    gather("Group", "Gap_Diff", 2:3)
  ggplot(data_to_plot) +
    geom_point(aes(clus_num, Gap_Diff, colour = Group))  +
    labs(x = "K", y = "Gap diff.") +
    theme(legend.title = element_blank())
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
# df is test_data_r$gap_diff, etc.
calc_K_est <- function(df){
  comp <- which(df$gap_diff > df$gap_threshold)
  if (sum(comp)> 0) {
    K_est <- min(comp) + 1
  } else {
    K_est <- 1
  }
  K_est
}

# Returns K identified with the conventional gap statistic procedure.
# df is all_data_r$clusGap[[1]]$gap[[1]], etc.
calc_K_gap <- function(df) {
  comp <- which.max(df$gap)
  df$clus_num[[comp]]
}
