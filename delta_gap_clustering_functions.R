
# Calculate the intra-cluster dispersions of data using mean Euclidean pairwise distances
# between cluster values. Returns a vector of dispersions (one element for each cluster) 
get_dispersion <- function(data_vec, cluster_labels, cluster_centers, d.power){
  # NB: use d.power = 2 to square Euclidean distances as in  Tibshirani et. al. 2001
  
  dispersions <- numeric(length(unique(cluster_labels))) # Initialise vector of dispersions
  
  for (k in 1:length(unique(cluster_labels))){
    # Get only the data in the relevant cluster (i.e. all the data with the same labels)
    data_in_cluster <- data_vec[cluster_labels == k]
    
    # Note: to follow the Tibshirani paper this code does NOT normalise by number of pairwise 
    # distances. Instead, it normalises by number of values in cluster. 
    # Note: the 'division' factor (2) shown in Tibshirani et. al equation is already implicity 
    # present: in Tibshirani et. al. the factor to divide by to convert pairwise distances to 
    # distances to cluster centres is not required here as only 'half' (the bottom triangle)
    # of pairwise distances are returned by the dist() function.
    dispersions[k] <- sum(dist(data_in_cluster, method="euclidean")^d.power)/
      length(data_in_cluster)
  }
  return(dispersions)
}

# Simulate datasets by drawing from uniform distributions on the interval [0, 1]. These 
# simulated datasets are NOT scaled so that the min=0 and the max=1
gen_unif_sims <- function(sims_per_dataset_size=1000, sim_dataset_sizes, k.max=8, d.power){
  
  clust_vec <- 1:k.max
  
  # Initialise storage arrays
  delta_gaps_store <- array(numeric(), 
                        dim=c(sims_per_dataset_size, 
                              k.max-1, 
                              length(sim_dataset_sizes)))
  
  dispersion_store <- array(numeric(), dim=c(sims_per_dataset_size, 
                                             k.max, 
                                             length(sim_dataset_sizes)))
  
  # Calculate dispersions for different dataset sizes separately
  for (n_data_pts in 1:length(sim_dataset_sizes)){
    dispersions_sim_data <- array(numeric(), dim=c(sims_per_dataset_size, k.max))
    
    # Generate simulated dataset
    for (sim in 1:sims_per_dataset_size){
      # Note: sample from uniform distribution is not set to scale between min=0, max=1
      sim_data <- scale_vec(runif(sim_dataset_sizes[n_data_pts]))
      
      for (k in clust_vec){
        kmeans_results <- kmeans(sim_data, centers=k, nstart=50, iter.max=100)
        dispersions_for_k <- get_dispersion(sim_data, kmeans_results$cluster,
                                                kmeans_results$centers, d.power=d.power)
        dispersions_sim_data[sim, k] <- log(mean(dispersions_for_k))
      }
    }
    
    # Subtract the mean (average dispersion) from each column to obtain the dispersion 'gaps'
    # for each dataset size.
    gaps_sim_data <- sweep(dispersions_sim_data, 2, colMeans(dispersions_sim_data), "-")
    # Subtract consecutive k dispersions to get the delta gaps (difference of consecutive gaps)
    delta_gaps_sim_data <- t(diff(t(gaps_sim_data)))
    # Store delta gaps and dispersions
    delta_gaps_store[,,n_data_pts] <- delta_gaps_sim_data
    dispersion_store[,,n_data_pts] <- dispersions_sim_data
  }
  # Save the file in the current working directory
  save(delta_gaps_store, dispersion_store, file="ModifiedGap/unif_sims_delta_gaps_and_dispersions.Rda")
}

# Generate simulated multimodal data according to specified paramenters and estimate 
# k_est for each generated dataset. Can take a *long* time to run, depending on the number
# factor combinations
get_mm_sim_k_ests <- function(n_mm_sims, k_vec, n_data_vec, 
                                  sep_sd_vec, threshold_fit_results, dispersion_fit_results){
  k.max <- max(k_vec)
  #  Initialise holding 4d array with named dimensions
  k_estimates_mm <- array(numeric(), dim=c(n_mm_sims, length(k_vec), 
                                           length(n_data_vec), length(sep_sd_vec)))
  dimnames(k_estimates_mm) <- list(sims=as.character(1:dim(k_estimates_mm)[[1]]),
                                   k=as.character(k_vec), n_data=as.character(n_data_vec),
                                   sd_sep=as.character(sep_sd_vec))
  
  # Evaluate k_est for each combination of factors
  for (s in 1:length(sep_sd_vec)){
    for (d in 1:length(n_data_vec)){
      # Use the fitted parameters to obtain thresholds and dispersions. Make up some data 
      # because the function expects a dataset and not a number of data points. Do it here
      # so we don't have to do it in the middle of the loop
      thresholds <- get_thresh_from_params(runif(n_data_vec[d]), k.max,
                                           threshold_fit_results$thresh_params[[1]])
      dispersions <- get_dispersions_from_params(runif(n_data_vec[d]), k.max,
                                                 dispersion_fit_results$dispersion_params[[1]])
      for (k in 1:length(k_vec)){
        for (i in 1:n_mm_sims){
          data <- get_multimodal_data_sample(n_data_vec[d], k_vec[k], sep_sd_vec[s])
          
          cluster_eval <- get_k_est(data, kmeans, nstart=50, iter.max=100, 
                                                  K.max=k.max, B=NULL, d.power=2,
                                                  thresholds=thresholds, dispersions=dispersions)
          k_estimates_mm[i, k, d, s] <- cluster_eval$k_est
        }
      }
      
    }
  }
  
  # Save the file in the specified working directory
  save(k_estimates_mm, file="ModifiedGap/k_estimates_mm.Rda")
  return(k_estimates_mm)
}

# Given a set of delta gaps across different k for simulated datasets, return the delta gaps 
# that constitute the cutoff threshold for significance
get_delta_gaps_thresh_vals <- function(delta_gaps_store, threshold_criterion=0.01){
  
  dims_delta_gaps <- dim(delta_gaps_store)
  n_sims <- dims_delta_gaps[1]
  n_cluster_pairs <- dims_delta_gaps[2]
  n_dataset_sizes <- dims_delta_gaps[3]
  
  # Initialise array that holds thresholds
  delta_gaps_thresh_vals <- array(numeric(), dim=c(n_dataset_sizes, n_cluster_pairs))
  
  for (dataset_size_ind in 1:n_dataset_sizes){
    # The number of delta gap values that are in the top 'threshold_criterion' percent
    n_exceeding_criterion <- ceiling(threshold_criterion*n_sims)
    # Sort each column (separate k) independently to easily find cutoff value
    delta_gaps_sorted <- apply(delta_gaps_store[,,dataset_size_ind], 2, sort)
    # Find cutoff delta gap value 
    delta_gaps_thresh_vals[dataset_size_ind,] <- delta_gaps_sorted[n_exceeding_criterion,]
  }
  
  delta_gaps_thresholds <- tibble(
    delta_gaps_thresh_vals = list(-delta_gaps_thresh_vals) # Switch to positive values
  )
  
  return(delta_gaps_thresholds)
}

# Calculate the average dispersion across all simulations for each k
get_mean_dispersions <- function(dispersion_store){
  k_dim <- 2
  pts_dim <- 3
  dims <- dim(dispersion_store)
  dispersions <- array(numeric(), dim=c(dims[pts_dim], dims[k_dim]))
  for (k in 1:dims[k_dim]){
    dispersions[,k] <- colMeans(dispersion_store[,k,])
  }
  return(dispersions)
}

# Fit delta_slope_thresh_vals with hyperbolic function
fit_thresholds <- function(delta_gaps_thresholds, dataset_sizes){
  dims <- dim(delta_gaps_thresholds$delta_gaps_thresh_vals[[1]])
  n_fit_parameters <- 3 # Parameters to be estimated for each model
  thresh_params <- array(numeric(), dim=c(dims[2], n_fit_parameters))
  fits <- array(numeric(), dim=dims)
  x <- dataset_sizes
  for (f in 1:dims[2]){
    y <- delta_gaps_thresholds$delta_gaps_thresh_vals[[1]][,f]
    # This is the hyperbolic model used in the MATLAB implementation
    # When d.power is 1 (not recommended), initialise a=3, b=0.8, c=0.04
    # When d.power is 2 (as per Tibshirani et. al), a=10, b=1, c=0.2
    # When d.power is 3 (not recommended), a = 15, b=1.5, c=1
    m <- nls(y ~ a/(x^b)+c, start=list(a=10, b=1, c=0.2))
    fits[ ,f] <- predict(m)
    thresh_params[f, ] <- coef(m)
  }
  
  fit_results <- tibble(
    thresh_params = list(thresh_params),
    threshold_fits = list(as_tibble(fits))
  )
  save(thresh_params, file="ModifiedGap/fitted_thresh_params.Rda")
  return(fit_results)
}  

# Fit dispersions with log function
fit_dispersions <- function(dispersions, dataset_sizes){
  dims <- dim(dispersions)
  dim_k <- 2
  n_fit_parameters <- 2 # number of parameters to be estimated for each model
  dispersion_params <- array(numeric(), dim=c(dims[dim_k], n_fit_parameters))
  dispersion_fits <- array(numeric(), dim=dims)
  x <- dataset_sizes
  for (f in 1:dims[dim_k]){
    y <- dispersions[,f]
    # Use a log model
    m <- nls(y ~ a*log(x)+b, start=list(a=1, b=1))
    dispersion_fits[ ,f] <- predict(m)
    dispersion_params[f, ] <- coef(m)
  }
  
  dispersion_fit_results <- tibble(
    dispersion_params = list(dispersion_params),
    dispersion_fits = list(as_tibble(dispersion_fits))
  )
  save(dispersion_params, file="ModifiedGap/fitted_dispersion_params.Rda")
  return(dispersion_fit_results)
}  

# Scaling function - currently configured to scale data between 0 and 1. Can also be set to
# scale data so that it has mean 0 and standard deviation 1 (note that this is ~50X slower)
scale_vec <- function(input_vector){
  input_vector <- input_vector[!is.na(input_vector)] # Remove NAs
  # Set the scale data between 0 and 1. Only accepts vectors, NOT matrices!!
  vec0 <- input_vector - min(input_vector)
  vec_out <- vec0 * (1/max(vec0))
  
  # vec_out <- scale(input_vector) # Uncomment this to switch scaling methods
  
  return(vec_out)
}

# Reconstruct thresholds from parameters
get_thresh_from_params <- function(data, K.max, thresh_params){
  thresholds <- numeric(K.max-1)
  for (k in 1:(K.max-1)){
    thresholds[k] <- thresh_params[k, 1]/(NROW(data)^thresh_params[k, 2]) + thresh_params[k, 3] 
  }
  return(thresholds)
}

# Use fitted parameters to reconstruct the reference dispersion if they are available,
get_dispersions_from_params <- function(data, K.max, dispersion_params){
  # Initialise the expected mean dispersion for each k 
  expected_dispersions <- numeric(K.max)
  # Reconstruct the expected mean dispersion from the supplied fitted parameters
  for (k in 1:K.max){
    expected_dispersions[k] <- dispersion_params[k, 1]*log(NROW(data)) + dispersion_params[k, 2]
  }
  return(expected_dispersions)
}


# Use the modfied gap statistic to find the estimated number of clusters in a dataset. If
# k_est is 1, there are no clusters. Note: needs threshold parameters, not thresholds themselves!
get_k_est <- function(input_data, FUNcluster, ..., K.max, B=NULL, d.power=2, 
                          thresholds, dispersions=NULL){
  
  # Can't use fitted dispersion paramenters AND also bootstrap. 
  stopifnot(is.null(B) | is.null(dispersions)) 
  
  # Initialise output
  out <- list()
  out$cluster <- array(numeric(), dim=c(NROW(input_data), K.max))
  
  # Ensure the data is scaled
  data <- scale_vec(input_data)
  
  # Get dispersions for the data
  data_dispersion <- array(numeric(), dim=c(K.max))

  # Note: FUNcluster is only designed to work with the 'kmeans' function as input.
  for (k in 1:K.max){
    kmeans_results <- FUNcluster(data, centers=k, ...)
    dispersions_for_k <- get_dispersion(data, kmeans_results$cluster,
                                        kmeans_results$centers, d.power=d.power)
    data_dispersion[k] <- log(mean(dispersions_for_k))
    out$cluster[,k] <- kmeans_results$cluster
  }
  
  # Initialise the expected mean dispersion for each k 
  expected_dispersions <- numeric(K.max)

  # For the reference dispersion, use supplied expected dispersions if they are available, otherwise
  # bootstrap by evaluating dispersions from B uniform distributions
  if (!is.null(dispersions)) { # Use supplied expected dispersions
    expected_dispersions <- dispersions
  }
  else if (!is.null(B)){ # Use bootstrapping instead (if the number of bootstraps is supplied). 
    # Construct expected mean dispersions by generating B reference datasets
    reference_dispersions <- array(numeric(), dim=c(K.max, B))
    for (ref in 1:B){
      ref_dist <- scale_vec(runif(length(data)))
      for (k in 1:K.max){
        kmeans_results <- FUNcluster(ref_dist, centers=k, ...)
        dispersions_for_k <- get_dispersion(ref_dist, kmeans_results$cluster,
                                                kmeans_results$centers, d.power=d.power)
        reference_dispersions[k, ref] <- log(mean(dispersions_for_k))
      }
    }
    expected_dispersions <- rowMeans(reference_dispersions) 
  }
  
  # Calculate delta gaps for each k pair
  gaps <- data_dispersion - expected_dispersions
  delta_gaps <- -diff(gaps)

  # Compare delta gap thresholds against the delta gaps to see if the data clustered for any k
  clustered_inds <- delta_gaps > thresholds
  
  k_est <- NA
  if (any(clustered_inds)){
    k_est <- which.max(delta_gaps-thresholds) + 1 # Add 1 to get from delta gap k to k
  } else {
    k_est <- 1
  }
  
  out$delta_gaps <- delta_gaps
  out$data_dispersion <- data_dispersion
  out$expected_dispersions <- expected_dispersions
  out$k_est <- k_est
  return(out)
}

# Synthesise multimodal dataset with specified characteristics. 
get_multimodal_data_sample <- function(n_data, n_modes, mode_sep_sd){
  mode_centres <- seq(from=0, by=mode_sep_sd, length.out=n_modes)
  data_in_modes <- rep(floor(n_data/n_modes), n_modes) # Initial allocation
  data_to_allocate_to_modes <- n_data %% n_modes # Calculate remaining unallocated data
  if (data_to_allocate_to_modes > 0){ # Is there unallocated data?
    data_in_modes[1:data_to_allocate_to_modes] <- data_in_modes[1:data_to_allocate_to_modes] + 1
  }
  # 'Collapse' data back into single 'pool'
  out <- purrr::reduce(purrr::map2(data_in_modes, mode_centres, rnorm, sd=1), c)
  return(out)
}





