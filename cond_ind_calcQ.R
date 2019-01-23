# Function to return a matrix of significant partial correlations.
# Significance is tested using a bootstrap method.
# Returns results as a data frame (Q_neurons).
# To use with corrplot, convert to a matrix:
# Q <- as.matrix(Q_neurons)
# colnames(Q) <- colnames(df)
# rownames(Q) <- colnames(df)


calcQ <- function(df) {
  Rho_neurons     <- cor(df)
  tol <- 0
  Q_neurons       <- corpcor::cor2pcor(Rho_neurons)
  Q_neurons       <- Q_neurons %>% as.data.frame
  
  ## bootstrap, resampling rows of data matrix
  M <- 1000
  index <- 1:nrow(df)
  ## index <- 1:length(unique(id))
  ## uniqueid <- unique(id)
  Q_neurons_star <- array(dim=c(ncol(df), ncol(df), M))
  tmp <- NULL
  for(i in 1:M){
    index_star <- base::sample(x=index, length(index), replace=TRUE)
    # tmp <- rbind(tmp, data.sc[index_star,]) # don't think we need this
    ## for(j in 1:length(uniqueid[index_star])){
    ##     tmp <- rbind(tmp, data.sc[which(id == uniqueid[index_star][j]),])
    ## }
    Rho_neurons_star <- cor(df[index_star,])
    Q_neurons_star[,,i] <- as.matrix(corpcor::cor2pcor(Rho_neurons_star))
    tmp <- NULL
  }
  
  Q_neurons_low <- apply(Q_neurons_star, c(1,2), quantile, 0.025)
  Q_neurons_upp <- apply(Q_neurons_star, c(1,2), quantile, 0.975)    
  CI <- Q_neurons_low * Q_neurons_upp
  CI[CI<0]  <- 0
  CI[CI!=0] <- 1
  CI <- as.data.frame(CI)
  Q_neurons <- CI*Q_neurons
  
  return(Q_neurons)
}