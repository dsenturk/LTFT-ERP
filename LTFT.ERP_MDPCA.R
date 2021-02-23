LTFT.ERP_MDPCA <- function(X.data,   # data.frame in wide format with (4 + m) labeled columns (described below) 
                                         # DATA.FRAME COLUMNS:
                                         # id: subject ID, i \in {1,...,N} (vector)
                                         # electrode: electrode ID, j \in {1,...,J} (vector) 
                                         # trial: experimental trial, s \in {s_min,...,S} (vector)
                                         # condition: experimental condition, l (Expected, Unexpected) (vector)
                                         # X(t_1): TFT power, X_{ijsl}(t_1), (vector)
                                         # ...
                                         # X(t_m): TFT power, X_{ijsl}(t_m), (vector)
                                     # and number of rows equal to the total number of TFT power vectors
                                     # across electrodes and subject-specific sets of non-missing 
                                     # trials (S_i) and respective conditions (L_{is}):
                                         # nrow = J x [sum(i = 1)^{N} sum(s \in S_i) (size of set L_{is})]
                           m,        # length of functional domain (scalar)
                           s,        # longitudinal domain grid (vector)
                           H,        # maximum number of leading eigenvectors (scalar)
                           k,        # maximum number of trials included in the sliding window (scalar)
                           s.set = s # subset of longitudinal grid included in marginal covariance estimation. Default: s. (vector)
                           ){
  #############################################################################
  ## Description: Function for performing step 3 of the LTFT-ERP algorithm in which a
  ##              data-driven dimension reduction of the TFT power vectors is employed
  ##              via multidimensional principal component analysis (MDPCA). Specifically, 
  ##              this function estimates the trial-specific covariance matrices, 
  ##              the marginal covariance matrix, the H leading eigenvectors, and 
  ##              the resulting MDPCA scores.
  ## Args:        (see above)
  ## Returns:     list()
  ##              X.Mean: Estimated TFT mean power vector. \bar{X} (vector, m x 1)
  ##              XX.s: Trial-specific covariance matrices. (list)
  ##              XX: Marginal covariance matrix. (matrix, m x m)
  ##              Phi: Estimated leading eigenvectors. (matrix, m x H)
  ##              Scores: data.table containing the estimated MDPCA scores for the H indicated  
  ##              leading eigenvectors with (m + H) columns c("id", "electrode", "trial",     
  ##              "condition", "scores1",..., "scoresH") and number of rows equal to nrow.
  ## MDPCA Outline:
  ##              1. Calculate average power vector
  ##              2. Calculate trial specific covariances
  ##              3. Calculate marginal covariance
  ##              4. Calculate MDPCA scores
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("data.table", "compiler", "dplyr", "coop")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  library(data.table)
  library(compiler)
  library(dplyr)
  library(coop)
  
  #############################################################################
  # 1. Calculate average power vector
  #############################################################################
  
  X.Mean <- colMeans(as.matrix(X.data[, -c(1:4)]))
  
  print("1. Calculate average power vector (completed)")
  
  #############################################################################
  # 2. Calculate trial specific covariances
  ############################################################################
  
  # Create sliding window indexes 
  S <- max(s); s_min <- min(s)
  lower_window <- ifelse(s.set < k/2 + s_min, s_min, # Lower bounds of sliding windows
                         ifelse(s.set > S - k/2, 2*s.set - S, s.set - k/2 + 1)) 
  upper_window <- ifelse(s.set < k/2 + s_min, 2*s.set - s_min, # Upper bounds of sliding windows
                         ifelse(s.set > S - k/2, S, s.set + k/2))
  
  # Calculate trial specific covariances 
  XX.s <- matrix(list(), nrow = length(s.set)) # Initialize list to store covariances
  for(x in 1:length(s.set)){
    trials <- c(lower_window[x]:upper_window[x]) # Set of trials in A_s
    XX.s[[x]] <- covar(X.data[trial %in% trials, -c(1:4)]) # Trial-specific covariance calculation, \hat{\Sigma_s} (m x m)
    names(XX.s)[x] <- paste("trial", s.set[x], sep = "_") # Label trial for calculated matrix
  }
  
  print("2. Trial-specific covariance calculation (completed)")
  
  #############################################################################
  # 3. Calculate marginal covariance
  ############################################################################
  
  # Calculate marginal covariance, \hat{\Sigma} (m x m)
  XX <- Reduce("+", XX.s) / length(XX.s) # Method of moments estimation 
  
  print("3. Calculate marginal covariance (completed)")
  
  #############################################################################
  # 4. Calculate MDPCA scores
  ############################################################################
  
  # Calculate leading eigenvectors
  XXE <- eigen(XX) # Eigendecomposition  
  Phi <- XXE$vectors[, 1:H] # Obtain matrix of leading eigenvectors, \phi_h (m x H)
  
  # Center the power vectors 
  Scores <- X.data[, c(1:4)] # Store index variables in data.frame
  X.data <- as.matrix(X.data[, -c(1:4)]) # Store power vectors in matrix (nrow x m)
  X.data <- scale(X.data, center = TRUE, scale = FALSE) # Mean-center TFT power vectors

  # Calculate the MDPCA scores
  for(h in 1:H){ 
    scores <- tcrossprod(X.data, t(Phi[, h])); gc() # Calculate MDPCA scores, Y^h_{ij(r)l(s)}
    Scores <- cbind(Scores, scores)
  }
  names(Scores)[5:(5 + H - 1)] <- paste("scores", 1:H, sep = "") # Name score variables
  Scores <- data.table(Scores) # Format as data.table
  
  print("4. Calculate MDPCA scores (completed)")
  
  # Store MDPCA results in list
  MDPCA = list(X.Mean = X.Mean, XX.s = XX.s, XX = XX, Phi = Phi, Scores = Scores)
  
  return(MDPCA)
}