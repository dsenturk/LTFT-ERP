LTFT.ERP_simulateTFT <- function(N.group,       # group (ASD/TD) sample size (scalar)
                                 SNR,           # targeted signal-to-noise ratio (scalar)
                                 freq,          # frequency band (Delta, Theta, Alpha, Beta, Gamma) (character)
                                 H,             # number of total leading eigenvectors, H (scalar)
                                 n.spline = 4,  # number of spline knots used in modeling longitudinal trends. Default n.spline = 4. (scalar)
                                 f_tot,         # total number of scale (frequency) parameters, F (scalar)
                                 d_tot,         # total number of translation (time, u) parameters, D (scalar)
                                 s,             # longitudinal domain grid (vector)
                                 MissProf,      # data.frame of missingness profiles of subjects
                                                # with 6 labeled columns (described below)
                                                   # DATA.FRAME COLUMNS:
                                                   # miss.id: missingness profile subject ID, i \in {1,...,N} (vector)
                                                   # condition: experimental condition, l (Expected, Unexpected) (vector)
                                                   # electrode: electrode ID, j \in {1,...,J} (vector)
                                                   # trial: experimental trial, s \in {s_{min},...,S} (vector)
                                                   # group: diagnostic group, (ASD, TD) (vector)
                                                   # region: electrode region ID, r \in {1,...,R} (vector)
                                                # and number of rows equal to the total number number of observations
                                                # across electrodes and subject-specific sets of non-missing trials (S_i) 
                                                # and respective conditions (L_{is}):
                                                   # nrow = J x [sum(i = 1)^{N} sum(s \in S_i) (length of set L_{is})]
                                 X.Mean,        # overall mean TFT power vector (vector, m x 1) 
                                                # (Note: m is the length of the functional domain)
                                 Phi,           # matrix of the H leading estimated eigenvectors (matrix, m x H)
                                 fixed.effects, # list() of length H such that each element contains the vector of estimated
                                                # fixed effects, B^h (40 x 1) for h = 1,...,H leading eigenvectors (list)
                                 random.effects # list() of length H such that each element contains a list() of length 2
                                                # of the estimated random effects for h = 1,...,H. For the hth element of the 
                                                # list() of length H, the first and second elements of the list() of length 2
                                                # are, respectively, the subject-specific, D^{1h} (matrix, (n.spline + 1) x (n.spline + 1)), 
                                                # and subject- and region-specific, D^{2h} (matrix, (n.spline + 1) x (n.spline + 1)),
                                                # random effects covariance matrices (list)
                                ){
  #############################################################################
  ## Description: Function for data generation of the simulated TFT power vectors from the delta frequency band as 
  ##              detailed in Section 3.1. Specifically, the H leading MDPCA score trajectories are simulated 
  ##              using the empirical mixed effects model estimates from our data analysis. Then, the 
  ##              (h = 1) simulated true mean and predicted MDPCA scores trajectories associated with the leading eigenvector (h = 1) are stored in a data.frame 
  ##              for future calculations of ME and PE. Missingness is then induced by removing a fraction of the 
  ##              TFT power vectors by sampling with replacement from the missingness profiles from subjects in 
  ##              our data. Using the simulated MDPCA score trajectories and the delta frequency band 
  ##              empirical estimates of the TFT mean power vector and the leading functional eigenvectors,
  ##              the simulated true TFT power vectors are reconstructed. After this, random noise vectors 
  ##              of length m = f_tot x d_tot are simulated idependently from a mean zero normal distribution 
  ##              whose variance is determined by the given SNR. The random noise vector is then wavelet
  ##              transformed to obtain the error TFT power vector. The error TFT power vector is 
  ##              added to the true TFT power vector to obtain the simulated TFT power vector for all 
  ##              observations and stored in a data.frame.
  ## Args:        (see above)
  ## Returns:     list()
  ##              Scores.Sim.pre: data.frame containing the true mean and predicted MDPCA score trajectories at h = 1 (needed for 
  ##              ME and PE calculations) with 9 labeled columns c("id", "electrode", "region", "trial", "condition", 
  ##              "group", "lobe", "fixed.scores1", "pred.scores1") described below 
  ##                  # DATA FRAME COLUMNS:
  ##                  # id: subject ID, i \in {1,...,2*N.group} (vector)
  ##                  # region: region, r \in {1,...,R} (vector)
  ##                  # electrode: electrode ID, j \in {1,...,J} (vector) 
  ##                  # condition: experimental condition, l (Expected, Unexpected) (vector)
  ##                  # group: diagnostic group (ASD = 0, TD = 1) (vector)
  ##                  # lobe: scalp section (Frontal = 0, Posterior = 1) (vector)
  ##                  # trial: experimental trial, s \in {s_min,...,S} (vector)
  ##                  # fixed.scores1: true estimated mean trajectory, E[Y^1_{ij(r)l}(s)] (vector)
  ##                  # pred.scores1: true estimated subject- and region-specific trajectory, E[Y^1_{ij(r)l}(s)|b^1_i, b^1_{ir}] (vector)
  ##              and number of rows equal to the number of observations across electrodes and subject-specific 
  ##              sets of non-missing trials (S_i) and respective conditions (L_{is}):
  ##                  # nrow2 = J x [sum(i = 1)^{2*N.group} sum(s \in S_i) (length of set L_{is})]
  ##              X.data: data.frame containing the simulated TFT power vectors in wide format with (4 + m) 
  ##              labeled columns c("id", "electrode", "trial", "condition", "X(t_1)",..., "X(t_m)")
  ##              where m is the length of the functional domain. Each row contains the (1 x m) X_{ijsl}
  ##              power vector stored in columns X(t_1),..., X(t_m) for a total of nrow2 rows.
  ## Simulation Outline:
  ##              1. Define model components
  ##              2. Simulate the H leading longtudinal MDPCA scores trajectories
  ##              3. Induce missingness from missingness profile 
  ##              4. Simulate the power vectors 
  ##              5. Add random noise to power vectors at specified SNR
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("data.table", "splines", "MASS", "compiler", "nlme", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  library(data.table)
  library(splines)
  library(MASS)
  library(compiler)
  library(nlme)
  library(dplyr)
  
  #############################################################################
  # 1. Define model components
  #############################################################################
  
  # Global variables 
  i <- c(1:(2*N.group)) # Subject index
  j <- unique(MissProf$electrode) # Eletrode index 
  r <- unique(MissProf$region) # Electrode region index
  l <- unique(MissProf$condition) # Condition index
  m_tot <- d_tot * f_tot # Length of functional domain (time x frequency, m = F x D)
  
  #############################################################################
  # 2. Simulate the H leading longtudinal MDPCA scores trajectories
  #############################################################################
  
  # Create cubic b-spline basis for trials
  spline <- as.data.frame(cbind(s, ns(s, df = n.spline))) # spline basis
  names(spline) <- c("trial", paste("spline", 1:n.spline, sep = ""))
  
  # Simulate the missingness profile 
  miss.id.td <- unique(MissProf[MissProf$group == "TD",]$miss.id) # Unique identifiers for TD
  miss.id.asd <- unique(MissProf[MissProf$group == "ASD",]$miss.id) # Unique identifiers for ASD
  miss.ids <- c(sample(miss.id.td, N.group, replace = TRUE), # Sample missingness profile identifiers with replacement
                sample(miss.id.asd, N.group, replace = TRUE)) 
  
  # Function simulating the longitudinal MDPCA scores using the empirical estimates from the mixed 
  # effects modeling from our motiviating study
  scores.sim <- function(h){
    
    # a. Define empirical model estimates 
    
    # Fixed effects, b^h (vector, 40 x 1)
    B <- fixed.effects[[h]]
    
    # Random effects
    D1 <- diag(sqrt(random.effects[[h]][[1]])) # Subject-specific, D^{1h} (matrix, (n.spline + 1) x (n.spline + 1))
    D2 <- diag(sqrt(random.effects[[h]][[2]])) # Subject- and region-specific, D^{2h} (matrix, (n.spline + 1) x (n.spline + 1))
    
    # b. Create data.frame to store indices and simulated MDPCA scores 
    
    # Create matrix to store calculated scores and index variables 
    Scores.Sim <- data.frame(id = rep(rep(i, each = length(j) * length(s)), length(l)), # Subject index, i
                             miss.id = rep(rep(miss.ids, each = length(j) * length(s)), length(l)), # Missingess profile index
                             electrode = rep(rep(j, each = length(s), times = length(i)), length(l)), # Electrode index, j
                             trial = rep(rep(s, length(i) * length(j)), 2), # Trials, s (longitudinal index)
                             condition = rep(c(0, 1), each = length(j) * length(s) * length(i))) # Condition, l (Expected = 0, Unexpected = 1)
    
    # Create electrode region index variable from missingness profile
    Scores.Sim <- inner_join(Scores.Sim, unique(MissProf[, c("electrode", "region")]))
    
    # Create subgroup variables used in mixed effects modeling
    Scores.Sim$lobe <- if_else(Scores.Sim$region %in% c("RF", "ZF", "LF"), # Create scalp section variable 
                               "Frontal", "Posterior")
    Scores.Sim$lobe <- ifelse(Scores.Sim$lobe == "Frontal", 0, 1) # (Frontal = 0, Posterior = 1)
    Scores.Sim$group <- if_else(Scores.Sim$miss.id %in% miss.id.td, "TD", "ASD") # Create diagnostic group variable 
    Scores.Sim$group <- ifelse(Scores.Sim$group == "ASD", 0, 1) # (ASD = 0, TD = 1)

    # Note: condition is both an index and subgroup variable
    
    # Merge cubic spline basis by trial variable 
    Scores.Sim <- inner_join(Scores.Sim, spline, by = "trial") 
    
    # c. Calculate true mean MDPCA scores trajectory
    
    # Create fixed effects design matrix 
    Q <- Scores.Sim %>% # Select key variables
      select(spline1, spline2, spline3, spline4, group, condition, lobe)
    Q <- cbind("(Inctercept)" = 1, Q) # Intercept variable
    Q <- Q %>% # Define interaction terms
      mutate(V9 = group * condition, V10 = group * lobe, V11 = condition * lobe,
             V12 = spline1 * group, V13 = spline2 * group, V14 = spline3 * group,
             V15 = spline4 * group, V16 = spline1 * condition, V17 = spline2 * condition, 
             V18 = spline3 * condition, V19 = spline4 * condition, V20 = spline1 * lobe, 
             V21 = spline2 * lobe, V22 = spline3 * lobe, V23 = spline4 * lobe,
             V24 = group * condition * lobe, V25 = spline1 * group * condition,
             V26 = spline2 * group * condition, V27 = spline3 * group * condition,
             V28 = spline4 * group * condition, V29 = spline1 * group * lobe, 
             V30 = spline2 * group * lobe, V31 = spline3 * group * lobe,
             V32 = spline4 * group * lobe, V33 = spline1 * condition * lobe,
             V34 = spline2 * condition * lobe, V35 = spline3 * condition * lobe,
             V36 = spline4 * condition * lobe, 
             V37 = spline1 * group * condition * lobe, 
             V38 = spline2 * group * condition * lobe,
             V39 = spline3 * group * condition * lobe, 
             V40 = spline4 * group * condition * lobe)
    Q <- as.matrix(Q) # Format as matrix 
    
    # Calculate true fixed effects estimates, E[Y^h_{ij(r)l}(s)] 
    fixed.scores <- tcrossprod(Q, t(as.matrix(B))) # (vector) 
    Scores.Sim$fixed.scores <- fixed.scores
    
    # d. Simulate random effects 
    
    # Simulate true subject-specific random effects, \b_{i}
    bi <- matrix(NA, nrow = length(i), ncol = (n.spline + 1)) # Create matrix to store random effects
    for(x in 1:length(i)){ # Simulate b_i from MVN 
      bi[x, ] <- c(rnorm(1, mean = 0, sd = D1[1]), 
                   rnorm(1, mean = 0, sd = D1[2]),
                   rnorm(1, mean = 0, sd = D1[3]), 
                   rnorm(1, mean = 0, sd = D1[4]),
                   rnorm(1, mean = 0, sd = D1[5]))
    }
    bi <- data.frame(id = i, miss.id = miss.ids, bi) # Format bi as data.frame
    names(bi) <- c("id", "miss.id", "bi.Intercept", "bi.spline1", "bi.spline2",
                   "bi.spline3", "bi.spline4") # Name columns
    
    # Simulate true subject- and region-specific random effects, \b_{ir}
    bir <- do.call(rbind, lapply(1:(length(i) * length(r)), function(x) # Simulate random effects from MVN
      {
      c(rnorm(1, mean = 0, sd = D2[1]), rnorm(1, mean = 0, sd = D2[2]),
        rnorm(1, mean = 0, sd = D2[3]), rnorm(1, mean = 0, sd = D2[4]),
        rnorm(1, mean = 0, sd = D2[5]))}))
    bir <- data.frame(id = rep(i, each = length(r)), 
                      miss.id = rep(miss.ids, each = length(r)), # Store random effects in data.frame
                      region = rep(r, length(i)), bir)
    bir$region <- factor(bir$region)
    names(bir) <- c("id", "miss.id", "region", "bir.Intercept", "bir.spline1",
                    "bir.spline2", "bir.spline3", "bir.spline4") # Name columns
    
    # Merge random effects with data
    Scores.Sim <- inner_join(Scores.Sim, bi, by = c("id", "miss.id")) # Merge subject-specific random effects
    Scores.Sim <- inner_join(Scores.Sim, bir, by = c("id", "miss.id", "region")) # Merge subject- and region-specific random effects
    
    # e. Calculate the true MDPCA score trajectories
    
    # Calculate true MDPCA score trajectories 
    Scores.Sim <- Scores.Sim %>% # Create variable to Scores.Sim via mutate
      mutate(scores = fixed.scores + bi.Intercept + # Simulated MDPCA score trajectories, Y^h_{ij(r)l}(s)  
               (bi.spline1 * spline1) + (bi.spline2 * spline2) + 
               (bi.spline3 * spline3) + (bi.spline4 * spline4) + bir.Intercept +
               (bir.spline1 * spline1) + (bir.spline2 * spline2) + 
               (bir.spline3 * spline3) + (bir.spline4 * spline4))
    
    # Format simulated hth leading MDPCA scores data.frame
    Scores.Sim <- Scores.Sim[, c("id", "miss.id", "region", "electrode", "condition", # Keep indicator and subgroup variables, and fixed and simulated scores
                                 "group", "lobe", "trial", "scores", 
                                 "fixed.scores")]
    names(Scores.Sim)[9] <- paste("scores", h, sep = "") # Name generated MDPCA score trajectory for model h
    names(Scores.Sim)[10] <- paste("fixed.scores", h, sep = "") # Name true mean MDPCA score trajectory from model h
    
    return(Scores.Sim)
  }

  # Simulate MDPCA scores for the H leading eigenvectors
  Scores.Sim <- list() # Create list to store simulated scores data.frames
  for(h in 1:H){
    Scores.Sim[[h]] <- scores.sim(h)
  }

  # Merge the simulated MDPCA scores data.frames into one data.table
  Scores.Sim <- Reduce(function(x, y)
    inner_join(x, y, by = c("id", "miss.id", "region", "electrode", "condition", 
                            "group", "lobe", "trial")), Scores.Sim)
  
  # Create data.table containing true MDPCA score trajectories for h = 1 (used for calculating ME and PE)
  Scores.Sim.pre <- Scores.Sim[, c("id", "electrode", "region", "trial", # Select variables to keep 
                                   "condition", "group", "lobe", 
                                   "fixed.scores1", "scores1")]
  Scores.Sim.pre$condition <- if_else(Scores.Sim.pre$condition == 0, "Expected", "Unexpected")

  # Format simulated MDPCA scores data 
  Scores.Sim <- Scores.Sim[, c("id", "miss.id", "electrode", "trial", # Select variables
                               "condition", "region", "group", "lobe", 
                               "scores1", "scores2", "scores3", "scores4", 
                               "scores5", "scores6")] 
  Scores.Sim <- data.table(Scores.Sim) # Format as data.table
  
  #############################################################################
  # 3. Induce missingness from missingness profile 
  #############################################################################
  
  # Format the missingess profile data
  MissProf$lobe <- ifelse(MissProf$region %in% c("LF", "RF", "ZF"), "Frontal", "Posterior") # Scalp section variable
  MissProf$lobe <- as.numeric(ifelse(MissProf$lobe == "Frontal", 0, 1))
  MissProf$condition <- as.numeric(ifelse(MissProf$condition == "Expected", 0, 1)) # Condition variable
  MissProf$group <- as.numeric(ifelse(MissProf$group == "ASD", 0, 1)) # Group variable
  
  # Down sample the simulated scores by missingness profile
  MissProf <- MissProf[MissProf$miss.id %in% miss.ids,] # Obtain missingness profile
  Scores.Sim <- inner_join(Scores.Sim, MissProf) # Merge simulated scores and missingness profile 
  Scores.Sim <- Scores.Sim[, -"miss.id"] # Remove unnecessary variable
  
  #############################################################################
  # 4. Reconstruct the true TFT power vectors
  #############################################################################
  
  # Separate scores and indicator variables
  score.columns <- ((ncol(Scores.Sim) - H + 1):ncol(Scores.Sim)) # Column indices for scores
  Scores.Sim.scores <- as.matrix(Scores.Sim[, ..score.columns], ncol = H) # Create matrix of scores (H columns)
  Scores.Sim <- Scores.Sim[, c("id", "electrode", "trial", "condition")] # Store indicator variables
  
  # Reconstruct the simulated true TFT power vectors 
  X.Signal <- t(vapply(1:nrow(Scores.Sim.scores), cmpfun(function(y){ # Reconstruct and store vectorized power signals, (matrix, m columns)
    X <- rowSums(Phi %*% diag(Scores.Sim.scores[y,])) + X.Mean # Calculate X_{ijsl}^{signal}
    return(X)}), numeric(m_tot)))
  rm(Scores.Sim.scores, score.columns) # Clean environment
  
  #############################################################################
  # 5. Add random noise to TFT power vectors at specified SNR
  #############################################################################
  
  # The following function was modified from the WaveletTransform function found in the 
  # WaveletComp package. For an in-depth explanation of details involved in the wavelet 
  # transformation utilized in this function, please refer to the WaveletComp documentation. 
  
  # a. Function for performing Morlet wavelet transformation
  WaveletTransform <- function (w,     # vectorized ERP waveform 
                                d_tot, # total number of translation parameters
                                f_tot, # total number of scale parameters
                                freq   # frequency band
  ){
    dt = 1/d_tot
    max_freq <- ifelse(freq == "Delta" | freq == "Theta", 12, # Upper bound of frequencies for specified frequency band
                       ifelse(freq == "Alpha", 16, ifelse(freq == "Beta", 32, 64)))
    # Note: max_freq is 12Hz for both Delta and Theta to be consistent with data analysis
    lowerPeriod <- 1/max_freq # The reciprocal of the maximum frequency (scale) parameter
    upperPeriod <- 1/0.5 # The reciprocal of the minimum frequency (scale) parameter
    series.length = length(w)
    pot2 = trunc(log2(series.length) + 0.5)
    pad.length = 2^(pot2 + 1) - series.length
    omega0 = 6 # Angular frequency
    fourier.factor = (2 * pi)/omega0
    min.scale = lowerPeriod/fourier.factor
    max.scale = upperPeriod/fourier.factor
    J = f_tot - 1
    p <- ceiling((J/log2(max.scale/min.scale)*2))/2 # Octaves per voice 
    dj <- 1/p # Frequency resolution, voices per octave 
    scales = min.scale * 2^((0:J) * dj)
    scales.length = length(scales)
    periods = fourier.factor * scales
    a_f = rev(1/periods)
    N = series.length + pad.length
    omega.k = 1:floor(N/2)
    omega.k = omega.k * (2 * pi)/(N * dt)
    omega.k = c(0, omega.k, -omega.k[floor((N - 1)/2):1])
    morlet.wavelet.transform = function(w) {
      xpad = c(w, rep(0, pad.length))
      fft.xpad = fft(xpad)
      wave = matrix(0, nrow = scales.length, ncol = N)
      wave = wave + (0+1i) * wave
      for (ind.scale in (1:scales.length)) {
        my.scale = scales[ind.scale]
        norm.factor = pi^(1/4) * sqrt(2 * my.scale/dt)
        expnt = -((my.scale * omega.k - omega0)^2/2) * (omega.k > 
                                                          0)
        daughter = norm.factor * exp(expnt)
        daughter = daughter * (omega.k > 0)
        wave[ind.scale, ] = fft(fft.xpad * daughter, inverse = TRUE)/N
      }
      wave = wave[, 1:series.length]
      return(wave)
    }
    Wave = morlet.wavelet.transform(w)
    return(Wave)
  }
  
  # b. Function to simulate and add the TFT noise vector to the true TFT power vector
  xnoise <- function(x){ # x = row index 
    sd.signal <- sd(X.Signal[x, ])
    c <- sqrt(sd.signal/(2 * SNR)) # SD of noise corrected for specified SNR
    noise <- rnorm(d_tot, mean = 0, sd = c) # Simulate random noise vector (m x 1)
    noise.p <- as.vector(t((Mod(WaveletTransform(noise, d_tot, f_tot, freq))^2)[f_tot:1, ])) # Vectorize the wavelet transform of random error vector, X^{noise}_{ijsl}
    x.sim <- as.numeric(t(X.Signal[x, ] + noise.p)) # Calculate simulated power vector, X_{ijsl} = X^{noise}_{ijsl} + X^{signal}_{ijsl}
    return(x.sim)
  }
  xnoise <- cmpfun(xnoise) # Pre-compile function
  
  # Add noise to simulate power signals
  X.Sim <- matrix(NA, nrow = nrow(X.Signal), ncol = m_tot)
  for(x in 1:nrow(X.Sim)){
    X.Sim[x,] <- xnoise(x)
  }
  rm(X.Signal); gc() # Clean environment
  
  # Format data of simulated power vectors 
  X.Sim <- data.table(cbind(Scores.Sim, X.Sim)) # Combine indicator variables and simulated power
  X.Sim$condition <- if_else(X.Sim$condition == 0, "Expected", "Unexpected")
  names(X.Sim) <- c("id", "electrode", "trial", "condition", paste("X(t_", 1:m_tot, ")", sep = "")) # Name columns
  rm(Scores.Sim); gc() # Clean environment
  
  # Return simulated data 
  Simulated.data <- list(Scores.Sim.pre = Scores.Sim.pre, X.Sim = X.Sim) 
  return(Simulated.data)
}