LTFT.ERP_MTFT.ERP_reshape <- function(data,       # data.frame in long format with six labeled columns (described below)
                                                    # DATA.FRAME COLUMNS:
                                                    # id: subject ID, i \in {1,...,N} (vector)
                                                    # electrode: electrode ID, j \in {1,...,J} (vector) 
                                                    # trial: experimental trial, s \in {s_min,...,S} (vector)
                                                    # condition: experimental condition, l (Expected, Unexpected) (vector)
                                                    # erp.time: time, u \in {1,...,U} (vector)
                                                    # value: microvoltage (W_{ijsl}(u)) (vector)
                                                  # and number of rows equal to the length of the vectorized ERP
                                                  # waveform observations across electrodes and subject-specific 
                                                  # sets of non-missing trials (S_i) and respective conditions (L_{is}):
                                                    # nrow = U x J x [sum(i = 1)^{N} sum(s \in S_i) (length of set L_{is})]
                                      freq,       # frequency band (Delta, Theta, Alpha, Beta, Gamma) (character)
                                      f_tot,      # total number of scale (frequency) parameters, F (scalar)
                                      d_tot       # total number of translation (time, u) parameters, D (scalar)
                                      ){
  ##############################################################################
  ## Description: Function for performing Step 1 and Step 2 of the LTFT-ERP algorithm
  ##              outlined in Section 2 of "A Study of Longitudinal Trends in Time-Frequency 
  ##              Transformations of EEG Data during a Learning Experiment" by Boland
  ##              et al. (2021). Specifically, this function filters the ERP waveforms 
  ##              for the specified frequency band, performs the wavelet transformations
  ##              to obtain the TFT power surfaces, and reshapes the subsequent 
  ##              TFT-surfaces into vectors.
  ## Args:        see above
  ## Returns:     list()
  ##              X.data: data.frame containing TFT power vectors in wide format with (4 + m) 
  ##              columns c("id", "electrode", "trial", "condition", "X(t_1)",..., "X(t_m)")
  ##              where m = f_tot x d_tot is the length of the functional domain. Each 
  ##              row contains the (1 x m) X_{ijsl} power vector stored in columns
  ##              X(t_1),..., X(t_m) for a total of nrow/U rows.
  ##              a_f: grid of scale (frequency) parameters (vector, F x 1)
  ##              b_d: grid of translation (time) parameters (vector, D x 1)
  ## MTFT.ERP_reshape Outline: 
  ##              1. Format data and filter by frequency band 
  ##              2. Perform wavelet transformation and reshape TFT surfaces
  ##              3. Output results
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("data.table", "signal", "dplyr", "compiler")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  library(data.table)
  library(signal)
  library(dplyr)
  library(compiler)
  
  #############################################################################
  # 1. Format data and filter by frequency band
  #############################################################################
  
  # Format data and sort
  data <- as.data.table(data) # Convert data to data.table format
  data <- data[order(id, electrode, condition, trial),] # Sort data.table
  
  # Define global variables
  i <- sort(unique(data$id)) # Subject ID index
  j <- sort(unique(data$channel)) # Electrode index
  s <- sort(unique(data$trial)) # Trial index
  l <- sort(unique(data$condition)) # Condition index
  u <- sort(unique(data$erp.time)) # ERP time index
  b_d <- u[round(seq(1, length(u), length.out = d_tot))] # Translation/time parameters, b_d
  m_tot <- d_tot * f_tot # Total grid points in functional domain, m
  
  # Filter DC Shift
  nyq <- length(u)/2 # Nyquist frequency: one-half of the sampling rate of ERP waveforms
  bfH <- butter(3, 1.25/nyq, type = 'high') # 3rd order high-pass butterworth filter at 1.25Hz
  data <- data[, valueH := signal::filter(bfH, value), ] # Filter ERP above 1.25Hz
  data <- data[, c("value"):= NULL];
  names(data)[6] <- "value"
  
  # Filter data by frequency band
  if(freq == "Delta"){ # Filter delta frequency band
    bfL4 <- butter(3, 4/nyq, type = 'low') # 3rd order low-pass butterworth filter at 4Hz
    data <- data[, valueL := signal::filter(bfL4, value),] # Filter ERP below 4Hz
    data <- data[, c("value"):= NULL];
    names(data)[6] <- "value"
  } else if (freq == "Theta") { # Filter theta frequency band
    bfH4 <- butter(3, 4/nyq, type = 'high') # 3rd order high-pass butterworth filter at 4Hz
    bfL8 <- butter(3, 8/nyq, type = 'low') # 3rd order low-pass butterworth filter at 8Hz
    data <- data[, valueH := signal::filter(bfH4, value),] # Filter ERP above 4Hz
    data <- data[, valueL := signal::filter(bfL8, valueH),] # Filter ERP below 8Hz
    data <- data[, c("value", "valueH"):= NULL];
    names(data)[6] <- "value"
  } else if (freq == "Alpha"){ # Filter alpha frequency band
    bfH8 <- butter(3, 8/nyq, type = 'high') # 3rd order high-pass butterworth filter at 8Hz
    bfL16 <- butter(3, 16/nyq, type = 'low') # 3rd order low-pass butterworth filter at 16Hz
    data <- data[, valueH := signal::filter(bfH8, value),] # Filter ERP above 8Hz
    data <- data[, valueL := signal::filter(bfL16, valueH),] # Filter ERP below 16Hz
    data <- data[, c("value", "valueH"):= NULL];
    names(data)[6] <- "value"
  } else if (freq == "Beta"){ # Filter beta frequency band
    bfH16 <- butter(3, 16/nyq, type = 'high') # 3rd order high-pass butterworth filter at 16Hz
    bfL32 <- butter(3, 32/nyq, type = 'low') # 3rd order low-pass butterworth filter at 32Hz
    data <- data[, valueH := signal::filter(bfH16, value),] # Filter ERP above 16Hz
    data <- data[, valueL := signal::filter(bfL32, valueH),] # Filter ERP below 32Hz
    data <- data[, c("value", "valueH"):= NULL];
    names(data)[6] <- "value"
  } else if (freq == "Gamma"){ # Filter gamma frequency band
    bfH32 <- butter(3, 32/nyq, type = 'high') # 3rd order high-pass butterworth filter at 32Hz
    data <- data[, valueH := signal::filter(bfH32, value),] # Filter ERP above 32Hz
    data <- data[, c("value"):= NULL];
    names(data)[6] <- "value"
  }

  #############################################################################
  # 2. Perform wavelet transformation and reshape TFT surfaces
  #############################################################################
  
  # Format the data from long to wide
  data <- dcast(data, ... ~ erp.time, value.var = "value") # Cast data from long to wide 
 
  # The following function was modified from the WaveletTransform function found in the 
  # WaveletComp package. For an in-depth explanation of details involved in the wavelet 
  # transformation utilized in this function, please refer to the WaveletComp documentation. 
  
  # Function for performing Morlet wavelet transformation
  WaveletTransform <- function (w,     # vectorized ERP waveform W_{ijsl}
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
    a_f = rev(1/periods) # Scale (frequency) parameters, a_f
    N = series.length + pad.length
    omega.k = 1:floor(N/2)
    omega.k = omega.k * (2 * pi)/(N * dt)
    omega.k = c(0, omega.k, -omega.k[floor((N - 1)/2):1])
    morlet.wavelet.transform = function(w) {
      w = (w - mean(w))/sd(w)
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
      wave = wave[, round(seq(1, series.length, length.out = d_tot))]
      return(wave)
    }
    Wave = morlet.wavelet.transform(w)
    return(list(Wave, a_f))
  }
  
  # Separate ERP waveforms and store in matrix
  W <- as.matrix(data[, -c("id", "electrode", "trial", "condition")]) # Store vectorized ERP waveform data in matrix (nrow x U)
  X.data <- data[, c("id", "electrode", "trial", "condition")] # Create data.frame to store index variables

  # Perform wavelet transformation on ERP waveforms
  X <- matrix(NA, nrow = nrow(W), ncol = m_tot) # Matrix to store TFT power vectors
  a_f <- rep(0, f_tot) # Vector to store scale (frequency) parameters
  for(g in 1:nrow(W)){
    wave <- WaveletTransform(W[g, ], d_tot, f_tot, freq) # Obtain wavelet coefficients, C_{ijsl}(a_f, b_d)
    if(g == 1){
      a_f <- wave[[2]] # Store scale (frequency) parameters
    }
    X[g,] <- as.vector(t((Mod(wave[[1]])^2)[f_tot:1, ])) # Calculate and vectorize power, X_{ijsl} 
  }
  rm(W, data)
  
  #############################################################################
  # 3. Output results
  #############################################################################
  
  # Create data.frame to store observations
  X.data <- data.frame(cbind(X.data, X)) # Combine index variables and TFT power vectors into data.frame
  names(X.data)[5:ncol(X.data)] <- paste("X(t_", 1:m_tot, ")", sep = "") # Name columns
  rm(X)
  
  # Create list to output results
  MTFT.ERP <- list(X.data = X.data, b_d = b_d, a_f = a_f)
  
  return(MTFT.ERP)
}
