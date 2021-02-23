LTFT.ERP_mixedEffectsModel <- function(Scores,           # data.frame of the estimated MDPCA scores from the H leading eigenvectors 
                                                         # with (7 + H) labeled columns (described below)
                                                              # DATA.FRAME COLUMNS:
                                                              # id: subject ID, i \in {1,...,N} (vector)
                                                              # electrode: electrode ID, j \in {1,...,J} (vector)
                                                              # trial: experimental trial, s \in {s_min,...,S} (vector)
                                                              # condition: experimental condition, l (Expected, Unexpected) (vector)
                                                              # region: electrode region ID, r \in {1,...,R}
                                                              # group: diagnostic group, (ASD, TD) (vector)
                                                              # lobe: scalp section, (Frontal, Posterior) (vector)
                                                              # scores1: MDPCA score for leading eigenvector h = 1, (Y^1_{ij(r)l}(s)) (vector)
                                                              # ...
                                                              # scoresH: MDPCA score for leading eigenvector h = H, (Y^H_{ij(r)l}(s)) (vector)
                                                          # and number of rows equal to the total number of observations across electrodes and 
                                                          # subject-specific sets of non-missing trials (S_i) and respective conditions (L_{is}):
                                                              # nrow = J x [sum(i = 1)^{N} sum(s \in S_i) (length of set L_{is})]
                                       H,                 # number of leading eigenvectors used (scalar)
                                       s,                 # longitudinal domain grid (vector)
                                       n.spline = 4,      # number of spline knots used in modeling longitudinal trends. Default n.spline = 4 (scalar)
                                       reControl = list() # control argument for lme() function found in the nlme package. Defaults to an empty list.  
                                                          # Note: Please refer to the lme() function documentation in the nlme package for more details
                                       ){ 

  #############################################################################
  ## Description: Function for performing step 4 of the LTFT-ERP algorithm in which the
  ##              longitudinal trends of the MDPCA scores are modeled via a linear mixed 
  ##              effects model.
  ## Args:        (see above)
  ## Returns:     list()
  ##              B: fixed effects estimates (list, length H)
  ##              bi: subject-specific random effects MDPCA score estimates, \hat{b_i^h} (list, length H)
  ##              bir: subject- and region-specific random effects MDPCA score estimates, \hat{b_{ir}^h} (list, length H)
  ##              D1: subject-specific random effects covariance matricx estimates, \hat{D^{1h}} (list, length H)
  ##              D2: subject- and region-specific random effects covariance matrix estimates, \hat{D^{2h}} (list, length H)
  ##              sigma2: random error variance, \hat{\sigma^2_h} (list, length H)
  ## mixedEffectsModel Outline:
  ##              1. Define Model Arguments
  ##              2. Fit mixed effects model 
  ##              3. Store and return model estimates
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("splines", "nlme", "data.table", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  library(splines)
  library(nlme)
  library(data.table)
  library(dplyr)
  
  #############################################################################
  # 1. Define Model Arguments
  ############################################################################# 
  
  # a. Format Scores data 
  
  # Format variables in Scores
  Scores$region <- factor(Scores$region) # Turn region variable into factor 
  Scores$group <- factor(Scores$group) # Turn group variable into factor 
  Scores$lobe <- factor(Scores$lobe) # Turn scalp section variable into factor 
  Scores$condition <- factor(Scores$condition) # Turn condition variable into factor
  
  # Create natural cubic B-spline basis for trial
  spline <- as.data.frame(cbind(s, ns(s, df = n.spline))) 
  names(spline) <- c("trial", paste("spline", 1:n.spline, sep = ""))
  splines <- names(spline)[-1] # Names of spline fixed effects variables
  Scores <- inner_join(Scores, spline, by = "trial") # Merge splines and scores data by trial
  
  # b. Define model arguments in lme function 
  
  # Note: The following code utilizes functions from the nlme package, particularly the lme function. 
  #       For a more detailed explanation of implementation and arguements of these functions, please '
  #       refer to the nlme package documentation.
  
  
  # Create list to define random effects matrices 
  r <- list(id=pdDiag(as.formula(paste("~", paste(splines, # Subject-level random effects, D^{1h} (matrix, (n.spline + 1) x (n.spline + 1))
                                                  collapse = " + "), sep = ""))),
            region=pdDiag(as.formula(paste("~", paste(splines, collapse = " + "), # Subject- and region-level random effects, D^{2h} (matrix, (n.spline + 1) x (n.spline + 1))
                                           sep = ""))))
  
  # Create design matrix for lme
  vars <- c("group", "condition", "lobe", "group*condition", "group*lobe", # Names of subgroup and subgroup interactions fixed effects variables 
            "condition*lobe", "group*condition*lobe") 
  n.vars <- length(vars) # Total number of subgroup fixed effects variables
  f <- lapply(1:H, function(h){ # Store the linear formula objects describing the fixed-effects, B^h (list)
    as.formula(paste(paste("scores", h, sep = ""), # Define linear formula object of "scoresh ~ fixed effects" 
                     paste(c(splines, vars,
                             apply(expand.grid(splines, vars), 1,
                                   paste, collapse="*")),
                           collapse = " + "), sep = " ~ "))})

  
  #############################################################################
  # 2. Fit mixed effects model 
  #############################################################################
  
  # Fit model
  fit <- lapply(1:H, function(h){ 
    fit <- lme(f[[h]], data=Scores, r, control=reControl, keep.data=FALSE)
    return(fit)
  })
  
  #############################################################################
  # 3. Store and return model estimates
  #############################################################################
  
  # Store fixed effects estimates
  B <- list()
  for(h in 1:H){
    B[[h]] <- fit[[h]]$coefficients$fixed # \hat{B^h} (vector, (n.var x (n.spline + 1)) x 1) 
  }
  
  # Store subject-specific random effects scores estimates
  bi <- list()
  for(h in 1:H){
    bi[[h]] <- fit[[h]]$coefficients$random$id # \hat{b_i^h} (matrix, length(i) x (n.spline + 1))
  }
  
  # Store subject- and region-specific random effects scores estimates
  bir <- list()
  for(h in 1:H){
    bir[[h]] <- fit[[h]]$coefficients$random$region # \hat{b_{ir}^h} (matrix, (length(i) * length(r)) x (n.spline + 1))
  }
  
  # Store subject-specific random effects covariance matrices estimates
  D1 <- list()
  for(h in 1:H){
    D1[[h]] <- diag(as.numeric(VarCorr(fit[[h]])[2:6, 1]))  # \hat{D^{1h}} (matrix, (n.spline + 1) x (n.spline + 1))
  }
  
  # Store subject- and region-specific random effects covariance matrices estimates
  D2 <- list()
  for(h in 1:H){
    D2[[h]] <- diag(as.numeric(VarCorr(fit[[h]])[8:12, 1])) # \hat{D^{2h}} (matrix, (n.spline + 1) x (n.spline + 1))
  }
  
  # Store random error variance estimates
  sigma2 <- list()
  for(h in 1:H){
    sigma2[[h]] <- as.numeric(VarCorr(fit[[h]])[13, 1]) # \hat{\sigma^2_h} (scalar)
  }
  
  # Note: In our empirical model, B is a (40 x 1) vector, and D1 and D2 are both (5 x 5) matrices
  
  # Create list to store model estimates 
  ME.model <- list(B = B, bi = bi, bir = bir, D1 = D1, D2 = D2, sigma2 = sigma2)
  
  return(ME.model)
}