#Regression between FT and Synapt chromatograms obtained for the same
#sample. Deviations from a straight curve might occur as different LC
#systems are involved. Therefore, a piecewise regression is performed
#including at most two knots. A hyperbolic curve was noted when
#aligning a Thermo Accela chromatogram with a Waters UPLC chromatogram.
#By looking at the result from multiple experiments, robust regression
#was more accurate than traditional least squares. The 'startpoint'
#allows to decide on the FT retention time from which regression
#should be started. 

RegressionPie_LCalign_SQL <- function(LCal, startpoint = 1, min_gap = 1) {
  # Extract numeric columns safely
  x <- as.numeric(LCal[, 2])  # FT retention times
  y <- as.numeric(LCal[, 4])  # Synapt retention times
  
  # Remove invalid rows
  valid_idx <- which(!is.na(x) & !is.na(y) & is.finite(x) & is.finite(y))
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  if(length(x) < 2 || length(y) < 2) {
    warning("Not enough valid data points for regression.")
    return(rep(NA, 6))
  }
  
  # Candidate knot positions: unique x above startpoint
  knots <- sort(unique(x[x >= startpoint]))
  n_knots <- length(knots)
  
  knottime1 <- c(0)
  knottime2 <- c(0)
  VarExp <- c(NA)
  
  # Initial simple linear regression
  out <- tryCatch(lm(y ~ x), error = function(e) NULL)
  VarExp[1] <- if(!is.null(out)) summary(out)$adj.r.squared else NA
  
  # One-knot piecewise regression
  for(k1 in knots) {
    t1 <- pmax(0, x - k1)
    out <- tryCatch(lm(y ~ x + t1), error = function(e) NULL)
    knottime1 <- c(knottime1, k1)
    knottime2 <- c(knottime2, 0)
    VarExp <- c(VarExp, if(!is.null(out)) summary(out)$adj.r.squared else NA)
  }
  
  # Two-knot piecewise regression (only consider knots with min_gap)
  if(n_knots >= 2) {
    for(i in 1:(n_knots-1)) {
      k1 <- knots[i]
      t1 <- pmax(0, x - k1)
      for(j in (i+1):n_knots) {
        k2 <- knots[j]
        if(k2 - k1 < min_gap) next  # skip too-close knots
        t2 <- pmax(0, x - k2)
        out <- tryCatch(lm(y ~ x + t1 + t2), error = function(e) NULL)
        knottime1 <- c(knottime1, k1)
        knottime2 <- c(knottime2, k2)
        VarExp <- c(VarExp, if(!is.null(out)) summary(out)$adj.r.squared else NA)
      }
    }
  }
  
  # Select the best model
  best_model <- data.frame(knottime1, knottime2, VarExp)
  best_model <- best_model[order(best_model$VarExp, decreasing = TRUE), ]
  
  t1 <- pmax(0, x - best_model$knottime1[1])
  t2 <- pmax(0, x - best_model$knottime2[1])
  
  library(MASS)
  exp.intercept <- exp.slope <- exp.t1 <- exp.t2 <- 0
  
  # Robust regression based on best model
  if(best_model$knottime1[1] == 0 && best_model$knottime2[1] == 0) {
    out <- tryCatch(rlm(y ~ x), error = function(e) NULL)
    if(!is.null(out)) { co <- summary(out)$coef; exp.intercept <- co[1,1]; exp.slope <- co[2,1] }
  } else if(best_model$knottime2[1] == 0) {
    out <- tryCatch(rlm(y ~ x + t1), error = function(e) NULL)
    if(!is.null(out)) { co <- summary(out)$coef; exp.intercept <- co[1,1]; exp.slope <- co[2,1]; exp.t1 <- co[3,1] }
  } else {
    out <- tryCatch(rlm(y ~ x + t1 + t2), error = function(e) NULL)
    if(!is.null(out)) { co <- summary(out)$coef; exp.intercept <- co[1,1]; exp.slope <- co[2,1]; exp.t1 <- co[3,1]; exp.t2 <- co[4,1] }
  }
  
  rg <- c(exp.intercept, exp.slope, best_model$knottime1[1], exp.t1,
          best_model$knottime2[1], exp.t2)
  return(rg)
}
