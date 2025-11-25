#Regression between FT and Synapt chromatograms obtained for the same
#sample. Deviations from a straight curve might occur as different LC
#systems are involved. Therefore, a piecewise regression is performed
#including at most two knots. A hyperbolic curve was noted when
#aligning a Thermo Accela chromatogram with a Waters UPLC chromatogram.
#By looking at the result from multiple experiments, robust regression
#was more accurate than traditional least squares. The 'startpoint'
#allows to decide on the FT retention time from which regression
#should be started. 

RegressionPie_LCalign_SQL <- function(LCal, startpoint = 1) {

  x <- as.numeric(LCal[, 2])  # FT retention times
  y <- as.numeric(LCal[, 4])  # Synapt retention times
  
  # Remove NA
  valid <- which(!is.na(x) & !is.na(y))
  x <- x[valid]
  y <- y[valid]
  
  # Prepare output vectors
  knottime1 <- c()
  knottime2 <- c()
  VarExp <- c()
  
  # Initial model: no knots
  out <- lm(y ~ x)
  knottime1 <- c(knottime1, 0)
  knottime2 <- c(knottime2, 0)
  VarExp <- c(VarExp, summary(out)$adj.r.squared)
  

  unique_x <- sort(unique(x[x >= startpoint]))
  
  # One-knot models
  for (k1 in unique_x) {
    t1 <- ifelse(x > k1, x - k1, 0)
    out <- lm(y ~ x + t1)
    knottime1 <- c(knottime1, k1)
    knottime2 <- c(knottime2, 0)
    VarExp <- c(VarExp, summary(out)$adj.r.squared)
  }
  
  # Two-knot models 
  if (length(unique_x) > 1) {
    for (i in seq_len(length(unique_x) - 1)) {
      k1 <- unique_x[i]
      t1 <- ifelse(x > k1, x - k1, 0)
      
      for (j in (i + 1):length(unique_x)) {
        k2 <- unique_x[j]
        t2 <- ifelse(x > k2, x - k2, 0)
        out <- lm(y ~ x + t1 + t2)
        
        knottime1 <- c(knottime1, k1)
        knottime2 <- c(knottime2, k2)
        VarExp <- c(VarExp, summary(out)$adj.r.squared)
      }
    }
  }
  
  # Select best model
  best.model <- data.frame(knottime1, knottime2, VarExp)
  best.model <- best.model[order(best.model$VarExp, decreasing = TRUE), ]
  
  k1_best <- best.model$knottime1[1]
  k2_best <- best.model$knottime2[1]
  
  t1 <- ifelse(x > k1_best, x - k1_best, 0)
  t2 <- ifelse(x > k2_best, x - k2_best, 0)

  library(MASS)
  
  if (k1_best == 0 && k2_best == 0) {
    out <- rlm(y ~ x)
    co <- summary(out)$coef
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- 0
    exp.t2 <- 0
    
  } else if (k2_best == 0) {
    out <- rlm(y ~ x + t1)
    co <- summary(out)$coef
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- co[3, 1]
    exp.t2 <- 0
    
  } else {
    out <- rlm(y ~ x + t1 + t2)
    co <- summary(out)$coef
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- co[3, 1]
    exp.t2 <- co[4, 1]
  }
  
  return(c(exp.intercept, exp.slope, k1_best, exp.t1, k2_best, exp.t2))
}
