# direc refers to elution order of substrate and product. It is specific for
# a certain conversion: 1, product elute earlier than substrate; 2, product
# elute later than substrate; 3, product might elute earlier or later, but not
# at the same time as the substrate.
# Prod.exp refers to the considered experiment: EXPID
# peakwidth is necessary to create a minimum retention time difference between
# substrate and product.
# min relates to the minimum intensity of the product ions.

conv.CSPP_SQL <- function(inp.x, mzdiff, direc,
                          peakwidth, mzerr, expid) {
  
  
  Prod.dat <- inp.x[order(inp.x$mass_measured), ]
  Sub.dat  <- Prod.dat
  
  cspp.df <- data.frame(
    COMPID.sub   = integer(),
    MZ.sub       = numeric(),
    IONS.sub     = integer(),
    COMPID.prod  = integer(),
    MZ.prod      = numeric(),
    IONS.prod    = integer(),
    COMMON_IONS  = integer(),
    DOT_IONS     = numeric(),
    COMMON_LOSS  = integer(),
    DOT_LOSS     = numeric(),
    FORW_IONS    = numeric(),
    REV_IONS     = numeric(),
    FORW_LOSS    = numeric(),
    REV_LOSS     = numeric(),
    stringsAsFactors = FALSE
  )
  
  tR.df <- data.frame(
    tR.sub  = numeric(),
    tR.prod = numeric(),
    stringsAsFactors = FALSE
  )
  
  i <- 1
  repeat {
    
    if (i > nrow(Sub.dat)) break
    
    mz.sub  <- Sub.dat$mass_measured[i]
    mz.prod <- mz.sub + mzdiff
    
    prd.low  <- mz.prod - mzerr
    prd.high <- mz.prod + mzerr
    
    j <- 1
    while (j <= nrow(Prod.dat)) {
      
      if (Prod.dat$mass_measured[j] < prd.low) {
        Prod.dat <- Prod.dat[-j, ]
        next
      }
      
      if (Prod.dat$mass_measured[j] <= prd.high) {
        
        out <- targMS2comp(
          Sub.dat$compound_id[i],
          Prod.dat$compound_id[j],
          subDB,
          AnalMS
        )
        
        rt.sub  <- Sub.dat$retention_time[i]
        rt.prod <- Prod.dat$retention_time[j]
        
        if (
          (direc == 1 && rt.prod < rt.sub - peakwidth) ||
          (direc == 2 && rt.prod > rt.sub + peakwidth) ||
          (direc == 3 && (rt.prod < rt.sub - peakwidth ||
                          rt.prod > rt.sub + peakwidth))
        ) {
          cspp.df <- rbind(cspp.df, out)
          tR.df   <- rbind(tR.df, c(rt.sub, rt.prod))
        }
        
        j <- j + 1
        next
      }
      
      if (Prod.dat$mass_measured[j] > prd.high) break
    }
    
    if (nrow(Prod.dat) == 0) break
    i <- i + 1
  }
  
  return(cspp.df)
}

