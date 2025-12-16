#' @title finding in-source fragments.
#' 
#' @description
#' finding whether a m/z feature is an in-source fragment (ISF) of another m/z
#' feature in a particular feature group. 
#'fg is the feature group number.
#'
#' @param XCMS the resulting `XcmsExperiment` object from the 
#'        `XCMSgrouping_SQL()` function.
#'        
#' @param fg set within the AutomFindingISF() function, represents the feature 
#'        group.
#'        
#' @param spd spectra Data of the `XcmsExperiment` object.
#' 
#' @param ms2.ft peaks data of the `XcmsExperiment` object.
#' 
#' @return An updated `XcmsExperiment` object with new annotation column 
#'         created inside FindingISF_SQL() function.
#'         
#'         
#' @import xcms
#'
#' @import MsExperiment
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom Spectra
#'
#' @importFrom dplyr
#'
#' @author Ahlam Mentag
#'
#' @export 

FindingISF_SQL <- function(fg, XCMS, spd, ms2.ft) {
  
  ## select feature group
  XCMS.sel <- XCMS[XCMS$CON.new == fg, ]
  XCMS.sel$feature_id <- rownames(XCMS.sel)
  
  isf <- rep(NA, nrow(XCMS.sel))
  
  ## MS2 spectra linked to these features
  idx <- which(spd$feature_id %in% XCMS.sel$feature_id)
  
  if (length(idx) > 0) {
    
    ## parent m/z
    parent <- round(XCMS.sel$mzmed, 0)
    
    ## collect product ions
    prod_ion <- c()
    for (i in idx) {
      pi <- ms2.ft[[i]][,1]  # m/z column
      par <- which(round(pi,0) %in% parent)
      if (length(par) > 0) pi <- pi[-par]
      prod_ion <- c(prod_ion, pi)
    }
    
    ## check ISF
    xcms.mz <- round(XCMS.sel$mzmed)
    isf[xcms.mz %in% round(prod_ion)] <- "ISF"
  }
  
  cbind(XCMS.sel, isf)
}

