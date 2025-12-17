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
  
  # Select features in this feature group
  XCMS.sel <- XCMS[XCMS$CON.new == fg, ]
  XCMS.sel$feature_id <- rownames(XCMS.sel)
  
  isf <- rep(NA_character_, nrow(XCMS.sel))
  
  ms2.fid <- spd$feature_id
  # Find indices of MS2 spectra corresponding to these features
  idx <- which(ms2.fid %in% XCMS.sel$feature_id)
  
  if (length(idx) > 0) {
    
    # Rounded parent m/z of features in this group
    parent <- round(XCMS.sel$mzmed, 0)
    
    # Collect all product ions across these MS2 spectra
    prod_ion <- c()
    for (i in idx) {
      pi <- ms2.ft[[i]][,1]           # m/z column
      pi_round <- round(pi, 0)        # round to integer
      # Remove parent m/z from product ions
      pi_round <- pi_round[!pi_round %in% parent]
      prod_ion <- c(prod_ion, pi_round)
    }
    
    # Rounded feature m/z for matching
    xcms_mz <- round(XCMS.sel$mzmed, 0)
    
    # Assign ISF where feature m/z is in product ions
    isf[xcms_mz %in% prod_ion] <- "ISF"
  }
  
  # Return the XCMS subset with isf column
  cbind(XCMS.sel, isf)
}


 

