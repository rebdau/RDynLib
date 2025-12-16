#' @title searching for combination and in-source fragments within the FGs 
#' using CID spectral data. 
#'
#' @description
#' The `XCMSgrouping2_SQL()` function searches for in-source fragments within 
#' the FGs using CID spectral data via the nested AutomFindingISF()
#' function In addition, the function searches for features that represent  
#' combinations of other features present in the proper and in the neighboring 
#' FGs via the SearchCombinat() function. The function takes the xcms resulting 
#' object from the XCMSgrouping() function as input and creates an upgraded 
#' xcms object.
#' 
#' @param x `XcmsExperiment` object resulting from the XCMSgrouping() function.
#' 
#' @param adjPeak `numeric()` number of adjacent feature groups in the 
#' chromatogram to search for possible combinations of features, 1 by default.
#' 
#' @param mode `character(1)` is the mode of the experiment, either neg or pos.
#' 
#' @return An update `XcmsExperiment` object with new and updated columns in 
#'         featureDefinitions .
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
XCMSgrouping2_SQL <- function(x, adjPeak = 1, mode = NULL) {
  
  ## Feature table
  dat <- featureDefinitions(x)
  dat$feature_id <- rownames(dat)
  
  ## Create names (same as before)
  dat$name <- .create_name(dat$mzmed, dat$rtmed)
  
  ## Spectra & MS2 peaks
  s <- spectra(x)
  ms2.ft <- peaksData(s)
  spd <- spectraData(s)
  
  ## Polarity
  if (is.null(mode)) {
    mode <- ifelse(spd$polarity[1] < 0, "neg", "pos")
  }
  
  ## ISF detection
  XCMS4 <- AutomFindingISF_SQL(dat, spd, ms2.ft)
  
  ## Downstream logic unchanged
  XCMS4 <- ConvertXCMS4(XCMS4)
  XCMS4 <- SearchCombinat(XCMS4, adjPeak, mode)
  
  ## Write back to featureDefinitions
  new_cols <- setdiff(colnames(XCMS4), colnames(featureDefinitions(x)))
  for (col in new_cols) {
    featureDefinitions(x)[[col]] <- XCMS4[[col]]
  }
  
  return(x)
}


#'helper function to a create a name for the features.
#'
.create_name <- function(mzmed, rtmed) {
  paste0(
    "M", round(mzmed, 0),
    "T", round(rtmed, 0)
  )
}
