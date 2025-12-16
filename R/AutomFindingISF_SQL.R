#' @title Performing a re-indexing and re-ordering of feature groups based  
#'        on retention time.
#'        
#' @description
#' The `AutomFindingISF_SQL()` iterates over feature groups derived from XCMS 
#' (CON.new) and detects in-source fragments (ISFs) within each group using
#'  MS/MS peak information.
#' 
#' @param XCMS the resulting `XcmsExperiment` object from the 
#'        `XCMSgrouping_SQL()` function.
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
AutomFindingISF_SQL <- function(XCMS, spd, ms2.ft) {
  
  XCMS3 <- data.frame()
  fgs <- sort(unique(XCMS$CON.new))
  
  for (fg in fgs) {
    XCMS.ISF <- FindingISF_SQL(fg, XCMS, spd, ms2.ft)
    XCMS3 <- rbind(XCMS3, XCMS.ISF)
  }
  
  XCMS3
}
