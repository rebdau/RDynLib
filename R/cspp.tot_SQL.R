# No rows representing the different COMPIDs of the experiment have to be
# added to compound_add.txt, this is done while running the function. In case
# no conversions are found for a particular COMPID, a row of NAs is added
# preceeded by the name of the COMPID, the function jumps then to the next
# COMPID. However, at the end of the COMPID list belonging to the experiment,
# some final expected rows might be missing whenever no conversions were
# found for the few last COMPIDs. This is not a problem as these last
# COMPIDs followed by NAs will be added upon running the next experiment.
# Note: Conversions for a particular COMPID are added to the row with row
# number equal to the COMPID.

cspp.tot_SQL <- function(sql_path, expid, mzerr = 0.015,
                         cspp = "cspp.txt", peakwidth = NULL) {
  
  data_con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(data_con), add = TRUE)
  
  if (is.null(peakwidth)) peakwidth <- 0.2
  
  ## get compound IDs for the experiment
  inp.x <- dbGetQuery(
    data_con,
    sprintf(
      "SELECT compound_id,
            mass_measured,
            retention_time
     FROM ms_compound
     WHERE expid = %d",
      expid
    )
  )
  
  
  ## get compound_add table
  comp_add <- dbGetQuery(
    data_con,
    "SELECT compound_id FROM compound_add"
  )
  
  ## read CSPP conversion table
  conv.table <- read.table(
    cspp,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  
  conver <- data.frame(
    conv.type  = as.character(conv.table[, 1]),
    conv.mz    = as.numeric(conv.table[, 3]),
    conv.direc = as.integer(conv.table[, 5]),
    conv.colmn = as.integer(conv.table[, 6]),
    stringsAsFactors = FALSE
  )
  
  
  for (k in seq_len(nrow(conver))) {
    
    mzdiff <- conver$conv.mz[k]
    direc  <- conver$conv.direc[k]
    conv.col <- conver$conv.colmn[k]
    
    cspp.df <- conv.CSPP(
      inp.x,
      mzdiff,
      direc,
      peakwidth,
      mzerr,
      expid
    )
    
    comp_add <- rank.cspp(cspp.df, conv.col, comp_add)
  }
  
  comp_add$compound_id <- inp.x$compound_id[seq_len(nrow(comp_add))]
  
  return(comp_add)
}


# cspp_add<-cspp.tot(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(cspp_add,"compound_add.txt",sep="\t",row.names=F)
