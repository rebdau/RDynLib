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
  
  if (nrow(inp.x) == 0) {
    stop("No compounds found for expid = ", expid)
  }

  comp_add <- dbGetQuery(
    data_con,
    "SELECT * FROM compound_add"
  )
  

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
    
    cspp.df <- conv.CSPP_SQL(
      inp.x,
      mzdiff    = conver$conv.mz[k],
      direc     = conver$conv.direc[k],
      peakwidth = peakwidth,
      mzerr     = mzerr,
      expid     = expid
    )
    
    comp_add <- rank.cspp(
      cspp.df,
      conver$conv.colmn[k],
      comp_add
    )
  }

  comp_add$compound_id <- inp.x$compound_id[
    seq_len(nrow(comp_add))
  ]
  

  cols_to_update <- setdiff(colnames(comp_add), "compound_id")
  
  dbBegin(data_con)
  
  for (i in seq_len(nrow(comp_add))) {
    
    cid <- comp_add$compound_id[i]
    
    for (col in cols_to_update) {
      
      dbExecute(
        data_con,
        sprintf(
          "UPDATE compound_add
           SET %s = ?
           WHERE compound_id = ?",
          col
        ),
        params = list(comp_add[i, col], cid)
      )
    }
  }
  
  dbCommit(data_con)
  

  invisible(comp_add)
}


# cspp_add<-cspp.tot(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(cspp_add,"compound_add.txt",sep="\t",row.names=F)
