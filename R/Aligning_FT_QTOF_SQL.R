Aligning_FT_QTOF_SQL <- function(
    FT_path, QTOF_path, expnr.ft, expnr.syn, 
    Assoc = NULL, err = 0.02, t.ini = 5,
    lc.err = 1, rng = 2, minIon = 0.6, startpoint = 1,
    save_assoc = FALSE 
) {
  
  library(DBI)
  library(RSQLite)
  
  stopifnot(length(FT_path) == 1, is.character(FT_path), file.exists(FT_path))
  stopifnot(length(QTOF_path) == 1, is.character(QTOF_path), file.exists(QTOF_path))
  

  if (is.character(Assoc)) {
    if (file.exists(Assoc)) {
      Assoc_df <- read.table(Assoc, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      assoc_file <- Assoc
    } else {
      Assoc_df <- data.frame(
        ref_compid = integer(0),
        target_compid = integer(0),
        ref_database = character(0),
        target_database = character(0),
        stringsAsFactors = FALSE
      )
      assoc_file <- NULL
    }
  } else if (is.data.frame(Assoc)) {
    Assoc_df <- Assoc
    assoc_file <- NULL
  } else {
    Assoc_df <- data.frame(
      ref_compid = integer(0),
      target_compid = integer(0),
      ref_database = character(0),
      target_database = character(0),
      stringsAsFactors = FALSE
    )
    assoc_file <- NULL
  }
  
  FT_con   <- DBI::dbConnect(RSQLite::SQLite(), FT_path)
  QTOF_con <- DBI::dbConnect(RSQLite::SQLite(), QTOF_path)
  
  #Detect polarity
  ft_mode  <- dbGetQuery(FT_con, "SELECT mode FROM experiment WHERE expid = ?", params = list(expnr.ft))$mode
  qtof_mode <- dbGetQuery(QTOF_con, "SELECT mode FROM experiment WHERE expid = ?", params = list(expnr.syn))$mode
  polarity_ft  <- ifelse(ft_mode == "neg", 0, 1)
  polarity_qtof <- ifelse(qtof_mode == "neg", 0, 1)
  
  #Generate LCal and remove outliers
  LCal <- Aligning_General_SQL(FT_con, QTOF_con, expnr.ft, expnr.syn, err, t.ini)
  LCal <- matchFTSyn_SQL(LCal, FT_con, QTOF_con, minIon = minIon)
  LCal <- RemoveOutliers(LCal, rng)
  
  #Regression
  rg <- RegressionPie_LCalign(LCal, startpoint)
  PlotPie_LCalign(LCal, rg)
  
  #Fill Assoc
  Assoc_df <- FillAssocFTnQTOFn_SQL(
    FT_con, QTOF_con, Assoc_df, expnr.ft, expnr.syn,
    cutoff = 1, rg = rg, lc.err = lc.err, err = err, minIon = minIon,
    polarity_ft = polarity_ft, polarity_qtof = polarity_qtof,
    FT_path = FT_path, QTOF_path = QTOF_path
  )
  

  dbDisconnect(FT_con)
  dbDisconnect(QTOF_con)
  
  #Save updated Assoc if requested 
  if (save_assoc && !is.null(assoc_file)) {
    write.table(Assoc_df, file = assoc_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  return(Assoc_df)
}
