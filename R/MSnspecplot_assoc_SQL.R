get_aligned_compids <- function(dbkey, db_name, assoc_path) {
  
  assoc <- read.table(assoc_path,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE)
  
  aligned <- assoc[
    assoc$ref_compid == dbkey &
      assoc$ref_database == db_name,
  ]
  
  if (nrow(aligned) == 0)
    return(NULL)
  
  return(aligned[, c("target_compid", "target_database")])
}

MSnspecplot_assoc_SQL <- function(sql_dir,
                            db_name,
                            dbkey,
                            assoc_path,
                            nr_col = 35,
                            nr_col2 = 6,
                            lc.err = 0.02,
                            mz.err = 0.001,
                            MS1 = NULL,
                            prcx = 0.6) {
  
  library(DBI)
  library(RSQLite)
  
  main_sql_path <- file.path(sql_dir, db_name)
  
  if (!file.exists(main_sql_path))
    stop("Main SQLite file not found.")
  
  #  Main DB connection
  con <- dbConnect(SQLite(), main_sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  #  CSPP + GNPS for main compoun
  cspp.res <- cspp.display_SQL(
    sql_path = main_sql_path,
    dbkey = dbkey,
    nr_col = nr_col
  )
  
  gnps.res <- gnps.display_SQL(
    sql_path = main_sql_path,
    dbkey = dbkey,
    nr_col = nr_col2
  )
  
  cat("\nAssociated CSPPs:\n")
  print(cspp.res)
  
  cat("\nAssociated GNPS conversions:\n")
  print(gnps.res)
  
  #  Plot main MSn
  MSnplot_SQL(
    sql_path = main_sql_path,
    dbkey = dbkey,
    prcx = prcx
  )
  
  #  Alignment handling
  aligned <- get_aligned_compids(
    dbkey = dbkey,
    db_name = db_name,
    assoc_path = assoc_path
  )
  
  if (!is.null(aligned)) {
    
    cat("\nAligned compounds in other experiments:\n")
    print(aligned)
    
    for (k in seq_len(nrow(aligned))) {
      
      target_db  <- aligned$target_database[k]
      target_key <- aligned$target_compid[k]
      
      target_sql_path <- file.path(sql_dir, target_db)
      
      if (!file.exists(target_sql_path)) {
        warning(paste("Database not found:", target_db))
        next
      }
      
      cat("\nPlotting aligned COMPID", target_key,
          "from", target_db, "\n")
      
      MSnplot_SQL(
        sql_path = target_sql_path,
        dbkey = target_key,
        prcx = prcx
      )
    }
  }
  
  return(list(
    cspp = cspp.res,
    gnps = gnps.res,
    aligned = aligned
  ))
}
