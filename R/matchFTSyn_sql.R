matchFTSyn_SQL <- function(LCal, FT_con, QTOF_con, minIon = 0.6) {
  
  library(DBI)
  library(RSQLite)
  
 
  subdb_ft  <- dbGetQuery(FT_con,  "SELECT retention_time, mass_measured, compound_id FROM ms_compound ORDER BY compound_id")
  subdb_syn <- dbGetQuery(QTOF_con, "SELECT retention_time, mass_measured, compound_id FROM ms_compound ORDER BY compound_id")
  
  ms2_ft_df <- dbGetQuery(FT_con, "
    SELECT c.compound_id,
           GROUP_CONCAT(DISTINCT CAST(ROUND(p.mz, 0) AS INT) ORDER BY CAST(ROUND(p.mz,0) AS INT)) AS peaklist
    FROM msms_spectrum c
    JOIN msms_spectrum_peak p ON c.spectrum_id = p.spectrum_id
    WHERE p.intensity > 0
    GROUP BY c.compound_id
    ORDER BY c.compound_id
  ")
  
  ms2_syn_df <- dbGetQuery(QTOF_con, "
    SELECT c.compound_id,
           GROUP_CONCAT(DISTINCT CAST(ROUND(p.mz, 0) AS INT) ORDER BY CAST(ROUND(p.mz,0) AS INT)) AS peaklist
    FROM msms_spectrum c
    JOIN msms_spectrum_peak p ON c.spectrum_id = p.spectrum_id
    WHERE p.intensity > 0
    GROUP BY c.compound_id
    ORDER BY c.compound_id
  ")
  
  # Align peak lists with compound_id
  ms2.ft  <- as.list(ms2_ft_df$peaklist[ match(subdb_ft[,3], ms2_ft_df$compound_id) ])
  ms2.syn <- as.list(ms2_syn_df$peaklist[ match(subdb_syn[,3], ms2_syn_df$compound_id) ])
  

  ms2.ft[ is.na(ms2.ft) ]   <- ""
  ms2.syn[ is.na(ms2.syn) ] <- ""
  
  i <- 1
  while(i <= nrow(LCal)) {
    
    ft.compid  <- LCal[i,1]
    syn.compid <- LCal[i,7]
    
    ft.compid_row  <- which(subdb_ft[,3] == ft.compid)
    syn.compid_row <- which(subdb_syn[,3] == syn.compid)
    
    # Skip if no match or empty MS2
    if (length(ft.compid_row) == 0 || length(syn.compid_row) == 0 ||
        ms2.ft[[ft.compid_row]] == "" || ms2.syn[[syn.compid_row]] == "") {
      LCal <- LCal[-i, ]
      next
    }
    
    ms2ion  <- as.integer(strsplit(ms2.ft[[ft.compid_row]], ",")[[1]])
    msmsion <- unique(as.integer(strsplit(ms2.syn[[syn.compid_row]], ",")[[1]]))
    
    same_ion <- sum(ms2ion %in% msmsion)
    
    if(same_ion / length(ms2ion) < minIon) {
      LCal <- LCal[-i, ]
      next
    }
    
    i <- i + 1
  }
  
  return(unique(LCal))
}
