#Remove false positive FTm/z,QTOFm/z peak pairs by checking
#whether a certain percentage of the MS2 spectral peaks
#can be traced in the MS/MS spectrum. By default, 60% of 
#the MS2 peaks should be traced in the MS/MS spectrum as
#the MS2 spectrum is much less crowded.

matchFTSyn_SQL <- function(LCal, FT_path, QTOF_path, minIon=NULL, chrg=NULL){
  
  if (is.null(minIon)) minIon = 0.6
  
  library(DBI)
  library(RSQLite)
  
  # Connect to SQL databases
  FT_con  <- dbConnect(SQLite(), FT_path)
  QTOF_con <- dbConnect(SQLite(), QTOF_path)
  
  # subdb_ft and subdb_syn must have compound_id in column 3 (same index as finlist)
  subdb_ft  <- dbGetQuery(FT_con,  "SELECT name, smiles, compound_id FROM ms_compound ORDER BY compound_id")
  subdb_syn <- dbGetQuery(QTOF_con, "SELECT name, smiles, compound_id FROM ms_compound ORDER BY compound_id")
  
  # Build MS2 peak lists EXACTLY like finlist (rounded integer, sorted, unique, comma string)
  ms2_ft_df <- dbGetQuery(FT_con, "
    SELECT c.compound_id,
           GROUP_CONCAT(DISTINCT CAST(ROUND(p.mz,0) AS INT) ORDER BY ROUND(p.mz,0)) AS peaklist
    FROM msms_spectrum c
    JOIN msms_spectrum_peak p ON c.spectrum_id = p.spectrum_id
    GROUP BY c.compound_id
    ORDER BY c.compound_id
  ")
  
  ms2_syn_df <- dbGetQuery(QTOF_con, "
    SELECT c.compound_id,
           GROUP_CONCAT(DISTINCT CAST(ROUND(p.mz,0) AS INT) ORDER BY ROUND(p.mz,0)) AS peaklist
    FROM msms_spectrum c
    JOIN msms_spectrum_peak p ON c.spectrum_id = p.spectrum_id
    GROUP BY c.compound_id
    ORDER BY c.compound_id
  ")
  
  dbDisconnect(FT_con)
  dbDisconnect(QTOF_con)
  
  # Align peak list vectors so indexes correspond to the compound_id in subdb tables
  ms2.ft  <- ms2_ft_df$peaklist[ match(subdb_ft[,3],  ms2_ft_df$compound_id) ]
  ms2.syn <- ms2_syn_df$peaklist[ match(subdb_syn[,3], ms2_syn_df$compound_id) ]
  
  # Replace missing spectra with "0" 
  ms2.ft[ is.na(ms2.ft) ]   <- "0"
  ms2.syn[ is.na(ms2.syn) ] <- "0"
  

  i = 1
  while(i <= dim(LCal)[1]){
    ft.compid = LCal[i,1]
    ft.compid_row <- which(as.integer(subdb_ft[,3]) %in% as.integer(ft.compid))
    syn.compid = LCal[i,7]
    syn.compid_row <- which(as.integer(subdb_syn[,3]) %in% as.integer(syn.compid))
    ms2ion  = as.integer(strsplit(ms2.ft[[as.integer(ft.compid_row)]], ",")[[1]])
    msmsion = unique(as.integer(strsplit(ms2.syn[[as.integer(syn.compid_row)]], ",")[[1]]))
    same_ion = which(ms2ion %in% msmsion)
    if(length(same_ion)/length(ms2ion) < minIon){
      LCal = LCal[-i, ]
      next
    }
    i = i + 1
  }
  
  return(unique(LCal))
}

