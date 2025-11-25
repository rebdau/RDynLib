matchFTSyn_SQL <- function(LCal, FT_con, QTOF_con, polarity_ft = 0, 
                           polarity_qtof = 0,minIon = 0.6, 
                           aggregated_Ft = FALSE, aggregated_QTOF = FALSE) {
  
  library(DBI)
  library(RSQLite)
  
  # Load FT and QTOF compound tables
  subdb_ft  <- dbGetQuery(FT_con,
                          "SELECT retention_time, mass_measured, compound_id
     FROM ms_compound ORDER BY compound_id")
  
  subdb_syn <- dbGetQuery(QTOF_con,
                          "SELECT retention_time, mass_measured, compound_id
     FROM ms_compound ORDER BY compound_id")
  
  
  # Retrieve MS2 peak lists
  # FTMS MS2 peaks
  ft_cols <- dbListFields(FT_con, "msms_spectrum")
  
  if ("spectrum_type" %in% ft_cols) {
    ft_filter <- if (!aggregated_Ft) {
      "(s.spectrum_type IS NULL OR s.spectrum_type = '' OR s.spectrum_type = 'not_assembled')"
    } else {
      "s.spectrum_type = 'assembled'"
    }
  } else {
    # no spectrum_type column → no filtering
    ft_filter <- "1=1"
  }
  
  
  ms2_ft_df <- dbGetQuery(FT_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_ft,ft_filter))
  
  
  # Round and split into lists per compound_id
  ms2_ft_list <- split(round(ms2_ft_df$mz), ms2_ft_df$compound_id)
  ms2_ft_list <- lapply(ms2_ft_list, unique)
  
  
  # QTOF MS2 peaks 
  qtof_cols <- dbListFields(QTOF_con, "msms_spectrum")
  
  if ("spectrum_type" %in% qtof_cols) {
    qtof_filter <- if (!aggregated_QTOF) {
      "(s.spectrum_type IS NULL OR s.spectrum_type = '' OR s.spectrum_type = 'intersect_single_energy')"
    } else {
      "s.spectrum_type = 'merged_MSMS_all_energies'"
    }
  } else {
    qtof_filter <- "1=1"
  }
  
  
  ms2_syn_df <- dbGetQuery(QTOF_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_qtof,qtof_filter ))
  
  
  ms2_syn_list <- split(round(ms2_syn_df$mz), ms2_syn_df$compound_id)
  ms2_syn_list <- lapply(ms2_syn_list, unique)
  
  # Align MS2 lists by compound ID order

  ms2.ft  <- ms2_ft_list[ as.character(subdb_ft$compound_id) ]
  ms2.syn <- ms2_syn_list[ as.character(subdb_syn$compound_id) ]
  
  # Replace NULL or NA lists with empty vectors
  ms2.ft[sapply(ms2.ft, is.null)]   <- list(integer(0))
  ms2.syn[sapply(ms2.syn, is.null)] <- list(integer(0))
  

  
  if (is.null(LCal) || nrow(LCal) == 0) {
    message("No candidate alignments found; returning empty LCal.")
    return(LCal)
  }
  
  i <- 1
  while (i <= nrow(LCal)) {
    
    ft.compid  <- LCal[i, 1]
    syn.compid <- LCal[i, 7]
    
    ft.row  <- which(subdb_ft$compound_id  == ft.compid)
    syn.row <- which(subdb_syn$compound_id == syn.compid)
    
    # Skip pair if compound_id not found
    if (length(ft.row) == 0 || length(syn.row) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    
    ms2ion  <- ms2.ft[[ft.row]]
    msmsion <- ms2.syn[[syn.row]]
    
    # Skip if MS2 list is empty 
    if (length(ms2ion) == 0 || length(msmsion) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    
    # Compute matching peak ratio
    same_ion <- sum(ms2ion %in% msmsion)
    
    if (same_ion / length(ms2ion) < minIon) {
      LCal <- LCal[-i, ]
      next
    }
    
    i <- i + 1
  }
  
  return(unique(LCal))
}
