FillAssocFTnQTOFn_SQL <- function(
    FT_con, QTOF_con, Assoc, expnr.ft, expnr.syn,
    cutoff, rg, lc.err, err, minIon = 0.6,
    polarity_ft = 0, polarity_qtof = 0,
    FT_path, QTOF_path
) {
  
  
  ft.exp <- dbGetQuery(FT_con, sprintf(
    "SELECT retention_time AS rt, mass_measured AS mz, compound_id, expid 
     FROM ms_compound WHERE expid = %d", expnr.ft))
  
  syn.exp <- dbGetQuery(QTOF_con, sprintf(
    "SELECT retention_time AS rt, mass_measured AS mz, compound_id, expid 
     FROM ms_compound WHERE expid = %d", expnr.syn))

  ms2_ft <- dbGetQuery(FT_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2 AND s.polarity = %d
    ORDER BY s.compound_id, p.mz", polarity_ft))
  ms2_ft_list <- split(round(ms2_ft$mz), ms2_ft$compound_id)
  
  ms2_qtof <- dbGetQuery(QTOF_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2 AND s.polarity = %d
    ORDER BY s.compound_id, p.mz", polarity_qtof))
  ms2_qtof_list <- split(round(ms2_qtof$mz), ms2_qtof$compound_id)
  
  #Loop over FT compounds
  for (i in seq_len(nrow(ft.exp))) {
    COMPID <- ft.exp$compound_id[i]
    x.tR   <- ft.exp$rt[i]
    if (x.tR < cutoff) next
    
    t1.tR <- ifelse(rg[3] != 0 & x.tR > rg[3], x.tR - rg[3], 0)
    t2.tR <- ifelse(rg[5] != 0 & x.tR > rg[5], x.tR - rg[5], 0)
    y.tR <- rg[1] + rg[2] * x.tR + rg[4] * t1.tR + rg[6] * t2.tR
    y.l  <- y.tR - lc.err
    y.h  <- y.tR + lc.err
    
    pres <- Find_cand_matches(COMPID, err, syn.exp, ft.exp)
    if (nrow(pres) == 0) next
    
    ms2ion <- ms2_ft_list[[as.character(COMPID)]]
    if (is.null(ms2ion)) next
    
    pres1 <- data.frame()
    w.diff <- numeric()
    
    for (w in seq_len(nrow(pres))) {
      if ((pres$rt[w] > y.l) & (pres$rt[w] < y.h)) {
        msmsion <- ms2_qtof_list[[as.character(pres$compound_id[w])]]
        if (is.null(msmsion)) next
        
        same_ion <- length(intersect(ms2ion, msmsion))
        if (same_ion / length(ms2ion) >= minIon) {
          pres1 <- rbind(pres1, pres[w, ])
          w.diff <- c(w.diff, abs(y.tR - pres$rt[w]))
        }
      }
    }
    
    if (nrow(pres1) == 0) next
    pres1 <- pres1[order(w.diff), ]
    
    new_row <- data.frame(
      ref_compid    = COMPID,
      target_compid = pres1$compound_id[1],
      ref_database  = basename(FT_path),
      target_database = basename(QTOF_path),
      stringsAsFactors = FALSE
    )
    
    Assoc <- rbind(Assoc, new_row)
  }
  
  return(Assoc)
}
