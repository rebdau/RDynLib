Find_cand_matches_SQL <- function(COMPID, err, syn.exp, ft.exp) {
  # Find the FT compound
  i <- which(ft.exp$compound_id == COMPID)
  if (length(i) == 0) return(data.frame(
    retention_time = numeric(0),
    mass_measured  = numeric(0),
    compound_id    = integer(0),
    expid          = integer(0)
  ))
  
  ft_mass <- ft.exp$mass_measured[i]
  lb <- ft_mass - err
  ub <- ft_mass + err
  
  # Order syn.exp by mass_measured
  syn.o <- syn.exp[order(syn.exp$mass_measured), ]
  
  # Filter candidates within tolerance
  pres <- syn.o[syn.o$mass_measured > lb & syn.o$mass_measured < ub, ]
  
  # Ensure output is always a data.frame with the right columns
  if (nrow(pres) == 0) {
    pres <- data.frame(
      retention_time = numeric(0),
      mass_measured  = numeric(0),
      compound_id    = integer(0),
      expid          = integer(0)
    )
  }
  
  return(pres)
}
