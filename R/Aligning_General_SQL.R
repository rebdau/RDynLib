Aligning_General_SQL<-function(FT_path, QTOF_path, expnr.ft, expnr.syn, err, t.ini){
  library(RSQLite)
  con_ft  <- dbConnect(SQLite(), FT_path)
  con_syn <- dbConnect(SQLite(), QTOF_path)
  
  syn_sql <- dbGetQuery(con_syn, "SELECT * FROM ms_compound")
  ft_sql  <- dbGetQuery(con_ft,  "SELECT * FROM ms_compound")
  
  # Reorder columns to match original format
  ft <- ft_sql[, c("retention_time", "mass_measured", "compound_id", "expid", "name", "smiles", "nodename")]
  syn <- syn_sql[, c("retention_time", "mass_measured", "compound_id", "expid", "name", "smiles", "nodename")]
  
  ft.exp<-ft[ft[,4]==expnr.ft,]
  syn.exp<-syn[syn[,4]==expnr.syn,]
  ft.exp.o<-ft.exp[order(ft.exp[,1]),]
  COMPID<-ft.exp.o[dim(ft.exp.o)[1]/2,3]
  reg<-dim(ft.exp.o)[1]
  if (reg%%2!=0) reg=reg-1
  ft.sh<-Make_loc_shortlist(COMPID,reg,ft.exp)
  ft.sh <- Make_loc_shortlist(COMPID, reg, ft.exp)
  
  #ESImode=NULL means that both chromatograms are obtained in negative
  #mode. Thus, if the ESImode is given as e.g. pos, the m/z values of
  #the corresponding chromatogram should be corrected to be able to
  #match with those of the key chromatogram.
  # ESImode adjustment
  if (!is.null(ESImode)) {
    syn.exp[,"m/z"] <- syn.exp[,"m/z"] - 2.01456
  }
  
  ft.sh <- Remov_empties(err, syn.exp, ft.sh, reg)
  
  LCal<-LCalign(err,t.ini,syn.exp,ft.sh)
  
  
  if (!is.null(ESImode)) {
    LCal[,5]=as.character(as.numeric(LCal[,5])+2.01456)
  }
  return(LCal)
}