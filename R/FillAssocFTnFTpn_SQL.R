FillAssocFTnFTpn_SQL(
  FTn_con, FTp_con, Assoc, FTn_expnr, FTp_expnr, cutoff = 1, rg = rg, lc.err, 
  err, minIon, polarity_ftn, polarity_ftp, FTn_path, FTp_path){
  
  
  # Load FT and QTOF compounds
  ftn.exp <- dbGetQuery(FTn_con, sprintf(
    "SELECT retention_time, mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", FTn_expnr))
  
  ftp.exp <- dbGetQuery(FTp_con, sprintf(
    "SELECT retention_time , mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", FTp_expnr))
  
  o<-order(ftn.exp$mass_measured)
  ftng.o<-ftn.exp[o,]
  o<-order(ftp.exp$mass_measured)
  ftps.o<-ftp.exp[o,]
  i=1
  repeat{
    mass.ftng<-ftng.o$mass_measured[i]
    time.ftng<-ftng.o$retention_time[i]
    COMPID.ftng<-ftng.o$compound_id[i]
    pres<-Find_pos_compound(mass.ftng,time.ftng,ftps.o,err,lc.err,rg)
    if(dim(pres)[1]==0){
      if(i==dim(ftng.o)[1])break
      i=i+1
      next
    }else{
      Assoc<-Sort_comp_matches(pres,Assoc,time.ftng,COMPID.ftng,rg)
    }
    if(i==dim(ftng.o)[1])break
    i=i+1
  }
  return(Assoc)
 
  
  
}