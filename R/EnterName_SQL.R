






EnterName_SQL <- function(data_path, CoMPID, SubDB){
 
  data_con   <- dbConnect(SQLite(), data_path)
  cmp_data  <- dbGetQuery(data_con, paste0("SELECT * FROM ms_compound where ",
                                       "compound_id = ", CoMPID))
  CmpCSV <- cmp_data[, c("retention_time", "mass_measured", "compound_id",
                       "expid", "name", "smiles", "nodename", "formula",
                       "ppm_deviation")]
  
  CmpCSV<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)
  # compound.csv might be a shrunk version of the initial database, just containing
  # the selected experiment, then, row number is not equal to the COMPID:
  dbkey<-which(as.integer(CmpCSV["compound_id"])%in%as.integer(dbkey))
  #ConcatName<-read.table("compound_name.txt",header=T,sep="\t",stringsAsFactors=F)
  print("Currently, this info for your dbkey is present:")
  print(paste("COMPID: ",CmpCSV[dbkey, "compound_id"],sep=""))
  print(paste("COMPNAME: ",CmpCSV[dbkey, "name"],sep=""))
  print(paste("FORMULA: ",CmpCSV[dbkey, "formula"],sep=""))
  print(paste("MASS_MEASURED: ",CmpCSV[dbkey, "mass_measured"],sep=""))
  print(paste("PPM_DEVIATION: ",CmpCSV[dbkey, "ppm_deviation"],sep=""))
  print(paste("RETENTION_TIME: ",CmpCSV[dbkey, "retention_time"],sep=""))
  print(paste("EXPID: ",CmpCSV[dbkey, "expid"],sep=""))
  print(paste("SMILES: ",CmpCSV[dbkey, "smiles"],sep=""))
  rm_sel=""
  while (rm_sel=="") rm_sel<-toupper(c(readline("Initials first and last name:")))
  if (nchar(rm_sel)!=2) {
    print("Should be 2 characters long.")
    return()
  }
  print("Info regarding structural confidence of the given name:")
  print("ID, identified via NMR or spiking of standard")
  print("AN, annotated, reasonably certain that the given name is correct")
  print("KA, characterized, highly certain that the main structure invoked by the name is ok,")
  print("    but e.g. moiety linkages or stereomer configuration might be wrong")
  print("PU, putative, uncertain - might as well be another structural isomer")
  rm_sel1=""
  while (rm_sel1=="") rm_sel1<-toupper(c(readline("Structural confidence level:")))
  if (nchar(rm_sel1)!=2) {
    print("Should be 2 characters long.")
    return()
  }
  rm_sel1<-paste(rm_sel1," ",sep="")
  #print("Give structural moieties separated by ' + ' as a concatenated name.")
  #rm_selQ<-readline("Concatenated name:")
  print("For the following name, replace comma by underscore and primes by the abbreviation pr")
  rm_sel2=""
  while (rm_sel2=="") rm_sel2<-tolower(c(readline("Name:")))
  rm_sel2<-paste(rm_sel2," ",sep="")
  print("Is it a buffer adduct? E.g. acetate or formate? If not, press enter.")
  rm_sel3<-substr(toupper(c(readline("Adduct:"))),1,4)
  if (rm_sel3!="") rm_sel3<-paste(rm_sel3," ",sep="")
  print("Is it an isotope, e.g. 13C, 34S, 37Cl, 13C34S? If not, press enter.")
  rm_sel4<-c(readline("Isotope:"))
  if (rm_sel4!="") rm_sel4<-paste(rm_sel4," ",sep="")
  print("Do you think it is an in-source fragment? If not, press enter, else say yes")
  rm_sel5<-c(readline("In-source fragment:"))
  if (rm_sel5!="") rm_sel5<-"ISF "
  print("Do you think it is an odd electron ion? If not, press enter, else say yes")
  rm_sel6<-c(readline("Odd electron ion:"))
  if (rm_sel6!="") rm_sel6<-"OE "
  print("Is it a complex between the neutral molecule and its ion, i.e. a 'dimer'? If not, press enter.")
  rm_sel8<-c(readline("Ion-neutral complex of molecule:"))
  if (rm_sel8!="") rm_sel8<-"DIM "
  print("Is it a complex between (a) neutral(s) and an ion of different species, i.e. a 'heteromer'? If not, press enter.")
  rm_sel9<-c(readline("Ion-neutral complex:"))
  if (rm_sel9!="") rm_sel9<-"ADD "
  
  rm_sel7<-c(readline("Is the name referring to part of the molecule? (y/n)"))
  if (rm_sel7=="n") {
    rmsl<-paste(rm_sel2,rm_sel6,rm_sel5,rm_sel3,rm_sel9,rm_sel8,rm_sel4,rm_sel1,rm_sel,sep="")
  } else {
    rmsl<-paste("[",rm_sel2,rm_sel6,rm_sel5,rm_sel3,rm_sel9,rm_sel8,rm_sel4,rm_sel1,rm_sel,"]",sep="")
  }
  print("If a formula should be entered, add it here. Else press enter.")
  formul<-gsub(" ","",toupper(c(readline("Formula:"))),fixed=T)
  print("If a SMILES structure should be entered, add it here. Else press enter.")
  smiles<-c(readline("SMILES:"))
  print("If the ppm deviation is known, add it here. Else press enter.")
  ppm<-c(readline("PPM:"))
  print("New data to be entered are:")
  print(paste("COMPNAME: ",rmsl,sep=""))
  print(paste("FORMULA: ",formul,sep=""))
  print(paste("PPM_DEVIATION: ",ppm,sep=""))
  print(paste("SMILES: ",smiles,sep=""))
  print("Are you satisfied? If not, press enter. Else say 'yes'")
  SATIS<-c(readline("Enter the data?"))
  
  #Ahlam : what is termen? the switch of directories? 
  #termen.txt was a list of random messages, in Flemish, telling the user that *again* he/she made a mistake and will have to try again.
  #it was meant as a joke, but in the mean time a necessary escape route if the entered information was not correct.
  #I replaced it here by a sinlge English message.
  if (SATIS!="yes"){
    message <- "Oops! Looks like you made a mistake. Don’t worry, you can try again by calling EnterName_SQL(data_path, CoMPID, SubDB) again."
    return(message)
  }else{
    setwd(TMP[[3]])	 
    CmpCSV[dbkey,3]<-rmsl
    if (formul!="") CmpCSV[dbkey, "formula"]<-formul
    if (ppm!="") CmpCSV[dbkey, "ppm_deviation"]<-ppm
    if (smiles!="") CmpCSV[dbkey, "smiles"]<-smiles # nog functie om te checken of SMILES string overeenkomt met chemische formula
    write.table(CmpCSV,"compound.csv",sep="\t",row.names=F)
    if(SubDB=="FTneg"){
      CmpCSV[dbkey, "name"]<-rmsl
      CmpCSV[dbkey, "smiles"]<-smiles
    } else if (SubDB=="FTpos") {
      CmpCSV[dbkey, "name"]<-rmsl
      CmpCSV[dbkey, "smiles"]<-smiles
    } else if (SubDB=="QTOFneg") {
      CmpCSV[dbkey, "name"]<-rmsl
      CmpCSV[dbkey, "smiles"]<-smiles
    } else if (SubDB=="QTOFpos") {
      CmpCSV[dbkey, "name"]<-rmsl
      CmpCSV[dbkey, "smiles"]<-smiles
    }
    #"database_FTMS_neg\CSV_add\compound_name.txt" was meant to contain conservative structures: an enumeration of substructures without specifying the links between them.
    #We can leave it out for now
    #setwd(TMP[[4]])
    #ConcatName[dbkey,2]<-rm_selQ
    #write.table(ConcatName,"compound_name.txt",sep="\t",row.names=F)
  }
}

# EnterName(base.dir,finlist,dbkey,SubDB)

