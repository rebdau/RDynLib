GetAllInfo<-function(base.dir,finlist,SubDB,dbkey){
	if(SubDB=="FTneg"){
	 chrg=NULL
	 if(dbkey>max(finlist[[1]][[1]][,3])) return("COMPID is not in the sub-database")
	 expnr.ft<-finlist[[1]][[1]][finlist[[1]][[1]][,3]==dbkey,4]
	 Assoc<-finlist[[4]]
	 Assoc_EXPID_FTn<-strsplit(finlist[[5]][,3],",")
	 for (i in 1:length(Assoc_EXPID_FTn)) {
	  syn.sel<-which(Assoc_EXPID_FTn[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 	#i reflects the row obtained by the former for loop:
	 syn.poss<-finlist[[5]][i,2]
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		QTOFneg experiment.")
	 }
		# if you are not satisfied with the proposed experiment number in the next
		# question, you can enter another number after e.g. checking the subDBexperiment.txt file:
	 expnr.syn<-as.integer(readline(paste("Which of these exp, i.e.",syn.poss,
			   ", is the QTOF experiment:",sep=" ")))
	 ft<-finlist[[1]][[1]]
	 syn<-finlist[[3]][[1]]
	 ft.exp<-ft[ft[,4]==expnr.ft,]
	 syn.exp<-syn[syn[,4]==expnr.syn,]
	} else if (SubDB=="QTOFneg") {
	 chrg=NULL
	 if(dbkey>max(finlist[[3]][[1]][,3])) return("COMPID is not in the sub-database")
	 expnr.syn<-finlist[[3]][[1]][finlist[[3]][[1]][,3]==dbkey,4]
	 Assoc<-finlist[[4]]
	 Assoc_EXPID_Synn<-strsplit(finlist[[5]][,2],",")
	 for (i in 1:length(Assoc_EXPID_Synn)) {
	  ft.sel<-which(Assoc_EXPID_Synn[[i]]%in%as.character(expnr.syn))
	  if (length(ft.sel)!=0) break
	 }
	 ft.poss<-finlist[[5]][i,3]
	 if (ft.poss=="NULL") {
	  return("The entered expnr.syn does not correspond with a 
	 		FTneg experiment.")
	 }
	 expnr.ft<-as.integer(readline(paste("Which of these exp, i.e.",ft.poss,
			   ", is the FT experiment:",sep=" ")))
	 ft<-finlist[[1]][[1]]
	 syn<-finlist[[3]][[1]]
	 ft.exp<-ft[ft[,4]==expnr.ft,]
	 syn.exp<-syn[syn[,4]==expnr.syn,]
	} else if (SubDB=="FTpos") {
	 chrg<-"pos"
	 if(dbkey>max(finlist[[2]][[1]][,3])) return("COMPID is not in the sub-database")
	 expnr.ft<-finlist[[2]][[1]][finlist[[2]][[1]][,3]==dbkey,4]
	 Assoc<-finlist[[10]]
	 Assoc_EXPID_FTp<-strsplit(finlist[[5]][,4],",")
	 for (i in 1:length(Assoc_EXPID_FTp)) {
	  syn.sel<-which(Assoc_EXPID_FTp[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 syn.poss<-finlist[[5]][i,1]
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		QTOFpos experiment.")
	 }
	 expnr.syn<-as.integer(readline(paste("Which of these exp, i.e.",syn.poss,
			   ", is the QTOF experiment:",sep=" ")))
	 ft<-finlist[[2]][[1]]
	 syn<-finlist[[9]][[1]]
	 ft.exp<-ft[ft[,4]==expnr.ft,]
	 syn.exp<-syn[syn[,4]==expnr.syn,]
	} else if (SubDB=="QTOFpos") {
	 	 chrg<-"pos"
	 if(dbkey>max(finlist[[9]][[1]][,3])) return("COMPID is not in the sub-database")
	 expnr.syn<-finlist[[9]][[1]][finlist[[9]][[1]][,3]==dbkey,4]
	 Assoc<-finlist[[10]]
	 Assoc_EXPID_Synp<-strsplit(finlist[[5]][,1],",")
	 for (i in 1:length(Assoc_EXPID_Synp)) {
	  ft.sel<-which(Assoc_EXPID_Synp[[i]]%in%as.character(expnr.syn))
	  if (length(ft.sel)!=0) break
	 }
	 ft.poss<-finlist[[5]][i,4]
	 if (ft.poss=="NULL") {
	  return("The entered expnr.syn does not correspond with a 
	 		FTpos experiment.")
	 }
	 expnr.ft<-as.integer(readline(paste("Which of these exp, i.e.",ft.poss,
			   ", is the FT experiment:",sep=" ")))
	 ft<-finlist[[2]][[1]]
	 syn<-finlist[[9]][[1]]
	 ft.exp<-ft[ft[,4]==expnr.ft,]
	 syn.exp<-syn[syn[,4]==expnr.syn,]
	}else{
	 return("Not a valid sub-database")
	}
	return(list(expnr.ft,expnr.syn,ft.exp,syn.exp,chrg))
}