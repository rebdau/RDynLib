#Remove false positive FTm/z,QTOFm/z peak pairs by checking
#whether a certain percentage of the MS2 spectral peaks
#can be traced in the MS/MS spectrum. By default, 60% of 
#the MS2 peaks should be traced in the MS/MS spectrum as
#the MS2 spectrum is much less crowded.

matchFTSyn<-function(LCal,finlist,minIon=NULL,chrg=NULL){
	if (is.null(minIon)) minIon=0.6
	if (is.null(chrg)) {
	 ms2.ft<-finlist[[1]][[3]]
	 ms2.syn<-finlist[[3]][[3]]
	 subdb_ft<-finlist[[1]][[1]]
	 subdb_syn<-finlist[[3]][[1]]
	} else {
	 ms2.ft<-finlist[[2]][[3]]
	 ms2.syn<-finlist[[9]][[3]]
	 subdb_ft<-finlist[[2]][[1]]
	 subdb_syn<-finlist[[9]][[1]]
	}
	i=1
	while(i<=dim(LCal)[1]){
	 ft.compid=LCal[i,1]
	 ft.compid_row<-which(as.integer(subdb_ft[,3])%in%as.integer(ft.compid))
	 syn.compid=LCal[i,7]
	 syn.compid_row<-which(as.integer(subdb_syn[,3])%in%as.integer(syn.compid))
	 ms2ion=as.integer(strsplit(ms2.ft[[as.integer(ft.compid_row)]],",")[[1]])
	 msmsion=unique(as.integer(strsplit(ms2.syn[[as.integer(syn.compid_row)]],",")[[1]]))
	 same_ion=which(ms2ion%in%msmsion)
	 if(length(same_ion)/length(ms2ion)<minIon){
	  LCal<-LCal[-i,]
	  next
	 }
	 i=i+1
	}
	return(unique(LCal))
}

