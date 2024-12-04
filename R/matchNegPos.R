#Remove false positive negFTm/z,posFTm/z peak pairs by 
#checking whether a certain percentage of the neg MS2  
#spectral peaks can be traced in the pos MS2 spectrum. 
#By default, 20% of the neg MS2 peaks are traced in the
#pos MS2 spectrum.

matchNegPos<-function(LCal,finlist,minIon=NULL,chrg=NULL){
	if (is.null(minIon)) minIon=0.2
	if (is.null(chrg)) {
	 ms2.ft<-finlist[[1]][[3]]
	 ms2.ftps<-finlist[[2]][[3]]
	} else {
	 ms2.ft<-finlist[[3]][[3]]
	 ms2.ftps<-finlist[[9]][[3]]	
	}
	i=1
	while(i<=dim(LCal)[1]){
	 ft.compid=LCal[i,1]
	 ftps.compid=LCal[i,7]
	 ms2ion=as.integer(strsplit(ms2.ft[[as.integer(ft.compid)]],",")[[1]])
	 posion=as.integer(strsplit(ms2.ftps[[as.integer(ftps.compid)]],",")[[1]])
	 msmsion=posion-2
	 same_ion=which(ms2ion%in%msmsion)
	 if(length(same_ion)/length(ms2ion)<minIon){
	  LCal<-LCal[-i,]
	  next
	 }
	 i=i+1
	}
	return(unique(LCal))
}
