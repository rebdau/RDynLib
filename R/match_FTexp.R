match_FTexp<-function(LCal,subDB,AnalMS,thr1=NULL,thr2=NULL,IntThres=NULL){
		# how many of the ions in spectrum 1 are found in 
		# spectrum 2: between 0 and 1
	if (is.null(thr1)) thr1=1
		# dot product of the ions in common
	if (is.null(thr2)) thr2=0.9
	if (AnalMS=="FT") {
	 if (is.null(IntThres)) IntThres<-100
	} else {
	 if (is.null(IntThres)) IntThres<-5
	}	
	i=1
	while(i<=dim(LCal)[1]){
	 dbkey1=as.integer(LCal[i,1])
	 dbkey2=as.integer(LCal[i,7])
	 SameMS2<-out2.targMS2comp(dbkey1,dbkey2,subDB,AnalMS,thr1,thr2,IntThres=IntThres)
	 if(is.na(SameMS2)){
	  LCal<-LCal[-i,]
	  next
	 }
	 i=i+1
	}
	return(unique(LCal))
}
