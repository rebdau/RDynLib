NodenameMass<-function(base.dir,finlist,SubDB,exp,err=NULL){
	if (is.null(err)) err=0.02
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	TMP.tmp<-TMP[[1]][[1]][TMP[[1]][[1]][,4]%in%exp,]
	nodnam<-readline("What is the NODENAME?")
	if(nodnam!=""){
	 return(TMP.tmp[TMP.tmp[,7]%in%nodnam,])
	}
	compnam<-readline("Give part of the trivial name:")
	if(compnam!=""){
	 sl<-c()
	 i=1
	 while(i<=dim(TMP.tmp)[1]){
	  if(grepl(compnam,TMP.tmp[i,5],fixed=TRUE)) sl<-append(sl,i)
	  i=i+1
	 }
	 return(TMP.tmp[sl,])
	}
	masss<-readline("What is the m/z value?")
	if(masss!=""){
	 lowmass=as.numeric(masss)-err
	 highmass=as.numeric(masss)+err
	 TMP.tmp<-TMP.tmp[order(TMP.tmp[,2]),]
	 sl<-c()
	 i=1
	 while(i<=dim(TMP.tmp)[1]){
	  if(TMP.tmp[i,2]<=lowmass){
	   i=i+1
	   next
	  } else if(TMP.tmp[i,2]<=highmass) {
	   sl<-append(sl,i)
	   i=i+1
	   next
	  } else if(TMP.tmp[i,2]>highmass) {
	   break
	  }
	 }
	 return(TMP.tmp[sl,])
	}
}

# NodenameMass(base.dir,finlist,SubDB,exp=197)




