dynamCIDplot<-function(dbkey1,dbkey2,base.dir,finlist,LCal,rg,pres1,pres2,
				z,v,TMP1,TMP2,TMP3,TMP4,SubDB,Assoc,Assoc2,prcx=NULL){
	if (is.null(prcx)) prcx<-0.6
	oldpar<-par(no.readonly=T)
	par(mfrow=c(3,3))
	par(mfg=c(1,1,3,3))
	MS2plot(dbkey=dbkey1,TMP=TMP1,prcx=prcx)	# MS2 plot of FT subDB
	par(mfg=c(2,1,3,3))
	MS2plot(dbkey=dbkey2,TMP=TMP2,prcx=prcx)	# MS/MS plot of QTOF subDB
	dbkey1_row<-which(as.integer(TMP1[[1]][[1]][,3])%in%as.integer(dbkey1))
	if(TMP1[[1]][[10]][[dbkey1_row]]!="NULL"){
	 msnle<-length(as.integer(unlist(strsplit(TMP1[[1]][[10]][[dbkey1_row]],split=","))))
	 if (msnle>2) msnle<-2
	 ms3i=1
	 while(ms3i<=msnle){
	  par(mfg=c(1,ms3i+1,3,3))
	  MS3specplot(dbkey=dbkey1,base.dir,finlist,SubDB=SubDB,TMP=TMP1,single=single,
				prcx=prcx,wh=ms3i)
	  ms3i=ms3i+1
	 }
	}
	dbalt1<-which(Assoc[,1]%in%dbkey1)
	dbkey3<-Assoc[dbalt1,2]
	dbalt2<-which(Assoc2[,1]%in%dbkey3)
	dbkey4<-Assoc2[dbalt2,3]
	if (length(dbkey4)!=0) {
	 if(!is.na(dbkey4)){
	  par(mfg=c(2,3,3,3))
	  MS2plot(dbkey=dbkey4,TMP=TMP4,prcx=prcx)	# MS/MS plot of QTOF subDB
	 }
	}
	if(length(dbkey3)!=0){
	 if(!is.na(dbkey3)){
	  par(mfg=c(3,1,3,3))
	  MS2plot(dbkey=dbkey3,TMP=TMP3,prcx=prcx)	# MS2 plot of FT subDB
	  dbkey3_row<-which(as.integer(TMP3[[1]][[1]][,3])%in%as.integer(dbkey3))
	  if(TMP3[[1]][[10]][[dbkey3_row]]!="NULL"){
	   msnle<-length(as.integer(unlist(strsplit(TMP3[[1]][[10]][[dbkey3_row]],split=","))))
	   if (msnle>2) msnle<-2
	   ms3i=1
	   while(ms3i<=msnle){
	    par(mfg=c(3,ms3i+1,3,3))
	    MS3specplot(dbkey=dbkey3,base.dir,finlist,SubDB=SubDB,TMP=TMP3,single=single,
	  			prcx=prcx,wh=ms3i)
	    ms3i=ms3i+1
	   }
	  }
	 }
	}
	par(mfg=c(2,2,3,3))
	Dynam_plotPie(LCal,rg,z,v,pres1,pres2)
}