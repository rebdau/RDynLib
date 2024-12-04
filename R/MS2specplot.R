MS2specplot<-function(dbkey,base.dir,finlist,SubDB,TMP=NULL,single=NULL,
				prcx=NULL){
	if(is.null(single)){
	 oldpar<-par(no.readonly=T)
	 par(mfrow=c(1,1))
	}
	if(is.null(TMP))TMP<-SelectSubDB(base.dir,finlist,SubDB)
	if(is.null(prcx))prcx=0.7
	MS2plot(dbkey,TMP,prcx)
}

# MS2specplot(dbkey=2459,base.dir,finlist,SubDB="FTneg")