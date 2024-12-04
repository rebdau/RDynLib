MS3specplot<-function(dbkey,base.dir,finlist,SubDB,TMP=NULL,single=NULL,
				prcx=NULL,wh=NULL){
	if(is.null(single)){
	 oldpar<-par(no.readonly=T)
	 par(mfrow=c(1,1))
	}
	if(is.null(TMP))TMP<-SelectSubDB(base.dir,finlist,SubDB)
	if(is.null(prcx))prcx=0.7
	if(is.null(wh))wh=1
	MS3plot(dbkey,TMP,prcx,wh)
}

# MS3specplot(dbkey=2459,base.dir,finlist,SubDB="FTneg",wh=2)