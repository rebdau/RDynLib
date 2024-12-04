MultExp_GNPS<-function(base.dir,finlist,SubDB,startexp,stopexp,inp.x=NULL,
				peakwidth=NULL,mzerr=NULL,min=NULL,adduct=NULL,
				thr1=NULL,thr2=NULL,thr3=NULL){
	if (is.null(peakwidth)) peakwidth<-0.2
	if (is.null(adduct)) adduct<-46.0055 
		# adduct=formic acid by default, acetic acid=60.0211 Da
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	AnalMS<-TMP[[2]]
	if (is.null(mzerr)) {
	 if (AnalMS=="FT") {
	  mzerr<-0.01
	  min<-100
	 } else if (AnalMS=="QTOF"){	
	  mzerr<-0.015
	  min<-5
	 }
	}
	if (is.null(thr1)) thr1=3
	if (is.null(thr2)) thr2=0.1
	if (is.null(thr3)) thr3=0.4
	if (is.null(inp.x)) {
	 setwd(TMP[[3]])
	 inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
      	# has problems when primes are used, replace them by pr
	}
	setwd(TMP[[4]])
	gnps_add<-read.table("gnps_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	i=startexp
	repeat{
	 Prod.exp=i
	 gnps_add<-conv.GNPS(base.dir,finlist,SubDB,Prod.exp,inp.x,peakwidth,
					mzerr,min,adduct,thr1,thr2,thr3,gnps_add)
	 if (i==stopexp) break
	 i=i+1
	}
	return(gnps_add)
}

# gnps_add<-MultExp_GNPS(base.dir,finlist,SubDB="FTneg",startexp=1,stopexp=3)
# write.table(gnps_add,"gnps_add.txt",sep="\t",row.names=F)