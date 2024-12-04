MultExp_CSPP<-function(base.dir,finlist,SubDB,startexp,stopexp,inp.x=NULL,
				peakwidth=NULL,mzerr=NULL){
	if (is.null(peakwidth)) peakwidth<-0.2
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	AnalMS<-TMP[[2]]
	if (is.null(mzerr)) {
	 if (AnalMS=="FT") {
	  mzerr<-0.01
	 } else if (AnalMS=="QTOF"){	
	  mzerr<-0.015
	 }
	}
	if (is.null(inp.x)) {
	 setwd(TMP[[3]])
	 inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
      	# has problems when primes are used, replace them by pr
	}
	setwd(TMP[[4]])
	comp_add<-read.table("compound_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	i=startexp
	repeat{
	 Prod.exp=i
	 comp_add<-cspp.tot(base.dir,finlist,SubDB,Prod.exp,inp.x,peakwidth,
					mzerr,comp_add)
	 if (i==stopexp) break
	 i=i+1
	}
	return(comp_add)
}

# comp_add<-MultExp_CSPP(base.dir,finlist,SubDB="FTneg",startexp=1,stopexp=4)
# write.table(comp_add,"compound_add.txt",sep="\t",row.names=F)