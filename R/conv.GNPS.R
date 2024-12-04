# No rows representing the different COMPIDs of the experiment have to be
# added to gnps_add.txt, this is done while running the function. In case
# no conversions are found for a particular COMPID, a row of NAs is added
# preceeded by the name of the COMPID, the function jumps then to the next
# COMPID. However, at the end of the COMPID list belonging to the experiment,
# some final expected rows might be missing whenever no conversions were
# found for the few last COMPIDs. This is not a problem as these last
# COMPIDs followed by NAs will be added upon running the next experiment.
# Note: Conversions for a particular COMPID are added to the row with row
# number equal to the COMPID.

conv.GNPS<-function(base.dir,finlist,SubDB,Prod.exp,inp.x=NULL,peakwidth=NULL,
				mzerr=NULL,min=NULL,adduct=NULL,thr1=NULL,
				thr2=NULL,thr3=NULL,gnps_add=NULL){
	if (is.null(peakwidth)) peakwidth<-0.2
	if (is.null(adduct)) adduct<-46.0055 
		# adduct=formic acid by default, acetic acid=60.0211 Da
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
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
	Prod.dat<-inp.x[inp.x[,9]==Prod.exp,]
		# sort on retention time
	Prod.dat<-Prod.dat[order(Prod.dat[,8]),]
	if (is.null(gnps_add)) {
	 setwd(TMP[[4]])		
	 gnps_add<-read.table("gnps_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	}
	i=1
	while(i<=dim(Prod.dat)[1]){
	 gnps.df<-data.frame(COMPID.sub=integer(),MZ.sub=numeric(),
				IONS.sub=integer(),COMPID.prod=integer(),
				MZ.prod=numeric(),IONS.prod=integer(),
				COMMON_IONS=integer(),DOT_IONS=numeric(),
				COMMON_LOSS=integer(),DOT_LOSS=numeric(),
				FORW_IONS=numeric(),REV_IONS=numeric(),
				FORW_LOSS=numeric(),REV_LOSS=numeric(),
				stringsAsFactors=FALSE)
	 sub.tR<-Prod.dat[i,8]+peakwidth
	 sub.mz<-Prod.dat[i,5]
	 sub.lowmz<-sub.mz-mzerr
	 sub.highmz<-sub.mz+1.0034+mzerr
	 forb.lowmz<-sub.mz+adduct-mzerr
	 forb.highmz<-sub.mz+adduct+1.0034+mzerr
	 j=i+1
	 while(j<=dim(Prod.dat)[1]){
		# selection on higher retention time means that GNPS 
		# conversions will only be found in one direction!
	  if(Prod.dat[j,8]>sub.tR){	
	   if((Prod.dat[j,5]<sub.lowmz)|(Prod.dat[j,5]>sub.highmz)){
	    if((Prod.dat[j,5]<forb.lowmz)|(Prod.dat[j,5]>forb.highmz)){
	     out<-targMS2comp(Prod.dat[i,1],Prod.dat[j,1],subDB,AnalMS,
					IntThres=min)
		# thresh1 is the average of the ions and losses in common
	     thresh1<-(out[1,7]+out[1,9])/2
		# thresh2 is the contribution of the average number of ions
		# losses to the total number of ions in the MS2 spectrum
		# with the lowest number of ions
	     thresh2<-thresh1/(min(c(out[1,3],out[1,6])))
		# thresh3 average dot product of the ions/losses
	     thresh3<-(out[1,8]+out[1,10])/2
	     if(is.na(thresh1)|is.na(thresh2)|is.na(thresh3)){
	      j=j+1
	      next
	     }
	     if((thresh1>thr1)&(thresh2>thr2)&(thresh3>thr3)){
		prnt1<-paste(out[1,1],out[1,2],out[1,4],out[1,5],sep=" ")
#	      print(prnt1)
	      gnps.df<-rbind(gnps.df,out)
	     }
	    }
	   }
	  }
	  j=j+1
	 }
	 if(dim(gnps.df)[1]!=0){
	  gnps_add<-rank.GNPS(gnps.df,gnps_add,Prod.dat,inp.x)
	 }
	 i=i+1
	}
	gnps_add[,1]<-inp.x[1:(dim(gnps_add)[1]),1]
	return(gnps_add)
}

# gnps_add<-conv.GNPS(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(gnps_add,"gnps_add.txt",sep="\t",row.names=F)
