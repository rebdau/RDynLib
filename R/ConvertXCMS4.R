ConvertXCMS4<-function(XCMS4){
	MeAb<-as.numeric(rep(NA,dim(XCMS4)[1]))
#	maxPkGrp<-max(XCMS4[,dim(XCMS4)[2]-4])
	Singlet<-which(as.integer(XCMS4[,dim(XCMS4)[2]-4])==1)
	MeAb[Singlet]<-XCMS4[Singlet,8]
	i=2
	while(length(which(XCMS4[,dim(XCMS4)[2]-4]%in%i))!=0){
	 PkGrp<-which(XCMS4[,dim(XCMS4)[2]-4]%in%i)
	 RtAve<-mean(XCMS4[PkGrp,8])
	 MeAb[PkGrp]<-RtAve
	 i=i+1
	}
	XCMS4<-data.frame(XCMS4,MeAb)
	XCMS4<-XCMS4[order(XCMS4[,dim(XCMS4)[2]]),]
	CorrPkGrp<-as.integer(rep(NA,dim(XCMS4)[1]))
	CorrPkGrp[1:2]<-c(1,2)
	i=2
	while(i<dim(XCMS4)[1]){
	 j=i+1
	 if((XCMS4[j,dim(XCMS4)[2]-5]==1)|
			(XCMS4[j,dim(XCMS4)[2]-5]!=XCMS4[i,dim(XCMS4)[2]-5])){
	  CorrPkGrp[j]<-CorrPkGrp[i]+1
	 }else{
	  CorrPkGrp[j]<-CorrPkGrp[i]
	 }
	 i=i+1
	}
	XCMS4<-data.frame(XCMS4,CorrPkGrp)
	return(XCMS4)
}

# XCMS4<-read.table("xcms4.txt",header=T,stringsAsFactor=F)