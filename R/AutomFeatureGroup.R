AutomFeatureGroup<-function(stp,y10,y20,fc){
	stp<-as.integer(stp)
	dat=read.table("xcms.tsv",header=T)
	o<-order(dat$rtmed)
	dat.o<-dat[o,]
	XCMS<-data.frame()
	i=1
	repeat{
	 j=i+1
	 strt<-stp[i]+1
	 dt<-dat.o[strt:stp[j],]
	 xcms2<-FeatureGroupXCMS(dt,y10,y20,fc)
	 XCMS<-rbind(XCMS,xcms2)	
	 if(j==length(stp))break
	 i=i+1
	}
	return(XCMS)
}
