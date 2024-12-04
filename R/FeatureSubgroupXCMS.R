# performs further grouping based on the median retention time
# of each current group

FeatureSubgroupXCMS<-function(XCMS,y10){
	XCMSfine<-data.frame()
	i=2
	repeat{
	 XCMSsub<-XCMS[XCMS[,dim(XCMS)[2]]==i,]
	 tRmed<-median(XCMSsub[,8])
	 tRmed.low<-tRmed-y10
	 tRmed.high<-tRmed+y10
	 subgrp<-as.integer(rep(NA,dim(XCMSsub)[1]))
	 j=1
	 while(j<=dim(XCMSsub)[1]){
	  if(XCMSsub[j,8]<tRmed.low){
	   subgrp[j]<-1
	  }else{
	   if(XCMSsub[j,8]>tRmed.high){
	    subgrp[j]<-3
	   }else{
	    subgrp[j]<-2
	   }
	  }
	  j=j+1
	 }
	 XCMSsub<-data.frame(XCMSsub,subgrp)
	 XCMSfine<-rbind(XCMSfine,XCMSsub)
	 if(i==max(XCMS[,dim(XCMS)[2]]))break
	 i=i+1
	}
	return(XCMSfine)
}
