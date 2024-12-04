# creates list of breaking rows to enable the subsequent peakgrouping of a
# xcms.tsv file. 'xcms.tsv" should be present in the respective experiment folder.
# tR (e.g. 1 sec) defines the retention time difference between two subsequent
# rows that is needed to define a breaking row

BreakRowsXCMS<-function(chosdir,tR){
	setwd(chosdir)
	dat=read.table("xcms.tsv",header=T)
	o<-order(dat$rtmed)
	dat.o<-dat[o,]
	mx<-dim(dat.o)[1]
	if(mx<=1000){
	 stp<-c(1,mx)
	}else{
	 stp<-c(0)
	 i=500
	 while(i+500<mx){
	  repeat{
	   j=i+1
	   ch.1<-dat.o[i,8]
	   ch.2<-dat.o[j,8]
	   ch=ch.2-ch.1
	   if(ch>=tR){					
	    stpm<-i
	    stp<-append(stp,stpm)
	    i=i+500
	    break
	   }else{
	    i=i+1
	    if(i==mx)break
	   }
	  }
	 }
	}
	stp<-append(stp,mx)
	return(stp)
}
