WithinExpAligned<-function(dbkey,TMP){
	setwd(TMP[[4]])
	AlExp<-read.table("RetT_alignment.txt",header=T,sep="\t",stringsAsFactors=F)
	COMPIDcollection<-as.integer(rep(NA,1000))
	if (!is.na(AlExp[AlExp[,1]==dbkey,2])) COMPIDcollection[1]<-as.integer(AlExp[AlExp[,1]==dbkey,2])
	sl<-which(AlExp[,2]==dbkey)
	if(length(sl)!=0){
	 for (i in 1:length(sl)) COMPIDcollection[1+i]<-as.integer(AlExp[sl[i],1])
 	}
	COMPIDcollection<-COMPIDcollection[!is.na(COMPIDcollection)]
	COMPIDcollection<-unique(COMPIDcollection[order(COMPIDcollection)])
	return(COMPIDcollection)
}