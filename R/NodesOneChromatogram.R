# Go to the directory that contains the .cdf file that has to be processed.

NodesOneChromatogram<-function(base.dir,finlist,SubDB,xcdf,exp.id,profstep=NULL,
					fwhm=NULL,snthresh=NULL,step=NULL,steps=NULL){
	if (is.null(profstep)) pr=0.1
	if (is.null(fwhm)) f=8
	if (is.null(snthresh)) sn=2
	if (is.null(step)) stp=0.01
	if (is.null(steps)) stps=2
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	setwd(TMP[[5]])
	Exp<-paste("Exp",exp.id,"_",sep="")
	chosdir<-paste(TMP[[5]],"/",list.files()[grep(Exp,list.files())],sep="")
	setwd(chosdir)
	library(xcms)
	xraw<-xcmsRaw(xcdf,profstep=pr,profmethod="bin")
	xpeaks<-findPeaks(xraw,fwhm=f,max=300,snthresh=sn,step=stp,steps=stps)
	rt_min<-xpeaks[,4]/60
	mzmed<-as.character(c())
	i=1
	while(i<=dim(xpeaks)[1]){	
	 mzmed[i]<-paste("M",round(xpeaks[i,1],0),"T",round(xpeaks[i,4],0),sep="")
	 i=i+1
	}
	options(stringsAsFactors=F)
	EmpCol<-as.character(rep("NULL",dim(xpeaks)[1]))
	nodes<-data.frame(rt_min,xpeaks[,1],mzmed,EmpCol,xpeaks[,7:8])
	colname<-c("rt(min)","mzmed","name","","into","intf")
	colnames(nodes)<-colname
	setwd(TMP[[5]])
	chosdir<-paste(TMP[[5]],"/",list.files()[grep(Exp,list.files())],sep="")
	setwd(chosdir)
	write.table(nodes,"nodes.txt",sep="\t",row.names=F)
}

# NodesOneChromatogram(base.dir,finlist,SubDB="FTneg",xcdf="popB_HSS3sl.cdf",exp.id=178)


