XCMSgrouping2<-function(base.dir,SubDB,exp,adjPeak=NULL,mode=NULL){
	if (is.null(adjPeak)) adjPeak=1
	if (is.null(mode)){
	 mode="neg"
	}else{
	 mode="pos"
	}
	finlist<-Get_All_DynLib_files(base.dir)
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	Exp<-paste("Exp",exp,"_",sep="")
	setwd(TMP[[5]])
	chosdir<-paste(TMP[[5]],"/",list.files()[grep(Exp,list.files())],sep="")
	setwd(chosdir)
	XCMS<-read.table("xcms2.txt",header=T,sep="\t",stringsAsFactors=F)
	ft<-TMP[[1]][[1]]
	ms2.ft<-TMP[[1]][[3]]
	XCMS4<-AutomFindingISF(XCMS,exp,ft,ms2.ft)
	XCMS4<-ConvertXCMS4(XCMS4)
	XCMS4<-SearchCombinat(XCMS4,adjPeak,mode)
	write.table(XCMS4,"xcms2.txt",sep="\t",row.names=F)	
}

# XCMSgrouping2(base.dir,SubDB="FTneg",exp=1)