XCMSgrouping<-function(base.dir,SubDB,exp,tR=NULL,y10=NULL,y20=NULL,y11=NULL){
	if (is.null(tR)) tR=1
	if (is.null(y10)) y10=1
	if (is.null(y20)) y20=0.8
	if (is.null(y11)) y11=46.0
	finlist<-Get_All_DynLib_files(base.dir)
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	rm(finlist)
	Exp<-paste("Exp",exp,"_",sep="")
	setwd(TMP[[5]])
	chosdir<-paste(TMP[[5]],"/",list.files()[grep(Exp,list.files())],sep="")
	stp<-BreakRowsXCMS(chosdir,tR)
	CBR<-CheckBreakRows(stp)
	if (CBR==0) break
	procfc<-LoadProc(exp,TMP)
	fc<-procfc[[2]]
	rm(TMP)
	XCMS<-AutomFeatureGroup(stp,y10,y20,fc)
	CON.new<-GeneralizeConXCMS(XCMS)
	XCMS[,dim(XCMS)[2]]<-CON.new
	o<-order(XCMS[,dim(XCMS)[2]])
	XCMS<-XCMS[o,]
	XCMS<-FeatureSubgroupXCMS(XCMS,y10)
	XCMS<-GeneralizeSubgroupXCMS(XCMS)
	XCMS<-FeatureAnnotXCMS(XCMS,y11,fc)
	write.table(XCMS,"xcms2.txt",sep="\t",row.names=F)
	nodes<-data.frame(XCMS[,8]/60,XCMS[,5],XCMS[,1],XCMS[,14:(length(fc)+13)])
	colnames(nodes)<-c("rt(min)","mzmed","name",colnames(XCMS[,14:
						(length(fc)+13)]))
	nodes<-nodes[order(nodes[,2]),]
	nodes<-data.frame(nodes[,1:3],rep(NA,dim(nodes)[1]),nodes[,4:dim(nodes)[2]])
	write.table(nodes,"nodes.txt",sep="\t",row.names=F)
}