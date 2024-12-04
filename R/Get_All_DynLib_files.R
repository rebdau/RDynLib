# e.g. base.dir<-"C:/kris/werk/databases/R scripts DynLib July2017"

Get_All_DynLib_files<-function(base.dir){
	ft.dir<-paste(base.dir,"/database_FTMS_neg/CSV",sep="")
	finlist_FTn<-Get_FT_DynLib_files(ft.dir)
	ftps.dir<-paste(base.dir,"/database_FTMS_pos/CSV",sep="")
	finlist_FTp<-Get_FT_DynLib_files(ftps.dir)
	syn.dir<-paste(base.dir,"/database_QTOF_neg/CSV",sep="")
	finlist_Synn<-Get_Synapt_DynLib_files(syn.dir)
	synps.dir<-paste(base.dir,"/database_QTOF_pos/CSV",sep="")
	finlist_Synp<-Get_Synapt_DynLib_files(synps.dir)
	DynLibassoc.dir<-paste(base.dir,"/DynLib subDB alignment",sep="")
	setwd(DynLibassoc.dir)
	Assoc_COMPID<-read.table("DynLibDBassociation.txt",header=T,sep="\t")
	Assoc_COMPID2<-read.table("DynLibDBassociation2.txt",header=T,sep="\t")
	Assoc_COMPID3<-read.table("DynLibDBassociation3.txt",header=T,sep="\t")
	Assoc_EXPID<-read.table("subDBexperiment.txt",header=T,sep="\t",
				stringsAsFactor=F)
	MS1feat<-read.table("MS1features.txt",header=T,sep="\t",
				stringsAsFactor=F)
	neut<-read.table("neutral_loss_candid.txt",header=T,sep="\t",
				stringsAsFactor=F)
	frag<-read.table("product_ion_candid.txt",header=T,sep="\t",
				stringsAsFactor=F)
	return(list(finlist_FTn,finlist_FTp,finlist_Synn,Assoc_COMPID,
				Assoc_EXPID,MS1feat,neut,frag,finlist_Synp,
				Assoc_COMPID2,Assoc_COMPID3))
}

# finlist<-Get_All_DynLib_files(base.dir)