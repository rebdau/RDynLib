# Can only be done when the same rows in the compound.csv and compound_add.txt
# files refer to the same COMPID

CSPPnameIncorp<-function(base.dir,finlist,SubDB,thr1=NULL,thr2=NULL){
	if (is.null(thr1)) thr1=0.9
	if (is.null(thr2)) thr2=5
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	setwd(TMP[[4]])
	comp_add<-read.table("compound_add.txt",header=T,sep="\t",
					stringsAsFactors=F)
	inp.x<-CSPPintoDynLib(inp.x,comp_add,thr1,thr2)
	msdet<-File_CSPPname(inp.x)
	msdet<-CSPPnameTransfer(msdet)
	inp.x[1:dim(msdet)[1],3]<-msdet[,2]
	return(inp.x)
}

# inp.x<-CSPPnameIncorp(base.dir,finlist,SubDB)