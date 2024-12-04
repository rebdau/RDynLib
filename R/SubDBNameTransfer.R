SubDBNameTransfer<-function(base.dir,finlist,assocfl=NULL){
	if(is.null(assocfl)){
	 SBDB<-as.character(c("FTneg","FTpos","QTOFneg"))
	 Names_lst<-Files_TransferNames(base.dir,finlist,SBDB)
	 transfertNam<-TransferNames(Names_lst,SBDB)
	}else if(assocfl==2){
	 SBDB<-as.character(c("FTpos","FTneg","QTOFpos"))
	 Names_lst<-Files_TransferNames(base.dir,finlist,SBDB,assocfl=2)
	 transfertNam<-TransferNames(Names_lst,SBDB)
	}else if(assocfl==3){
	 SBDB<-as.character(c("QTOFneg","QTOFpos","FTneg"))
	 Names_lst<-Files_TransferNames(base.dir,finlist,SBDB,assocfl=3)
	 transfertNam<-TransferNames(Names_lst,SBDB)
	}
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[1])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	inp.x[1:dim(transfertNam[[1]][[1]])[1],3]<-transfertNam[[1]][[1]][,2]
	inp.x[1:dim(transfertNam[[1]][[1]])[1],10]<-transfertNam[[1]][[1]][,3]
	inp.x[1:dim(transfertNam[[1]][[1]])[1],13]<-transfertNam[[1]][[1]][,4]
	write.table(inp.x,"compound.csv",sep="\t",row.names=F)
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[2])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	inp.x[1:dim(transfertNam[[1]][[2]])[1],3]<-transfertNam[[1]][[2]][,2]
	inp.x[1:dim(transfertNam[[1]][[2]])[1],10]<-transfertNam[[1]][[2]][,3]
	inp.x[1:dim(transfertNam[[1]][[2]])[1],13]<-transfertNam[[1]][[2]][,4]
	write.table(inp.x,"compound.csv",sep="\t",row.names=F)
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[3])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	inp.x[1:dim(transfertNam[[1]][[3]])[1],3]<-transfertNam[[1]][[3]][,2]
	inp.x[1:dim(transfertNam[[1]][[3]])[1],10]<-transfertNam[[1]][[3]][,3]
	inp.x[1:dim(transfertNam[[1]][[3]])[1],13]<-transfertNam[[1]][[3]][,4]
	write.table(inp.x,"compound.csv",sep="\t",row.names=F)
	return(transfertNam[[2]])
}

# transfertNam<-SubDBNameTransfer(base.dir,finlist)
