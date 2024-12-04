Files_TransferNames<-function(base.dir,finlist,SBDB,assocfl=NULL){
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[1])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	ftng<-data.frame(inp.x[,1],inp.x[,3],inp.x[,10],inp.x[,13])
	ftng[,2]<-as.character(ftng[,2])
	ftng[,3]<-as.character(ftng[,3])
	ftng[,4]<-as.character(ftng[,4])
	rm(inp.x)
	colnames(ftng)<-c("COMPID","COMPNAME","SUBSID","SMILES")
	j=2
	while (j<=dim(ftng)[2]) {
	 i=1
	 while (i<=dim(ftng)[1]) {
	  if (is.na(ftng[i,2])) ftng[i,2]<-"NULL"
	  i=i+1
	 }
	 j=j+1
	}
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[2])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	ftps<-data.frame(inp.x[,1],inp.x[,3],inp.x[,10],inp.x[,13])
	ftps[,2]<-as.character(ftps[,2])
	ftps[,3]<-as.character(ftps[,3])
	ftps[,4]<-as.character(ftps[,4])
	colnames(ftps)<-c("COMPID","COMPNAME","SUBSID","SMILES")
	rm(inp.x)
	j=2
	while (j<=dim(ftps)[2]) {
	 i=1
	 while (i<=dim(ftps)[1]) {
	  if (is.na(ftps[i,2])) ftps[i,2]<-"NULL"
	  i=i+1
	 }
	 j=j+1
	}
	TMP<-SelectSubDB(base.dir,finlist,SubDB=SBDB[3])
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
	syn<-data.frame(inp.x[,1],inp.x[,3],inp.x[,10],inp.x[,13])
	syn[,2]<-as.character(syn[,2])
	syn[,3]<-as.character(syn[,3])
	syn[,4]<-as.character(syn[,4])
	colnames(syn)<-c("COMPID","COMPNAME","SUBSID","SMILES")
	rm(inp.x)
	j=2
	while (j<=dim(syn)[2]) {
	 i=1
	 while (i<=dim(syn)[1]) {
	  if (is.na(syn[i,2])) syn[i,2]<-"NULL"
	  i=i+1
	 }
	 j=j+1
	}
	if(is.null(assocfl)){
	 Assoc<-finlist[[4]]
	}else if(assocfl==2){
	 Assoc<-finlist[[10]]
	}else if(assocfl==3){
	 Assoc<-finlist[[11]]
	}
	return(list(ftng,ftps,syn,Assoc)) # ftng, ftps and syn are the subDB corresponding to the SBDB entries
}

# Names_lst<-Files_TransferNames(base.dir,finlist)
