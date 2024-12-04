EnteringSmiles<-function(base.dir,finlist,SubDB){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	setwd(TMP[[3]])
	comp<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)
	library(rcdk)
	wdnew<-paste(TMP[[6]],"/mol files",sep="")
	setwd(wdnew)
	struc<-c(list.files())
	mols<-load.molecules(molfiles=struc[1:length(struc)]) # Note that the full path to the files should be provided. 
	i=1
	while (i<=length(mols)) {
	 compid<-as.integer(unlist(strsplit(struc[i],"[.]"))[1])
	 compid_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(compid))
	 comp[compid_row,13]<-get.smiles(mols[[i]])
	 i=i+1
	}
	write.table(comp,"compound.csv",sep="\t",row.names=F)
}

# EnteringSmiles(base.dir,finlist,SubDB="FTneg")