Smiles_Fillup<-function(base.dir,finlist,SubDB){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)
	i=1
	while (i<=dim(inp.x)[1]) {
	 if (inp.x[i,13]=="NULL") {
	  inp.x[i,13]<-"ZZZZ"
	 } else if (is.na(inp.x[i,13])) {
	  inp.x[i,13]<-"ZZZZ"
	 }
	 i=i+1
	}
	inp.x<-inp.x[order(inp.x[,3],inp.x[,13]),]
	i=2
	while (i<=dim(inp.x)[1]) {
	 j=i-1
	 if ((inp.x[i,3]==inp.x[j,3])&(inp.x[j,13]!="ZZZZ")) {
	  inp.x[i,13]<-inp.x[j,13]
	 }
	 i=i+1
	}
	inp.x<-inp.x[order(inp.x[,1]),]
	i=1
	while (i<=dim(inp.x)[1]) {
	 if (inp.x[i,13]=="ZZZZ") inp.x[i,13]<-"NULL"
	 i=i+1
	}
	write.table(inp.x,"compound.csv",sep="\t",row.names=F)
}

# Smiles_Fillup(base.dir,finlist,SubDB)