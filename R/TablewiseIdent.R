TablewiseIdent<-function(base.dir,finlist,SubDB,nOc=NULL){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	AnalMS<-TMP[[2]]
	setwd(TMP[[3]])
	inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
      	# has problems when primes are used, replace them by pr
	if (is.null(nOc)) nOc<-dim(inp.x)[1]
	reren<-MS2simOfIdentMasses(inp.x,finlist,subDB,AnalMS,nOc)
	for (i in 1:dim(reren)[1]){
	 known<-reren[i,3]
	 unknown<-reren[i,1]
	 inp.x[unknown,3]<-inp.x[known,3]	# COMPNAME
	 inp.x[unknown,4]<-inp.x[known,4]	# FORMULA
	 inp.x[unknown,13]<-inp.x[known,13] # SMILES
	 if (substr(reren[i,4],1,4)=="CSPP") {	# SUBSID
	  inp.x[unknown,10]<-paste("0",as.character(known),sep=" ") 
		#'0' because the unknown got the name from a CSPP present
		#in the same subDB
	 }
	}
	return(inp.x)
}

# write.table(inp.x,"compound2.csv",quote=F,sep="\t",row.names=F)