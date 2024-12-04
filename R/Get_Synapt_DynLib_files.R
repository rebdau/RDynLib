Get_Synapt_DynLib_files<-function(syn.dir){
 	setwd(syn.dir)
 	comp<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)

#tR, m/z, COMPID, EXPID, COMPNAME and SMILES

 	syn<-data.frame(round(as.numeric(comp[,8]),2),
				round(as.numeric(comp[,5]),3),
				as.integer(comp[,1]),
				as.integer(comp[,9]),
				as.character(comp[,3]),
				as.character(comp[,13]),
				as.character(comp[,2]),
				stringsAsFactors=F)
	colnames(syn)<-c("tR","m/z","COMPID","EXPID","COMPNAME","SMILES","NODENAME")
	rm(comp)
 	inp<-scan("MS2spectra.csv",what=list("","","","",""),sep="\t")

#again COMPID, but this is associated with MS2spectrum file 
#so that "syn" can be sorted if necessary

 	compid.ms2<-as.integer(inp[[1]][-1])
 	ms2peaklist<-as.list(inp[[4]][-1])
 	ms2intensitylist<-as.list(inp[[5]][-1])
 	rm(inp)
 	fin.list<-list(syn,compid.ms2,ms2peaklist,ms2intensitylist)
 	return(fin.list)
}

