Get_FT_DynLib_files<-function(ft.dir){
	setwd(ft.dir) 
 	comp<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)

#extract just the tR, m/z, COMPID, EXPID, COMPNAME and SMILES columns

	ft<-data.frame(round(as.numeric(comp[,8]),2),
				round(as.numeric(comp[,5]),5),
				as.integer(comp[,1]),
				as.integer(comp[,9]),
				as.character(comp[,3]),
				as.character(comp[,13]),
				as.character(comp[,2]),
				stringsAsFactors=F)
	colnames(ft)<-c("tR","m/z","COMPID","EXPID","COMPNAME","SMILES","NODENAME")			
	rm(comp)
 	inp<-scan("MS2spectra.csv",what=list("","","","","",""),sep="\t")

#again COMPID, but this is associated with MS2spectrum file so that
#"ft" can be sorted if necessary

 	compid.ms2<-as.integer(inp[[1]][-1])
 	ms2peaklist<-as.list(inp[[4]][-1])
 	ms2intensitylist<-as.list(inp[[5]][-1])
	ms3idlist<-as.list(inp[[6]][-1])
 	rm(inp)
 	inp<-scan("MS3spectra.csv",what=list("","","","","",""),sep="\t")

#again COMPID, but this is associated with MS3spectrum file.

	compid.ms3<-as.integer(inp[[1]][-1])
 	ms3id<-as.integer(inp[[2]][-1])
 	ms3peaklist<-as.list(inp[[4]][-1])
 	ms3intensitylist<-as.list(inp[[5]][-1])

#MS2ions from which MS3 spectra are derived

 	parent_ionms3.num<-as.numeric(inp[[3]][-1])
 	rm(inp)
 	fin.list<-list(ft,compid.ms2,ms2peaklist,ms2intensitylist,
	 compid.ms3,ms3id,ms3peaklist,ms3intensitylist,parent_ionms3.num,ms3idlist)
 	return(fin.list)
}
