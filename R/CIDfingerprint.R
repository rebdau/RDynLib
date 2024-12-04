# with MS2series(), check first the quality of the CID spectra for which you
# want to extract the common mass spectral fingerprint. Only common ions are
# returned, not the intensities. Product ions are searched as integers.

CIDfingerprint<-function(nrCOMPID,base.dir,finlist,SubDB){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	nrdec1<-0
#	nrdec1<-ifelse(TMP[[2]]=="FT",0,1)	
#	nrdec2<-ifelse(TMP[[2]]=="FT",3,2)
	oldpar<-par(no.readonly=T)
	nrCOMPID_row<-which(subDB[[1]][,3]%in%as.integer(nrCOMPID))
	common_ion<-round(as.numeric(unlist(strsplit(
				subDB[[3]][[nrCOMPID_row[1]]],split=","))),nrdec1)
	i=1
	for(i in 1:length(nrCOMPID_row)){
	 dbkey<-nrCOMPID_row[i]
	 prod_ion.int<-round(as.numeric(unlist(strsplit(
				subDB[[3]][[dbkey]],split=","))),nrdec1)
	 common_ion<-intersect(common_ion,prod_ion.int)
	}
	return(common_ion)
}