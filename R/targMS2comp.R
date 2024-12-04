targMS2comp<-function(dbkey1,dbkey2,subDB,AnalMS,IntThres=NULL){
	# in the initial script, dbkey1 and dbkey2 referred to the row in the subdatabase (subdb)
	# as this rownumber was the same as the COMPID. However, to allow working with subdatabases
	# that are shrunk to a few experiments, dbkey1 and dbkey2 should really refer to the rows
	# corresponding with the intended COMPIDs:
	dbkey1<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey1))
	dbkey2<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey2))
	if (AnalMS=="FT") {
	 if (is.null(IntThres)) IntThres<-100
	} else {
	 if (is.null(IntThres)) IntThres<-5
	}
	pP<-subDB[[1]][dbkey1,2]	#p, precursor MS2 spectrum; P, parent ion 
	pMS2int<-as.numeric(unlist(strsplit(subDB[[4]] #int, intensity
					[[dbkey1]],split=",")))
	pMS2rint<-as.integer(pMS2int/max(pMS2int)*100)	#rint, relative intensity
	pMS2i<-as.integer(unlist(strsplit(subDB[[3]] # i = ions
				[[dbkey1]],split=",")))			
	pMS2n<-as.integer(round(pP)-pMS2i)	# n = neutrals
	cP<-subDB[[1]][dbkey2,2]	#c, compound MS2 spectrum
	cMS2int<-as.numeric(unlist(strsplit(subDB[[4]][[dbkey2]],split=",")))
	cMS2rint<-as.integer(cMS2int/max(cMS2int)*100)
	cMS2i<-as.integer(unlist(strsplit(subDB[[3]][[dbkey2]],split=",")))
	cMS2n<-as.integer(round(cP)-cMS2i)
	ac<-data.frame(cbind(pMS2i,pMS2rint,pMS2int))
	ac<-ac[!ac[,3]<IntThres,]
	pnr<-dim(ac)[1]	#nr, number of ions with intensity > IntThres
	bd<-data.frame(cbind(cMS2i,cMS2rint,cMS2int))
	bd<-bd[!bd[,3]<IntThres,]
	cnr<-dim(bd)[1]
	DotIons<-CommonDotProd(ac,bd)
	ac<-data.frame(cbind(pMS2n,pMS2rint,pMS2int))
	ac<-ac[!ac[,3]<IntThres,]
	bd<-data.frame(cbind(cMS2n,cMS2rint,cMS2int))
	bd<-bd[!bd[,3]<IntThres,]
	DotNeutrals<-CommonDotProd(ac,bd)
	dbkey1=as.integer(subDB[[1]][dbkey1,3])
	dbkey2=as.integer(subDB[[1]][dbkey2,3])
	fin<-data.frame(as.integer(dbkey1),pP,as.integer(pnr),as.integer(dbkey2),
		cP,as.integer(cnr),as.integer(DotIons[[1]]),round(DotIons[[2]],2),
		as.integer(DotNeutrals[[1]]),round(DotNeutrals[[2]],2),
		round(DotIons[[1]]/pnr,2),round(DotIons[[1]]/cnr,2),
		round(DotNeutrals[[1]]/pnr,2),round(DotNeutrals[[1]]/cnr,2))
	colnames(fin)<-c("prec","prec_parent","prec_nr_ions","cmp","cmp_parent",
			"cmp_nr_ions","common_ions","dot_ions","common_losses",
			"dot_losses","forward_ions","reverse_ions","forward_losses",
			"reverse_losses")
	return(fin)
}

