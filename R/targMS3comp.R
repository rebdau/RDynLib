targMS3comp<-function(dbkey1,dbkey2,subDB,AnalMS,IntThres=NULL){
	dbkey1<-which(subDB[[6]]==dbkey1)
	dbkey2<-which(subDB[[6]]==dbkey2)
	if (AnalMS=="FT") {
	 if (is.null(IntThres)) IntThres<-100
	} else {
	 if (is.null(IntThres)) IntThres<-5
	}
	pP<-subDB[[9]][[dbkey1]]	#p, precursor MS3 spectrum; P, parent ion 
	pMS3int<-as.numeric(unlist(strsplit(subDB[[8]] #int, intensity
					[[dbkey1]],split=",")))
	pMS3rint<-as.integer(pMS3int/max(pMS3int)*100)	#rint, relative intensity
	pMS3i<-as.integer(unlist(strsplit(subDB[[7]] # i = ions
				[[dbkey1]],split=",")))			
	pMS3n<-as.integer(round(pP)-pMS3i)	# n = neutrals
	cP<-subDB[[9]][[dbkey2]]	#c, compound MS2 spectrum
	cMS3int<-as.numeric(unlist(strsplit(subDB[[8]][[dbkey2]],split=",")))
	cMS3rint<-as.integer(cMS3int/max(cMS3int)*100)
	cMS3i<-as.integer(unlist(strsplit(subDB[[7]][[dbkey2]],split=",")))
	cMS3n<-as.integer(round(cP)-cMS3i)
	ac<-data.frame(cbind(pMS3i,pMS3rint,pMS3int))
	ac<-ac[!ac[,3]<IntThres,]
	pnr<-dim(ac)[1]	#nr, number of ions with intensity > IntThres
	bd<-data.frame(cbind(cMS3i,cMS3rint,cMS3int))
	bd<-bd[!bd[,3]<IntThres,]
	cnr<-dim(bd)[1]
	DotIons<-CommonDotProd(ac,bd)
	ac<-data.frame(cbind(pMS3n,pMS3rint,pMS3int))
	ac<-ac[!ac[,3]<IntThres,]
	bd<-data.frame(cbind(cMS3n,cMS3rint,cMS3int))
	bd<-bd[!bd[,3]<IntThres,]
	DotNeutrals<-CommonDotProd(ac,bd)
	dbkey1=subDB[[6]][[dbkey1]]
	dbkey2=subDB[[6]][[dbkey2]]
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

