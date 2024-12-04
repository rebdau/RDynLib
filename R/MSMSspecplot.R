MSMSspecplot<-function(dbkey,base.dir,finlist,SubDB,err=NULL,minum=NULL,
				nr_col=NULL,comp_add=NULL,nr_col2=NULL,gnps_add=NULL){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	if (is.null(err)) err=0.015
	if (is.null(minum)) minum=2
	if (is.null(nr_col)) nr_col<-35
	if (is.null(nr_col2)) nr_col2<-6
	if (is.null(comp_add)) {
	 setwd(TMP[[4]])
	 comp_add<-read.table("compound_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	}
	if (is.null(gnps_add)) {
	 setwd(TMP[[4]])
	 gnps_add<-read.table("gnps_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	}
	if(SubDB=="QTOFneg"){
	 ass<-which(as.integer(finlist[[4]][,3])%in%as.integer(dbkey))
	 if(length(ass)!=0){
	  writeLines("\nNegative ion MSn spectra:")
	  print(as.integer(finlist[[4]][ass,1]))
	  writeLines("\nPositive ion MSn spectra:")
	  print(as.integer(finlist[[4]][ass,2]))
	 }
	}else{
	 ass<-which(as.integer(finlist[[11]][,2])%in%as.integer(dbkey))
	 if(length(ass)!=0){
	  writeLines("\nNegative ion MS/MS spectra:")
	  print(as.integer(finlist[[11]][ass,1]))
	 }
	}
	cspp.res<-cspp.display(dbkey,base.dir,finlist,SubDB,nr_col=nr_col,
					comp_add=comp_add)
	writeLines("\nAssociated CSPPs:")
	print(cspp.res)
	gnps.res<-gnps.display(dbkey,base.dir,finlist,SubDB,nr_col2=nr_col2,
					gnps_add=gnps_add)
	writeLines("\nAssociated GNPS conversions:")
	print(gnps.res)
	oldpar<-par(no.readonly=T)
	MSMSplot(dbkey,base.dir,finlist,SubDB,err,minum,oldpar,nl="nl")
	par(oldpar)
	return(list(cspp.res,gnps.res))
}


# MSMSspecplot(25120,base.dir,finlist,SubDB="QTOFneg")