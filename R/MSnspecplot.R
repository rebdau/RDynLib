# Yields MSn spectra of a selected COMPID and the MS3 spectra of COMPIDs from other experiments
# within the same sub-database that represent the same ion. In case MS1 is not NULL,
# an MS1 spectrum is displayed instead of the latter MS3 spectra. In addition, CSPPs and GNPS-like
# mass differences with other COMPIDs from the same experiment are given. 

MSnspecplot<-function(dbkey,base.dir,finlist,SubDB,nr_col=NULL,
				comp_add=NULL,lc.err=NULL,mz.err=NULL,MS1=NULL,
				nr_col2=NULL,gnps_add=NULL,prcx=NULL){
	# Settings:
	if (is.null(lc.err)) lc.err=0.02
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	if (is.null(mz.err)){
	 if (TMP[[2]]=="FT"){
	  mz.err=0.001
	 } else if (TMP[[2]]=="QTOF"){
	  mz.err=0.02
	 } else {
	  print("Unknown MS Analyzer, set SubDB equal to FTneg,
					 FTpos or QTOFneg")
	 }
	}
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
	# Checking for associated COMPIDs from other sub-databases
	if(SubDB=="FTneg"){
	 ass<-which(as.integer(finlist[[4]][,1])%in%as.integer(dbkey))
	 if(length(ass)!=0){
	  writeLines("\nNegative ion MSMS spectra:")
	  print(as.integer(finlist[[4]][ass,3]))
	  writeLines("\nPositive ion MSn spectra:")
	  print(as.integer(finlist[[4]][ass,2]))
	 }
	}else{
	 ass<-which(as.integer(finlist[[10]][,1])%in%as.integer(dbkey))
	 if(length(ass)!=0){
	  writeLines("\nPositive ion MSMS spectra:")
	  print(as.integer(finlist[[10]][ass,3]))
	  writeLines("\nNegative ion MSn spectra:")
	  print(as.integer(finlist[[10]][ass,2]))
	 }
	}
	# Generation of local CSPP and GNPS like mass differences
	cspp.res<-cspp.display(dbkey,base.dir,finlist,
					SubDB,nr_col,comp_add)
	writeLines("\nAssociated CSPPs:")
	print(cspp.res)
	gnps.res<-gnps.display(dbkey,base.dir,finlist,SubDB,
					nr_col2,gnps_add)
	writeLines("\nAssociated GNPS conversions:")
	print(gnps.res)
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	# MS1 and MSn data plotting for COMPID
	if (!is.null(MS1)) {
	 if(is.null(prcx))prcx=0.6
	 if(subDB[[10]][[dbkey_row]]!="NULL"){
	  ms3id.int<-as.integer(unlist(strsplit(subDB[[10]][[dbkey_row]],split=",")))
	  msnle<-length(ms3id.int)+2
	 }else{
	  msnle<-2
	 }
	 row.nr<-floor(sqrt(msnle))
	 col.nr<-ceiling(msnle/row.nr)
	 oldpar<-par(no.readonly=T)
	 par(mfrow=as.integer(c(row.nr,col.nr)))
	 single<-"single"
		# when single is NULL, the plot specifications of the nested function are used rather
		# than those demanded in this function.
	 MS1peaks<-MS1specplot(dbkey,base.dir,finlist,SubDB,TMP=TMP,lc.err=lc.err,mz.err=mz.err,
					single=single,prcx=prcx)
	 MSnplot(dbkey,base.dir,finlist,SubDB,TMP=TMP,single=single,prcx=prcx)
	 writeLines("\nMS1 peaks:")
	 print(MS1peaks)
	 return(list(MS1peaks,cspp.res,gnps.res))
	}else{
	 AllCompid<-WithinExpAligned(dbkey,TMP)
	 # MSn data of COMPID and MS3 spectra of other COMPIDs representing the same ion:
	 if(length(AllCompid)!=0){
	  if(is.null(prcx))prcx=0.6
	  print(AllCompid)
	  ans<-readline("Enter additional COMPIDs separated by a comma for which you want to see the MSn spectra:")	 
	  ans<-as.integer(unlist(strsplit(ans,split=",")))
	  ans<-append(dbkey,ans)
			# on the previous question, only COMPIDs that are mentioned in the line above should be
			# added; not the selected dbkey itself! The latter will be anyway attached
	  msnle=1
	  i=1
	  while(i<=length(ans)){
	   ans_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(ans[i]))
	   if(subDB[[10]][[ans_row]]!="NULL"){
	    ms3id.int<-as.integer(unlist(strsplit(subDB[[10]][[ans_row]],split=",")))
	    msnle<-msnle+length(ms3id.int)
	   }
	   i=i+1
	  }
	  row.nr<-floor(sqrt(msnle))
	  col.nr<-ceiling(msnle/row.nr)
	  oldpar<-par(no.readonly=T)
	  par(mfrow=as.integer(c(row.nr,col.nr)))
	  single<-"single"
	  MSnplot(dbkey,base.dir,finlist,SubDB,TMP=TMP,single=single,prcx=prcx)
	  i=2
	  while(i<=length(ans)){
	   ans_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(ans[i]))
	   if(subDB[[10]][[ans_row]]!="NULL"){
	    ms3id.int<-as.integer(unlist(strsplit(subDB[[10]][[ans_row]],split=",")))
	    ms3i=1
	    while(ms3i<=length(ms3id.int)){
	     MS3specplot(ans[i],base.dir,finlist,SubDB,TMP=TMP,single=single,prcx=prcx,wh=ms3i)
	     ms3i=ms3i+1
	    }
	   }
	   i=i+1
	  }
	  return(list(cspp.res,gnps.res))
	 # Just the MSn data of the COMPID:
	 }else{
	  MSnplot(dbkey,base.dir,finlist,SubDB,TMP=TMP)
	  return(list(cspp.res,gnps.res))
	 }
	}
}

# MSndata<-MSnspecplot(dbkey,base.dir,finlist,SubDB,MS1="pos")
# MSndata<-MSnspecplot(dbkey,base.dir,finlist,SubDB)