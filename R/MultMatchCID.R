MultMatchCID<-function(base.dir,finlist,SubDB,dbkey,minIons=NULL){
	if (is.null(minIons)) minIons=4
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	AnalMS<-TMP[[2]]
	MS2list<-subDB[[3]]	
		# MS2list is a list of all subDB COMPIDs, with each list entry 
		# containing the CID product ions
	unk<-as.integer(unlist(strsplit(unlist(MS2list[dbkey_row]),split=",")))
		# unk comprises the CID product ions of the COMPID of interest
	if(length(unk)<1) return("Not enough product ions")
	if(!(minIons<length(unk))) return("minIons should be reduced")
	msnSect<-MultMatch(MS2list,unk,minIons)
		# msnSect contains very small subsets of the sub-database
		# (each containing e.g. 4 to 5 COMPIDs) for which the number of
		# product ions in common with those in the CID spectrum of the
		# COMPID of interest is above a preset threshold. These subsets
		# should be further refined to contain only one CID spectrum per
		# subset. This is done by the subMultMatch() function.
	ENTRIES<-subMultMatch(msnSect,unk,MS2list,minIons)
	fin<-data.frame(as.integer(),as.numeric(),as.integer(),as.integer(),
				as.numeric(),as.integer(),as.numeric())
	compname.cha<-as.character()
	rettime.num<-as.numeric()
	expid.int<-as.integer()
	i=1
		# for all retained COMPIDs from the sub-database that each have a
		# minimum number of product ions in common with those in the CID
		# spectrum of the COMPID of interest, the dot product is computed
		# and added to a final table. 
	while (i<=length(ENTRIES)) {
	 ENTRIES_compid<-subDB[[1]][ENTRIES[i],3]
	 out<-targMS2comp(dbkey1=dbkey,dbkey2=ENTRIES_compid,subDB,AnalMS)
	 compname.cha<-append(compname.cha,subDB[[1]][ENTRIES[i],5])
	 rettime.num<-append(rettime.num,subDB[[1]][ENTRIES[i],1])
	 expid.int<-append(expid.int,subDB[[1]][ENTRIES[i],4])
	 fin<-rbind(fin,out[1,c(4:10)])
	 i=i+1
	}
	fin<-data.frame(fin,compname.cha,rettime.num,expid.int)
	colnames(fin)<-c("COMPID","m/z","#FragIon","ComIons","DotIons",
			"ComLoss","DotLoss","COMPNAME","tR","EXPID")
	fin<-fin[order(fin$ComIons*fin$DotIons/fin[,3],decreasing=T),]
	writeLines(paste("The CID spectrum of compound",dbkey,
				"with m/z",subDB[[1]][dbkey_row,2],
				"matches to:",split=""))
	return(fin)
}

# MultMatchCID(base.dir,finlist,SubDB="FTneg",dbkey=1,minIons=3)