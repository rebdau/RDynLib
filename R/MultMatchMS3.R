MultMatchMS3<-function(base.dir,finlist,SubDB,dbkey,minIons=NULL){
	if (is.null(minIons)) minIons=4
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	AnalMS<-TMP[[2]]
	MS3list<-subDB[[7]]
	dbms.ms3<-as.numeric()
	if (subDB[[10]][[dbkey_row]]=="NULL") {
	 return("No MS3 spectra available for the selected COMPID.")
	} else {
	 dbms<-which(subDB[[5]]==dbkey)
	 ms3id.int<-subDB[[6]][dbms]
	 dbms.ms3<-subDB[[9]][dbms]
#	 dbms<-as.integer(unlist(strsplit(unlist(subDB[[10]][[dbkey]]),
#					split=",")))
#	 for(i in 1:length(dbms)){
#	  dbms.ms3[i]<-subDB[[9]][[dbms[i]]]
#	 }
	}
	fin.lst<-list()
	mess.lst<-list()
	i=1
	while (i<=length(dbms)) {
	 unk<-as.integer(unlist(strsplit(unlist(MS3list[dbms[i]]),split=",")))
	 if(length(unk)<1){
	  print(paste("Not enough second-order product ions for m/z",
				dbms.ms3[i],split=""))
	  i=i+1
	  next
	 }
	 if(!(minIons<length(unk))){
	  print(paste("minIons is bigger or equal to the number of second-order product ions for m/z",
				dbms.ms3[i],split=""))
	  i=i+1
	  next
	 }
	 msnSect<-MultMatch(MS3list,unk,minIons)
	 ENTRIES<-subMultMatch(msnSect,unk,MS3list,minIons)
	 fin<-data.frame(as.integer(),as.numeric(),as.integer(),as.integer(),
				as.numeric(),as.integer(),as.numeric())
	 compname.cha<-as.character()
	 rettime.num<-as.numeric()
	 compid.int<-as.integer()
	 expid.int<-as.integer()
	 j=1
	 while (j<=length(ENTRIES)) {
	  ENTRIES_ms3id<-subDB[[6]][ENTRIES[j]]
	  out<-targMS3comp(dbkey1=ms3id.int[i],dbkey2=ENTRIES_ms3id,subDB,AnalMS)
	  rannr<-which(subDB[[1]][,3]%in%(subDB[[5]][ENTRIES[j]]))
#	  minMS2=1
#	  maxMS2=length(subDB[[10]])
#	  while(minMS2<maxMS2){
#	   rannr<-sample(minMS2:maxMS2,1)
#	   if(subDB[[10]][[rannr]]=="NULL") {
#	    next
#	   }else{
#	    v<-as.integer(unlist(strsplit(unlist(subDB[[10]][[rannr]]),split=",")))
#	    min.v<-min(v)
#	    max.v<-max(v)
#	    if(length(which(v%in%ENTRIES[j]))!=0){
#	     break	     
#	    }else if(min.v>ENTRIES[j]){
#	     maxMS2=rannr
#	    }else if(max.v<ENTRIES[j]){
#	     minMS2=rannr
#	    }
#	    next
#	   }
#	  }
	  compname.cha<-append(compname.cha,subDB[[1]][rannr,5])
	  rettime.num<-append(rettime.num,subDB[[1]][rannr,1])
	  compid.int<-append(compid.int,subDB[[1]][rannr,3])
	  expid.int<-append(expid.int,subDB[[1]][rannr,4])
	  fin<-rbind(fin,out[1,c(4:10)])
	  j=j+1
	 }
	 fin<-data.frame(fin,compname.cha,rettime.num,compid.int,expid.int)
	 colnames(fin)<-c("COMPID_MS3","m/z","#FragIon","ComIons","DotIons",
			"ComLoss","DotLoss","COMPNAME","tR","COMPID","EXPID")
	 fin<-fin[order(fin$ComIons*fin$DotIons/fin[,3],decreasing=T),]
	 fin.lst[[i]]<-fin
	 mess.lst[[i]]<-paste("The product ion at m/z",dbms.ms3[i],
				"of the MS2 spectrum of compound",dbkey,
				"with m/z",subDB[[1]][dbkey_row,2],
				"has an MS3 spectrum that matches:",split="")
	 i=i+1
	}
	LST<-list()
	for(i in 1:length(fin.lst)){
	 j=2*i
	 LST[[j-1]]<-mess.lst[[i]]
	 LST[[j]]<-fin.lst[[i]]
	}
	return(LST)
}

# MultMatchMS3(base.dir,finlist,SubDB,dbkey)