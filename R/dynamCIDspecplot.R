dynamCIDspecplot<-function(dbkey,base.dir,finlist,SubDB,err=NULL,t.ini=NULL,
					lc.err=NULL,rng=NULL,startpoint=NULL,cutoff=NULL,
					minIon=NULL,prcx=NULL) {
	if (is.null(prcx)) prcx=0.6
	if (is.null(err)) err=0.02 
	if (is.null(t.ini)) t.ini=5
	if (is.null(lc.err)) lc.err=1
	if (is.null(rng)) rng=2
	if (is.null(startpoint)) startpoint=1
	if (is.null(cutoff)) cutoff=1
	if (is.null(minIon)) minIon=0.6
	Getinf<-GetAllInfo(base.dir,finlist,SubDB,dbkey)
	expnr.ft<-as.integer(Getinf[[1]])
	expnr.syn<-as.integer(Getinf[[2]])
	ft.exp<-Getinf[[3]]
	syn.exp<-Getinf[[4]]
	chrg<-Getinf[[5]]
	LCalNM<-paste(SubDB,expnr.ft,"-",expnr.syn,".txt",sep="")
	setwd(paste(base.dir,"/DynLib subDB alignment/LCal",sep=""))
	if (LCalNM%in%list.files()) {
	 LCal<-read.table(LCalNM,header=T,sep="\t",stringsAsFactors=F)
	}else{
		# Below, the order of the LCal columns is always the same whether
		# you entered a dbkey from a QTOF subDB or from a FT subDB. First
		# column is the COMPID (FT), then the tR (FT), the m/z (FT), the
		# tR (QTOF), the m/z (QTOF), the tR difference and the COMPID (QTOF)
	 LCal<-Aligning_General(expnr.ft,expnr.syn,err,t.ini,finlist,chrg=chrg)
	 write.table(LCal,LCalNM,sep="\t",row.names=F)
	}
	LCal<-matchFTSyn(LCal,finlist,minIon=minIon,chrg=chrg)
	LCal<-RemoveOutliers(LCal,rng)
	rg<-RegressionPie_LCalign(LCal,startpoint)
	pres<-Find_isomers(dbkey,SubDB,syn.exp,ft.exp,err=err)
	pres1<-pres[[1]][order(pres[[1]][,1]),]
	print(pres1)
	pres2<-pres[[2]][order(pres[[2]][,1]),]
	print(pres2)
	if( (dim(pres1)[1]==0)|(dim(pres2)[1]==0) ) return("No aligned COMPIDs between FT and QTOF subDB")
	sdb<-SubDBcombin(finlist,SubDB)
	TMP1<-SelectSubDB(base.dir,finlist,SubDB=sdb[[1]])
	TMP2<-SelectSubDB(base.dir,finlist,SubDB=sdb[[2]])
	TMP3<-SelectSubDB(base.dir,finlist,SubDB=sdb[[3]])
	TMP4<-SelectSubDB(base.dir,finlist,SubDB=sdb[[4]])
	z=1
	v=1
	repeat{
	 COMPID_FT<-pres1[z,3]
	 COMPID_QTOF<-pres2[v,3]
	 print(paste("Checking FT COMPID ",COMPID_FT," and QTOF COMPID ",COMPID_QTOF,sep=""))
	 dynamCIDplot(dbkey1=COMPID_FT,dbkey2=COMPID_QTOF,base.dir,finlist,LCal,rg,pres1,pres2,
				z,v,TMP1,TMP2,TMP3,TMP4,SubDB,Assoc=sdb[[5]],Assoc2=sdb[[6]],prcx=prcx)
	 print("Enter 21 if you want to move to another FT isomer")
	 print("Enter 12 if you want to move to another QTOF isomer")
#	 print("Enter 3 if you want to make an association between the current FT and QTOF isomers")
	 print("Enter 0 if you want to quit the function")
	 verif<-as.integer(readline("What do you want to do:"))
	 if (verif==12) {
	  v=v+1
	 } else if (verif==21) {
	  z=z+1
#	 } else if (verif==3) {
#	  if((SubDB=="FTneg")|(SubDB=="QTOFneg")){
#	   print(finlist[[4]][which(finlist[[4]][,1]%in%COMPID_FT),])
#	   chg_ass<-readline("Enter new association?(y/n):")
#	   if (chg_ass=="y") finlist[[4]][which(finlist[[4]][,1]%in%COMPID_FT),3]<-COMPID_QTOF
#	  }else if ((SubDB=="FTpos")|(SubDB=="QTOFpos")){
#	   print(finlist[[10]][which(finlist[[10]][,1]%in%COMPID_FT),])
#	   chg_ass<-readline("Enter new association?(y/n):")
#	   if (chg_ass=="y") finlist[[10]][which(finlist[[10]][,1]%in%COMPID_FT),3]<-COMPID_QTOF
#	  }else{
#	   print("SubDB is wrongly defined, cannot write association.")
#	  }
	 } else if (verif==0){
	  break
	 }
	 dev.off()
	 if(z>dim(pres1)[1])z=1
	 if(v>dim(pres2)[1])v=1
	}
	dev.off()
}

# dbkey=9941
# SubDB="FTneg"
# dynamCIDspecplot(dbkey,base.dir,finlist,SubDB,err=NULL,t.ini=NULL,lc.err=NULL,rng=NULL,
#                        startpoint=NULL,cutoff=NULL,minIon=NULL,prcx=NULL)