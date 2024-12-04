MS1plot<-function(dbkey,exp,subDB,proc,fc,finlist,lc.err,mz.err,prcx){
	dynlb<-subDB[[1]][subDB[[1]][,4]==exp,]
	nodename<-as.character(dynlb[dynlb[,3]==dbkey,7])
#	rettime<-dynlb[dynlb[,3]==dbkey,1]
#	lowRt<-rettime-lc.err
#	highRt<-rettime+lc.err
	proc<-proc[order(proc[,1]),]
	sel<-which(proc[,3]%in%nodename)
	rettime<-proc[sel,1]
	lowRt<-rettime-lc.err
	highRt<-rettime+lc.err
	mzMS1<-as.numeric()
	tRMS1<-as.numeric()
	namMS1<-as.character()
	abunMS1<-as.numeric()
	i=sel
	while(proc[i,1]>lowRt){
	 mz.MS1<-proc[i,2]
	 tR.MS1<-proc[i,1]
	 nam.MS1<-proc[i,3]
	 abun.MS1<-mean(as.numeric(t(proc[i,5:dim(proc)[2]])),na.rm=T)
 	 mzMS1<-append(mzMS1,mz.MS1)
	 tRMS1<-append(tRMS1,tR.MS1)
	 namMS1<-append(namMS1,nam.MS1)
	 abunMS1<-append(abunMS1,abun.MS1)
	 i=i-1
	}
	i=sel+1
	while(proc[i,1]<highRt){
	 mz.MS1<-proc[i,2]
	 tR.MS1<-proc[i,1]
	 nam.MS1<-proc[i,3]
	 abun.MS1<-mean(as.numeric(t(proc[i,5:dim(proc)[2]])),na.rm=T)
 	 mzMS1<-append(mzMS1,mz.MS1)
	 tRMS1<-append(tRMS1,tR.MS1)
	 namMS1<-append(namMS1,nam.MS1)
	 abunMS1<-append(abunMS1,abun.MS1)
	 i=i+1
	}
	relabun<-round(100*abunMS1/abunMS1[1],digits=1)
	mzdiff<-round(mzMS1-mzMS1[1],digits=4)	
	compID<-as.integer(rep(NA,length(mzMS1)))
	i=1
	while(i<=length(mzMS1)){
	 sl<-which(dynlb[,7]%in%namMS1[i])
	 if(length(sl)!=0){
	  compID[i]<-dynlb[sl,3]
	 }
	 i=i+1
	}
	FeatName<-as.character(rep(NA,length(mzMS1)))
	i=1
	while(i<=length(mzMS1)){
	 lowthr<-abs(mzdiff[i])-mz.err
	 highthr<-abs(mzdiff[i])+mz.err
	 sl1<-which(finlist[[6]]>lowthr)
	 sl2<-which(finlist[[6]]<highthr)
	 sl<-intersect(sl1,sl2)
	 if (length(sl)==1) {
	  FeatName[i]<-finlist[[6]][sl,2]
	 }
	 i=i+1
	}
	massMS1<-as.integer(round(mzMS1+1.0073,2)*100)
	ion_neut<-outer(massMS1,massMS1,"+")
	ion_neut_ind<-which(ion_neut%in%massMS1)
	ion_neut_ind_row<-as.integer(rep(NA,length(ion_neut_ind)))
	ion_neut_ind_column<-as.integer(rep(NA,length(ion_neut_ind)))
	INcomp<-as.integer(rep(NA,length(ion_neut_ind)))
	for (i in 1:length(ion_neut_ind)){
	 ion_neut_ind_row[i]<-ceiling(ion_neut_ind[i]/length(massMS1))
	 ion_neut_ind_column[i]<-ifelse(ion_neut_ind[i]%%length(massMS1)==0,
					length(massMS1),ion_neut_ind[i]%%length(massMS1))
	 INcomp[i]<-which(massMS1==ion_neut[ion_neut_ind_row[i],
				ion_neut_ind_column[i]])
	}
	Combin<-as.character(rep(NA,length(mzMS1)))
	Combin[INcomp]<-paste(round(mzMS1[ion_neut_ind_row],2),"/",
				round(mzMS1[ion_neut_ind_column],2),split="")
	if(all(is.na(fc))){
	 dMS1<-data.frame(mzMS1,mzdiff,relabun,FeatName,compID,Combin)
	 colnames(dMS1)<-c("m/z","m/z diff","Rel.Int.","Type","COMPID","Combin")
	} else {
	 corrcoef<-MS1corr(fc,proc,namMS1)
	 dMS1<-data.frame(mzMS1,mzdiff,relabun,FeatName,compID,Combin,corrcoef)
	 colnames(dMS1)<-c("m/z","m/z diff","Rel.Int.","Type","COMPID","Combin",
					"R2")
	}
	adj_intens<-max(relabun)*1.1
	adj_prod_min<-min(mzMS1)*0.9
	adj_prod_max<-max(mzMS1)*1.1
	par(cex=prcx)
	plot(mzMS1,relabun,type="h",xlab="m/z",ylab="relative intensity",
	 main=paste("MS1: m/z",mzMS1[1],sep=" "),xlim=c(adj_prod_min,adj_prod_max),
			ylim=c(0,adj_intens))
	text(mzMS1,relabun,round(mzMS1,digits=3),pos=3)
	text(mzMS1,relabun,round(mzdiff,digits=3),pos=3,offset=1.5,col=2)
	return(dMS1)
}

