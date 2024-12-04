# not perfect: XCMS4 is ordered by the peak groups, but not fully by retention time!!!
# This is solved by first performing ConvertXCMS4.R!


SearchCombinat<-function(XCMS4,adjPeak,mode){
	if (mode=="neg") {
	 prtn=1.0073
	}else{
	 prtn=-1.0073
	}
	PeakGrp_col<-dim(XCMS4)[2]
	EndGrp<-1+dim(XCMS4)[1]
	StartGrp<-1+adjPeak
	massMS1<-as.integer(round(XCMS4[,5]+prtn,2)*100)
	Comb<-as.character(rep("NULL",dim(XCMS4)[1]))
	XCMS4<-data.frame(XCMS4,massMS1,Comb)
	XCMS4[,dim(XCMS4)[2]]<-as.character(XCMS4[,dim(XCMS4)[2]])
	i=StartGrp
	while(i<=EndGrp){
	 subXCMS<-XCMS4[XCMS4[,PeakGrp_col]%in%c((i-adjPeak):(i+adjPeak)),]
	 PkGrpOfInt<-which(subXCMS$CorrPkGrp%in%i)
	 ion_neut<-outer(subXCMS$massMS1,subXCMS$massMS1,"+")
	 ion_neut_ind<-which(ion_neut%in%subXCMS$massMS1)
	 ion_neut_ind_row<-as.integer(rep(NA,length(ion_neut_ind)))
	 ion_neut_ind_column<-as.integer(rep(NA,length(ion_neut_ind)))
	 INcomp<-as.integer(rep(NA,length(ion_neut_ind)))
	 for (j in 1:length(ion_neut_ind)){
	  ion_neut_ind_row[j]<-ceiling(ion_neut_ind[j]/length(subXCMS$massMS1))
	  ion_neut_ind_column[j]<-ifelse(ion_neut_ind[j]%%length(subXCMS$massMS1)==0,
				dim(subXCMS)[1],ion_neut_ind[j]%%length(subXCMS$massMS1))
	  INcomp[j]<-which(subXCMS$massMS1==ion_neut[ion_neut_ind_row[j],
	 			ion_neut_ind_column[j]])
	 }
	 Combin<-as.character(rep(NA,length(subXCMS[,5])))
	 Combin[INcomp]<-paste(round(subXCMS[,5][ion_neut_ind_row],2),"/",
				round(subXCMS[,5][ion_neut_ind_column],2),split="")
	 CombSel1<-which(!is.na(Combin))
	 CombSel<-intersect(CombSel1,PkGrpOfInt)
	 nameSel<-subXCMS[CombSel,1]
	 XCMS4[which(XCMS4[,1]%in%nameSel),dim(XCMS4)[2]]<-Combin[CombSel]
	 i=i+1
	}
	return(XCMS4)
}