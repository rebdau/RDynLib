Aligning_FT_QTOF<-function(expnr.ft,finlist,err=NULL,t.ini=NULL,lc.err=NULL,
				rng=NULL,startpoint=NULL,cutoff=NULL,minIon=NULL,
				Assoc=NULL){
	if (is.null(Assoc)) {
	 max_expnr.ft<-max(as.integer(finlist[[1]][[1]][,4]))
	 if (expnr.ft>max_expnr.ft){
	  return("The entered expnr.ft is not in the FTneg DynLib database.")
	 }
	} else {
	 max_expnr.ft<-max(as.integer(finlist[[2]][[1]][,4]))
	 if (expnr.ft>max_expnr.ft){
	  return("The entered expnr.ft is not in the FTpos DynLib database.")
	 }	
	}
	if (is.null(minIon)) minIon=0.6
	if (is.null(err)) err=0.02 #Based on QTOF m/z accuracy
	if (is.null(t.ini)) t.ini=5
		#lc.err represents the error on the corrected QTOF retention time,
		#i.e. after applying the regression model.
	if (is.null(lc.err)) lc.err=1
		#rng is the number of standard deviations that the retention time
		#difference can deviate from the mean retention time difference:
	if (is.null(rng)) rng=2
	if (is.null(startpoint)) startpoint=1
		#Beyond which retention time should FT peaks be connected to QTOF
		#peaks in the DynLib database:
	if (is.null(cutoff)) cutoff=1
	if (is.null(Assoc)) {
	 chrg=NULL
	 Assoc<-finlist[[4]]
	 Assoc_EXPID_FTn<-strsplit(finlist[[5]][,3],",")
	 for (i in 1:length(Assoc_EXPID_FTn)) {
	  syn.sel<-which(Assoc_EXPID_FTn[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 	#i reflects the row obtained by the former for loop:
	 syn.poss<-finlist[[5]][i,2]
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		QTOFneg experiment.")
	 }
	} else {
	 chrg<-"pos"
	 Assoc<-finlist[[10]]
	 Assoc_EXPID_FTn<-strsplit(finlist[[5]][,4],",")
	 for (i in 1:length(Assoc_EXPID_FTn)) {
	  syn.sel<-which(Assoc_EXPID_FTn[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 	#i reflects the row obtained by the former for loop:
	 syn.poss<-finlist[[5]][i,1]
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		QTOFpos experiment.")
	 }	
	}
	expnr.syn<-as.integer(readline(paste("Which of these exp, i.e.",syn.poss,
			", is the QTOF experiment:",sep=" ")))
	LCal<-Aligning_General(expnr.ft,expnr.syn,err,t.ini,finlist,chrg=chrg)
	LCal<-matchFTSyn(LCal,finlist,minIon=minIon,chrg=chrg) 
	LCal<-RemoveOutliers(LCal,rng)
	rg<-RegressionPie_LCalign(LCal,startpoint)
	PlotPie_LCalign(LCal,rg)
	contQues<-as.integer(readline("Press 1 if you want to connect
			FT and Synapt peaks or 0 if you want to exit:"))
	if(contQues==1){
	 Assoc<-FillAssocFTnQTOFn(finlist,Assoc,expnr.ft,expnr.syn,cutoff,rg,lc.err,err,minIon=minIon,chrg=chrg)
	}
	return(Assoc)
}