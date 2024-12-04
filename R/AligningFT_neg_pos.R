AligningFT_neg_pos<-function(expnr.ft,finlist,err=NULL,t.ini=NULL,lc.err=NULL,
				rng=NULL,minIon=NULL,Assoc=NULL){
	if (is.null(Assoc)) {			
	 max_expnr.ft<-max(as.integer(finlist[[1]][[1]][,4]))
	 if (expnr.ft>max_expnr.ft){
	  return("The entered expnr.ft is not in the FTneg DynLib database.")
	 }
	 if (is.null(err)) err=0.01 #Based on FT m/z accuracies
	} else {
	 max_expnr.ft<-max(as.integer(finlist[[1]][[1]][,4]))
	 if (expnr.ft>max_expnr.ft){
	  return("The entered expnr.ft is not in the QTOFneg DynLib database.")
	 }
	 if (is.null(err)) err=0.02 #Based on QTOF m/z accuracies
	}
	if (is.null(minIon)) minIon=0.2
	if (is.null(t.ini)) t.ini=5
		#Neg and pos should be run close in time on same instrument.
		#peak width is 12 sec or 0.2 min, lc.err is half a peak width
		#larger by default:
	if (is.null(lc.err)) lc.err=0.3
		#rng is the number of standard deviations that the retention time
		#difference can deviate from the mean retention time difference:
	if (is.null(rng)) rng=2
	if (is.null(Assoc)) {
	 chrg=NULL
	 Assoc<-finlist[[4]]
	 Assoc_EXPID_FTn<-strsplit(finlist[[5]][,3],",")
	 for (i in 1:length(Assoc_EXPID_FTn)) {
	  syn.sel<-which(Assoc_EXPID_FTn[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 	#i reflects the row obtained by the former for loop:
	 syn.poss<-finlist[[5]][i,4] 
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		FTpos experiment.")
	 }
	} else {
	 chrg="pos"
	 Assoc<-finlist[[11]]
	 Assoc_EXPID_FTn<-strsplit(finlist[[5]][,2],",")
	 for (i in 1:length(Assoc_EXPID_FTn)) {
	  syn.sel<-which(Assoc_EXPID_FTn[[i]]%in%as.character(expnr.ft))
	  if (length(syn.sel)!=0) break
	 }
	 	#i reflects the row obtained by the former for loop:
	 syn.poss<-finlist[[5]][i,1] 
	 if (syn.poss=="NULL") {
	  return("The entered expnr.ft does not correspond with a 
	 		FTpos experiment.")
	 }	
	}
	expnr.syn<-as.integer(readline(paste("Which of these exp, i.e.",syn.poss,
			", is the FTpos experiment:",sep=" ")))
	LCal<-Aligning_General(expnr.ft,expnr.syn,err,t.ini,finlist,ESImode="pos",chrg=chrg)
	LCal<-matchNegPos(LCal,finlist,minIon=minIon,chrg=chrg)
	LCal<-RemoveOutliers(LCal,rng)
	rg<-RegressionLin_LCalign(LCal)
	PlotLin_LCalign(LCal,rg)
	contQues<-as.integer(readline("Press 1 if you want to connect
			FTneg and FTpos peaks or 0 if you want to exit:"))
	if(contQues==1){
	 Assoc<-FillAssocFTnFTp(finlist,Assoc,rg,err,lc.err,expnr.ft,expnr.ftps=expnr.syn,chrg=chrg)
	}
	return(Assoc)

}