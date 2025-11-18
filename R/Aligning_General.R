Aligning_General<-function(expnr.ft,expnr.syn,err,t.ini,finlist,ESImode=NULL,chrg=NULL){
		#If chrg is not NULL, then the QTOFpos database is involved.
	if (is.null(chrg)) {
	 ft<-finlist[[1]][[1]]  #ftmsneg
	 if (is.null(ESImode)){
	  syn<-finlist[[3]][[1]] #qtofneg
	 } else {
	  syn<-finlist[[2]][[1]]  #ftmspos
	 }
	} else {
	 if (is.null(ESImode)){
	  ft<-finlist[[2]][[1]]  #ftmspos
	  syn<-finlist[[9]][[1]]  #qtofpos
	 } else {
	  ft<-finlist[[3]][[1]] #qtofneg
	  syn<-finlist[[9]][[1]] #qtofpos
	 }
	}
	ft.exp<-ft[ft[,4]==expnr.ft,]
	syn.exp<-syn[syn[,4]==expnr.syn,]
	ft.exp.o<-ft.exp[order(ft.exp[,1]),]
	COMPID<-ft.exp.o[dim(ft.exp.o)[1]/2,3]
	reg<-dim(ft.exp.o)[1]
	if (reg%%2!=0) reg=reg-1
	ft.sh<-Make_loc_shortlist(COMPID,reg,ft.exp)
	ft.sh <- Make_loc_shortlist(COMPID, reg, ft.exp)
	
		#ESImode=NULL means that both chromatograms are obtained in negative
		#mode. Thus, if the ESImode is given as e.g. pos, the m/z values of
		#the corresponding chromatogram should be corrected to be able to
		#match with those of the key chromatogram.
	# ESImode adjustment
	if (!is.null(ESImode)) {
	  #m/z is it the mass_measured of the feature?
	  syn.exp[,"m/z"] <- syn.exp[,"m/z"] - 2.01456
	}
	
	ft.sh <- Remov_empties(err, syn.exp, ft.sh, reg)
	
	LCal<-LCalign(err,t.ini,syn.exp,ft.sh)
	
	
	if (!is.null(ESImode)) {
	 LCal[,5]=as.character(as.numeric(LCal[,5])+2.01456)
	}
	return(LCal)
}