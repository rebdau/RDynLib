Aligning_Exp<-function(expnr1,expnr2,base.dir,finlist,SubDB,err=NULL,t.ini=NULL,
				rng=NULL,thr1=NULL,thr2=NULL,IntThres=NULL,startpoint=NULL){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB=TMP[[1]]
	AnalMS=TMP[[2]]
		# error of m/z value
	if (is.null(err)) err<-ifelse(AnalMS=="FT",0.01,0.02)
		# how many neighboring peaks around the considered peak should be taken
	if (is.null(t.ini)) t.ini=5
		# number of SD that retention time difference can be outside the mean
	if (is.null(rng)) rng=2
		# retention time from which to start the regression
	if (is.null(startpoint)) startpoint=1
	if (is.null(thr1)) thr1=1
	if (is.null(thr2)) thr2=0.9
	if (AnalMS=="FT") {
	 if (is.null(IntThres)) IntThres<-100
	} else {
	 if (is.null(IntThres)) IntThres<-5
	}	
	ft<-subDB[[1]]
	int.ft<-subDB[[4]]
	ms2.ft<-subDB[[3]]
	ft.exp<-ft[ft[,4]==expnr1,]
	ftps.exp<-ft[ft[,4]==expnr2,]
	ft.exp.o<-ft.exp[order(ft.exp[,1]),]
	COMPID<-ft.exp.o[dim(ft.exp.o)[1]/2,3]
	reg<-dim(ft.exp.o)[1]
	if (reg%%2!=0) reg=reg-1
	ft.sh<-Make_loc_shortlist(COMPID,reg,ft.exp)
	ft.sh<-Remov_empties(err,ftps.exp,ft.sh,reg)
	LCal<-LCalign(err,t.ini,ftps.exp,ft.sh)
	ftps<-ft
	int.ftps<-int.ft
	ms2.ftps<-ms2.ft
	LCal<-match_FTexp(LCal,subDB,AnalMS,thr1=thr1,thr2=thr2,IntThres=IntThres)
		# only pairs are retained in which all product ions of one MS2 spectrum
		# are present in the other MS2 spectrum (thr1=1) with a dot product > 0.9
		# (thr2). Furthermore, both MS2 spectra should contain at least 2 ions. 
	LCal.lst<-plotAlignment(LCal,rng,startpoint)
		# Up to this stage, a set of peaks in the reference chromatogram has been
		# paired to the same peaks in the second chromatogram and a piecewise
		# retention time model has been built. 
	LCal.lst<-LimKeys(LCal.lst)
	Update_ExpAlignment(TMP,LCal.lst,expnr1,expnr2)
	ChkNames(LCal.lst,TMP)
}

# Aligning_Exp(expnr1=24,expnr2=25,base.dir,finlist,SubDB="FTneg")