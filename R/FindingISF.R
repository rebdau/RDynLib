# finding whether a m/z feature is an in-source fragment (ISF) of another m/z
# feature in a particular feature group. 
# fg and exp are the feature group number and the experiment (EXPID).

FindingISF<-function(fg,exp,ft,ms2.ft,XCMS){
	ftep<-ft[ft[,4]==exp,]
	fgcol=dim(XCMS)[2]-3
	XCMS.sel<-XCMS[XCMS[,fgcol]==fg,]
	XCMS.sel[,1]<-as.character(XCMS.sel[,1])
	isf<-as.character(rep(NA,length(XCMS.sel[,1])))
	xd<-which(ftep[,7]%in%XCMS.sel[,1])
	if(length(xd)>0){
	 print(xd)
	 parent<-round(ftep[xd,2],0)
	 COMPID<-ftep[xd,3]
	 prod_ion<-c()
	 for(i in 1:length(COMPID)){
	  prod_ion.int<-as.integer(unlist(strsplit(ms2.ft[[COMPID[i]]],split=",")))
	  par<-which(prod_ion.int%in%parent[i])
	  if(length(par)>0){
	   prod_ion.int<-prod_ion.int[-par]
	  }
	  prod_ion<-append(prod_ion,prod_ion.int)
	 }
	 xcms.mz<-as.integer(round(XCMS.sel[,5]))
	 xd2<-which(xcms.mz%in%prod_ion)
	 isf[xd2]<-'ISF'
	}
	XCMS.ISF<-data.frame(XCMS.sel,isf)
	return(XCMS.ISF)
}

