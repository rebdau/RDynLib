FillAssocFTnFTp<-function(finlist,Assoc,rg,err,lc.err,expnr.ft,expnr.ftps,chrg=NULL){
	if (is.null(chrg)) {
	 ft<-finlist[[1]][[1]]
	 ftps<-finlist[[2]][[1]]
	} else {
	 ft<-finlist[[3]][[1]]
	 ftps<-finlist[[9]][[1]]
	}
	ftng.exp<-ft[ft[,4]==expnr.ft,]
	ftps.exp<-ftps[ftps[,4]==expnr.ftps,]
	print(paste("FTneg COMPIDs from",min(ftng.exp[,3]),"to",max(ftng.exp[,3]),
			sep=" "))
	print(paste("FTpos COMPIDs from",min(ftps.exp[,3]),"to",max(ftps.exp[,3]),
			sep=" "))
	o<-order(ftng.exp[,2])
	ftng.o<-ftng.exp[o,]
	o<-order(ftps.exp[,2])
	ftps.o<-ftps.exp[o,]
	i=1
	repeat{
	 mass.ftng<-ftng.o[i,2]
	 time.ftng<-ftng.o[i,1]
	 COMPID.ftng<-ftng.o[i,3]
	 pres<-Find_pos_compound(mass.ftng,time.ftng,ftps.o,err,lc.err,rg)
	 if(dim(pres)[1]==0){
	  if(i==dim(ftng.o)[1])break
	  i=i+1
	  next
       }else{
	  Assoc<-Sort_comp_matches(pres,Assoc,time.ftng,COMPID.ftng,rg)
	 }
	 if(i==dim(ftng.o)[1])break
	 i=i+1
	}
	return(Assoc)
}