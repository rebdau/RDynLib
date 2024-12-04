FillAssocFTnQTOFn<-function(finlist,Assoc,expnr.ft,expnr.syn,cutoff,rg,lc.err,
					err,minIon=NULL,chrg=NULL){
	if (is.null(minIon)) minIon=0.6
	if (is.null(chrg)) {
	 ft<-finlist[[1]][[1]]
	 syn<-finlist[[3]][[1]]
	} else {
	 ft<-finlist[[2]][[1]]
	 syn<-finlist[[9]][[1]] 
	}
	ft.exp<-ft[ft[,4]==expnr.ft,]
	print(paste("FT COMPIDs from",min(ft.exp[,3]),"to",max(ft.exp[,3]),
			sep=" "))
	syn.exp<-syn[syn[,4]==expnr.syn,]
	print(paste("QTOF COMPIDs from",min(syn.exp[,3]),"to",
			max(syn.exp[,3]),sep=" "))
	if (is.null(chrg)) {
	 ms2.ft<-finlist[[1]][[3]]
	 ms2.syn<-finlist[[3]][[3]]
	} else {
	 ms2.ft<-finlist[[2]][[3]]
	 ms2.syn<-finlist[[9]][[3]]	
	}
	v=1
	while(v<=dim(ft.exp)[1]){
	 COMPID=ft.exp[v,3]
	 x.tR=ft.exp[v,1]
	 if(x.tR<cutoff){
	  v=v+1
	  next
	 }
	 if(rg[3]!=0){
	  t1.tR=ifelse(x.tR>rg[3],x.tR-rg[3],0)
	 }else{
	  t1.tR=0
	 }
	 if(rg[5]!=0){
	  t2.tR=ifelse(x.tR>rg[5],x.tR-rg[5],0)
	 }else{
	  t2.tR=0
	 }
	 y.tR=rg[1]+rg[2]*x.tR+rg[4]*t1.tR+rg[6]*t2.tR
	 y.l=y.tR-lc.err
	 y.h=y.tR+lc.err
	 pres<-Find_cand_matches(COMPID,err,syn.exp,ft.exp)
	 if(dim(pres)[1]==0){
	  v=v+1
	  next
	 }
	 ms2ion=as.integer(strsplit(ms2.ft[[COMPID]],",")[[1]])
	 pres1<-data.frame()
	 w.diff=c()
	 w=1
	 while(w<=dim(pres)[1]){
	  if((pres[w,1]>y.l)&(pres[w,1]<y.h)){
	   msmsion=unique(as.integer(strsplit(ms2.syn[[pres[w,3]]],",")[[1]]))
	   same_ion=which(ms2ion%in%msmsion)
	   if(length(same_ion)/length(ms2ion)>minIon){
	    pres1<-rbind(pres1,pres[w,])
	    w.diff<-append(w.diff,abs(y.tR-pres[w,1]))
	   }
	  } 
	  w=w+1
	 }
	 if(dim(pres1)[1]!=0){
	  pres1<-data.frame(pres1,w.diff)
	  pres1<-pres1[order(pres1[,5]),]
	  Assoc[COMPID,3]<-pres1[1,3]
	 }
	 v=v+1
	}
	return(Assoc)
}