Find_cand_matches<-function(COMPID,err,syn.exp,ft.exp){
	o<-order(syn.exp[,2])
	syn.o<-syn.exp[o,]
	i<-which(ft.exp[,3]==COMPID)
	lb=ft.exp[i,2]-err
	ub=ft.exp[i,2]+err
	pres<-array(dim=c(0,4))
	j=1
	while (j<=dim(syn.o)[1]) {
	 if (syn.o[j,2]>lb&syn.o[j,2]<ub){
	   int<-syn.o[j,1:4]
	   pres<-rbind(pres,int)
	 }
	 if (syn.o[j,2]>=ub) break
	 j=j+1
	}
	return(pres)
}

