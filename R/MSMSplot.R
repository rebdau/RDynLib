MSMSplot<-function(dbkey,base.dir,finlist,SubDB,err=NULL,minum=NULL,
				oldpar=NULL,nl=NULL){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	if (is.null(err)) err=0.015
	if (is.null(minum)) minum=2
	prdion<-finlist[[8]]
	neutloss<-finlist[[7]]
	prod_ion.num<-round(as.numeric(unlist(strsplit(subDB[[3]][[dbkey_row]],split=","))),2)
	writeLines("\nCandidate product ions:")
	ProdIonMatch(prod_ion.num,prdion,err)
	intens_ion.num<-as.numeric(unlist(strsplit(subDB[[4]][[dbkey_row]],split=",")))
	minum_sel<-which(intens_ion.num>=minum)
	neutprod.num<-as.numeric(round(as.numeric(subDB[[1]][dbkey_row,2])-prod_ion.num,2))
	writeLines("\nCandidate neutral losses:")
	NeutLossMatch(neutprod.num,neutloss,err)
	writeLines("\nComplementary product ions:")
	ComplemIons(prod_ion.num,intens_ion.num,neutprod.num,err)
	if (is.null(oldpar)) oldpar<-par(no.readonly=T)
	par(mfrow=c(1,1))
	par(cex=0.6)
	adj_intens<-max(intens_ion.num)*1.1
	adj_prod_min<-min(prod_ion.num)*0.9
	adj_prod_max<-max(prod_ion.num)*1.1
	plot(prod_ion.num,intens_ion.num,type="h",xlab="m/z",ylab="ion intensity",
	 main=round(as.numeric(subDB[[1]][dbkey_row,2]),2),xlim=c(adj_prod_min,adj_prod_max),ylim=c(0,adj_intens))
	text(prod_ion.num[minum_sel],intens_ion.num[minum_sel],prod_ion.num[minum_sel],pos=3)
	if (!is.null(nl)) {
	 text(prod_ion.num[minum_sel],intens_ion.num[minum_sel],neutprod.num[minum_sel],pos=3,offset=1.5,col=2)
	 writeLines("\nIf you enter 0 for the next question, you exit the function.")
	 i=2.5
	 repeat{
	  chosen_ion<-as.numeric(readline("From which product ion do you want to see the neutral losses?"))
	  if(chosen_ion==0)break
	  sl<-which(prod_ion.num%in%chosen_ion)
	  neutprod.num<-as.numeric(round(prod_ion.num[sl]-prod_ion.num,2))
	  NeutLossMatch(neutprod.num,neutloss,err)
	  text(prod_ion.num[minum_sel],intens_ion.num[minum_sel],neutprod.num[minum_sel],pos=3,offset=i,col=3)
	  i=i+1
	 }
	}
	if (is.null(oldpar)) par(oldpar)
}

# MSMSplot(25120,base.dir,finlist,SubDB="QTOFneg")