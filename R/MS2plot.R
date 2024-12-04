MS2plot<-function(dbkey,TMP,prcx){
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	par(cex=prcx)
	prod_ion.int<-round(as.numeric(unlist(strsplit
				(subDB[[3]][[dbkey_row]],split=","))),0)
	intens_ion.num<-as.numeric(unlist(strsplit
				(subDB[[4]][[dbkey_row]],split=",")))
	adj_intens<-max(intens_ion.num)*1.1
	adj_prod_min<-min(prod_ion.int)*0.9
	adj_prod_max<-max(prod_ion.int)*1.1
	plot(prod_ion.int,intens_ion.num,type="h",xlab="m/z",
		ylab="ion intensity",main=paste("MS2: m/z ",round(as.numeric
		(subDB[[1]][dbkey_row,2]),3),sep=""),
		xlim=c(adj_prod_min,adj_prod_max),
		ylim=c(0,adj_intens))
	text(prod_ion.int,intens_ion.num,prod_ion.int,pos=3)
}