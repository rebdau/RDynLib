MS3plot<-function(dbkey,TMP,prcx,wh){
	subDB<-TMP[[1]]
#	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	par(cex=prcx)
	ms3id.row<-which(subDB[[5]]==dbkey)
	ms3id.int<-subDB[[6]][ms3id.row]
#	ms3id.int<-as.integer(unlist(strsplit(
#				subDB[[10]][[dbkey_row]],split=",")))
	prod_ion.int<-as.integer(unlist(strsplit(
				subDB[[7]][[ms3id.row[wh]]],split=",")))
	intens_ion.num<-as.numeric(unlist(strsplit(
				subDB[[8]][[ms3id.row[wh]]],split=",")))
	adj_intens<-max(intens_ion.num)*1.1
	adj_prod_min<-min(prod_ion.int)*0.9
	adj_prod_max<-max(prod_ion.int)*1.1
	plot(prod_ion.int,intens_ion.num,type="h",
				xlab="m/z",ylab="ion intensity",
	 			main=paste("MS3: m/z ",round(subDB[[9]][ms3id.row[wh]],0),sep=""),
				xlim=c(adj_prod_min,
				adj_prod_max),ylim=c(0,adj_intens))
	text(prod_ion.int,intens_ion.num,prod_ion.int,pos=3)
}