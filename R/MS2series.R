MS2series<-function(nrCOMPID,base.dir,finlist,SubDB){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	nrdec1<-ifelse(TMP[[2]]=="FT",0,2)	
	nrdec2<-ifelse(TMP[[2]]=="FT",3,2)
	oldpar<-par(no.readonly=T)
	nrCOMPID<-as.integer(nrCOMPID)
	row.nr<-floor(sqrt(length(nrCOMPID)))
	col.nr<-ceiling(length(nrCOMPID)/row.nr)
	par(mfrow=c(row.nr,col.nr))
	par(cex=0.5)
	i=1
	for(i in 1:length(nrCOMPID)){
	 dbkey<-nrCOMPID[i]
	 dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	 prod_ion.int<-round(as.numeric(unlist(strsplit(
				subDB[[3]][[dbkey_row]],split=","))),nrdec1)
	 intens_ion.num<-as.numeric(unlist(strsplit(
				subDB[[4]][[dbkey_row]],split=",")))
	 adj_intens<-max(intens_ion.num)*1.1
	 adj_prod_min<-min(prod_ion.int)*0.9
	 adj_prod_max<-max(prod_ion.int)*1.1
	 plot(prod_ion.int,intens_ion.num,type="h",
			xlab="m/z",ylab="ion intensity",
	  		main=paste(nrCOMPID[i],"_",round(as.numeric
			(subDB[[1]][dbkey_row,2]),nrdec2),
			collapse=""),xlim=c(adj_prod_min,adj_prod_max),
			ylim=c(0,adj_intens))
	 text(prod_ion.int,intens_ion.num,prod_ion.int,pos=3)
	}
	par(oldpar)
}

# nrCOMPID<-c(101,102,103,104,105)
# MS2series(nrCOMPID,base.dir,finlist,SubDB="FTneg")
