allGNPS<-function(row.start,row.end,nr_col,gnps_add,thr1=NULL,thr2=NULL,thr3=NULL){
	if (is.null(thr1)) thr1=3
	if (is.null(thr2)) thr2=0.1
	if (is.null(thr3)) thr3=0.1
	gnps.res<-data.frame(compid.sub=as.integer(),compid.prod=as.integer(),mass.diff=as.numeric(),
			ions.prod=as.integer(),ave.common=as.numeric(),ave.dot=as.numeric(),
			stringsAsFactors=FALSE)
	i=row.start
	k=1
	while(i<=row.end){
	 j=2
	 repeat{
	  compid.sub<-as.integer(gnps_add[i,1])
	  compid.prod<-as.integer(strsplit(strsplit(gnps_add[i,j],"!!")[[1]][2],"!")[[1]][4])
	  if(is.na(compid.prod)){
	   if(j==nr_col)break
	   j=j+1
	   next
	  }else{
	   gnps.res[k,1]<-compid.sub
	   gnps.res[k,2]<-compid.prod
	   gnps.res[k,3]<-as.numeric(strsplit(gnps_add[i,j],"!!")[[1]][3])
	   gnps.res[k,4]<-as.integer(strsplit(strsplit(gnps_add[i,j],"!!")[[1]][2],"!")[[1]][1])
	   gnps.res[k,5]<-as.numeric((strsplit(strsplit(gnps_add[i,j],"!!")[[1]][2],"!")[[1]][2]))
	   gnps.res[k,6]<-as.numeric((strsplit(strsplit(gnps_add[i,j],"!!")[[1]][2],"!")[[1]][3]))
	  }
	  if(j==nr_col)break
	  j=j+1
	  k=k+1
	 }
	 i=i+1
	}
	gnps.res<-gnps.res[(gnps.res[,4]>thr1)&(gnps.res[,5]>thr2)&(gnps.res[,6]>thr3),]
	gnps.res<-gnps.res[order(gnps.res[,1]),]
	gnps.res<-gnps.res[!is.na(gnps.res[,1]),]
	return(gnps.res)
}
