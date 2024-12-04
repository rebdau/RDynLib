allCSPP<-function(row.start,row.end,nr_col,comp_add,thr1=NULL,thr2=NULL,thr3=NULL){
	if (is.null(thr1)) thr1=3
	if (is.null(thr2)) thr2=0.1
	if (is.null(thr3)) thr3=0.1
	cspp.res<-data.frame(compid.sub=integer(),compid.prod=integer(),conv.type=character(),
			ions.prod=as.integer(),ave.common=as.numeric(),ave.dot=as.numeric(),
			stringsAsFactors=FALSE)
	i=row.start
	k=1
	while(i<=row.end){
	 j=3
	 repeat{
	  compid.sub<-as.integer(comp_add[i,1])
	  compid.prod<-as.integer(strsplit(comp_add[i,j],"!!")[[1]][3])
	  if(is.na(compid.prod)){
	   if(j==nr_col)break
	   j=j+1
	   next
	  }else{
	   cspp.res[k,1]<-compid.sub
	   cspp.res[k,2]<-compid.prod
	   cspp.res[k,3]<-colnames(comp_add)[j]
	   cspp.res[k,4]<-as.integer((strsplit(strsplit(comp_add[i,j],"!!")[[1]][2],"!")[[1]][1]))
	   cspp.res[k,5]<-as.numeric((strsplit(strsplit(comp_add[i,j],"!!")[[1]][2],"!")[[1]][2]))
	   cspp.res[k,6]<-as.numeric((strsplit(strsplit(comp_add[i,j],"!!")[[1]][2],"!")[[1]][3]))
	  }
	  if(j==nr_col)break
	  j=j+1
	  k=k+1
	 }
	 i=i+1
	}
	cspp.res<-cspp.res[(cspp.res[,4]>thr1)&(cspp.res[,5]>thr2)&(cspp.res[,6]>thr3),]
	cspp.res<-cspp.res[order(cspp.res[,1]),]
	cspp.res<-cspp.res[!is.na(cspp.res[,1]),]
	return(cspp.res)
}
