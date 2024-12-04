ChkNames<-function(LCal.lst,TMP){
	i=1
	while(i<=dim(LCal.lst[[1]])[1]){
	 sel1<-which(TMP[[1]][[1]][,3]%in%LCal.lst[[1]][i,1])
	 sel2<-which(TMP[[1]][[1]][,3]%in%LCal.lst[[1]][i,7])
	 if(TMP[[1]][[1]][sel1,5]!=TMP[[1]][[1]][sel1,5]) print(paste(sel1," and ",sel2," names are ",TMP[[1]][[1]][sel1,5]," and ",TMP[[1]][[1]][sel1,5],sep=""))
	 i=i+1
	}
}