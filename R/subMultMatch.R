subMultMatch<-function(msnSect,unk,MSnlist,minIons){
	ENTRYretained<-as.integer()
	for (i in 1:dim(msnSect)[1]) {
	 ENTRYretained<-c(ENTRYretained,seq(msnSect[i,1],msnSect[i,2],1))
	}
	ENTRYretained<-as.integer(sort(unique(ENTRYretained)))
	i=1
	repeat{
	 Seq_ms2<-as.integer(unlist(strsplit(unlist(MSnlist[ENTRYretained[i]]),split=",")))
	 CommonIons<-intersect(unk,Seq_ms2)
	 if(length(CommonIons)>minIons){
	  if (i==length(ENTRYretained)) break
	  i=i+1
	  next
	 }else{
	  ENTRYretained<-ENTRYretained[-i]
	  if (i>length(ENTRYretained)) break
	  next
	 }
	}
	return(ENTRYretained)
}