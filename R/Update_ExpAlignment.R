Update_ExpAlignment<-function(TMP,LCal.lst,expnr1,expnr2){
	setwd(TMP[[4]])
	ExpA<-read.table("ExpAlign_Formula.txt",header=T,sep="\t",stringsAsFactor=F)
	RetT<-read.table("RetT_alignment.txt",header=T,sep="\t",stringsAsFactor=F)
	row.sel<-which(ExpA[,1]%in%expnr1)
	if(is.na(ExpA[row.sel,2])){
	 ExpA[row.sel,2]<-as.character(expnr2)
	 ExpA[row.sel,3]<-as.character(LCal.lst[[2]][1])
	 ExpA[row.sel,4]<-as.character(LCal.lst[[2]][2])
	 ExpA[row.sel,5]<-as.character(LCal.lst[[2]][3])
	 ExpA[row.sel,6]<-as.character(LCal.lst[[2]][4])
	 ExpA[row.sel,7]<-as.character(LCal.lst[[2]][5])
	 ExpA[row.sel,8]<-as.character(LCal.lst[[2]][6])
	}else{
	 ExpA[row.sel,2]<-paste(ExpA[row.sel,2],as.character(expnr2),sep=",")
	 ExpA[row.sel,3]<-paste(ExpA[row.sel,3],as.character(LCal.lst[[2]][1]),sep=",")
	 ExpA[row.sel,4]<-paste(ExpA[row.sel,4],as.character(LCal.lst[[2]][2]),sep=",")
	 ExpA[row.sel,5]<-paste(ExpA[row.sel,5],as.character(LCal.lst[[2]][3]),sep=",")
	 ExpA[row.sel,6]<-paste(ExpA[row.sel,6],as.character(LCal.lst[[2]][4]),sep=",")
	 ExpA[row.sel,7]<-paste(ExpA[row.sel,7],as.character(LCal.lst[[2]][5]),sep=",")
	 ExpA[row.sel,8]<-paste(ExpA[row.sel,8],as.character(LCal.lst[[2]][6]),sep=",")
	}
	write.table(ExpA,"ExpAlign_Formula.txt",sep="\t",row.names=F)	
	i=1
	while(i<=dim(LCal.lst[[1]])[1]){
	 sel<-which(RetT[,1]%in%LCal.lst[[1]][i,1])
	 if(is.na(RetT[sel,2])){
	  RetT[sel,2]<-as.character(LCal.lst[[1]][i,7])
	 }else{
	  RetT[sel,2]<-paste(RetT[sel,2],as.character(LCal.lst[[1]][i,7]),sep="$")
	 }
#	 sel<-which(RetT[,1]%in%LCal.lst[[1]][i,7])
#	 if(is.na(RetT[sel,2])){
#	  RetT[sel,2]<-as.character(LCal.lst[[1]][i,1])
#	 }else{
#	  RetT[sel,2]<-paste(RetT[sel,2],as.character(LCal.lst[[1]][i,1]),sep=",")
#	 }
	 i=i+1
	}
	write.table(RetT,"RetT_alignment.txt",sep="\t",row.names=F)
}