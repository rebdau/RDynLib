# Yields a list of 'proc' and 'fc', the former representing the 'nodes.txt'
# and the latter including the sample grouping. If this function is used as
# part of XCMSgrouping.R, then proc is empty as a 'nodes.txt' file has to be
# created. If this function is part of MSnspecplot.R, then the 'nodes.txt'
# file is explicitly needed to create the MS1 spectrum.


LoadProc<-function(exp,TMP){
	Exp<-paste("Exp",exp,sep="")
	setwd(TMP[[5]])
	chosdir<-paste(TMP[[5]],"/",list.files()[grep(Exp,list.files())],sep="")
	setwd(chosdir)
	nodes.pres<-which(list.files()%in%"nodes.txt")
	if (length(nodes.pres)!=0) {
	 proc<-read.table("nodes.txt",header=T,sep="\t",stringsAsFactor=F)
	} else {
	 print(paste("No nodes.txt file present in folder of experiment",
				exp,sep=" "))
	 proc<-"NA"
	}
	fc.pres<-which(list.files()%in%"fc.txt")
	if (length(fc.pres)!=0) {
	 fc<-read.table("fc.txt",header=T,sep="\t")[,2]
	} else {
	 print(paste("No fc.txt file present in folder of experiment",
				exp,sep=" "))
	 fc<-as.integer(c(NA))
	}
	return(list(proc,fc))
}
