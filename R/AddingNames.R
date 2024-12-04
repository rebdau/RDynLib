AddingNames<-function(base.dir,SubDB){
	if (SubDB=="FTneg") {
	 subdb_dir<-paste(base.dir,"/database_FTMS_neg",sep="")
	 csv_dir<-paste(base.dir,"/database_FTMS_neg/CSV",sep="")
	} else if (SubDB=="FTpos") {
	 subdb_dir<-paste(base.dir,"/database_FTMS_pos",sep="")
	 csv_dir<-paste(base.dir,"/database_FTMS_pos/CSV",sep="")
	} else if (SubDB=="QTOFneg") {
	 subdb_dir<-paste(base.dir,"/database_QTOF_neg",sep="")
	 csv_dir<-paste(base.dir,"/database_QTOF_neg/CSV",sep="")
	} else if (SubDB=="QTOFpos") {
	 subdb_dir<-paste(base.dir,"/database_QTOF_pos",sep="")
	 csv_dir<-paste(base.dir,"/database_QTOF_pos/CSV",sep="")
	} else {
	 print('For the subDB, you can choose between FTneg, FTpos, QTOFneg and QTOFpos')
	}
	setwd(subdb_dir)
	mat<-read.table("matches.txt",header=F,fill=T,stringsAsFactors=F)   
		# There should be 9 columns
	setwd(csv_dir)
	com<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)   
		# There should be 16 columns
	com[,3]<-as.character(com[,3])
	nam2<-"!"
	i=1
	ln<-dim(mat)[1]
	repeat{
	 if (mat[i,6]==""){
	  if (i==ln) break
	  i=i+1
	  next
	 }
	 com.ind<-as.integer(mat[i,2])
	 string<-mat[i,6]
	 ind<-as.integer(gsub("[^0-9]","",string))
	 nam<-com[ind,3]
	 glo.c<-mat[i,7]
	 glo.c1<-as.integer(gsub("[^0-9]","",glo.c))
	 glo.s<-mat[i,9]
	 glo.s1<-unlist(strsplit(as.character(glo.s),"="))
	 glo.s2<-round(as.numeric(glo.s1[2]),2)
	 if (nam=="NULL") {
	  if (i==ln) break
	  i=i+1
	  next
	 }else{
	  com[com.ind,3]<-paste(nam2,glo.c1,nam2,glo.s2,nam2,nam,sep="")
	 }
	 if (i==ln) break
	 i=i+1
	}
	write.table(com,"compound.csv",sep="\t",row.names=F)
}
