SelectSubDB<-function(base.dir,finlist,SubDB){
	if (SubDB=="FTneg") {
	 csv_dir<-paste(base.dir,"/database_FTMS_neg/CSV",sep="")
	 csvadd_dir<-paste(base.dir,"/database_FTMS_neg/CSV_add",sep="")
	 chosendir<-paste(base.dir,"/database_FTMS_neg/Experimenten",sep="")
	 strucdir<-paste(base.dir,"/database_FTMS_neg/structuren",sep="")
	 finlist.tmp=finlist[[1]]
	 AnalMS="FT"
	} else if (SubDB=="FTpos") {
	 csv_dir<-paste(base.dir,"/database_FTMS_pos/CSV",sep="")
	 csvadd_dir<-paste(base.dir,"/database_FTMS_pos/CSV_add",sep="")
	 chosendir<-paste(base.dir,"/database_FTMS_pos/Experimenten",sep="")
	 strucdir<-paste(base.dir,"/database_FTMS_pos/structuren",sep="")
	 finlist.tmp=finlist[[2]]
	 AnalMS="FT"
	} else if (SubDB=="QTOFneg") {
	 csv_dir<-paste(base.dir,"/database_QTOF_neg/CSV",sep="")
	 csvadd_dir<-paste(base.dir,"/database_QTOF_neg/CSV_add",sep="")
	 chosendir<-paste(base.dir,"/database_QTOF_neg/Experimenten",sep="")
	 strucdir<-paste(base.dir,"/database_QTOF_neg/structuren",sep="")
	 finlist.tmp=finlist[[3]]
	 AnalMS="QTOF"
	} else if (SubDB=="QTOFpos") {
	 csv_dir<-paste(base.dir,"/database_QTOF_pos/CSV",sep="")
	 csvadd_dir<-paste(base.dir,"/database_QTOF_pos/CSV_add",sep="")
	 chosendir<-paste(base.dir,"/database_QTOF_pos/Experimenten",sep="")
	 strucdir<-paste(base.dir,"/database_QTOF_pos/structuren",sep="")
	 finlist.tmp=finlist[[9]]
	 AnalMS="QTOF" 
	} else {
	 print('For the subDB, you can choose between FTneg, FTpos, QTOFneg and QTOFpos')
	}
	return(list(finlist.tmp,AnalMS,csv_dir,csvadd_dir,chosendir,strucdir))
}

# TMP<-SelectSubDB(base.dir,finlist,SubDB="FTneg")