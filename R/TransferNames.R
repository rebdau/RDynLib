TransferNames<-function(Names_lst,SBDB){
	Assoc<-Names_lst[[4]]			# DynLibDBassociation file / finlist[[4]], [[10]] or [[11]] 
	No.name=as.character(c("NULL","!"))
	Mult.1<-c()
	Mult.2<-c()
	Mult.3<-c()
	Mult.4<-c()
	i=1
	repeat{
	 Assoc_line<-as.integer(Assoc[i,])	# e.g.: [1] 11026    NA 19636    NA    NA    NA
	 if (sum(is.na(Assoc_line))==5) {	# if all entries, except the first, for Assoc row are NA, skip
	  if (i==dim(Assoc)[1]) break
	  i=i+1
	  next
	 }else{
	  sel<-which(!is.na(Assoc_line))	# sel is Assoc column(s) for which the cell(s) is(are) not NA, e.g. for the example Assoc_line, sel is: [1] 1 3
	  name.cha<-Name_cha(sel,Names_lst,Assoc,i)	# returns names for Assoc[i,sel]-containing COMPID, here COMPID 11026 and 19636 of the FTneg and QTOFneg subDB are proanthocyanidines
	  HasName<-Has_Name(name.cha,No.name,sel)		# returns the 'sel' entries for which the names are not NULL or !; here: [1] 1 3
	  if (sum(HasName)==0) {		# go to next Assoc entry if all names are either NULL or starting with ! for the considered Assoc_line
	   if (i==dim(Assoc)[1]) break
	   i=i+1
	   next
	  }
	 }
	 full.lst<-Real_Name(HasName,sel,Assoc_line,Names_lst,SBDB)
	 Names_lst<-full.lst[[1]]
	 if (length(full.lst[[2]])!=1) {
	  Mult.1<-append(Mult.1,full.lst[[2]][1])
	  Mult.2<-append(Mult.2,full.lst[[2]][2])
	  Mult.3<-append(Mult.3,full.lst[[2]][3])
	  Mult.4<-append(Mult.4,full.lst[[2]][4])
	 }
	 if (i==dim(Assoc)[1]) break
	 i=i+1
	}
	print(paste(length(Mult.1),class(Mult.2),class(Mult.3),class(Mult.4),sep=" "))
	Mult.frame<-data.frame(as.integer(Mult.1),Mult.2,Mult.3,Mult.4)
	if(dim(Mult.frame)[1]!=0){
	 Mult.frame[,2]<-as.character(Mult.frame[,2])
	 Mult.frame[,3]<-as.character(Mult.frame[,3])
	 Mult.frame[,4]<-as.character(Mult.frame[,4])
	 colnames(Mult.frame)<-c("COMPID",SBDB[1],SBDB[2],SBDB[3])
	}else{
	 Mult.frame<-"empty"
	}
	return(list(Names_lst,Mult.frame))
}

