SubDBcombin<-function(finlist,SubDB){
	if((SubDB=="FTneg")|(SubDB=="QTOFneg")){
	 sdb1="FTneg"
	 sdb2="QTOFneg"
	 sdb3="FTpos"
	 sdb4="QTOFpos"
	 Assoc=finlist[[4]]
	 Assoc2=finlist[[10]]
	}else if((SubDB=="FTpos")|(SubDB=="QTOFpos")){
	 sdb1="FTpos"
	 sdb2="QTOFpos"
	 sdb3="FTneg"
	 sdb4="QTOFneg"
	 Assoc=finlist[[10]]
	 Assoc2=finlist[[4]]
	}else{
	 print("No valuable subDB")
	}
	return(list(sdb1,sdb2,sdb3,sdb4,Assoc,Assoc2))
}