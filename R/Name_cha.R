Name_cha<-function(sel,Names_lst,Assoc,i){
	name.cha<-as.character()
	for (j in sel) {				# sel contains at least two entries, i.e. the first and one of the other entries of Assoc_line
	 name<-Names_lst[[j]][Assoc[i,j],2]	# Assoc[i,j] refers to the COMPID or row of 
							# the compound.csv-derived data.frame, the second
							# column of which contains the COMPNAME 
	 # sel refers to a sub-database that is entered into Names_lst 
	 # in the same order as its column number in the Assoc file
	 name.cha<-append(name.cha,name)
	}
#	print(name.cha)
	return(name.cha)
}
