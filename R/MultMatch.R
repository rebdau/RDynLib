MultMatch<-function(MSnlist,unk,minIons){
	end_ms2<-length(MSnlist)
	ms2_sect<-data.frame(COMPID.start=integer(),COMPID.stop=integer(),IonCount=integer())
	step=end_ms2-1
	i=1
	j=i+step
	Seq_ms2<-as.integer(unlist(strsplit(unlist(MSnlist[i:j]),split=",")))
		# Seq_ms2 contains all product ions in the considered part of the sub-database
	CommonIons<-intersect(unk,Seq_ms2)
		# Determine the number of ions in common between the CID spectrum of the candidate
		# COMPID and those present in the considered part of the sub-database and check (below)
		# whether they surpass a preset threshold. If so, add a first row to ms2_sect containing
		# the borders of the considered part of the sub-database by referring to the smallest
		# COMPID and the largest COMPID in the sub-database part and add the number of common ions
		# > ms2_sect
		#   i     j length.CommonIons.
		# 1 1 57292                  5
 	if(length(CommonIons)>minIons){
 	 out<-data.frame(i,j,length(CommonIons))
 	 ms2_sect<-rbind(ms2_sect,out)
	} 
	k=1
		# Below, a row of ms2_sect is taken and checked whether the number of common ions (third
		# column) surpasses the preset threshold. If so, the sub-database part bordered by the
		# lowest COMPID (first column) and the largest COMPID (third column) is split into two
		# halves.
	while(ms2_sect[k,3]>minIons){
	 i=ms2_sect[k,1]		# lowest COMPID in the considered sub-database part
	 end_ms2=ms2_sect[k,2]
	 if ((end_ms2-i)<=0) {
	  k=k+1
	  next
	 }
	 step=as.integer(((end_ms2-i)/2)+0.5) 
		# step to reach the COMPID in the middle of the considered sub-database
	 j=i+step			# j is the highest COMPID of the first half of the split sub-database part 
	 repeat{			
		# for each half of the sub-database part, determine the ions in common with the
		# those in the CID spectrum of the COMPID of interest
	  Seq_ms2<-as.integer(unlist(strsplit(unlist(MSnlist[i:j]),split=",")))
		# Seq_ms2 contains all product ions in the considered half of the sub-database part
	  CommonIons<-intersect(unk,Seq_ms2)
		# Determine the number of ions in common between the CID spectrum of the candidate
		# COMPID and those present in the considered half of the sub-database part and check (below)
		# whether they surpass a preset threshold. If so, add a row to ms2_sect containing
		# the borders of the considered half of the sub-database part by referring to the smallest
		# COMPID and the largest COMPID in the sub-database half and add the number of common ions.
		# Example of the full repeat loop for a particular row in ms2_sect:
		# > ms2_sect
		#       i     j length.CommonIons.
		# 1     1 57292                  5
		# 2     1 28647                  5
		# 3 28648 57292                  5
	  if(length(CommonIons)>minIons){
	   out<-data.frame(i,j,length(CommonIons))
	   ms2_sect<-rbind(ms2_sect,out)
	  } 
	  i=j+1
	  j=j+step
	  if(j>end_ms2){
	   if(i<end_ms2){
	    j=end_ms2
	   }else{
	    break
	   }
	  }
	 }
#	 print(paste(i,j,k,split=" "))	#!!!!!
#	 print(ms2_sect)				#!!!!!
	 if ((k<dim(ms2_sect)[1])&(all(ms2_sect[k,]==ms2_sect[dim(ms2_sect)[1],]))) {
	  ms2_sect<-ms2_sect[-k,]
	  k=k+1
	  if (k>dim(ms2_sect)[1]) break
	  next
	 }
	 ms2_sect<-ms2_sect[-k,]
		# In the command line above, the first row of ms2_sect is removed as this contained the
		# info of the initial sub-database part that is now, however, replaced by the info
		# of the two halves of the initial sub-database part.
	 if (k>dim(ms2_sect)[1]) break	
	}
	return(ms2_sect)
}

# msnSect<-MultMatch(MSnlist,unk)

