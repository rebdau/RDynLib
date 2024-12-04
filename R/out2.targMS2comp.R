out2.targMS2comp<-function(dbkey1,dbkey2,subDB,AnalMS,thr1,thr2,IntThres=NULL){
	if (AnalMS=="FT") {
	 if (is.null(IntThres)) IntThres<-100
	} else {
	 if (is.null(IntThres)) IntThres<-5
	}
	out<-targMS2comp(dbkey1,dbkey2,subDB,AnalMS,IntThres=IntThres)
	outcm<-((out$forward_ions>=thr1)|(out$reverse_ions>=thr1))&(out$prec_nr_ions>2)&(out$cmp_nr_ions>2)
	if (is.na(outcm)){
	 return(NA)
	}else{
	 if (outcm==TRUE){
	  if (is.nan(out$dot_ions)){
	   return(0)
	  }else{
	   if (out$dot_ions>thr2){
	    return(out$dot_ions)
	   }else{
	    return(NA)
	   }
	  }
	 }else{
	  return(NA)
	 }
	}
}
