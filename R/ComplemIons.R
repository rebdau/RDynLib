ComplemIons<-function(prod_ion.num,intens_ion.num,neutprod.num,err){
	comp.ion<-c(rep(NA,4))
	i=1
	repeat{
	 complem_ion<-neutprod.num[i]-1.008
	 complem_ion_low<-complem_ion-err
	 complem_ion_high<-complem_ion+err
	 j=1
	 repeat{
	  if((prod_ion.num[j]>=complem_ion_low)&(prod_ion.num[j]<=complem_ion_high)){
	   cmp.tmp<-c(prod_ion.num[i],intens_ion.num[i],prod_ion.num[j],intens_ion.num[j])
	   comp.ion<-rbind(comp.ion,cmp.tmp)
	  }
	  if(j==length(prod_ion.num))break
	  j=j+1
	 }
	 if(i==length(neutprod.num))break
	 i=i+1
	}
	if(is.logical(comp.ion))return()
	comp.ion<-comp.ion[-1,]
	retain_rows<-dim(comp.ion)[1]/2
	colnames(comp.ion)<-c("Ion 1","Intensity 1","Ion 2","Intensity 2")
	comp.ion<-comp.ion[1:retain_rows,]
	print(comp.ion)
}

