MS2simOfIdentMasses<-function(inp.x,finlist,subDB,AnalMS,nOc){	# finlist does not appear anymore in this function and its nested functions! Could be removed!
	sel<-c(1:nOc)
	inp2<-inp.x[sel,]							
	inp2<-inp2[order(inp2[,5]),]
	inp2[,5]<-round(inp2[,5],digits=2)
	res1<-as.integer()
	res2<-as.character()
	res3<-as.integer()
	res4<-as.character()
	res5<-as.numeric()
	i=1
		# 1. repeat over all particular masses
	while (i<dim(inp2)[1]){							
	 j=i+1
	print(i)
	 if(j>=dim(inp2)[1])break
		# 2. repeat across all peaks for a particular mass
	 repeat{
		# 3. check when mass of j becomes different than that of i.
		#As long as it is not different, you immediately jump to the end
		#whereby j is increased by 1									
	  if (inp2[i,5]!=inp2[j,5]){						
		#maken van array om op te vullen met similariteitswaarden 
		#gebaseerd op 'out.targMS2comp'
	   sim.mat<-array(NA,dim=c(j-i,j-i))
	   jk=j-1
	   dbkey.vec<-inp2[i:jk,1]
		#contains COMPIDs for all entries having the same mass
	   for (k in 1:dim(sim.mat)[1]){
		#computes the dot.products for all pairwise combinations 
		#of the set of compounds with the same mass
	    for (l in 1:dim(sim.mat)[1]){
	     sim.mat[k,l]<-out.targMS2comp(dbkey.vec[k],dbkey.vec[l],subDB,AnalMS)
	    }
	   }
	   sim.tri<-sim.mat[lower.tri(sim.mat)]
		#FUNCTIE om uit de similariteitsmatrix de juiste namen toe te kennen.
	   if (length(sim.tri[(!is.na(sim.tri))&sim.tri>0.8])==0){			
			# 4. als er geen similariteit is, ga naar volgende massa i
	    i=j
	    break
	   }else{											
			# 4. als er similariteit is, check het
	    sim.mat.ind<-which(sim.mat>0.8,arr.ind=T)					
		#all those with dot products in the sim.mat matrix
	    sim.mat.ind<-sim.mat.ind[sim.mat.ind[,2]<sim.mat.ind[,1],]		
		#select those above the diagonal in the sim.mat matrix
	    if (is.matrix(sim.mat.ind)){							
		#if there is only one line in the sim.mat.ind matrix, then it will
		#not be recognized as a matrix, you have to correct for that
	     sim.mat.ind<-array(sim.mat.ind,dim=dim(sim.mat.ind))
	    }else{
	     sim.mat.ind<-array(sim.mat.ind,dim=c(1,2))
	    }
	    sim.i=1
	    repeat{											
			# 5. check similariteit voor elke j
	     if (inp.x[dbkey.vec[sim.mat.ind[sim.i,1]],9]==
	                        inp.x[dbkey.vec[sim.mat.ind[sim.i,2]],9]){	
			# 6. als beide tot zelfde experiment horen, ga nr volgende j
	      if (sim.i==dim(sim.mat.ind)[1]) break
	      sim.i=sim.i+1
	      next
	     }else{											
			# 6. vervang de naam voor j of i
	      grep1.one<-grep("!",inp.x[dbkey.vec[sim.mat.ind[sim.i,1]],3])     
		#dbkeyvec[sim.mat.ind[sim.i,1]] geeft de COMPID, dus hier wordt
		#de naam bekeken van de geselecteerde entries met zelfde massa 
		#en hoge match
	      grep1.two<-grep("NULL",inp.x[dbkey.vec[sim.mat.ind[sim.i,1]],3])  
		#grep(pattern,name)
	      grep1<-as.integer(c(grep1.one,grep1.two))
	      grep2.one<-grep("!",inp.x[dbkey.vec[sim.mat.ind[sim.i,2]],3])
		#check both matched COMPID whether they have ! or NULL in their name
	      grep2.two<-grep("NULL",inp.x[dbkey.vec[sim.mat.ind[sim.i,2]],3])
	      grep2<-as.integer(c(grep2.one,grep2.two))
	      if ((1%in%grep1)&(1%in%grep2)){
	       if (sim.i==dim(sim.mat.ind)[1]) break
	       sim.i=sim.i+1
	       next
	      }else{
	       if ((1%in%grep1)|(1%in%grep2)){
	        if (1%in%grep1){
	         res1<-append(res1,as.integer(dbkey.vec[sim.mat.ind[sim.i,1]]))
	         res2<-append(res2,as.character(inp.x[dbkey.vec[sim.mat.ind[sim.i,1]],3]))
	         res3<-append(res3,as.integer(dbkey.vec[sim.mat.ind[sim.i,2]]))
	         res4<-append(res4,as.character(inp.x[dbkey.vec[sim.mat.ind[sim.i,2]],3]))
	         res5<-append(res5,as.numeric(sim.mat[sim.mat.ind[sim.i,1],
					sim.mat.ind[sim.i,2]]))
	         if (sim.i==dim(sim.mat.ind)[1]) break
	         sim.i=sim.i+1
	         next
	        }else{
	         res1<-append(res1,as.integer(dbkey.vec[sim.mat.ind[sim.i,2]]))
	         res2<-append(res2,as.character(inp.x[dbkey.vec[sim.mat.ind[sim.i,2]],3]))
	         res3<-append(res3,as.integer(dbkey.vec[sim.mat.ind[sim.i,1]]))
	         res4<-append(res4,as.character(inp.x[dbkey.vec[sim.mat.ind[sim.i,1]],3]))
	         res5<-append(res5,as.numeric(sim.mat[sim.mat.ind[sim.i,1],
					sim.mat.ind[sim.i,2]]))
	         if (sim.i==dim(sim.mat.ind)[1]) break
	         sim.i=sim.i+1
	         next
	        }
	       }else{
	        if (sim.i==dim(sim.mat.ind)[1]) break
	        sim.i=sim.i+1
	        next
	       }
	      }
	     }	# 6.
	     if (sim.i==dim(sim.mat.ind)[1]) break
	     sim.i=sim.i+1
	    }		# 5.
	    i=j
	    break
	   }		# 4. hier werd gecheckt of er al of niet similariteiten waren
	  }else{	# 3. if mass j is not different from i, increase j
	   j=j+1
	  }		# 3.
	  if(j>=dim(inp2)[1])break
	 }		# 2.
	}		# 1.
	re<-data.frame(res1,res2,res3,res4,res5)
	colnames(re)<-c("COMPID1","NAME1","COMPID2","NAME2","DOTPRODUCT")
	res<-re[order(re[,1],-re[,5]),]
	res.ind<-duplicated(res[,1])
	ress<-res[!res.ind,]											
	return(ress)
}

