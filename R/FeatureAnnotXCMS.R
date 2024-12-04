# Annotating peaks in peak groups. 
# 13C isotopes, buffer adduct and ion/neutral complex annotation.
# Works only on nominal mass
# fc is a vector of sample groups
# y11 is the mass of the buffer adduct 
# (formic acid: 46.0 Da, acetic acid: 60.0 Da)

FeatureAnnotXCMS<-function(XCMS,y11,fc){
	buf<-y11						
	dt=XCMS
	anno<-as.character(c(rep("ann",length(dt[,1]))))
	anno.1<-as.character(c(rep("ann",length(dt[,1]))))
	anno.2<-as.character(c(rep("ann",length(dt[,1]))))
	mx<-max(dt$CON.new)		# highest feature group number
	fc3<-max(fc)			# number of sample groups
	lfc<-length(fc)+13		# all columns of xcms.tsv file
	findmn<-fc3+lfc+4			# all columns of xcms.tsv + 
			#additional columns of feature grouping + 1??
	f=2
	repeat{
		f=f+1				# starts at feature group 3??
		dt.s<-dt[dt$CON.new==f,]
		o<-order(dt.s$mzmed)
		dt.s<-dt.s[o,]
		if (length(dt.s[,1])<2){
			if (f>mx) break
			next
		}
		mng<-c()
		i=1
		repeat{		# bereken gemiddelde waarde van de pieken
			ir<-as.numeric(dt.s[i,14:lfc])
			mn<-mean(ir,na.rm=T)
			mng<-rbind(mng,mn)
			if (i==length(dt.s[,1])) break
			i=i+1
		}
		dt.s<-cbind(dt.s,mng)	# hier wordt een kolom met de 
				#gemiddelde waarde van de pieken 
				#over de verschillende chromatogrammen toegevoegd
		z3=outer(dt.s$mzmed,dt.s$mzmed,"-")
		z3=z3[lower.tri(z3)]
		z3<-round(z3,digits=1)	# verschillen worden herleid tot 
				#nominale waarden, als je afrondt op een 
				#aantal decimalen, doe het dan ook verder bij (A)
		k=length(dt.s$mzmed)
		i=1
		rw=as.integer(c())
		while (i<k){
			rw.i<-seq(i+1,k,1)
			rw<-append(rw,rw.i)
			i=i+1
		}
		cw=as.integer(c())
		i=1
		ki<-k-1
		while (i<k){
			cw.i<-rep(i,ki)
			cw<-append(cw,cw.i)
			i=i+1	
			ki=ki-1
		}
		mtr<-data.frame(cbind(cw,rw,z3))	# coordinaten voor pieken
				#terug te vinden die de massaverschillen opleveren
		u=0
		s=length(mtr[,3])
		if (s==0){
			if (f>mx) break
			next
		}
		repeat{
			u=u+1
			if (round(mtr[u,3],1)==round(dt.s[mtr[u,1],5]+1,1)) {  
					# hier wordt gekeken of het massaverschil 
					#gelijk is aan dat van de neutrale component 
					#(vandaar +1), voer decimalen in als je meer 
					#precies wilt zijn (A)
				nmp<-dt.s[mtr[u,2],1]
				anno.2[which(dt[,1]==nmp)]="dimer"
#			}else{
#				if (u==s) break
#				next
			}
			if (u==s) break
		}
		s1=length(mtr[z3==1.0,3])	# (A) indien nodig, voer decimalen 
					#in voor meer accuraatheid
		s2=length(mtr[z3==buf,3])	# (A) definieer buf meer accuraat 
					#indien nodig
		if (s1==0) {
			if (s2==0) {
				if (f>mx) break
				next
			} else {
				anno.1<-SearchAdduct(mtr,buf,dt.s,dt,s2,anno.1)
			}	
		} else {
			if (s2==0) {
				anno<-SearchC13(mtr,dt.s,dt,s1,anno,findmn)
			} else {
				anno<-SearchC13(mtr,dt.s,dt,s1,anno,findmn)
				anno.1<-SearchAdduct(mtr,buf,dt.s,dt,s2,anno.1)
			}
		}
		if (f>mx) break
	}
	xcms3<-cbind(dt,anno,anno.1,anno.2)
	return(xcms3)
}

