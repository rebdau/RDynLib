ViewStruc2<-function(base.dir,finlist,SubDB,dbkey){
	library(rcdk)
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	if(subDB[[1]][dbkey,6]=="NULL"){
	 print("No structure available")
	}else{
	 smile<-as.character(c(subDB[[1]][dbkey,6]))
	 mol<-parse.smiles(smile)
	 view.molecule.2d(mol[[1]])
	}
}

# ViewStruc2(base.dir,finlist,SubDB="FTneg",dbkey)