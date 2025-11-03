# Hard-coded base path (adjust as needed)
# e.g. base.dir <- "D:/RDynLib/inst"


Get_All_DynLib_files <- function(base.dir) {
  
  # Define sub-directories
  ft.dir    <- file.path(base.dir, "database_FTMS_neg", "CSV")
  ftps.dir  <- file.path(base.dir, "database_FTMS_pos", "CSV")
  syn.dir   <- file.path(base.dir, "database_QTOF_neg", "CSV")
  synps.dir <- file.path(base.dir, "database_QTOF_pos", "CSV")
  assoc.dir <- file.path(base.dir, "DynLib subDB alignment")
  
  # Load data using helper functions
  finlist_FTn  <- Get_FT_DynLib_files(ft.dir)
  finlist_FTp  <- Get_FT_DynLib_files(ftps.dir)
  finlist_Synn <- Get_Synapt_DynLib_files(syn.dir)
  finlist_Synp <- Get_Synapt_DynLib_files(synps.dir)
  
  # load txt files in association directory
  Assoc_COMPID  <- read.table(file.path(assoc.dir, "DynLibDBassociation.txt"), header = TRUE, sep = "\t")
  Assoc_COMPID2 <- read.table(file.path(assoc.dir, "DynLibDBassociation2.txt"), header = TRUE, sep = "\t")
  Assoc_COMPID3 <- read.table(file.path(assoc.dir, "DynLibDBassociation3.txt"), header = TRUE, sep = "\t")
  Assoc_EXPID   <- read.table(file.path(assoc.dir, "subDBexperiment.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  MS1feat       <- read.table(file.path(assoc.dir, "MS1features.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  neut          <- read.table(file.path(assoc.dir, "neutral_loss_candid.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  frag          <- read.table(file.path(assoc.dir, "product_ion_candid.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Return everything in one list
  return(list(
    finlist_FTn = finlist_FTn,
    finlist_FTp = finlist_FTp,
    finlist_Synn = finlist_Synn,
    finlist_Synp = finlist_Synp,
    Assoc_COMPID = Assoc_COMPID,
    Assoc_COMPID2 = Assoc_COMPID2,
    Assoc_COMPID3 = Assoc_COMPID3,
    Assoc_EXPID = Assoc_EXPID,
    MS1feat = MS1feat,
    neut = neut,
    frag = frag
  ))
}
  
# finlist<-Get_All_DynLib_files(base.dir)