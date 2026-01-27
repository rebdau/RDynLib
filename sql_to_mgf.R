sql_to_mgf <- function(sqlite_path,
                       expid,
                       spectrum_type = NULL,
                       output_file) {
  
  # Packages
  library(Spectra)
  library(CompoundDb)
  library(MsBackendMgf)
  
  # GNPS helper functions
  source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")
  

  ## Load database as Spectra
  cdb <- CompDb(sqlite_path)
  
  cdb <- addJoinDefinition(
    cdb,
    table_a  = "ms_compound",
    table_b  = "feature_matrix",
    column_a = "nodename",
    column_b = "nodename"
  )
  
  sps <- Spectra(cdb)
  
  ## Filtering
  sps <- sps[sps$expid %in% expid]
  
  if (!is.null(spectrum_type)) {
    sps <- sps[sps$spectrum_type %in% spectrum_type]
  }
  
  if (length(sps) == 0)
    stop("No spectra left after filtering")
  
  # Restrict variables to GNPS relevant ones
  sps <- selectSpectraVariables(
    sps,
    c("precursorMz", "rtime", "nodename",  "dataOrigin")
  )
  

  ## Map nodename to feature_id (GNPS requirement)
  sps$feature_id <- sps$nodename
  
  # Drop nodename if not needed anymore
  sps <- selectSpectraVariables(
    sps,
    c("precursorMz", "rtime", "feature_id",  "dataOrigin")
  )
  
  ## Format spectra for GNPS 
  sps_gnps <- formatSpectraForGNPS(sps)
  
  ## Export MGF using MsBackendMgf
  export(
    sps_gnps,
    backend = MsBackendMgf(),
    file    = output_file
  )
  
  invisible(TRUE)
}


sql_to_mgf(
  sqlite_path   = "C:/Users/ament/Desktop/Github/RDynLib/data/DYL_FLAX_FTMSneg.sqlite",
  expid         = 266,
  spectrum_type = "not_assembled",
  output_file   = "C:/Users/ament/Desktop/Github/RDynLib/data/DYL_FLAX_FTMSneg.mgf"
)



