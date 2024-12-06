library(Spectra)
library(dplyr)
library(MSnbase)


mzxml_files <- list.files(path = "D:/amiens/Ahlam/XCMS/flaxseed/S19", pattern = "\\.mzXML$", full.names = TRUE)
spectra <- readMSData(mzxml_files, mode = "onDisk")


library(Spectra)

library(purrr)

extract_compound_data_optimized <- function(spectra) {
  file_names <- processingData(spectra)@files
  file_names_basename <- basename(file_names)
  
  data_list <- map_dfr(seq_along(spectra), function(i) {
    precursor_mz <- if (msLevel(spectra[[i]]) > 1) precursorMz(spectra[[i]]) else NA
    file_name <- file_names_basename[i]
    
    tibble(
      compound_id = rep(i, length(mz(spectra[[i]]))),
      mz = mz(spectra[[i]]),
      intensity = intensity(spectra[[i]]),
      precursorMz = precursor_mz,
      msLevel = msLevel(spectra[[i]]),
      scan_number = scanIndex(spectra[[i]]),
      fileName = file_name
    )
  })
  
  return(data_list)
}
spectra_subset <- spectra[1:100]  
compound_data_subset <- extract_compound_data_optimized(spectra_subset)

head(compound_data_subset)

metadata <- data.frame(
  name = c("dataset_name", "creator", "source", "source_version", "organism", "url", "source_date"),
  value = c("Flaxseed Study", "Ahlam", "Flaxseed MS Analysis", "1.0", "Unknown", "http://example.com", Sys.Date()),
  stringsAsFactors = FALSE
)

msms_spectra <- tibble(
  compound_id = compound_data_subset$compound_id,  
  mz = compound_data_subset$mz,                     
  intensity = compound_data_subset$intensity         
)
compound_data_subset <- compound_data_subset %>%
  mutate(
    name = NA_character_,
    inchi = NA_character_,
    inchikey = NA_character_,
    formula = NA_character_,
    exactmass = NA_real_,
    synonyms = list(NA)  
  )
msms_spectra <- msms_spectra %>%
  mutate(spectrum_id = paste0("spectrum_", 1:nrow(msms_spectra)))

msms_spectra$compound_id <- as.character(msms_spectra$compound_id)
compound_data_subset$compound_id <- as.character(compound_data_subset$compound_id)
msms_spectra$mz <- lapply(msms_spectra$mz, function(x) list(x))
msms_spectra$intensity <- lapply(msms_spectra$intensity, function(x) list(x))
msms_spectra <- msms_spectra[!is.na(msms_spectra$mz) & !is.na(msms_spectra$intensity), ]

str(msms_spectra)
createCompDb(x = compound_data_subset, 
             metadata = metadata, 
             msms_spectra = msms_spectra, 
             path = ".", 
             dbFile = "my_compound_db_ahl.sqlite")

library(DBI)
library(RSQLite)

con <- dbConnect(RSQLite::SQLite(), "my_compound_db_ahl.sqlite")
dbListTables(con)
dbReadTable(con, "ms_compound")
dbReadTable(con, "synonym")
spectra_object <- Spectra(msms_spectra)
compound_db <- CompDb(con)
insertSpectra(compound_db, spectra_object)
dbReadTable(con, "msms_spectrum")
dbReadTable(con, "msms_spectrum_peak")


