library(dplyr)
library(readr)
library(purrr)

read_csv_safely <- function(file_path, sep = "\t", quote = "\"", header = TRUE, col_types = NULL) {
  if (!file.exists(file_path)) {
    stop(paste("Le fichier n'existe pas :", file_path))
  }
  
  tryCatch({
    data <- read_delim(file_path, delim = sep, quote = quote, col_names = header, col_types = col_types)
    return(as.data.frame(data))  
  }, error = function(e) {
    message("Erreur lors de la lecture avec read_delim : ", e$message)
    message("Essai avec read.csv...")
    data <- read.csv(file_path, sep = sep, quote = quote, header = header, stringsAsFactors = FALSE)
    return(as.data.frame(data))  
  })
}

parse_spectral_data <- function(data) {
  data %>%
    mutate(
      mz = map(mz, ~ as.numeric(unlist(strsplit(as.character(.x), ","), use.names = FALSE))),
      intensity = map(intensity, ~ as.numeric(unlist(strsplit(as.character(.x), ","), use.names = FALSE)))
    )
}


compound_path <- "D:/These/RDynLib/RDynLib/database_FTMS_neg/CSV/compound.csv"
ms2_spectra_path <- "D:/These/RDynLib/RDynLib/database_FTMS_neg/CSV/MS2spectra.csv"
ms3_spectra_path <- "D:/These/RDynLib/RDynLib/database_FTMS_neg/CSV/MS3spectra.csv"
experiment_path <- "D:/These/RDynLib/RDynLib/database_FTMS_neg/CSV/experiment.csv"

col_types_spectra <- cols(
  MS2PEAKLIST = col_character(),
  MS2INTENSITYLIST = col_character(),
  MS3PEAKLIST = col_character(),
  MS3INTENSITYLIST = col_character()
)




compound <- read_csv_safely(compound_path)
ms2_spectra <- read_csv_safely(ms2_spectra_path, sep = "\t", col_types = col_types_spectra)
ms3_spectra <- read_csv_safely(ms3_spectra_path, sep = "\t", col_types = col_types_spectra)
experiment <- read_csv_safely(experiment_path, sep = "\t")

ms2_spectra <- ms2_spectra %>%
  select(MS2ID, PARENT_ION, MS2PEAKLIST, MS2INTENSITYLIST) %>%
  rename(COMPID = MS2ID, PrecursorMz = PARENT_ION, mz = MS2PEAKLIST, intensity = MS2INTENSITYLIST) %>%
  mutate(msLevel = 2)

ms3_spectra <- ms3_spectra %>%
  select(MS3ID, PARENT_ION, MS3PEAKLIST, MS3INTENSITYLIST) %>%
  rename(COMPID = MS3ID, PrecursorMz = PARENT_ION, mz = MS3PEAKLIST, intensity = MS3INTENSITYLIST) %>%
  mutate(msLevel = 3)

ms2_spectra <- parse_spectral_data(ms2_spectra)
ms3_spectra <- parse_spectral_data(ms3_spectra)

spectra <- bind_rows(ms2_spectra, ms3_spectra)

final_table <- compound %>%
  left_join(spectra, by = "COMPID")

final_table_Exp <- experiment %>%
  left_join(final_table, by = "EXPID")

print(head(final_table))
print(tail(final_table_Exp))

print(head(final_table_Exp$mz))
print(head(final_table_Exp$intensity))

library(dplyr)

final_table_Exp <- final_table_Exp %>%
  select(-...16, -...17)


colnames(final_table_Exp)[colnames(final_table_Exp) == "MODE"] <- "polarity"
colnames(final_table_Exp)[colnames(final_table_Exp) == "NODENAME"] <- "feature_name"
colnames(final_table_Exp)[colnames(final_table_Exp) == "COMPNAME"] <- "name"
colnames(final_table_Exp)[colnames(final_table_Exp) == "CE"] <- "collisionEnergy"
colnames(final_table_Exp)[colnames(final_table_Exp) == "RETENTION_TIME"] <- "rtime"



missing_in_filtered_db <- setdiff(colnames(final_table_Exp), colnames(filtered_db))


for (col in missing_in_filtered_db) {
  filtered_db[[col]] <- NA
}

missing_in_final_table_Exp <- setdiff(colnames(filtered_db), colnames(final_table_Exp))

for (col in missing_in_final_table_Exp) {
  final_table_Exp[[col]] <- NA
}

colnames(filtered_db)
colnames(final_table_Exp)
print(tail(final_table_Exp))
print(tail(filtered_db,4),width = Inf)

dynlib_version <- dynlib_version[, !colnames(dynlib_version) %in% "rtime"]

colnames(dynlib_version)[colnames(dynlib_version) == "RETENTION_TIME"] <- "rtime"
saveRDS(final_table_Exp, "D:/These/Ahlam/XCMS/flaxseed/S_all2/dynlib_versionr2.rds")
dynlib_version <- readRDS("D:/These/Ahlam/XCMS/flaxseed/S_all2/dynlib_versionr2.rds")
dynlib_version
tail(dynlib_version)
dynlib_version$mz |>head()
filtered_db$polarity <- ifelse(filtered_db$polarity == 0, "neg", filtered_db$polarity)
saveRDS(filtered_db, "D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db_neg.rds")
filtered_db <- readRDS("D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db_neg.rds")
filtered_db
filtered_trees$polarity <- ifelse(filtered_trees$polarity == 0, "neg", filtered_trees$polarity)
saveRDS(filtered_trees, "D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_trees_neg.rds")

selected_trees$polarity <- ifelse(selected_trees$polarity == 0, "neg", selected_trees$polarity)
saveRDS(selected_trees, "D:/These/Ahlam/XCMS/flaxseed/S_all2/selected_trees_neg.rds")



# Afficher le résultat
print(dynlib_version)



