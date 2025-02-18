
library(MSnbase)
library(MsExperiment)
library(Spectra)
library(dplyr)
library(tidyr)
library(BiocParallel)

register(SerialParam())

data_path <- "D:/These/Ahlam/XCMS/flaxseed/S_all"
fls <- dir(data_path, pattern = "\\.mzXML$", full.names = TRUE, recursive = TRUE)
cat("Chargement des fichiers mzXML...\n")
mse <- readMsExperiment(fls)

mse <- filterSpectra(mse, filterEmptySpectra)
cwp <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10)
mse_find <- findChromPeaks(mse, param = cwp, BPPARAM = SerialParam())
msrefine <- refineChromPeaks(mse_find, MergeNeighboringPeaksParam())
'''
mfp <- MatchedFilterParam(fwhm = 15, max = 300, snthresh = 2,steps = 2, binSize = 0.01, mzdiff = 0.01)
mse_mf <- findChromPeaks(mse, param = mfp, BPPARAM = SerialParam())
'''
sample_info <- data.frame(
  sampleNames = basename(fileNames(msrefine)),
  group = sub("^(S\\d+R\\d+).*", "\\1", basename(fileNames(msrefine)))
)

sampleData(msrefine) <- DataFrame(sample_info)

pdp <- PeakDensityParam(
  sampleGroups = sampleData(msrefine)$group,
  bw = 10,
  minFraction = 0.5
)
mse_grp1 <- groupChromPeaks(msrefine, param = pdp)

pyp <- PeakGroupsParam(
  minFraction = 0.9,
  extraPeaks = 1,
  smooth = "loess",
  span = 0.2,
  family = "gaussian"
)
mse_rt_corr <- adjustRtime(mse_grp1, param = pyp)

pdp_corr <- PeakDensityParam(
  sampleGroups = sampleData(mse_rt_corr)$group,
  bw = 10,
  minFraction = 0.5
)
rf_grp_corr <- groupChromPeaks(mse_rt_corr, param = pdp_corr)

if (inherits(mf_grp_corr, "XcmsExperiment")) {
  mse_grp_corr_xcmsn <- as(rf_grp_corr, "XCMSnExp") 
  
  
  fill_param <- FillChromPeaksParam(
    expandMz = 0,
    expandRt = 0,
    ppm = 10
  )
  
  xdata_filled <- fillChromPeaks(mse_grp_corr_xcmsn, param = fill_param)
  
  
  head(xdata_filled)
} else {
  stop("L'objet n'est pas un XcmsExperiment.")
}


corrected_features <- featureDefinitions(xdata_filled)


final_feature_matrix <- data.frame(
  `rt(min)` = corrected_features$rtmed / 60, 
  mzmed = corrected_features$mzmed,         
  feature_name = paste0(                     
    "M", round(corrected_features$mzmed, 0), 
    "T", round(corrected_features$rtmed, 0)
  )
)


feature_indices <- featureValues(xdata_filled, value = "index")
final_feature_matrix <- cbind(final_feature_matrix, feature_indices)


library(dplyr)
final_feature_matrix <- final_feature_matrix %>%
  group_by(feature_name) %>%
  mutate(compound_id = paste0("DYNLIB", sprintf("%08d", cur_group_id()))) %>%
  ungroup()


output_path <- "D:/These/Ahlam/XCMS/flaxseed/S_all/feature_matrix_rf.txt"
write.table(
  final_feature_matrix,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


head(final_feature_matrix)


library(dplyr)
final_feature_matrix <- final_feature_matrix %>%
  group_by(feature_name) %>%
  mutate(compound_id = paste0("DYNLIB", sprintf("%08d", cur_group_id()))) %>%
  ungroup()


output_path <- "D:/These/Ahlam/XCMS/flaxseed/S_all/feature_matrix_mse.txt"
write.table(
  final_feature_matrix,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


head(final_feature_matrix)




library(MSnbase)
library(MsExperiment)
library(Spectra)
library(dplyr)
library(tidyr)
library(BiocParallel)


register(SerialParam())

data_path <- "D:/These/Ahlam/XCMS/flaxseed/S_all"
fls <- dir(data_path, pattern = "\\.mzXML$", full.names = TRUE, recursive = TRUE)
cat("Chargement des fichiers mzXML...\n")
mse <- readMsExperiment(fls)

process_mse_data <- function(mse) {
  
  mse <- filterSpectra(mse, filterEmptySpectra)
  cwp <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10)
  mse_find <- findChromPeaks(mse, param = cwp, BPPARAM = SerialParam())
  msrefine <- refineChromPeaks(mse_find, MergeNeighboringPeaksParam())
 
  sample_info <- data.frame(
    sampleNames = basename(fileNames(msrefine)),
    group = sub("^(S\\d+R\\d+).*", "\\1", basename(fileNames(msrefine)))
  )
  
  sampleData(msrefine) <- DataFrame(sample_info)
  
  pdp <- PeakDensityParam(
    sampleGroups = sampleData(msrefine)$group,
    bw = 10,
    minFraction = 0.5
  )
  mse_grp1 <- groupChromPeaks(msrefine, param = pdp)
  
  pyp <- PeakGroupsParam(
    minFraction = 0.9,
    extraPeaks = 1,
    smooth = "loess",
    span = 0.2,
    family = "gaussian"
  )
  mse_rt_corr <- adjustRtime(mse_grp1, param = pyp)
  
  pdp_corr <- PeakDensityParam(
    sampleGroups = sampleData(mse_rt_corr)$group,
    bw = 10,
    minFraction = 0.5
  )
  rf_grp_corr <- groupChromPeaks(mse_rt_corr, param = pdp_corr)
  
  
  return(rf_grp_corr)
}

library(dplyr)

library(MSnbase)
library(dplyr)



extract_spectra_data <- function(spectra_obj, ms_level, intensity_threshold = 0) {
  
  spectra_metadata <- spectraData(spectra_obj) %>%
    as.data.frame()
  
  
  spectra_data <- spectra_metadata %>%
    mutate(
      ms_level = ms_level,
      mz_list = sapply(seq_along(spectra_obj), function(i) {
        mz_values <- mz(spectra_obj)[[i]]
        intensity_values <- intensity(spectra_obj)[[i]]
        if (length(mz_values) > 0) {
          filtered_indices <- which(intensity_values > intensity_threshold)
          if (length(filtered_indices) > 0) {
            paste(mz_values[filtered_indices], collapse = ", ")
          } else {
            NA_character_
          }
        } else {
          NA_character_
        }
      }),
      intensity_list = sapply(seq_along(spectra_obj), function(i) {
        intensity_values <- intensity(spectra_obj)[[i]]
        if (length(intensity_values) > 0) {
          filtered_indices <- which(intensity_values > intensity_threshold)
          if (length(filtered_indices) > 0) {
            paste(intensity_values[filtered_indices], collapse = ", ")
          } else {
            NA_character_
          }
        } else {
          NA_character_
        }
      })
    )
  
  
  spectra_data$mz_list <- as.character(spectra_data$mz_list)
  spectra_data$intensity_list <- as.character(spectra_data$intensity_list)
  
  cat("Extracted spectra data for MS level", ms_level, "with intensity threshold", intensity_threshold, ":", nrow(spectra_data), "rows\n")
  return(spectra_data)
}

process_single_feature <- function(feature_index, feature_matrix, rf_grp_corr, mse, intensity_threshold = 0) {
  if (feature_index > nrow(feature_matrix) || feature_index < 1) {
    warning(paste("Feature index", feature_index, "is out of range. Skipping."))
    return(NULL)
  }
  
  feature <- feature_matrix[feature_index, ]
  cat("\nProcessing feature index:", feature_index, "\n")
  print(feature)
  
  ms2 <- featureSpectra(rf_grp_corr, msLevel = 2, features = feature_index)
  cat("Number of MS2 spectra found:", length(ms2), "\n")
  
  if (length(ms2) == 0) {
    warning(paste("No MS2 spectra found for feature index:", feature_index))
    return(NULL)
  }
  
  feature_data <- list()
  
  for (origin in unique(spectraData(ms2)$dataOrigin)) {
    cat("Processing data origin:", origin, "\n")
    
    ms2_subset <- ms2[spectraData(ms2)$dataOrigin == origin]
    cat("Number of MS2 spectra for this origin:", length(ms2_subset), "\n")
    
    ms3 <- filterMsLevel(spectra(mse), msLevel = 3)
    ms3_filtered <- ms3[spectraData(ms3)$dataOrigin == origin & spectraData(ms3)$precScanNum %in% spectraData(ms2_subset)$scanIndex]
    cat("Number of MS3 spectra filtered by origin:", length(ms3_filtered), "\n")
    
    ms4 <- filterMsLevel(spectra(mse), msLevel = 4)
    ms4_filtered <- ms4[spectraData(ms4)$dataOrigin == origin & spectraData(ms4)$precScanNum %in% spectraData(ms3_filtered)$scanIndex]
    cat("Number of MS4 spectra filtered by origin:", length(ms4_filtered), "\n")
    
    
    ms2_data <- extract_spectra_data(ms2_subset, ms_level = 2, intensity_threshold = intensity_threshold)
    ms3_data <- extract_spectra_data(ms3_filtered, ms_level = 3, intensity_threshold = intensity_threshold)
    ms4_data <- extract_spectra_data(ms4_filtered, ms_level = 4, intensity_threshold = intensity_threshold)
    
    feature_data[[origin]] <- bind_rows(ms2_data, ms3_data, ms4_data)
  }
  
  final_feature_data <- bind_rows(feature_data, .id = "file_origin")
  cat("Final combined data for feature index", feature_index, "has", nrow(final_feature_data), "rows\n")
  
  final_feature_data <- final_feature_data %>%
    mutate(
      feature_name = feature$feature_name,
      compound_id = feature$compound_id
    )
  
  return(final_feature_data)
}


assign_msntree_id <- function(ms_data) {
  scan_index <- ms_data$scanIndex
  prec_scan_num <- ms_data$precScanNum
  ms_level <- ms_data$msLevel
  data_origin <- ms_data$dataOrigin
  
  tree_id <- rep(NA, length(scan_index))
  global_counter <- 0
  
  for (origin in unique(data_origin)) {
    origin_indices <- which(data_origin == origin)
    origin_scan_index <- scan_index[origin_indices]
    origin_prec_scan_num <- prec_scan_num[origin_indices]
    origin_ms_level <- ms_level[origin_indices]
    
    origin_tree_id <- rep(NA, length(origin_scan_index))
    
    for (i in seq_along(origin_scan_index)) {
      if (origin_ms_level[i] == 4) {
        
        parent_ms3 <- which(origin_scan_index == origin_prec_scan_num[i] & origin_ms_level == 3)
        if (length(parent_ms3) > 0) {
          parent_ms2 <- which(origin_scan_index == origin_prec_scan_num[parent_ms3] & origin_ms_level == 2)
          if (length(parent_ms2) > 0) {
            origin_tree_id[i] <- origin_tree_id[parent_ms2]
          }
        }
      } else if (origin_ms_level[i] == 3) {
        
        parent_ms2 <- which(origin_scan_index == origin_prec_scan_num[i] & origin_ms_level == 2)
        if (length(parent_ms2) > 0) {
          origin_tree_id[i] <- origin_tree_id[parent_ms2]
        }
      } else if (origin_ms_level[i] == 2) {
        global_counter <- global_counter + 1
        origin_tree_id[i] <- global_counter
      }
    }
    tree_id[origin_indices] <- origin_tree_id
  }
  
  return(tree_id)
}

create_database <- function(feature_matrix, mse, intensity_threshold) {
  cat("Création de la base de données...\n")
  library(xcms)
  library(dplyr)
  
  rf_grp_corr <- process_mse_data(mse)
  
  db <- bind_rows(lapply(seq_len(nrow(feature_matrix)), function(i) {
    process_single_feature(i, feature_matrix, rf_grp_corr, mse, intensity_threshold)
  }))
  
  if ("file_origin" %in% colnames(db) && "ms_level" %in% colnames(db) && "dataStorage" %in% colnames(db)) {
    db <- db %>% select(-file_origin, -ms_level, -dataStorage)
  }
  
  cat("Assigning MSn tree IDs...\n")
  db <- db %>%
    mutate(
      MSntreeID = assign_msntree_id(db)
    )
  
  if ("mz_list" %in% colnames(db) && "intensity_list" %in% colnames(db)) {
    db <- db %>%
      mutate(
        mz = sapply(mz_list, function(x) ifelse(length(x) > 0, paste(x, collapse = ", "), NA_character_)),
        intensity = sapply(intensity_list, function(x) ifelse(length(x) > 0, paste(x, collapse = ", "), NA_character_))
      ) %>%
      select(-mz_list, -intensity_list)
  }
  
  db <- db %>%
    group_by(scanIndex, dataOrigin) %>%
    summarise(
      peak_id = paste(unique(peak_id), collapse = ", "),
      across(-peak_id, ~ unique(.)[1]), 
      .groups = "drop"
    )
  
  db <- db %>%
    mutate(n_peaks = sapply(strsplit(peak_id, ", "), length))
  
  db <- db %>%
    mutate(
      ExpID = 1,
      name = NA_character_,
      inchi = NA_character_,
      inchikey = NA_character_,
      formula = NA_character_,
      exactmass = NA_real_,
      Smiles = NA_character_,
      ISOTOPE_RATIO = NA_character_,
      DRIFT_TIME = NA_real_,
      COMPOSITION = NA_character_,
      synonyms = list(NA)
    )
  

  cat("Base de données crée avec succès!\n")
  
  return(db)
}







add_spectra_to_database <- function(feature_matrix, existing_db, mse, ppm_tolerance = 20, rt_tolerance = 60, exp_id = NULL) {
  
  
  mse_grp_corr <- process_mse_data(mse)
  
  new_db <- bind_rows(lapply(seq_len(nrow(feature_matrix)), function(i) {
    process_single_feature(i, feature_matrix, mse_grp_corr, mse,intensity_threshold )
  }))
  
  
  if ("file_origin" %in% colnames(new_db) && "ms_level" %in% colnames(new_db) && "dataStorage" %in% colnames(new_db)) {
    new_db <- new_db %>%
      select(-file_origin, -ms_level, -dataStorage)
  }
  
  
  if ("mz_list" %in% colnames(new_db) && "intensity_list" %in% colnames(new_db)) {
    new_db <- new_db %>%
      mutate(
        mz = sapply(mz_list, function(x) ifelse(length(x) > 0, paste(x, collapse = ", "), NA_character_)),
        intensity = sapply(intensity_list, function(x) ifelse(length(x) > 0, paste(x, collapse = ", "), NA_character_))
      ) %>%
      select(-mz_list, -intensity_list)
  }
  
  
  if (!is.null(exp_id)) {
    new_db <- new_db %>%
      mutate(ExpID = exp_id)
  } else {
    
    next_exp_id <- if ("ExpID" %in% names(existing_db)) {
      max(existing_db$ExpID, na.rm = TRUE) + 1
    } else {
      1
    }
    
    
    new_db <- new_db %>%
      mutate(ExpID = next_exp_id)
  }
  
  
  new_db <- new_db %>%
    mutate(
      MSntreeID = assign_msntree_id(new_db)
    )
  
  
  new_db <- new_db %>%
    mutate(
      MSntreeID = paste(MSntreeID, paste0("Exp", ExpID), sep = "_")
    )
  
  
  existing_db$MSntreeID <- as.character(existing_db$MSntreeID)
  new_db$MSntreeID <- as.character(new_db$MSntreeID)
  
  
  updated_db <- bind_rows(existing_db, new_db)
  
  return(updated_db)
}







cat("Chargement de la matrice des features...\n")
feature_matrix <- read.table("D:/These/Ahlam/XCMS/flaxseed/S_all/feature_matrix_rf.txt", sep = "\t", header = TRUE)
feature_matrix$feature_name <- as.character(feature_matrix$feature_name)
feature_matrix$compound_id <- as.character(feature_matrix$compound_id)


filtered_db <- create_database(feature_matrix, mse,intensity_threshold=100)

head(filtered_db)
saveRDS(filtered_db, "D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db.rds")
cat("Le fichier a été enregistré avec succès sous 'filtered_db'.\n")
filtered_db <- readRDS("D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db.rds")
filtered_db
tail(filtered_db)



print(head(filtered_db, 10), width = Inf)

filtered_db <- filtered_db %>%
  arrange(feature_name, MSntreeID, scanIndex) %>%
  group_by(feature_name, MSntreeID) %>%
  mutate(
    unique_tree_id = cumsum(is.na(precScanNum) | !precScanNum %in% lag(scanIndex, default = 0)) 
  ) %>%
  ungroup()

# Étape 2: Calcul du dernier niveau MS pour chaque arbre unique
last_level_per_tree <- filtered_db %>%
  group_by(feature_name, MSntreeID, unique_tree_id) %>%
  summarise(
    last_ms_level = max(msLevel, na.rm = TRUE),  
    .groups = "drop"
  )

# Étape 3: Comptage des arbres fragmentés à chaque niveau pour chaque feature
ms_level_counts <- last_level_per_tree %>%
  group_by(feature_name) %>%
  summarise(
    MS2 = sum(last_ms_level == 2),  
    MS3 = sum(last_ms_level == 3), 
    MS4 = sum(last_ms_level == 4),  
    .groups = "drop"
  )


ms_level_counts_summary <- ms_level_counts %>%
  count(MS2, MS3, MS4) %>%
  rename(feature_number = n)


print(ms_level_counts_summary)
tail(ms_level_counts_summary)




total_features <- sum(ms_level_counts_summary$feature_number)


ms2_counts <- ms_level_counts_summary %>%
  count(MS2, wt = feature_number) %>%  
  mutate(percentage = n / total_features * 100)


ms3_counts <- ms_level_counts_summary %>%
  count(MS3, wt = feature_number) %>%
  mutate(percentage = n / total_features * 100)


ms4_counts <- ms_level_counts_summary %>%
  count(MS4, wt = feature_number) %>%
  mutate(percentage = n / total_features * 100)

sum(ms2_counts$percentage) 
sum(ms3_counts$percentage) 
sum(ms4_counts$percentage) 

#

library(dplyr)
library(ggplot2)
library(tidyr)
ggplot(ms2_counts, aes(x = MS2, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  labs(
    title = "Distribution de MS2 (Pourcentage de features)",
    x = "Nombre d'arbres fragmentés à MS2",
    y = "Pourcentage de features"
  ) +
  ylim(0, 50) +
  xlim(0, 50) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),   
    axis.title.x = element_text(size = 14,face = "bold"),               
    axis.title.y = element_text(size = 14, face = "bold")               
  )


ggplot(ms3_counts, aes(x = MS3, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  labs(
    title = "Distribution de MS3 (Pourcentage de features)",
    x = "Nombre d'arbres fragmentés à MS3",
    y = "Pourcentage de features"
  ) +
  ylim(0, 30) +  
  xlim(0, 20) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),   
    axis.title.x = element_text(size = 14,face = "bold"),               
    axis.title.y = element_text(size = 14, face = "bold")               
  )



ggplot(ms4_counts, aes(x = MS4, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightcoral", color = "black") +
  labs(
    title = "MS4 Distribution (Percentage of Features)",
    x = "Number of Fragmented Trees at MS4",
    y = "Percentage of Features"
  ) +
  ylim(0, 30) + 
  xlim(0, 20) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),   
    axis.title.x = element_text(size = 14,face = "bold"),               
    axis.title.y = element_text(size = 14, face = "bold")             
  )






filtered_db<- filtered_db %>%
  arrange(feature_name, MSntreeID, scanIndex) %>% 
  group_by(feature_name, MSntreeID) %>%
  mutate(
    unique_tree_id = cumsum(is.na(precScanNum) | !precScanNum %in% lag(scanIndex, default = 0))
  ) %>%
  ungroup()


last_level_per_tree <- filtered_db %>%
  group_by(feature_name, MSntreeID, unique_tree_id) %>%
  summarise(
    last_ms_level = max(msLevel, na.rm = TRUE), 
    .groups = "drop"
  )


tree_counts_per_level <- last_level_per_tree %>%
  group_by(last_ms_level) %>%
  summarise(
    tree_count = n(),  # Compter les arbres uniques
    .groups = "drop"
  ) %>%
  arrange(last_ms_level)

# Afficher les résultats
print(tree_counts_per_level)




# Étape 5: Ajustement des données pour calculer les pourcentages et préparer les graphiques
ms_level_counts_summary <- ms_level_counts %>%
  
  mutate(total = MS2 + MS3 + MS4) %>%
  pivot_longer(cols = MS2:MS4, names_to = "ms_level", values_to = "count") %>%
  filter(count > 0) %>%  
  group_by(ms_level) %>%
  summarise(
    total_count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(percentage = total_count / sum(total_count) * 100)



ggplot(ms_level_counts_summary, aes(x = "", y = percentage, fill = ms_level)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +  
  theme_void() +  
  theme(legend.title = element_blank()) +  
  labs(title = "Distribution of Fragmentation Levels (MS2, MS3, MS4)") +
  scale_fill_manual(values = c("MS4" = "red", "MS2" = "blue", "MS3" = "green")) + 
  geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold")  
















