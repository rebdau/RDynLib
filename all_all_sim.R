library(Spectra)


# Load required data
filtered_db <- readRDS("C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db_neg.rds")
selected_trees_neg <- readRDS("C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/selected_trees_neg.rds")
dynlib_version <- readRDS("C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/FTMS_QTOF.rds")

# Créer la nouvelle colonne CE_ramp en fonction de la valeur de collisionEnergy
dynlib_version$CE_ramp <- ifelse(dynlib_version$collisionEnergy == "ramp", TRUE, FALSE)
selected_trees_neg$CE_ramp <- ifelse(selected_trees_neg$collisionEnergy == "ramp", TRUE, FALSE)


# Convert selected_trees_neg to a dataframe
st_sps <- as.data.frame(selected_trees_neg)
st_sps$msLevel <- as.integer(st_sps$msLevel)
st_sps$rtime <- as.numeric(st_sps$rtime)
st_sps$acquisitionNum <- as.integer(st_sps$acquisitionNum)
st_sps$scanIndex <- as.integer(st_sps$scanIndex)
st_sps$dataOrigin <- as.character(st_sps$dataOrigin)
st_sps$polarity <- as.integer(ifelse(st_sps$polarity == "Neg", 0, 1))
st_sps$precScanNum <- as.integer(st_sps$precScanNum)
st_sps$precursorMz <- as.numeric(st_sps$precursorMz)
st_sps$precursorIntensity <- as.numeric(st_sps$precursorIntensity)
st_sps$precursorCharge <- as.integer(st_sps$precursorCharge)
st_sps$isolationWindowLowerMz <- as.numeric(st_sps$isolationWindowLowerMz)
st_sps$isolationWindowTargetMz <- as.numeric(st_sps$isolationWindowTargetMz)
st_sps$isolationWindowUpperMz <- as.numeric(st_sps$isolationWindowUpperMz)

.process_peak <- function(x, split = ",") {
  x <- as.character(x)
  
  x <- gsub("intensity=", "", x)
  x <- gsub("[^0-9.,-]", "", x)
  
  if (x == "" || is.na(x) || x == "NA") {
    return(NA_real_)
  }
  
  valeurs <- unlist(strsplit(x, split = "[, ]+"))
  
  valeurs <- as.numeric(valeurs)
  
  if (all(is.na(valeurs))) {
    return(NA_real_)
  }
  
  return(valeurs)
}


st_sps$mz <- as.character(st_sps$mz)
st_sps$intensity <- as.character(st_sps$intensity)

st_sps$mz <- lapply(st_sps$mz, .process_peak)
st_sps$intensity <- lapply(st_sps$intensity, .process_peak)
st_sps <- Spectra(data = st_sps)
head(st_sps$intensity, 10)
head(st_sps$mz, 10)

#dynlib_version to spectra

dy_sps <- as.data.frame(dynlib_version)
dy_sps$msLevel <- as.integer(dy_sps$msLevel)
dy_sps$rtime <- as.numeric(dy_sps$rtime)
dy_sps$acquisitionNum <- as.integer(dy_sps$acquisitionNum)
dy_sps$scanIndex <- as.integer(dy_sps$scanIndex)
dy_sps$dataOrigin <- as.character(dy_sps$dataOrigin)
dy_sps$polarity <- as.integer(ifelse(dy_sps$polarity == "neg", 0, 1))
dy_sps$precScanNum <- as.integer(dy_sps$precScanNum)
dy_sps$precursorMz <- as.numeric(dy_sps$PrecursorMz)
dy_sps$precursorIntensity <- as.numeric(dy_sps$precursorIntensity)
dy_sps$precursorCharge <- as.integer(dy_sps$precursorCharge)
dy_sps$isolationWindowLowerMz <- as.numeric(dy_sps$isolationWindowLowerMz)
dy_sps$isolationWindowTargetMz <- as.numeric(dy_sps$isolationWindowTargetMz)
dy_sps$isolationWindowUpperMz <- as.numeric(dy_sps$isolationWindowUpperMz)
dy_sps$collisionEnergy <- as.numeric(dy_sps$collisionEnergy)

.process_peak <- function(x, split = ",") {
  x <- as.character(x)
  x <- gsub("[^0-9.,-]", "", x)
  
  if (length(x) > 1) {
    return(as.numeric(unlist(strsplit(x, split = split, fixed = TRUE))))
  }
  
  if (is.na(x) || x == "NA") {
    return(NA_real_)
  }
  
  if (is.numeric(x)) {
    return(x)
  }
  
  valeurs <- unlist(strsplit(x, split = split, fixed = TRUE))
  return(as.numeric(valeurs))
}


dy_sps$mz <- as.character(dy_sps$mz)
dy_sps$intensity <- as.character(dy_sps$intensity)

dy_sps$mz <- lapply(dy_sps$mz, .process_peak)
dy_sps$intensity <- lapply(dy_sps$intensity, .process_peak)
dy_sps <- Spectra(data = dy_sps)
head(dy_sps$intensity, 10)

head(dy_sps$mz, 10)



library(Spectra)
library(dplyr)
library(S4Vectors)
library(IRanges)
library(tidyverse)

neutral_loss <- function(mz_values, precursorMz) {
  round(precursorMz - mz_values)
}

round_perl <- function(number) {
  return(floor(number + 0.5))
}

remove_duplicates <- function(mz, intensity) {
  mz <- round_perl(mz)
  unique_mz <- unique(mz)
  max_intensity <- sapply(unique_mz, function(m) max(intensity[mz == m], na.rm = TRUE))
  return(list(mz = unique_mz, intensity = max_intensity))
}

symmetric_dotproduct_combined <- function(x, y, xPrecursorMz, yPrecursorMz, n = 3, m = 0.6, ...) {
  
  mz1 <- as.numeric(unlist(x[, 1L]))
  intensity1 <- as.numeric(unlist(x[, 2L]))
  mz2 <- as.numeric(unlist(y[, 1L]))
  intensity2 <- as.numeric(unlist(y[, 2L]))
  
  if (length(mz1) == 0 || length(mz2) == 0) 
    return(list(common_peaks = 0, similarity_score = 0))
  
  cleaned1 <- remove_duplicates(mz1, intensity1)
  cleaned2 <- remove_duplicates(mz2, intensity2)
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  
  Fd_denom1 <- sum(intensity1^2)
  Fd_denom2 <- sum(intensity2^2)
  
  common_mz <- intersect(mz1, mz2)
  Fd_nom <- 0
  remaining_idx1 <- rep(TRUE, length(mz1))
  remaining_idx2 <- rep(TRUE, length(mz2))
  
  if (length(common_mz) > 0) {
    idx1 <- match(common_mz, mz1)
    idx2 <- match(common_mz, mz2)
    
    valid_indices <- !is.na(idx1) & !is.na(idx2)
    idx1 <- idx1[valid_indices]
    idx2 <- idx2[valid_indices]
    
    if (length(idx1) > 0 && length(idx2) > 0) {
      Fd_nom <- sum(intensity1[idx1] * intensity2[idx2])
    }
    
    remaining_idx1[idx1] <- FALSE
    remaining_idx2[idx2] <- FALSE
  }
  
  # Calcul des neutral losses sur les pics restants
  nl1 <- neutral_loss(mz1[remaining_idx1], xPrecursorMz)
  nl2 <- neutral_loss(mz2[remaining_idx2], yPrecursorMz)
  
  common_nl <- intersect(nl1, nl2)
  
  if (length(common_nl) > 0) {
    for (nl in common_nl) {
      matches1 <- which(nl1 == nl)
      matches2 <- which(nl2 == nl)
      
      if (length(matches1) > 0 && length(matches2) > 0) {
        Fd_nom <- Fd_nom + sum(intensity1[remaining_idx1][matches1]) * 
          sum(intensity2[remaining_idx2][matches2])
      }
    }
  }
  
  FD <- if (Fd_nom > 0) (Fd_nom^2) / (Fd_denom1 * Fd_denom2) else 0
  
  total_common <- length(common_mz) + length(common_nl)
  
  return(FD)
}

count_common_peaks <- function(x, y, xPrecursorMz, yPrecursorMz, ...) {
  # Extract m/z and intensity values from the spectra
  mz1 <- mz(x)[[1]]  # Extract the numeric vector from the list
  intensity1 <- intensity(x)[[1]]
  mz2 <- mz(y)[[1]]
  intensity2 <- intensity(y)[[1]]
  
  # Check for empty spectra
  if (length(mz1) == 0 || length(mz2) == 0)
    return(0)
  
  # Remove duplicates
  cleaned1 <- remove_duplicates(mz1, intensity1)
  cleaned2 <- remove_duplicates(mz2, intensity2)
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  # Find common m/z values
  common_mz <- intersect(mz1, mz2)
  remaining_idx1 <- rep(TRUE, length(mz1))
  remaining_idx2 <- rep(TRUE, length(mz2))
  
  if (length(common_mz) > 0) {
    idx1 <- match(common_mz, mz1)
    idx2 <- match(common_mz, mz2)
    
    valid_indices <- !is.na(idx1) & !is.na(idx2)
    idx1 <- idx1[valid_indices]
    idx2 <- idx2[valid_indices]
    
    remaining_idx1[idx1] <- FALSE
    remaining_idx2[idx2] <- FALSE
  }
  
  common_peaks <- length(common_mz)
  
  # Calculate neutral losses only on remaining peaks
  nl1 <- neutral_loss(mz1[remaining_idx1], xPrecursorMz)
  nl2 <- neutral_loss(mz2[remaining_idx2], yPrecursorMz)
  
  common_nl <- intersect(nl1, nl2)
  
  # Add common neutral losses to common peaks count
  common_peaks <- common_peaks + length(common_nl)
  
  return(common_peaks)
}





compare_and_add_spectra_all_levels <- function(st_sps, filtered_dy_sps, threshold = 0.8) {
  ensure_best_match_column <- function(spectra_obj) {
    metadata <- spectraData(spectra_obj)
    if (!"best_match" %in% colnames(metadata)) {
      metadata$best_match <- NA_character_
      spectraData(spectra_obj) <- metadata
    }
    return(spectra_obj)
  }
  
  st_sps <- ensure_best_match_column(st_sps)
  filtered_dy_sps <- ensure_best_match_column(filtered_dy_sps)
  
  # Filtrage des spectres MS selon les critères de machine et polarité
  filtered_dy_sps_ms <- filtered_dy_sps[grepl("FTMS", spectraData(filtered_dy_sps)$MACHINE) & 
                                          spectraData(filtered_dy_sps)$polarity == 0]
  
  metadata_filtered <- spectraData(filtered_dy_sps_ms)
  
  # Nettoyage et conversion de l'encodage des noms
  clean_names <- function(names) {
    valid_names <- iconv(names, from = "latin1", to = "UTF-8", sub = "byte")  # ou "ASCII" si nécessaire
    return(valid_names)
  }
  
  metadata_filtered$name <- clean_names(metadata_filtered$name)
  
  valid_names <- !is.na(metadata_filtered$name) & 
    metadata_filtered$name != "NULL" & 
    metadata_filtered$name != "" & 
    !grepl("^!", metadata_filtered$name)
  
  filtered_dy_sps_ms <- filtered_dy_sps_ms[valid_names]
  
  
  if (length(st_sps) > 0 && length(filtered_dy_sps_ms) > 0) {
    similarity_matrix <- compareSpectra(
      st_sps,
      filtered_dy_sps_ms,
      MAPFUN = joinPeaksNone,
      FUN = function(x, y, xPrecursorMz, yPrecursorMz, ...) {
        symmetric_dotproduct_combined(x, y, xPrecursorMz, yPrecursorMz)
      }
    )
    similarity_matrix <- as.matrix(similarity_matrix)
    print("Dimensions de la matrice de similarité :")
    print(dim(similarity_matrix))
  } else {
    stop("Aucun spectre trouvé dans les ensembles de données.")
  }
  
  metadata <- spectraData(st_sps)
  total_spectra <- length(st_sps)
  
  for (i in seq_along(st_sps)) {
    cmp_precursorMz <- precursorMz(st_sps)[i]
    cmp_mslevel <- msLevel(st_sps)[i]
    
    if (cmp_mslevel %in% c(3, 4)) {
      cmp_precursorMz <- round_perl(cmp_precursorMz)
    }
    
    best_match <- NA_character_
    
    if (length(filtered_dy_sps_ms) > 0) {
      best_score <- 0
      best_common_peaks <- 0
      best_name <- NA_character_
      best_compid <- NA_character_
      
      for (j in seq_along(filtered_dy_sps_ms)) {
        ref_precursorMz <- precursorMz(filtered_dy_sps_ms)[j]
        ref_name <- spectraData(filtered_dy_sps_ms)$name[j]
        ref_compid <- spectraData(filtered_dy_sps_ms)$COMPID[j]
        
        if (is.na(ref_name) || ref_name == "NULL" || ref_name == "" || grepl("^!", ref_name)) {
          next
        }
        
        ref_mslevel <- msLevel(filtered_dy_sps_ms)[j]
        if (ref_mslevel %in% c(3, 4)) {
          ref_precursorMz <- round_perl(ref_precursorMz)
        }
        
        if (is.na(cmp_precursorMz) || is.na(ref_precursorMz)) next
        if (abs(cmp_precursorMz - ref_precursorMz) >= 0.005) next  
        
        similarity_score <- similarity_matrix[i, j]
        
        if (similarity_score >= threshold) {
          best_score <- similarity_score
          best_name <- ref_name
          best_compid <- ifelse(!is.na(ref_compid), ref_compid, "NA")
          
          best_common_peaks <- count_common_peaks(
            st_sps[i], filtered_dy_sps_ms[j],
            cmp_precursorMz, ref_precursorMz
          )
        }
      }
      
      if (!is.na(best_name) && best_name != "NULL") {
        best_match <- paste("!", best_common_peaks, "!", best_score, "!", best_name, "!", best_compid, "!", cmp_mslevel, sep=" ")
      }
    }
    
    if (is.na(best_match)) {
      best_match <- NA_character_
    }
    
    metadata$best_match[i] <- best_match
    cat(sprintf("Traitement du spectre %d sur %d\n", i, total_spectra))
  }
  
  spectraData(st_sps) <- metadata
  expid_filtered <- spectraData(filtered_dy_sps)$ExpID
  expid_st <- spectraData(st_sps)$ExpID
  
  max_expid_filtered <- suppressWarnings(max(as.numeric(expid_filtered), na.rm = TRUE))
  if (is.infinite(max_expid_filtered)) max_expid_filtered <- 0
  
  unique_old_ids <- unique(expid_st)
  offset_ids <- seq_along(unique_old_ids) + max_expid_filtered
  id_mapping <- setNames(offset_ids, unique_old_ids)
  
  new_expid_st <- as.character(id_mapping[as.character(expid_st)])
  spectraData(st_sps)$ExpID <- new_expid_st
  
 
  combined_spectra <- c(filtered_dy_sps, st_sps)
  
  return(combined_spectra)

}



results_all_levels <- compare_and_add_spectra_all_levels(st_sps, dy_sps, threshold = 0.8)




saveRDS(results_all_levels, file = "C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/similarity_all_levels_adjusted.rds")

adj <- readRDS("C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/similarity_all_levels_adjusted.rds")
spectraData(adj)
bm<-adj[!is.na(adj$best_match), ]
print(bm$best_match)

cat("Number of matches found:", sum(!is.na(adj$best_match)), "\n")

#rounding scores
# Fonction pour arrondir les scores dans les chaînes best_match
arrondir_scores <- function(x) {
  # Split par "!" sans espace, puis nettoyer les espaces autour
  parts <- trimws(unlist(strsplit(x, "!")))
  parts <- parts[parts != ""]  # Retirer les éléments vides
  
  if (length(parts) >= 3) {
    score <- as.numeric(parts[2])
    parts[2] <- sprintf("%.3f", round(score, 3))
    return(paste("!", paste(parts, collapse = " ! "), sep = " "))
  } else {
    return(x)
  }
}



# Appliquer à tous les best_match dans l'objet adj
adj$best_match <- sapply(adj$best_match, arrondir_scores)


#adjusting Expnum where we have best match we increment






spectra <- sim
mz_values <- mz(spectra)
intensity_values <- intensity(spectra)

spectra_data_df <- as.data.frame(spectraData(sim))


combined_data <- data.frame()

for (i in 1:length(mz_values)) {
  metadata_row <- spectra_data_df[i, ]
  
  mz_vals <- mz_values[[i]]
  intensity_vals <- intensity_values[[i]]
  
  metadata_repeated <- metadata_row[rep(1, length(mz_vals)), ]
  
  temp_df <- data.frame(
    mz = mz_vals,
    intensity = intensity_vals
  )
  
  temp_df <- cbind(metadata_repeated, temp_df)
  
  combined_data <- rbind(combined_data, temp_df)
}

head(combined_data)








n_existing <- length(dy_sps)
n_new <- length(st_sps)
combined_metadata <- spectraData(sim)


n_existing <- length(dy_sps)
n_new <- length(st_sps)


combined_metadata$ExpID <- NA_integer_
combined_metadata$ExpID[1:n_existing] <- 1:n_existing

combined_metadata$ExpID[(n_existing + 1):(n_existing + n_new)] <- (n_existing + 1):(n_existing + n_new)

spectraData(sim) <- combined_metadata









sim <- readRDS("C:/Users/amentag/Desktop/these/These/Ahlam/XCMS/flaxseed/S_all2/similarity_all_levels.rds")
spectraData(sim)
bm<-sim[!is.na(sim$best_match), ]
print(bm$best_match)

cat("Number of matches found:", sum(!is.na(sim$best_match)), "\n")
