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

#' @description
#'
#' Function to calculate the symmetric dot product between two spectra. The
#' function also:
#'
#' - *cleans* the spectra keeping only a single fragment peak for groups of
#'   peaks with a m/z smaller 0.5
#'
#' @param x `numeric` 2-column `matrix` with the m/z and intensity values
#'
#' @param y `numeric` 2-column `matrix` with the m/z and intensity values
#'
#' @param xPrecursorMz `numeric(1)` with the precursor m/z of the spectrum `x`
#'
#' @param yPrecursorMz `numeric(1)` with the precursor m/z of the spectrum `y`
#'
#' @param n
#'
#' @param m
#'
#' @author Ahlam Mentag
symmetric_dotproduct_combined <- function(x, y, n = 3, m = 0.6, ...) {
  if (!nrow(x) || !nrow(y))
    return(list(common_peaks = 0, similarity_score = 0))
  
  mz1 <- x[, 1L]
  intensity1 <- x[, 2L]
  mz2 <- y[, 1L]
  intensity2 <- y[, 2L]
  
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  
  denom1 <- sum(intensity1^2)
  denom2 <- sum(intensity2^2)
  
  common_mz <- intersect(mz1, mz2)
  
  nom <- 0
  if (length(common_mz) > 0) {
    idx1 <- match(common_mz, mz1)
    idx2 <- match(common_mz, mz2)
    nom <- sum(intensity1[idx1] * intensity2[idx2])
  }
  
  similarity <- if (nom > 0) (nom^2) / (denom1 * denom2) else 0
  
  return(list(common_peaks = length(common_mz), similarity_score = similarity))
}



#' @description
#'
#' Matches/maps peaks between the two provided peak matrices.
#'
#' @return
#'
#' a list with numeric matrices containing only the matched peaks between
#' the two input spectra.
#'
#' @author Ahlam Mentag
dynlib_map <- function(x, y, xPrecursorMz, yPrecursorMz, ...) {
  if (!nrow(x) || !nrow(y))
    return(list(x = matrix(numeric(), ncol = 2, nrow = 0),
                y = matrix(numeric(), ncol = 2, nrow = 0)))
  
  cleaned1 <- remove_duplicates(x[, 1L], x[, 2L])
  cleaned2 <- remove_duplicates(y[, 1L], y[, 2L])
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  # Matching exact m/z
  common_mz <- intersect(mz1, mz2)
  matched1 <- cbind(mz1[mz1 %in% common_mz], intensity1[mz1 %in% common_mz])
  matched2 <- cbind(mz2[mz2 %in% common_mz], intensity2[mz2 %in% common_mz])
  
  # Remaining peaks
  remaining_idx1 <- !(mz1 %in% common_mz)
  remaining_idx2 <- !(mz2 %in% common_mz)
  
  nl1 <- neutral_loss(mz1[remaining_idx1], xPrecursorMz)
  nl2 <- neutral_loss(mz2[remaining_idx2], yPrecursorMz)
  common_nl <- intersect(nl1, nl2)
  
  for (nl in common_nl) {
    idx_nl1 <- which(nl1 == nl)
    idx_nl2 <- which(nl2 == nl)
    
    matched_nl1 <- cbind(mz1[remaining_idx1][idx_nl1],
                         intensity1[remaining_idx1][idx_nl1])
    matched_nl2 <- cbind(mz2[remaining_idx2][idx_nl2],
                         intensity2[remaining_idx2][idx_nl2])
    
    matched1 <- rbind(matched1, matched_nl1)
    matched2 <- rbind(matched2, matched_nl2)
  }
  
  return(list(x = matched1, y = matched2))
}




dynlib_fun <- function(x, y, ...) {
  
  peaks_x <- peaks(x)
  peaks_y <- peaks(y)
  

  precursor_x <- precursorMz(x)
  precursor_y <- precursorMz(y)
  

  matched <- dynlib_map(peaks_x, peaks_y, precursor_x, precursor_y)
  
 
  result <- symmetric_dotproduct_combined(matched$x, matched$y)
  
 
  return(result$similarity_score)
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





library(BiocParallel)
library(Spectra)

final_sim <- function(st_sps, filtered_dy_sps, threshold = 0.8, ppm = 5, tolerance = 0.005) {
  library(MetaboAnnotation)
  
  param <- CompareSpectraParam(
    ppm = ppm,
    tolerance = tolerance,
    requirePrecursor = TRUE,
    FUN = dynlib_fun,
    threshold = threshold
  )
  
  matches <- matchSpectra(
    query = st_sps,
    target = filtered_dy_sps,
    param = param
  )
  
  df <- matchedData(matches)
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  
  return(df)
}
