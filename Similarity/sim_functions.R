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
symmetric_dotproduct_combined <- function(x, y, xPrecursorMz,
                                          yPrecursorMz, n = 3, m = 0.6, ...) {
  if (!nrow(x) || !nrow(y))
    return(list(common_peaks = 0, similarity_score = 0))

  cleaned1 <- remove_duplicates(x[, 1L], x[, 2L])
  cleaned2 <- remove_duplicates(z[, 1L], y[, 2L])

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


#' @description
#'
#' Matches/maps peaks between the two provided peak matrices.
#'
#' @return
#'
#' a `list` with `numeric` matrices containing only the matched peaks between
#' the two input spectra.
#'
#' @author Ahlam Mentag
dynlib_map <- function(x, y, xPrecursorMz, yPrecursorMz, ...) {
    if (!nrow(x) || !nrow(y))
        return(list(x = matrix(numeric(), ncol = 3, nrow = 0),
                    y = matrix(numeric(), ncol = 3, nrow = 0)))

    cleaned1 <- remove_duplicates(x[, 1L], x[, 2L])
    cleaned2 <- remove_duplicates(z[, 1L], y[, 2L])

    mz1 <- cleaned1$mz
    mz2 <- cleaned2$mz
    ## Continue from here - keep only the rows in `cleaned1` or `cleaned2` that
    ## have matching m/z between the two.
    common_mz <- intersect(mz1, mz2)
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

dynlib_fun <- function() {
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

compare_with_metaboannotation <- function(st_sps, filtered_dy_sps, threshold = 0.8, ppm = 5, tolerance = 0.005) {
  library(MetaboAnnotation)

  param <- CompareSpectraParam(
    ppm = ppm,
    tolerance = tolerance,
    requirePrecursor = TRUE,
    FUN = symmetric_dotproduct_combined,
    threshold = threshold
  )

  matches <- matchSpectra(
    query = st_sps,
    target = filtered_dy_sps,
    param = param
  )

  df <- as.data.frame(matches)
  df <- df[df$score >= threshold, ]

  return(df)
}
