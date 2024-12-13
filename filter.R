
library(MSnbase)
library(MsExperiment)
library(xcms)
data_path <- "D:/These/Ahlam/XCMS/flaxseed/S19"
fls <- dir(data_path, pattern = "\\.mzXML$", full.names = TRUE, recursive = TRUE)
pd <- data.frame(
  file = basename(fls),
  sample = c("S1", "S9", "S3"),
  injection_index = c(1, 2, 3),
  group = rep("flaxseed", 3)
)
mse <- readMsExperiment(fls, sampleData = pd)

spectra_data <- spectra(mse)
spectra_columns <- spectraVariables(spectra_data)

spectra_info <- data.frame(
  lapply(spectra_columns, function(col) {
    tryCatch({
      spectra_col_data <- get(col)(spectra_data)
      return(spectra_col_data)
    }, error = function(e) {
      message(paste("Erreur avec la colonne", col, ":", e$message))
      return(NA) 
    })
  })
)

names(spectra_info) <- spectra_columns
head(spectra_info)
write.table(
  spectra_info,
  file = "D:/These/Ahlam/XCMS/flaxseed/S19/spectra_data_comp3.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)




cwp <- CentWaveParam(ppm = 25, peakwidth = c(5, 20), snthresh = 10)
mse_find <- findChromPeaks(mse, param = cwp)
sample_info <- sampleData(mse_find)
#sample_info <- phenoData(mse_find)
pdp <- PeakDensityParam(
  sampleGroups = sample_info$group,
  bw = 10,
  minFraction = 0.5
)
mse_grp1 <- groupChromPeaks(mse_find, param = pdp)

pyp <- PeakGroupsParam(
  minFraction = 0.9,
  extraPeaks = 1,
  smooth = "loess",
  span = 0.2,
  family = "gaussian",
  peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
  subset = integer(),
  subsetAdjust = c("average", "previous")
)
mse_rt_corr <- adjustRtime(mse_grp1,param=pyp)
pdp_corr <- PeakDensityParam(
  sampleGroups = sampleData(mse_rt_corr)$group,  
  bw = 10,                                      
  minFraction = 0.5                            
)
mse_grp_corr <- groupChromPeaks(mse_rt_corr, param = pdp_corr)
corrected_features <- featureDefinitions(mse_grp_corr)
plotAdjustedRtime(mse_rt_corr, adjustedRtime = F)
fill_param <- FillChromPeaksParam()
xdata <- as(mse_grp_corr, "XCMSnExp")
xdata_filled <- fillChromPeaks(xdata, param = fill_param)
rt_raw <- rtime(xdata_filled, adjusted = FALSE)
feature_matrix <- featureValues(
  xdata_filled,
  value = "into",
  filled = TRUE
)

corrected_features <- featureDefinitions(xdata_filled)
final_feature_matrix <- data.frame(
  `rt(min)` = corrected_features$rtmed / 60,
  mzmed = corrected_features$mzmed,
  name = rownames(corrected_features),
  feature_matrix
)


rt_seconds <- final_feature_matrix$'rt.min.'*60

final_feature_matrix$name <- paste0("M",round(final_feature_matrix$mzmed,0),"T",round(final_feature_matrix$'rt.min.'*60,0))
head(final_feature_matrix)
write.table(
  final_feature_matrix,
  file = "D:/These/Ahlam/XCMS/flaxseed/S19/corrected_feature_matrix3_comp.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
head(final_feature_matrix)

spectra_data  <- read.table("D:/These/Ahlam/XCMS/flaxseed/S19/spectra_data_comp3.txt", sep = "\t", header = TRUE)
feature_matrix <- read.table("D:/These/Ahlam/XCMS/flaxseed/S19/corrected_feature_matrix3_comp.txt", sep = "\t", header = TRUE)





#new
library(dplyr)
library(MSnbase)
library(MsExperiment)
library(xcms)
library(MSnbase)
library(MsCoreUtils)
library(Spectra)



str(feature_matrix)


feature_matrix$rt.min. <- as.numeric(feature_matrix$rt.min.)
feature_matrix$mzmed <- as.numeric(feature_matrix$mzmed)

rt_min_values <- feature_matrix$rt.min. * 60  
mz_values <- feature_matrix$mzmed

filtered_mse <- mse

do_filter <- function(rt_min, mz, data) {
  rt_range <- c(rt_min - 1 , rt_min +1 )  
  mz_range <- c(mz- 0.0001 , mz+ 0.0001 )      
  
  filtered <- filterMz(data, mz = mz_range)
  filtered <- filterRt(filtered, rt = rt_range)
  
  return(filtered)
}


filtered_mse_list <- lapply(
  1:length(rt_min_values),
  function(i) do_filter(rt_min_values[i], mz_values[i], filtered_mse)
)


spectra_data <- lapply(filtered_mse_list, spectra)


str(spectra_data)

spectra_df <- do.call(rbind, lapply(spectra_data, spectraData))



str(spectra_df)
head(spectra_df)

nrow(spectra_df)

nrow(spectra_info)   


#new


library(dplyr)





spectra_df_clean <- spectra_df[!is.na(spectra_df$precursorIntensity) & !is.na(spectra_df$msLevel), ]


filtered_spectra_info <- data.frame()


unique_spectrum_ids <- unique(spectra_df_clean$spectrumId)

for (spectrum_id in unique_spectrum_ids) {
  
  spectrum_candidates <- spectra_df_clean[spectra_df_clean$spectrumId == spectrum_id, ]
  
  
  if (nrow(spectrum_candidates) > 0) {
    best_spectrum <- spectrum_candidates[which.max(spectrum_candidates$precursorIntensity), , drop = TRUE]
    filtered_spectra_info <- rbind(filtered_spectra_info, best_spectrum)
  }
}


if (nrow(filtered_spectra_info) == 0) {
  cat("Aucune correspondance trouvée avec les critères donnés.\n")
} else {
  print(head(filtered_spectra_info))
}


write.table(
  filtered_spectra_info,
  file = "filtered_spectra_info_alternative.txt",  
  sep = "\t",  
  row.names = FALSE,  
  quote = FALSE  
)



head(spectra_filtered)


nrow(spectra_df)
nrow(filtered_spectra_info)
nrow(spectra_info)


