
library(MSnbase)
library(MsExperiment)
library(Spectra)
library(dplyr)
library(tidyr)
library(BiocParallel)



register(SerialParam())


filtered_trees <- function(db){
  
  feat.block <- split(db, f = db$feature_name)
  selection <- function(feature.block) {
    
    max.number <- max(table(feature.block$MSntreeID))
    FTmax <- names(table(feature.block$MSntreeID))[table(feature.block$MSntreeID) == max.number]
    feat.block.sel1 <- feature.block[feature.block$MSntreeID %in% as.numeric(FTmax),]
    
    
    feat.block.sel1.MS2 <- feat.block.sel1[feat.block.sel1$msLevel == 2,]
    index <- order(feat.block.sel1.MS2$precursorIntensity, decreasing = TRUE)[1]
    selected.FT <- feat.block.sel1.MS2$MSntreeID[index]
    feature.block[feature.block$MSntreeID == selected.FT,]
  }
  
  
  selected.trees.list <- lapply(feat.block, FUN = selection)
  
  selected.trees <- do.call(rbind, selected.trees.list)
  
  selected_trees_result <- selected.trees[, c("feature_name", "MSntreeID")]
  
  db <- db %>%
    filter(MSntreeID %in% selected_trees_result$MSntreeID)
  
  cat("Base de données filtrée avec succès!\n")
  return(db)
}


global_db <- readRDS("D:/These/Ahlam/XCMS/flaxseed/S_all2/filtered_db.rds")

selected_trees <- filtered_db(global_db)

saveRDS(selected.trees, "D:/These/Ahlam/XCMS/flaxseed/S_all2/selected_trees.rds")
selected_trees <- readRDS("D:/These/Ahlam/XCMS/flaxseed/S_all2/selected_trees.rds")
selected_trees








