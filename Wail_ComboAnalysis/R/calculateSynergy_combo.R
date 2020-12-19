require(PharmacoGx)
require(magicaxis)
library(abind)
library(robustbase)
library(Biobase)
library(synergyfinder)

source("R/synergySurf.R")
source("R/fixCellLinesNames.R")
source("R/SynVolume.r")
source("R/AnalyzeCombo_as_function.R")


listOfCombos <- list()
listOfCombos_stat <- list()

############################################################
############################################################
############################################################

# parse sensitivity raw files

# DP_Fluva
Combo <- "FLUVA_DP"


x1 <- c(0, 1.25, 2.5, 5, 10, 20)
x2 <- c(0.0,	.078,	.156,	.3125,	.625,	1.25,	2.5,	5,	10,	20)

files.ls <- dir("data/rawData_combo/FluvaDP/Processed")
files.ls <- files.ls[grep("txt",files.ls)]

files.ls <- paste("data/rawData_combo/FluvaDP/Processed/",files.ls,sep = "")
sensitivityData <- getSensitivityComboMatrix(viabilityFiles = files.ls,x1 = x1,x2 =x2)

rownames(sensitivityData) <- gsub(" DP/FLUVA","",rownames(sensitivityData))

sensitivityData <- sensitivityData[-which(rownames(sensitivityData)=="HCT116" | rownames(sensitivityData)=="HCT116_1" | rownames(sensitivityData)=="HCT116_2" | rownames(sensitivityData)=="HCT116_3"
                                          | rownames(sensitivityData)=="DU145" | rownames(sensitivityData)=="DU145_1" | rownames(sensitivityData)=="DU145_2" | rownames(sensitivityData)=="DU145_3"),,]


sensitivityData <- sensitivityData[order(rownames(sensitivityData)),,]

#saveRDS(sensitivityData,"sensitivityData_FLUVA_DP_23Sep19.rda")



# DP_Fluva
Combo <- "FLUVA_DP"

#sensitivityData <- readRDS(paste("sensitivityData_",Combo,"_23Sep19.rda",sep = ""))

# c("PairIndex","cell_line","Drug1","Drug2","Conc1","Conc2","Response","ConcUnit")

DrugR <- strsplit(Combo,split = "_")[[1]][2]
DrugC <- strsplit(Combo,split = "_")[[1]][1]

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

BlissList <- list()
LoeweList <- list()
HSAList <- list()
ZIPList <- list()

for (i in seq_along(1:dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  mat <- sensitivityData[i,,]*100
  mat <- mat[order(as.numeric(rownames(mat))),]
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  meta <- data.frame(drug.col = DrugC, drug.row = DrugR, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1"=mat), drug.pairs = meta)
  
  bliss_score <- NA; bliss_mat <- NA; bliss <- NA
  loewe_score <- NA; loewe_mat <- NA; loewe <- NA
  hsa_score <- NA; hsa_mat <- NA; hsa <- NA
  zip_score <- NA; zip_mat <- NA; zip <- NA
  
  bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss",adjusted = F)
  }, error = function(err) { print("Bliss error") ; bliss <- NA
  })
  if(!is.na(bliss_exe)){
    bliss_score <- round(median(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    bliss_mat <- bliss$scores[[1]]/100.0
    BlissList[[i]] <- bliss_exe
  } else {
    bliss_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(bliss_mat) <- c(combo_row_concentration)
    colnames(bliss_mat) <- c(combo_col_concentration)
    BlissList[[i]] <- NULL
  }
  
  loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe",adjusted = F)
  }, error = function(err) { print("Loewe error") ; loewe <- NA
  })
  if(!is.na(loewe_exe)){
    loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    loewe_mat <- loewe$scores[[1]]/100.0
    LoeweList[[i]] <- loewe_exe
  } else {
    loewe_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(loewe_mat) <- c(combo_row_concentration)
    colnames(loewe_mat) <- c(combo_col_concentration)
    LoeweList[[i]] <- NULL
  }
  
  hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA",adjusted = F)
  }, error = function(err) { print("HSA error"); hsa <- NA
  })
  if(!is.na(hsa_exe)){
    hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    hsa_mat <- hsa$scores[[1]]/100.0
    HSAList[[i]] <- hsa_exe
  } else {
    hsa_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(hsa_mat) <- c(combo_row_concentration)
    colnames(hsa_mat) <- c(combo_col_concentration)
    HSAList[[i]] <- NULL
  } 
  
  zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP",adjusted = F) 
  }, error = function(err) { print("ZIP error"); zip <- NA 
  })
  if(!is.na(zip_exe)){
    zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)			
    zip_mat <- zip$scores[[1]]/100.0
    ZIPList[[i]] <- zip_exe
  } else {
    zip_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(zip_mat) <- c(combo_row_concentration)
    colnames(zip_mat) <- c(combo_col_concentration)
    ZIPList[[i]] <- NULL
  } 
  
  result <- c(cellID, DrugR, DrugC, bliss_score, loewe_score, hsa_score, zip_score)
  
  outfile_stats <- rbind(outfile_stats, result)
  
  
  mat_out <- data.frame(
    name = rep(cellID, length(c(mat))),
    drugA = rep(DrugR, length(c(mat))),
    drugB = rep(DrugC, length(c(mat))),
    concA = as.numeric(rep(colnames(mat), each=nrow(mat))),
    concB = as.numeric(rep(rownames(mat), ncol(mat))),
    experiments = as.numeric(c(mat))
  )
  
  raw_matrix <- rbind(raw_matrix, mat_out)
  
  bliss_out <- data.frame(
    name = rep(cellID, length(c(bliss_mat))),
    drugA = rep(DrugR, length(c(bliss_mat))),
    drugB = rep(DrugC, length(c(bliss_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(bliss_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(bliss_mat))),
    bliss = as.numeric(c(bliss_mat))
  )
  bliss_matrix <- rbind(bliss_matrix, bliss_out)
  
  loewe_out <- data.frame(
    name = rep(cellID, length(c(loewe_mat))),
    drugA = rep(DrugR, length(c(loewe_mat))),
    drugB = rep(DrugC, length(c(loewe_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(loewe_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(loewe_mat))),
    loewe = as.numeric(c(loewe_mat))
  )
  loewe_matrix <- rbind(loewe_matrix, loewe_out)
  
  hsa_out <- data.frame(
    name = rep(cellID, length(c(hsa_mat))),
    drugA = rep(DrugR, length(c(hsa_mat))),
    drugB = rep(DrugC, length(c(hsa_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(hsa_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(hsa_mat))),
    hsa = as.numeric(c(hsa_mat))
  )
  hsa_matrix <- rbind(hsa_matrix, hsa_out)
  
  zip_out <- data.frame(
    name = rep(cellID, length(c(zip_mat))),
    drugA = rep(DrugR, length(c(zip_mat))),
    drugB = rep(DrugC, length(c(zip_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(zip_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(zip_mat))),
    zip = as.numeric(c(zip_mat))
  )
  zip_matrix <- rbind(zip_matrix, zip_out)
}

names(BlissList) <- rownames(sensitivityData)
names(LoeweList) <- rownames(sensitivityData)
names(HSAList) <- rownames(sensitivityData)
names(ZIPList) <- rownames(sensitivityData)

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")

mat_total <- as.matrix(outfile_stats[,c("Bliss","Loewe","HSA","ZIP")])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]

mat_total <- as.data.frame(mat_total,stringsAsFactors=F)


listOfCombos[[Combo]] <- list("Bliss"=BlissList,"Loewe"=LoeweList,"HSA"=HSAList,"ZIP"=ZIPList)
listOfCombos_stat[[Combo]] <- mat_total

############################################################
############################################################
############################################################

# NFV_Fluva
Combo <- "FLUVA_NFV"

x1 <- c(0,0.625, 1.25, 2.5, 5, 10)
x2 <- c(0.0,	.078,	.156,	.3125,	.625,	1.25,	2.5,	5,	10,	20)

files.ls <- dir("data/rawData_combo/FluvaNFV/processed/")
files.ls <- files.ls[grep("txt",files.ls)]

files.ls <- paste("data/rawData_combo/FluvaNFV/processed/",files.ls,sep = "")
sensitivityData <- getSensitivityComboMatrix(viabilityFiles = files.ls,x1 = x1,x2 =x2)

sensitivityData <- sensitivityData[-which(rownames(sensitivityData)=="HCT116" | rownames(sensitivityData)=="HCT116_1" | rownames(sensitivityData)=="HCT116_2" | rownames(sensitivityData)=="HCT116_3"
                                          | rownames(sensitivityData)=="DU145" | rownames(sensitivityData)=="DU145_1" | rownames(sensitivityData)=="DU145_2" | rownames(sensitivityData)=="DU145_3"),,]


sensitivityData <- sensitivityData[order(rownames(sensitivityData)),,]

#saveRDS(sensitivityData,"sensitivityData_FLUVA_NFV_23Sep19.rda")


#sensitivityData <- readRDS(paste("sensitivityData_",Combo,"_23Sep19.rda",sep = ""))

DrugR <- strsplit(Combo,split = "_")[[1]][2]
DrugC <- strsplit(Combo,split = "_")[[1]][1]

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

BlissList <- list()
LoeweList <- list()
HSAList <- list()
ZIPList <- list()

for (i in seq_along(1:dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  mat <- sensitivityData[i,,]*100
  mat <- mat[order(as.numeric(rownames(mat))),]
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  meta <- data.frame(drug.col = DrugC, drug.row = DrugR, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1"=mat), drug.pairs = meta)
  
  bliss_score <- NA; bliss_mat <- NA; bliss <- NA
  loewe_score <- NA; loewe_mat <- NA; loewe <- NA
  hsa_score <- NA; hsa_mat <- NA; hsa <- NA
  zip_score <- NA; zip_mat <- NA; zip <- NA
  
  bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss",adjusted = F)
  }, error = function(err) { print("Bliss error") ; bliss <- NA
  })
  if(!is.na(bliss_exe)){
    bliss_score <- round(median(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    bliss_mat <- bliss$scores[[1]]/100.0
    BlissList[[i]] <- bliss_exe
  } else {
    bliss_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(bliss_mat) <- c(combo_row_concentration)
    colnames(bliss_mat) <- c(combo_col_concentration)
    BlissList[[i]] <- NULL
  }
  
  loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe",adjusted = F)
  }, error = function(err) { print("Loewe error") ; loewe <- NA
  })
  if(!is.na(loewe_exe)){
    loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    loewe_mat <- loewe$scores[[1]]/100.0
    LoeweList[[i]] <- loewe_exe
  } else {
    loewe_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(loewe_mat) <- c(combo_row_concentration)
    colnames(loewe_mat) <- c(combo_col_concentration)
    LoeweList[[i]] <- NULL
  }
  
  hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA",adjusted = F)
  }, error = function(err) { print("HSA error"); hsa <- NA
  })
  if(!is.na(hsa_exe)){
    hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    hsa_mat <- hsa$scores[[1]]/100.0
    HSAList[[i]] <- hsa_exe
  } else {
    hsa_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(hsa_mat) <- c(combo_row_concentration)
    colnames(hsa_mat) <- c(combo_col_concentration)
    HSAList[[i]] <- NULL
  } 
  
  zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP",adjusted = F) 
  }, error = function(err) { print("ZIP error"); zip <- NA 
  })
  if(!is.na(zip_exe)){
    zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)			
    zip_mat <- zip$scores[[1]]/100.0
    ZIPList[[i]] <- zip_exe
  } else {
    zip_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(zip_mat) <- c(combo_row_concentration)
    colnames(zip_mat) <- c(combo_col_concentration)
    ZIPList[[i]] <- NULL
  } 
  
  result <- c(cellID, DrugR, DrugC, bliss_score, loewe_score, hsa_score, zip_score)
  
  outfile_stats <- rbind(outfile_stats, result)
  
  
  mat_out <- data.frame(
    name = rep(cellID, length(c(mat))),
    drugA = rep(DrugR, length(c(mat))),
    drugB = rep(DrugC, length(c(mat))),
    concA = as.numeric(rep(colnames(mat), each=nrow(mat))),
    concB = as.numeric(rep(rownames(mat), ncol(mat))),
    experiments = as.numeric(c(mat))
  )
  
  raw_matrix <- rbind(raw_matrix, mat_out)
  
  bliss_out <- data.frame(
    name = rep(cellID, length(c(bliss_mat))),
    drugA = rep(DrugR, length(c(bliss_mat))),
    drugB = rep(DrugC, length(c(bliss_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(bliss_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(bliss_mat))),
    bliss = as.numeric(c(bliss_mat))
  )
  bliss_matrix <- rbind(bliss_matrix, bliss_out)
  
  loewe_out <- data.frame(
    name = rep(cellID, length(c(loewe_mat))),
    drugA = rep(DrugR, length(c(loewe_mat))),
    drugB = rep(DrugC, length(c(loewe_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(loewe_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(loewe_mat))),
    loewe = as.numeric(c(loewe_mat))
  )
  loewe_matrix <- rbind(loewe_matrix, loewe_out)
  
  hsa_out <- data.frame(
    name = rep(cellID, length(c(hsa_mat))),
    drugA = rep(DrugR, length(c(hsa_mat))),
    drugB = rep(DrugC, length(c(hsa_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(hsa_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(hsa_mat))),
    hsa = as.numeric(c(hsa_mat))
  )
  hsa_matrix <- rbind(hsa_matrix, hsa_out)
  
  zip_out <- data.frame(
    name = rep(cellID, length(c(zip_mat))),
    drugA = rep(DrugR, length(c(zip_mat))),
    drugB = rep(DrugC, length(c(zip_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(zip_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(zip_mat))),
    zip = as.numeric(c(zip_mat))
  )
  zip_matrix <- rbind(zip_matrix, zip_out)
}

names(BlissList) <- rownames(sensitivityData)
names(LoeweList) <- rownames(sensitivityData)
names(HSAList) <- rownames(sensitivityData)
names(ZIPList) <- rownames(sensitivityData)

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")

mat_total <- as.matrix(outfile_stats[,c("Bliss","Loewe","HSA","ZIP")])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]

mat_total <- as.data.frame(mat_total,stringsAsFactors=F)

#hist(mat_total[,"ZIP"])

#PlotSynergy(BlissList$JIMT1_1, type = "all")



listOfCombos[[Combo]] <- list("Bliss"=BlissList,"Loewe"=LoeweList,"HSA"=HSAList,"ZIP"=ZIPList)
listOfCombos_stat[[Combo]] <- mat_total


############################################################
############################################################
############################################################

# HNK_Fluva
Combo <- "FLUVA_HNK"


x1 <- c(0, 1.25, 2.5, 5, 10, 20)
x2 <- c(0.0,	.078,	.156,	.3125,	.625,	1.25,	2.5,	5,	10,	20)

files.ls <- dir("data/rawData_combo/FluvaHNK/processed/")
files.ls <- files.ls[grep("txt",files.ls)]

files.ls <- paste("data/rawData_combo/FluvaHNK/processed/",files.ls,sep = "")
sensitivityData <- getSensitivityComboMatrix(viabilityFiles = files.ls,x1 = x1,x2 =x2)

sensitivityData <- sensitivityData[-which(rownames(sensitivityData)=="HCT116" | rownames(sensitivityData)=="HCT116_1" | rownames(sensitivityData)=="HCT116_2" | rownames(sensitivityData)=="HCT116_3"
                                          | rownames(sensitivityData)=="DU145" | rownames(sensitivityData)=="DU145_2" | rownames(sensitivityData)=="DU145_3"),,]


sensitivityData <- sensitivityData[order(rownames(sensitivityData)),,]

#saveRDS(sensitivityData,paste("sensitivityData_",Combo,"_23Sep19.rda",sep = ""))




#sensitivityData <- readRDS(paste("sensitivityData_",Combo,"_23Sep19.rda",sep = ""))

DrugR <- strsplit(Combo,split = "_")[[1]][2]
DrugC <- strsplit(Combo,split = "_")[[1]][1]

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

BlissList <- list()
LoeweList <- list()
HSAList <- list()
ZIPList <- list()

for (i in seq_along(1:dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  mat <- sensitivityData[i,,]*100
  mat <- mat[order(as.numeric(rownames(mat))),]
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  meta <- data.frame(drug.col = DrugC, drug.row = DrugR, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1"=mat), drug.pairs = meta)
  
  bliss_score <- NA; bliss_mat <- NA; bliss <- NA
  loewe_score <- NA; loewe_mat <- NA; loewe <- NA
  hsa_score <- NA; hsa_mat <- NA; hsa <- NA
  zip_score <- NA; zip_mat <- NA; zip <- NA
  
  bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss",adjusted = F)
  }, error = function(err) { print("Bliss error") ; bliss <- NA
  })
  if(!is.na(bliss_exe)){
    bliss_score <- round(median(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    bliss_mat <- bliss$scores[[1]]/100.0
    BlissList[[i]] <- bliss_exe
  } else {
    bliss_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(bliss_mat) <- c(combo_row_concentration)
    colnames(bliss_mat) <- c(combo_col_concentration)
    BlissList[[i]] <- NULL
  }
  
  loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe",adjusted = F)
  }, error = function(err) { print("Loewe error") ; loewe <- NA
  })
  if(!is.na(loewe_exe)){
    loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    loewe_mat <- loewe$scores[[1]]/100.0
    LoeweList[[i]] <- loewe_exe
  } else {
    loewe_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(loewe_mat) <- c(combo_row_concentration)
    colnames(loewe_mat) <- c(combo_col_concentration)
    LoeweList[[i]] <- NULL
  }
  
  hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA",adjusted = F)
  }, error = function(err) { print("HSA error"); hsa <- NA
  })
  if(!is.na(hsa_exe)){
    hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
    hsa_mat <- hsa$scores[[1]]/100.0
    HSAList[[i]] <- hsa_exe
  } else {
    hsa_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(hsa_mat) <- c(combo_row_concentration)
    colnames(hsa_mat) <- c(combo_col_concentration)
    HSAList[[i]] <- NULL
  } 
  
  zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP",adjusted = F) 
  }, error = function(err) { print("ZIP error"); zip <- NA 
  })
  if(!is.na(zip_exe)){
    zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)			
    zip_mat <- zip$scores[[1]]/100.0
    ZIPList[[i]] <- zip_exe
  } else {
    zip_mat <- matrix(rep(NA, length(c(mat))), ncol=ncol(mat), nrow=nrow(mat))
    rownames(zip_mat) <- c(combo_row_concentration)
    colnames(zip_mat) <- c(combo_col_concentration)
    ZIPList[[i]] <- NULL
  } 
  
  result <- c(cellID, DrugR, DrugC, bliss_score, loewe_score, hsa_score, zip_score)
  
  outfile_stats <- rbind(outfile_stats, result)
  
  
  mat_out <- data.frame(
    name = rep(cellID, length(c(mat))),
    drugA = rep(DrugR, length(c(mat))),
    drugB = rep(DrugC, length(c(mat))),
    concA = as.numeric(rep(colnames(mat), each=nrow(mat))),
    concB = as.numeric(rep(rownames(mat), ncol(mat))),
    experiments = as.numeric(c(mat))
  )
  
  raw_matrix <- rbind(raw_matrix, mat_out)
  
  bliss_out <- data.frame(
    name = rep(cellID, length(c(bliss_mat))),
    drugA = rep(DrugR, length(c(bliss_mat))),
    drugB = rep(DrugC, length(c(bliss_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(bliss_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(bliss_mat))),
    bliss = as.numeric(c(bliss_mat))
  )
  bliss_matrix <- rbind(bliss_matrix, bliss_out)
  
  loewe_out <- data.frame(
    name = rep(cellID, length(c(loewe_mat))),
    drugA = rep(DrugR, length(c(loewe_mat))),
    drugB = rep(DrugC, length(c(loewe_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(loewe_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(loewe_mat))),
    loewe = as.numeric(c(loewe_mat))
  )
  loewe_matrix <- rbind(loewe_matrix, loewe_out)
  
  hsa_out <- data.frame(
    name = rep(cellID, length(c(hsa_mat))),
    drugA = rep(DrugR, length(c(hsa_mat))),
    drugB = rep(DrugC, length(c(hsa_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(hsa_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(hsa_mat))),
    hsa = as.numeric(c(hsa_mat))
  )
  hsa_matrix <- rbind(hsa_matrix, hsa_out)
  
  zip_out <- data.frame(
    name = rep(cellID, length(c(zip_mat))),
    drugA = rep(DrugR, length(c(zip_mat))),
    drugB = rep(DrugC, length(c(zip_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(zip_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(zip_mat))),
    zip = as.numeric(c(zip_mat))
  )
  zip_matrix <- rbind(zip_matrix, zip_out)
}

names(BlissList) <- rownames(sensitivityData)
names(LoeweList) <- rownames(sensitivityData)
names(HSAList) <- rownames(sensitivityData)
names(ZIPList) <- rownames(sensitivityData)

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")

mat_total <- as.matrix(outfile_stats[,c("Bliss","Loewe","HSA","ZIP")])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]

mat_total <- as.data.frame(mat_total,stringsAsFactors=F)



listOfCombos[[Combo]] <- list("Bliss"=BlissList,"Loewe"=LoeweList,"HSA"=HSAList,"ZIP"=ZIPList)
listOfCombos_stat[[Combo]] <- mat_total


save(listOfCombos,listOfCombos_stat,file = "data/SynergyStat_Statin_combo_23Sep19.RData")
