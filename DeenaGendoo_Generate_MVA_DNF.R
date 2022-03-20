# Deena M.A. Gendoo
# Created on January 10, 2017
# Generate the MVA-DNF matrix, using set of MVA-defined genes (genes identified by Dr. Linda Penn's lab)
########################################################################
########################################################################
library(PharmacoGx) 
library(apcluster)
library(rcdk)
library(fingerprint)
library(annotate)
library(org.Hs.eg.db)
library(SNFtool)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

mva<-read.csv("Data/MVA_genes.csv")
mvagenes<-as.character(sort(mva$Gene.name))

source('FunctionsBank.R')
load("Data/PrePros.RData")
pertData<-pertData[(which(rownames(pertData) %in% mva$Gene.name)),] 
dim(pertData)
rownames(pertData)

# Reduce All Matrices to lowest common set of drugs 
# Get 237 drugs now in the reduced sets
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                     sort(colnames(pertData))))
strcData<- strcData[commonDrugs]
sensData <- sensData[,commonDrugs]
pertData <- pertData[,commonDrugs]

#Sanity Checks
if (ncol(sensData) != ncol(pertData)) stop(sprintf("error!"))
if (ncol(sensData) !=length(strcData)) stop(sprintf("error!"))
if (all(colnames(pertData) != colnames(sensData))) stop(sprintf("error!"))
if (all(colnames(pertData) != names(strcData))) stop(sprintf("error!"))

# Correlations of Drugs Closest to DP, individual layers
  InterestDrug<-"DIPYRIDAMOLE"
  # Correlation for Sensivity (Pearson)
  sensCor <- cor(sensData, method = "pearson", use = "pairwise.complete.obs")
  sensCor <- apply(sensCor, 1, function(x) ifelse(is.na(x),0,x))
  sensCorDP<-sensCor[InterestDrug,]
  # Correlation for Perturbation (Perason)
  pertCor <- cor(pertData, method = "pearson", use = "pairwise.complete.obs")
  pertCorDP<-pertCor[InterestDrug,]
  ## Correlation for Structure (Tanimoto metric)
  fpSim <- fingerprint::fp.sim.matrix(strcData, method = "tanimoto")
  rownames(fpSim) <- names(strcData)
  colnames(fpSim) <- names(strcData)
  strcCorDP<-fpSim[InterestDrug,]

## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)
integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)

# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)

# Save an RData Object with all matrices: MVA-DNF and single layer taxonomies
save(integrtStrctSensPert,strcAffMat,sensAffMat,pertAffMat, file="MVA_DNF.RData")







