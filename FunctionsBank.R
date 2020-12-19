# LIST OF FUNCTIONS
# Functions developed by Gendoo, Ghoraie, and El-Hachem for DNF publication (Cancer Research 2017)

##### constPert #####

constPerturbationLayer <- function(pertDat) {
  
  # Correlation for Perturbation
  pertCor <- cor(pertDat, method = "pearson", use = "pairwise.complete.obs")
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  pertAff <- SNFtool::affinityMatrix(1-pertCor, 20, 0.5)
  
  return(pertAff)
  
}

##### constSens #####

constSensitivityLayer <- function(sensDat) {
  
  # Correlation for Sensivity
  sensCor <- cor(sensDat, method = "pearson", use = "pairwise.complete.obs")
  ## if NA remaining in cor matrix, replace with 0s, not very clean but no other choices for now
  sensCor <- apply(sensCor, 1, function(x) ifelse(is.na(x),0,x))
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  sensAff <- SNFtool::affinityMatrix(1-sensCor, 20, 0.5)
  
  return(sensAff)
  
}

##### Merger  #####

integrateStrctSensPert <- function(sensAff, strcAff, pertAff) {
  
  integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff))
  colnames(integration) <- rownames(integration) <- colnames(strcAff)
  
  return(integration)
  
}


##### constStruct #####

constStructureLayer <-  function(targFps) {
  ## Correlation for Structure (Tanimoto metric)
  fpSim <- fingerprint::fp.sim.matrix(targFps, method = "tanimoto")
  rownames(fpSim) <- names(targFps)
  colnames(fpSim) <- names(targFps)
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  fpAff <- SNFtool::affinityMatrix(1-fpSim, 20, 0.5)
  return(fpAff)
}

