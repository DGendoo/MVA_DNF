synergySurf <- function(x1, #m1-by-1 matrix of concentrations of drug 1 tested in monotherapy
                        y1, #m1-by-1 matrix of viabilities corresponding to concentrations x1
                        x2, #m2-by-1 matrix of concentrations of drug 2 tested in monotherapy
                        y2, #m2-by-1 matrix of viabilities corresponding to concentrations x2
                        combination_concentrations1, #m3-by-1 matrix of concentrations of drug 1 tested in combination
                        combination_concentrations2, #m4-by-1 matrix of concentrations of drug 2 tested in combination
                        combination_viability, #m3-by-m4 matrix of viabilities of each combination of one of the m1 concentrations
                        #of drug 1 with one of the m2 concentrations of drug 2
                        hill_fit1, #list of 3 Hill fit parameters (HS, E_inf, EC50) for monotherapy with drug 1
                        hill_fit2, #list of 3 Hill fit parameters (HS, E_inf, EC50) for monotherapy with drug 2
                        method = c("Bliss", "Loewe", "HighestSingleEffect", "ResponseAdditivity"),
                        reference_viability, #viability of reference concentration for Loewe Additivity
                        conc_as_log = FALSE, #are concentrations in log?
                        viability_as_pct = FALSE, #are viabilities in pct?
                        fold_change = FALSE, #should synergy scores be reported as fold changes? See details.
                        type = TRUE, retRef = F) { #subtract observed from added? else divide observed over added [fold_change=false]
  
  method <- match.arg(method)
  
  if(!missing(hill_fit1) && is.list(hill_fit1)) hill_fit1 <- unlist(hill_fit1)
  if(!missing(hill_fit2) && is.list(hill_fit2)) hill_fit2 <- unlist(hill_fit2)
  
  if(missing(y1)){
    if (0 %in% as.numeric(dimnames(combination_viability)[[1]])){
      y1 <- combination_viability["0",]
      y1 <- y1[2:length(y1)]
      combination_viability <- combination_viability[-grep(rownames(combination_viability), pattern="^0$"),]
    }
  }
  if(missing(y2)){
    if (0 %in% as.numeric(dimnames(combination_viability)[[2]])){
      y2 <- combination_viability[,"0"]
      combination_viability <- combination_viability[,-grep(colnames(combination_viability), pattern="^0$")]
    }
  }
  #CONVERT VIABILITY TO DECIMAL AND CONCENTRATION TO LOG
  if (viability_as_pct) {
    combination_viability <- combination_viability / 100
    if (!missing(y1)) {
      y1 <- y1 / 100
    }
    if (!missing(y2)) {
      y2 <- y2 / 100
    }
    if (!missing(hill_fit1)) {
      hill_fit1[2] <- hill_fit1[2] / 100
    }
    if (!missing(hill_fit2)) {
      hill_fit2[2] <- hill_fit2[2] / 100
    }
    if (!missing(reference_viability)) {
      reference_viability <- reference_viability / 100
    }
  }
  
  
  #SET MISSING INPUTS TO DEFAULT VALUES
  if (missing(combination_concentrations1)) {
    if(is.null(dimnames(combination_viability)[[1]])){
      if(!missing(x1)){
        combination_concentrations1 <- x1
      } else {
        stop("No concentrations provided for the combination_viability for drug 1 in either rownames of single drug concentration variables")
      }
    } else {
      combination_concentrations1 <- as.numeric(dimnames(combination_viability)[[2]])
    }
  }
  
  if (missing(x1)){
    x1 <- combination_concentrations1
  }
  
  if (missing(combination_concentrations2)) {
    if(is.null(dimnames(combination_viability)[[2]])){
      if(!missing(x2)){
        combination_concentrations2 <- x2
      } else {
        stop("No concentrations provided for the combination_viability for drug 2 in either rownames of single drug concentration variables")
      }
    } else {
      combination_concentrations2 <- as.numeric(dimnames(combination_viability)[[1]])
      
    }
  }
  if (missing(x2)){
    x2 <- combination_concentrations2
  }
  if (!conc_as_log) {
    x1 <- log10(x1)
    x2 <- log10(x2)
    if (!missing(combination_concentrations1)) {
      combination_concentrations1 <- log10(combination_concentrations1)
    }
    if (!missing(combination_concentrations2)) {
      combination_concentrations2 <- log10(combination_concentrations2)
    }
    if (!missing(hill_fit1)) {
      hill_fit1[3] <- log10(hill_fit1[3])
    }
    if (!missing(hill_fit2)) {
      hill_fit2[3] <- log10(hill_fit2[3])
    }
  }
  
  
  
  if (missing(hill_fit1) && all(x1 == combination_concentrations1) ) {
    y1 <- PharmacoGx:::.Hill(combination_concentrations1, unlist(logLogisticRegression(conc = x1,
                                                                                       viability = y1,
                                                                                       conc_as_log = TRUE,
                                                                                       viability_as_pct = FALSE)))
  } else {
    if (missing(y1) || any(x1 != combination_concentrations1)) {
      if(missing(hill_fit1)){
        stop("If the combination concentrations for drug 1 provided don't match the monotherapy concentrations and no monotherapy viabilities are passed in, or if monotherapy viabilities are not included in the combination matrix, 
             hill_fit1 parameters must be provided to estimate the viabilities.")
      } 
      y1 <- PharmacoGx:::.Hill(combination_concentrations1, hill_fit1)
      }
  }
  if (missing(hill_fit2) && all(x2 == combination_concentrations2)) {
    y2 <- PharmacoGx:::.Hill(combination_concentrations2, unlist(logLogisticRegression(conc = x2,
                                                                                       viability = y2,
                                                                                       conc_as_log = TRUE,
                                                                                       viability_as_pct = FALSE)))
  } else {
    if (missing(y2) || any(x2 != combination_concentrations2)) {
      if(missing(hill_fit2)){
        stop("If the combination concentrations for drug 2 provided don't match the monotherapy concentrations and no monotherapy viabilities are passed in, or if monotherapy viabilities are not included in the combination matrix, 
             hill_fit2 parameters must be provided to estimate the viabilities.")
      } 
      y2 <- PharmacoGx:::.Hill(combination_concentrations2, hill_fit2)
      } 
  }
  
  #CALCULATE PREDICTED VIABILITY SURFACE
  if (method == "Bliss") {
    predictedSurface <- y1 %*% t(y2)
  } else if (method == "Loewe") {
    if (missing(reference_viability)) {
      stop("No reference viability provided for Loewe Additivity calculation.")
    } else {
      m1 <- lm((y1 - 1) ~ (10 ^ x1) - 1)[[1]][[2]]
      m2 <- lm((y2 - 1) ~ (10 ^ x2) - 1)[[1]][[2]]
      r1 <- (reference_viability - 1) / m1
      r2 <- (reference_viability - 1) / m2
      predictedSurface <- matrix(NA, nrow = NROW(y1), ncol = NROW(y2))
      for (i in 1:nrow(y1)) {
        for (j in 1:nrow(y2)) {
          predictedSurface[i, j] <- 1 - (y1[i] / r1 + y2[j] / r2) * (1 - reference_viability)
        }
      }
    }
  } else if (method == "HighestSingleEffect") {
    predictedSurface <- matrix(NA, nrow = NROW(y1), ncol = NROW(y2))
    for (i in 1:nrow(predictedSurface)) {
      for (j in 1:ncol(predictedSurface)) {
        predictedSurface[i, j] <- min(y1[i], y2[j])
      }
    }
  } else {
    predictedSurface <- matrix(NA, nrow = NROW(y1), ncol = NROW(y2))
    for (i in 1:nrow(predictedSurface)) {
      for (j in 1:ncol(predictedSurface)) {
        predictedSurface[i, j] <- y1[i] + y2[j] - 1
      }
    }
  }
  if(fold_change){
    synergySurface <- (predictedSurface - combination_viability) / (1 - predictedSurface) * sign(1 - predictedSurface)
  } else {
    if(type){
      synergySurface <- t(predictedSurface) - combination_viability 
    }else{
      synergySurface <- (t(predictedSurface)) / (combination_viability)
    }
  }
  if(viability_as_pct){
    synergySurface <- synergySurface * 100
  }
  if(retRef){
    return(t(predictedSurface))
  }
    return(synergySurface)
  } 