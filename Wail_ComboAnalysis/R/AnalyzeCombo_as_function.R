

getSensitivityComboMatrix <- function(viabilityFiles,x1,x2){
  
  require(abind)
  require(robustbase)
  sensitivityData <- lapply(viabilityFiles, function(x, x1,x2){
    
    print(x)
    xx <- parseComboViabilityFile(x, x1,x2)
    
    return(xx)
    
  }, x1=x1,x2=x2)
  
  sensitivityData <- abind( sensitivityData, along=1)
  
  xx <- unlist(strsplit(gsub("_[1-9]$", "", rownames(sensitivityData)), split=' (?=[^ ]+$)', perl=TRUE))
  
  #cells_names <- unlist(lapply(strsplit(rownames(sensitivityData)," "),"[[",1))
  #cells_names <- cells_names[cells_names!="1"]
  #cells_names <- cells_names[cells_names!="2"] 
  #cells_names <- cells_names[cells_names!="3"]
  #cells_names <- cells_names[cells_names!="4"]
  #cells_names <- unique(cells_names)
  
  
#  sensitivityData_Final <- array(NA, dim=c(length(cells_names),dim(sensitivityData)[2], dim(sensitivityData)[3]), dimnames=list(cells_names, rownames(sensitivityData[1,,]), colnames(sensitivityData[1,,])))
  
#  for(cell_tmp  in cells_names){
#    ibx <- grep(paste("^",cell_tmp,sep = ""),rownames(sensitivityData))

#      if(length(ibx)>1){
#        sensitivityData_Final[cell_tmp,,] <- apply(sensitivityData[ibx,,],MARGIN = 3,FUN = function(x){
#          return(colMedians(x))
#      })
#    }else{
#      sensitivityData_Final[cell_tmp,,] <- sensitivityData[ibx,,]
#    }
    
    
#  }
  
#  sensitivityData <- sensitivityData_Final
  
  
  
  for( i in 1:dim(sensitivityData)[1]){
    
#    sensitivityData[i,,] <- sensitivityData[i,,] - 0.045
    control <- sensitivityData[i,,][dim(sensitivityData[i,,])[1],1]
    sensitivityData[i,,] <- sensitivityData[i,,]/control
    
  }
  
  

  return(sensitivityData)
}




parseComboViabilityFile <- function(fileCombo, conc1, conc2){
  # conc1 <- x1
  #  conc2 <- x2
  
  #  fileCombo <- viabilityFiles[1]
  
  conc1 <- rev(conc1)
  viability <- read.csv(fileCombo, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  
  ##this function assumes in the text file produced by the machine always there is field called "Plate:" before plate names  
  indicies <- grep("Plate", viability$V1)
  plates <- viability$V2[indicies]
  
  data.blocks <- table(is.na(viability[indicies[1] + 3, ]))["FALSE"] %/% 10#ncol(viability) %/% 10
  
  raw.sensitivity.rows <- NULL
  for(plate in plates) {
    
    xx <- unlist(strsplit(gsub("[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))
    # xx <- unlist(strsplit(gsub("[1-9]$", "", xx[1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    
    drugs <- unlist(strsplit(xx[length(xx)], split='/'))
    
    if(data.blocks > 1){
      raw.sensitivity.rows <- c(raw.sensitivity.rows, do.call(c, lapply(1:data.blocks, function(x){paste(paste(cell, drugs, sep="_"), x, sep="_")})))
    }else{
      raw.sensitivity.rows <- c(raw.sensitivity.rows,cell)
    }
  }
  
  if(length(raw.sensitivity.rows)!=length(unique(raw.sensitivity.rows))){
    stop("Cell lines names are not unique")
  }
  
  
  raw.sensitivity <- array(NA, dim=c(length(raw.sensitivity.rows),length(conc1), length(conc2)), dimnames=list(raw.sensitivity.rows, conc1, conc2))
  
  for(index in indicies) {
    start <- 3
    
    plate <- viability$V2[index]
    
    xx <- unlist(strsplit(gsub("_[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))
    #    xx <- unlist(strsplit(gsub("_[1-9]$", "", xx[1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    drugs <- paste(xx[(length(xx))], sep=" ")
    
    
    
    for(i in 1:length(conc1)){
      raw.sensitivity[cell,i,] <- as.numeric(do.call(c,viability[(index + i + 2), (start + 1):(start + 10)]))
      if(!is.na(viability[(index + i + 2),start])){
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - viability[(index + i + 2),start]
      }else{
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - 0.045
      }
    }
    
    
  }
  
  return(raw.sensitivity)
  
}




ourHeatmap <- function(mat, title="", coldivs=200, min_of_scale=-1, max_of_scale=1, drug1, drug2){
  
  require(ggplot2)
  require(reshape2)
  Synergy <- "Synergy"
  mat.m <- melt(mat, varnames = c( "Celline", drug1,drug2), value.name=c("Synergy"))
  mat.m <- mat.m[,c(drug1,drug2,"Celline","Synergy")]
  colscheme = colorRampPalette(c('#fc8d59', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'))
  mat.m$Synergy <- pmin(mat.m$Synergy,max_of_scale)
  mat.m$Synergy <- pmax(mat.m$Synergy,min_of_scale)
  
  x1 <- unique(mat.m[,1])
  x2 <- unique(mat.m[,2])
  
  
  cols <- c('#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')
#  cols <- c( '#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4')
  
  p <- ggplot(mat.m, aes_string(drug2,drug1 , fill=Synergy)) + geom_tile() + scale_x_continuous(breaks=mat.m[,2], trans="log10",expand = c(0,0), labels = function(x){sprintf("%3g",x)})  + scale_y_continuous(breaks=mat.m[,1], trans = "log10",expand=c(0,0))  + xlab(paste(drug2,"(µM)")) + ylab(paste(drug1,"(µM)")) + scale_fill_gradientn(limits = c(-1,1),colours=rev(cols)) + theme(axis.text = element_text(size=12),title=element_text(size=16), axis.text.x=element_text(size=12,hjust=1,vjust=0.5,angle=90), strip.text.x = element_text(size=14), legend.text=element_text(size=14, hjust=1)) + facet_wrap(~Celline, ncol=8) + ggtitle(paste(title))
  #  p <- ggplot(mat.m, aes_string(drug2,drug1 , fill=Synergy)) + geom_tile() + scale_x_continuous(breaks=mat.m[,2], trans="log10",expand = c(0,0), labels = sprintf("%3g",x))  + scale_y_continuous(breaks=mat.m[,1], trans = "log10",expand=c(0,0))  + xlab(paste(drug2,"(µM)")) + ylab(paste(drug1,"(µM)")) + scale_fill_gradientn(colours=rev(c( '#d73027','#fc8d59','#fee08b','#d9ef8b','#91cf60','#1a9850'))) + theme(axis.text = element_text(size=12),title=element_text(size=16), axis.text.x=element_text(size=12,hjust=1,vjust=0.5,angle=90), strip.text.x = element_text(size=14), legend.text=element_text(size=14, hjust=1)) + facet_wrap(~Celline, ncol=8) + ggtitle(paste(title))
  
  return(p)
}



ourHeatmap_CFI945 <- function(mat, title="", coldivs=200, min_of_scale=-1, max_of_scale=1, drug1, drug2){
  
  require(ggplot2)
  require(reshape2)
  Synergy <- "Synergy"
  mat.m <- melt(mat, varnames = c( "Celline", drug1,drug2), value.name=c("Synergy"))
  mat.m <- mat.m[,c(drug1,drug2,"Celline","Synergy")]
  colscheme = colorRampPalette(c('#fc8d59', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'))
  mat.m$Synergy <- pmin(mat.m$Synergy,max_of_scale)
  mat.m$Synergy <- pmax(mat.m$Synergy,min_of_scale)
  
  t1 <- unique(mat.m[,1])
  t2 <- (mat.m[,2])
  
  mat.m$CFI945[which(mat.m$CFI945==10000)] <- 6666.6666
  names(t2) <- mat.m$CFI945
  #  p <- ggplot(mat.m, aes_string(drug2,drug1 , fill=Synergy)) + geom_tile() + scale_x_continuous(breaks=mat.m[,2], trans="log10",expand = c(0,0), labels = function(x){sprintf("%3g",x)})  + scale_y_continuous(breaks=mat.m[,1], trans = "log10",expand=c(0,0))  + xlab(paste(drug2,"(µM)")) + ylab(paste(drug1,"(µM)")) + scale_fill_gradientn(colours=rev(c( '#d73027','#fc8d59','#fee08b','#d9ef8b','#91cf60','#1a9850'))) + theme(axis.text = element_text(size=12),title=element_text(size=16), axis.text.x=element_text(size=12,hjust=1,vjust=0.5,angle=90), strip.text.x = element_text(size=14), legend.text=element_text(size=14, hjust=1)) + facet_wrap(~Celline, ncol=8) + ggtitle(paste(title))
  p <- ggplot(mat.m, aes_string(drug2,drug1 , fill=Synergy)) + geom_tile() + scale_x_continuous(breaks=as.numeric(names(t2)), trans="log10",expand = c(0,0), labels = sprintf("%3g",as.numeric((t2))))  + scale_y_continuous(breaks=mat.m[,1], trans = "log10",expand=c(0,0))  + xlab(paste(drug2,"(nM)")) + ylab(paste(drug1,"(nM)")) + scale_fill_gradientn(limits = c(min_of_scale,max_of_scale),colours=rev(c( '#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4'))) + theme(axis.text = element_text(size=12),title=element_text(size=16), axis.text.x=element_text(size=12,hjust=1,vjust=0.5,angle=90), strip.text.x = element_text(size=14), legend.text=element_text(size=14, hjust=1)) + facet_wrap(~Celline, ncol=8) + ggtitle(paste(title))
  
  return(p)
}

