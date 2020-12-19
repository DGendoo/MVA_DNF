
SynElem <- function(i,j,DoseRow, DoseCol, DiffMat){
  SynElem <- 0.5*(DoseRow[(i+1)]-DoseRow[(i-1)])*0.5*(DoseCol[(j+1)]-DoseCol[(j-1)])*DiffMat[i,j]
  return(SynElem)
}


SynVolume <- function(DoseRow, DoseCol, DiffMat){
  
  DoseRow <- log10(DoseRow)
  DoseCol <- log10(DoseCol)
  SynVol <- 0
  SynVol <- SynVol+0.5*(DoseRow[2]-DoseRow[1])*0.5*(DoseCol[2]-DoseCol[1])*DiffMat[1,1]
  SynVol <- SynVol+0.5*(DoseRow[length(DoseRow)]-DoseRow[(length(DoseRow)-1)])*0.5*(
    DoseCol[length(DoseCol)]-DoseCol[(length(DoseCol)-1)])*DiffMat[length(DoseRow),length(DoseCol)]
  SynVol <- SynVol+0.5*(DoseRow[length(DoseRow)]-DoseRow[(length(DoseRow)-1)])*0.5*(DoseCol[2]-DoseCol[1])*DiffMat[length(DoseRow),1]
  SynVol <- SynVol+0.5*(DoseRow[2]-DoseRow[1])*0.5*(DoseCol[length(DoseCol)]-DoseCol[(length(DoseCol)-1)])*DiffMat[1,length(DoseCol)]

  for(i in 2:(length(DoseRow)-1)){
    for(j in 2:(length(DoseCol)-1)){
      SynVol <- SynVol+SynElem(i,j,DoseRow,DoseCol,DiffMat)
    }
  }
  
  return(SynVol)
}




