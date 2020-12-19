fixCellLinesNames <- function(x, ref){



cell_all <- read.csv(ref, na.strings=c("", " ", "NA"))


cell_unannot <- unique(x)

#cell_annot <- gsub(x=cell_unannot, pattern="^[0-9]+-", replacement="")

cell_annot <- toupper(cell_unannot)

cell_unannot <- cell_unannot[order(cell_annot)]
cell_annot <- cell_annot[order(cell_annot)]

cellMatch <- data.frame("annotated"=cell_annot, "unannotated"=cell_unannot)


badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


closeMatches <- lapply(cellMatch[,"annotated"], function(x){
  x <- as.character(x)
  myx <- as.character(cell_all$unique.cellid[which(x == toupper(cell_all[,"unique.cellid"]))])
  if(length(myx)==0){
    if(grepl(pattern="^MDA", x)){ ## different way of writing these cell lines 
      print(x)
      if(is.na(strsplit(x, split="B")[[1]][2])){
        x <- paste(strsplit(x, split="A")[[1]][1], "A", "MB",  strsplit(x, split="A")[[1]][2])
      }else{
        x <- paste(strsplit(x, split="A")[[1]][1], "A", "MB",  strsplit(x, split="B")[[1]][2])
      }
    }
    if(grepl(pattern="^UAC[1-9]+", x)){ #common typo
      print(x)
      x <- paste(strsplit(x, split="C")[[1]][1], "CC", strsplit(x, split="C")[[1]][2])
    }
    myx <- grep(pattern=toupper(gsub(badchars, "",x)), x=toupper(gsub(badchars, "",cell_all$unique.cellid)))
    myx <- as.character(cell_all$unique.cellid[myx])
  }
  if(x=="MFN223"){
    myx <- "MFM-223"
  } ## typo in data
  if(x=="AU655"){
    myx <- "AU565"
    print(x)
  } ## typo in data
  if(x=="HCL70"){
    myx <- "HCC70"
  } ## typo in data
  if(x=="MPE600"){
    myx <- "600MPE"
  } ## typo in data 
  if(x=="OCUB-1" | x=="OCUB1"){
    myx <- "OCUB-M"
  } ## Closest match (is child of cell line)
  return(ifelse(length(myx)>0,myx,NA))
})
sum(!sapply(closeMatches,is.na))
### Fix OCUB HCL70? 
########### I still have problems matching the cell lines HCL 70 and OCUB-1.
########### I assume HCL70 is HCC70 
########### MPE600 -> 600MPE
cellMatch[,"closeMatches"] <- unlist(closeMatches)
#cellMatch[2,"closeMatches"] <- "AU565"
warning("All matches are now OK, but this may change if other cells are added in the future")
cellMatch$annotated <- trimws(as.character(cellMatch$annotated))
cellMatch$unannotated <- trimws(as.character(cellMatch$unannotated))
cellMatch$closeMatches <- trimws(as.character(cellMatch$closeMatches))

cellMatch$annotated[!is.na(cellMatch$closeMatches)] <- na.omit(cellMatch$closeMatches)
x1 <- cellMatch$annotated[match(trimws(x), cellMatch$unannotated)]

return(x1)
}
