# Deena M.A. Gendoo
# Accompanying code for publication release, Dec-2020

# Version Control: 
# Code developed by DG, updated October-31-2018
# Updated in Sept-2019 for focus on top drug hits

# GOAL: Look at the expression of MVA genes in the drug perturbation signatures
# Do this top drug hits from MVA-DNF algorithm

###############################################################################################
###############################################################################################
library(snowfall)
library(piano)
library(GSA)
library(GSVA)
library(org.Hs.eg.db)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(org.Hs.eg.db)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
colsch2<-rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))


#####################################################
# MVA GENES OF INTEREST
#####################################################
#Genes of the Mavelonate pathway that were used in DNF
MVAGenes<-c("ACAT2","ACLY","FDFT1","HMGCR","HMGCS1","INSIG1")

#####################################################
# DRUGS OF INTEREST
#####################################################
#Drugs of Interest: Top Drug Hits based on MVA-DNF 
# List is composed of DP + statins + Top 10 drugs 
load("Data/MVA_DNF_DrugSigs.RData")

# Match those drugs against the 414 drugs from the general DNFs
DrugsOfInterest<-read.csv("Data/MVA_DNF_DrugList.csv")

#Intentionally plotted as the full heatmap without fluvastatin or DP
colnames(MVA_DNF_DrugSigs)<-DrugsOfInterest$DrugName[match(colnames(MVA_DNF_DrugSigs),DrugsOfInterest$DrugID)]

Top25<-MVA_DNF_DrugSigs[MVAGenes,]
Top23_MVA<-MVA_DNF_DrugSigs[MVAGenes,]
Top23_MVA<-Top23_MVA[,-c(24,25)] 


Top19_MVA<-Top23_MVA[,-c(20:23)] 

par(mfrow=c(1,2),mar=c(1,6,6,6))
pdf("DrugPertSignatures_MVAgenes.pdf",width = 10,height = 10)
res<-heatmap.2(main = "Drug Pert Sigs - MVA genes",Top19_MVA,col = colsch2,
               trace="none",cexRow = 1,cexCol = 0.6,Rowv = T,Colv = T)
genelist<-rev(res$rowInd)
heatmap.2(main = "Drug Pert Sigs - MVA genes",Top25[genelist,],col = colsch2,
          trace="none",cexRow = 1,cexCol = 0.6,Rowv = F,Colv = F)
dev.off()