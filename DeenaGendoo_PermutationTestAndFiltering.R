# Deena M.A. Gendoo
# Created March-13-2018
# Calculate the significance of getting a particular drug ranking in MVA-DNF
################################################################################
################################################################################
Pathway="MVA"

source('FunctionsBank.R')

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
library(pheatmap)

InterestDrug<-"DIPYRIDAMOLE" #dipyridamole,fluvastatin,methotrexate
NumHits<-238 #Number of top hits to tabulate agains the query drug. Use all drugs here

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

Iterations<-999
DrugHitMatrix<-matrix(data = NA,nrow = 238,ncol = 1)

################################################
load("Data/MVA_DNF.RData")
W<-integrtStrctSensPert
InterestDrug.SNF<-W[InterestDrug,]
DP_Hits_MVA<-names(sort(InterestDrug.SNF,decreasing = T))

set.seed(12345)

#Load the Permutation Matrix Used in the publication
load("Data/DrugHitMatrix.RData") 
if(!exists("DrugHitMatrix"))
{
  #Iterative repeat
  for(count in 1:Iterations)
  {
    message("Iteration: ",count)
  
    load("Data/PrePros.RData") #loading data layers
    
    RandomGenes<-sample(rownames(pertData),6)
    pertData<-pertData[(which(rownames(pertData) %in% RandomGenes)),] 
  
    commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                         sort(colnames(pertData))))
    strcData<- strcData[commonDrugs]
    sensData <- sensData[,commonDrugs]
    pertData <- pertData[,commonDrugs]
    
    strcAffMat <- constStructureLayer(strcData)
    sensAffMat <- constSensitivityLayer(sensData)
    pertAffMat <- constPerturbationLayer(pertData)
    integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)
    
    InterestDrug.SNF<-integrtStrctSensPert[InterestDrug,]
    TOPS<-names(sort(InterestDrug.SNF,decreasing = T))
    
    DrugHitMatrix<-cbind(DrugHitMatrix,TOPS)
  }
  
  DrugHitMatrix<-DrugHitMatrix[,-1]
  colnames(DrugHitMatrix)<-c(paste("Round_",rep(1:ncol(DrugHitMatrix)),sep = ""))
  
  DrugHitMatrix<-cbind("MVA"=DP_Hits_MVA,DrugHitMatrix)
  
  save(DrugHitMatrix,file="DrugHitMatrix.RData")
  
}


DrugList<-DrugHitMatrix[,1]

Matrix_RANK<-data.frame(DrugList)  
for(sample in 1:ncol(DrugHitMatrix))
{
  TempRankHuman<-DrugHitMatrix[,sample]
  TempRankHuman<-TempRankHuman[which(TempRankHuman %in% DrugList)]
  Matrix_RANK[,(colnames(DrugHitMatrix)[sample])]<-match(factor(DrugList),factor(TempRankHuman))
}
rownames(Matrix_RANK)<-DrugList
Matrix_RANK<-Matrix_RANK[,-1]


RanksZScore<-RanksPVal<-RanksPVal2<-NULL

for(counter in 1:nrow(Matrix_RANK))
{
  DrugRanksForSelectedDrug<-as.numeric(Matrix_RANK[counter,])
  
  pop_sd<-sd(DrugRanksForSelectedDrug)
  pop_mean <- mean(DrugRanksForSelectedDrug)
  zscore <- (counter - pop_mean) / pop_sd
  pvalue <- pnorm(counter, pop_mean, pop_sd)

  RanksZScore<-c(RanksZScore,zscore)
  RanksPVal<-c(RanksPVal,pvalue)
}

FullMatrix<-(cbind("MVA_Rank"=Matrix_RANK[,1],RanksZScore,RanksPVal))
rownames(FullMatrix)<-rownames(Matrix_RANK)

#FINAL OUTPUT
FullMatrix<-(FullMatrix[(FullMatrix[,3]<0.05 & FullMatrix[,2]<(-1.8)),])

write.csv(FullMatrix,file="DrugHit_StatScores.csv")
