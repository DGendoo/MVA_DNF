# Deena M.A. Gendoo

# Created by DG Sept-2017, updated on Sept-2019
# Show contributions of each layer in MVA-DNF drug hits
##########################################################################################
##########################################################################################

library(pheatmap)
library(fmsb)

TOPDrugs<-c("DIPYRIDAMOLE","SELUMETINIB","NELFINAVIR","MITOXANTRONE","DOXORUBICIN","HONOKIOL","CLOTRIMAZOLE",
"SULFATHIAZOLE","VEMURAFENIB","CHROMOMYCINA3","BACCATINIII","NOSCAPINE",
"METHOTREXATE","CADMIUMCHLORIDE","RHAMNETIN","PENTAMIDINE","EMODIN","FLUOROURACIL",
"TRYPTOPHAN","ALVOCIDIB")    

InterestDrug<-"DIPYRIDAMOLE" #dipyridamole,fluvastatin,methotrexate
NumHits<-19 #Number of top hits to tabulate agains the query dru. Use ALL DRUGS 

Matrices<-c("integrtStrctSensPert","strcAffMat","sensAffMat","pertAffMat")
MatricesNames<-c("DNF","Structure","Sensitivity","Perturbation")

load("Data/MVA_DNF.RData")

  ################################################
  # Get Affinity Matrix weights for the top-ranked drugs 
  ################################################
  WeightList<-matrix(nrow = (NumHits+1),ncol = 4) 
  colnames(WeightList)<-c("DNF","Structure","Sensitivity","Perturbation")
  
  W<-integrtStrctSensPert
  grep(InterestDrug,rownames(W),perl = F,value = T) #Sanity check
  InterestDrug.SNF<-W[InterestDrug,]
  TOPS<-as.matrix(InterestDrug.SNF[TOPDrugs])
  rownames(WeightList)<-rownames(TOPS)
  WeightList[,1]<-TOPS[,1]
  
  for (layercount in 2:length(Matrices)) 
  {
    W<- get(Matrices[layercount])
    InterestDrug.SNF<-W[InterestDrug,]
    
    for(drugcount in 1:nrow(WeightList)) 
    {
      WeightList[drugcount,layercount]<-InterestDrug.SNF[as.character(rownames(WeightList)[drugcount])]
    }
  }
  
  WeightSub<-t(WeightList[,2:4])
  WeightSub2<-prop.table(WeightSub, margin=2) #Converts the weight contributions it to a % of the total
  
  SPECIFIC<-WeightSub2
  
  rm(list = c("integrtStrctSensPert","strcAffMat","sensAffMat","pertAffMat"))
  
SpecificNew<-SPECIFIC[,TOPDrugs]
SpecificNew <- cbind(SpecificNew) #Add NA for the legend margin

#################################################
# Create a radar plot (aka star plot) showing the layer contributions
##################################################
SpecificNew=rbind(rep(1,ncol(SpecificNew)) , rep(0,ncol(SpecificNew)) , SpecificNew)
SpecificNew<-as.data.frame(SpecificNew)

write.csv(SpecificNew[-c(1:2),],file="ComparisonLayerWeights_Percentage.csv")

SpecificNew<-SpecificNew[,-1]

pdf(file="ComparisonOfLayerWeights.pdf",width = 20,height = 10)
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
radarchart( SpecificNew  , axistype=1 , title="MVA-DNF",
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.2,5), cglwd=0.8,
            #custom labels
            vlcex=0.9 
)
legend(x=1, y=1, legend = rownames(SpecificNew[-c(1,2),]), bty = "n", pch=20 , 
       col=colors_in , text.col = "grey", cex=0.9, pt.cex=2)
dev.off()


