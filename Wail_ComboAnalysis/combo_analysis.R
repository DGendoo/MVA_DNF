require(PharmacoGx)
require(magicaxis)
library(abind)
library(robustbase)
library(Biobase)
library(synergyfinder)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(reshape)
library(snowfall)
library(GSA)
library(piano)
library(scales)
library(ggrepel)
library(dplyr)
library(tidyr)
library(BoutrosLab.plotting.general)
library(UpSetR)




##################
####Important#####
##################
# set the working directory to "FluvaDP_repo-main" folder
setwd("Wail_ComboAnalysis")



# set whether you want to recalculate the synergy values from raw data
PerformSynergyCalcFromRaw = F




if(PerformSynergyCalcFromRaw){
  # code used to calculate the synergy values from raw data
  source("R/calculateSynergy_combo.R")
}else{
  load("data/SynergyStat_Statin_combo_23Sep19.RData",verbose = T)
}

# code used to standardize cell lines names across experiments 
source("R/fixCellLinesNames.R")


#######################################
#######################################
#######################################
# Calculate Fluvastatin mono effect from the combo assays

listOfCombos_MonoFluva <- lapply(names(listOfCombos), function(combo){
  
  monoResults <- do.call(rbind,lapply(names(listOfCombos[[combo]]$Bliss),function(sample){
    
    listOfCombos[[combo]]$Bliss[[sample]]$dose.response.mats$`1`[1,]
  }))
  
  rownames(monoResults) <- names(listOfCombos[[combo]]$Bliss)
  colnames(monoResults) <- names(listOfCombos[[combo]]$Bliss[[1]]$dose.response.mats$`1`[1,])
  
  uniq_samples <- unique(unlist(lapply(strsplit(rownames(monoResults),"_"),"[[",1)))
  monoResults_final <- do.call(rbind,lapply(uniq_samples, function(x){
    ibx <- grep(x,x = rownames(monoResults))
    return(colMedians(monoResults[ibx,,drop=F]))
  }))
  
  rownames(monoResults_final) <- uniq_samples
  
  rownames(monoResults_final) <- fixCellLinesNames(rownames(monoResults_final),"data/CL_data/cell_annotation_all.csv")
  
  fluvaMono <- lapply(1:dim(monoResults_final)[1],function(x){
    
    ic50 <- computeIC50(concentration = as.numeric(colnames(monoResults_final)[-1]),viability = monoResults_final[x,-1],viability_as_pct = T)
    aac <- computeAUC(concentration = as.numeric(colnames(monoResults_final)[-1]),viability = monoResults_final[x,-1],viability_as_pct = T)/100
    return(c(ic50,aac))
  })
  
  fluvaMono <- do.call(rbind,fluvaMono)
  
  rownames(fluvaMono) <- rownames(monoResults_final)
  colnames(fluvaMono) <- c("IC50","AAC")
  
  return(fluvaMono)
  
})

names(listOfCombos_MonoFluva) <- names(listOfCombos)

common <- Reduce(intersect,list(rownames(listOfCombos_MonoFluva$FLUVA_DP),rownames(listOfCombos_MonoFluva$FLUVA_NFV),rownames(listOfCombos_MonoFluva$FLUVA_HNK)))
Fluva <- colMeans(rbind(listOfCombos_MonoFluva$FLUVA_DP[common,"AAC"],listOfCombos_MonoFluva$FLUVA_NFV[common,"AAC"],listOfCombos_MonoFluva$FLUVA_HNK[common,"AAC"]))
FluvaMonocommonAvg_IC50 <- colMeans(rbind(listOfCombos_MonoFluva$FLUVA_DP[common,"IC50"],listOfCombos_MonoFluva$FLUVA_NFV[common,"IC50"],listOfCombos_MonoFluva$FLUVA_HNK[common,"IC50"]))


#######################################
#######################################
#######################################
# Prepare the data for Fig 4A and Fig S6B


# SMOD2 breast cancer subtypes for the cell lines used in this paper
subtypes <- read.csv("data/CL_data/subtypes.txt",header = T,row.names = 1,stringsAsFactors = F,sep = "\t")



#summarize combo stat using average Bliss index across concentrations
listOfCombos_stat_summarized = list()

for (i in 1:length(listOfCombos_stat)) {
  combo <- names(listOfCombos_stat)[i]
  stat <- listOfCombos_stat[[combo]]
  
  sampleIDs <- rownames(stat)
  
  sampleIDs <- unique(unlist(lapply(strsplit(sampleIDs,"_"),"[[",1)))
  stat_summarized <- matrix(nrow = length(sampleIDs),ncol = ncol(stat),dimnames = list(sampleIDs,colnames(stat)))
  for (sample in sampleIDs) {
    ibx <- grep(sample,rownames(stat))
    stat_summarized[sample,] <- colMeans(stat[ibx,],na.rm = T)
    
    if(startsWith(sample,"X")){
      rownames(stat_summarized)[grep(sample,rownames(stat_summarized))] <- paste("MDAMB",gsub("X","",sample),sep = "")
    }
    
  }
  
  final <- fixCellLinesNames(rownames(stat_summarized),"data/CL_data/cell_annotation_all.csv")
  cbind(rownames(stat_summarized),final)
  rownames(stat_summarized) <- final
  listOfCombos_stat_summarized[[combo]] <- stat_summarized
  
}


unionCell <- Reduce(union,list(rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[1]]])
                               ,rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[2]]])
                               ,rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[3]]])))



BlissMat_summarized_all <- lapply(listOfCombos_stat_summarized, function(x){
  
  tmp <- c(x[,"Bliss"],rep(NA,length(setdiff(unionCell,names(x[,"Bliss"])))))
  names(tmp) <- c(rownames(x),(setdiff(unionCell,names(x[,"Bliss"]))))
  return(tmp[order(names(tmp))])
})

common <- Reduce(union,list(names(BlissMat_summarized_all[[1]]),names(BlissMat_summarized_all[[2]]),names(BlissMat_summarized_all[[3]])))

BlissMat_summarized_all <- do.call(rbind,lapply(BlissMat_summarized_all, function(x){
  return(x[common])
}))


commonCells <- Reduce(intersect,list(rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[1]]])
                                     ,rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[2]]])
                                     ,rownames(listOfCombos_stat_summarized[[names(listOfCombos_stat_summarized)[3]]])))

BlissMat_summarized <- do.call(rbind,lapply(listOfCombos_stat_summarized, function(x){
  return(x[commonCells,"Bliss"])
}))

commonMono <- intersect(names(Fluva),colnames(BlissMat_summarized_all))

Fluva.sub <- Fluva[commonMono]
Fluva.sub <- sort(Fluva.sub)
order <- names(sort(BlissMat_summarized["FLUVA_DP",]))
column_ha = HeatmapAnnotation(FLUVA.AAC = anno_barplot(Fluva.sub))
subtypes_final <- subtypes[colnames(BlissMat_summarized),"SCMOD2"]
names(subtypes_final) <- colnames(BlissMat_summarized)
column_ha = HeatmapAnnotation(Subtype = subtypes_final[order], col = list(Subtype=c("Basal"="#4daf4a","HER2"="#377eb8","LumB"="#984ea3","LumA"="#e78ac3")))


##############
### Fig 4A ###
##############

pdf("Plots/Fig4_A.pdf",width = 15,height = 2.5)
Heatmap(BlissMat_summarized_all[,order],cluster_columns = F,cluster_rows = F,top_annotation = column_ha
        ,row_split = c("A","B","C"), row_title = NULL,na_col = "black"
        ,col = colorRamp2(c(-1,0,1),(c("red","white","blue")))
        ,heatmap_legend_param = list(
          title = "Bliss", at = c(-1, 1), 
          labels = c("Synergy", "Antagonism"))
)
dev.off()

BlissMat_summarized_all_final <- data.frame(t(BlissMat_summarized),"Subtype"=subtypes[colnames(BlissMat_summarized),"SCMOD2"],stringsAsFactors = F)
BlissMat_summarized_all_final$cellID <- rownames(BlissMat_summarized_all_final)

BlissMat_summarized_all_final <- melt(BlissMat_summarized_all_final, id=c("cellID","Subtype")) 
BlissMat_summarized_all_final$Subtype <- as.factor(BlissMat_summarized_all_final$Subtype)
BlissMat_summarized_all_final$variable <- as.factor(BlissMat_summarized_all_final$variable)

colnames(BlissMat_summarized_all_final) <- c("cellID",   "Subtype",  "DrugCombo", "Synergy")

my_comparisons = list( c("Basal", "LumB"), c("LumB", "HER2"), c("Basal", "HER2") )



Subtype_col=c("Basal"="#4daf4a","HER2"="#377eb8","LumB"="#984ea3","LumA"="#e78ac3")




##############
### Fig S6B ##
##############

pdf("Plots/FigS6B.pdf",width = 11,height = 6)
ggplot(BlissMat_summarized_all_final,mapping = aes(x=Subtype,y=Synergy,fill=Subtype)) +
  geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0,lty=2,col="red") +
  geom_jitter(width = 0.2) +
  facet_wrap(~DrugCombo) + 
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.25, 0.20, 0.15))+
  stat_compare_means(label.y = 0.4) +
  scale_fill_manual(values=Subtype_col) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),strip.text.x =element_text(size=14,face="bold"))

dev.off()


#######################################
#######################################
#######################################

RNAseq <- readRDS("data/CL_data/cell_lines_expression_matrix.rds")

commonSamples <- intersect(rownames(RNAseq),colnames(BlissMat_summarized))

gene_mappings <- readRDS("data/CL_data/genes_ids_mappings.rds")

geneAssociations_cor_Mono <- apply(RNAseq[commonSamples,], 2, function(x){
  a <- cor.test(x,Fluva[commonSamples])
  return(c(a$estimate,a$p.value))
})

geneAssociations_cor_Mono <- t(geneAssociations_cor_Mono)
geneAssociations_cor_Mono <- cbind(geneAssociations_cor_Mono,p.adjust(geneAssociations_cor_Mono[,2],method = "fdr"))
colnames(geneAssociations_cor_Mono) <- c("estimate","pval","fdr")


geneAssociations_cor_Mono <- geneAssociations_cor_Mono[order(geneAssociations_cor_Mono[,"fdr"]),]

geneAssociations_cor_Mono_final <- data.frame(geneAssociations_cor_Mono,"Symbol"=gene_mappings[rownames(geneAssociations_cor_Mono),"Symbol"])
genesSymbols <- gene_mappings[rownames(geneAssociations_cor_Mono),"Symbol"]
ibx <- which(duplicated(genesSymbols))


gene_stat_levels <-  -geneAssociations_cor_Mono[-ibx,"estimate"]
names(gene_stat_levels) <- genesSymbols[-ibx]

gene_stat_levels <- gene_stat_levels[!is.na(gene_stat_levels)]

set.seed(3425)


nPerm <- 10000
gsc1 <- loadGSC("data/pathwaysInfo/h.all.v6.2.symbols.gmt")


gseaRes <- piano::runGSA(geneLevelStats = gene_stat_levels,gsc = gsc1,adjMethod = "none",nPerm = nPerm,geneSetStat = "fgsea"
                         ,ncpus = 4)

gseaResSummary <- piano::GSAsummaryTable(gseaRes)
if(length(colnames(gseaResSummary))<8){
  if("p (dist.dir.up)" %in% colnames(gseaResSummary) & !("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
    gseaResSummary[,"p (dist.dir.dn)"] <- NA
  }else if(!("p (dist.dir.up)" %in% colnames(gseaResSummary)) & ("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
    gseaResSummary[,"p (dist.dir.up)"] <- NA
  }
}
gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(nPerm+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr") # fdr correction
gseares1 <- gseares1[,-(4:5)]

gseares1 <- gseares1[order(gseares1$FDR,na.last = T),]
gseares1 <- gseares1[order(gseares1$FDR,gseares1$pval,-abs(gseares1$`Stat (dist.dir)`),na.last = T),]

gseares1Sig <- gseares1[gseares1$FDR<=0.05,]
gseares1Sig <- gseares1Sig[order(abs(gseares1Sig$`Stat (dist.dir)`),decreasing=T),]
gseares1Sig_Fluva_mono <- gseares1Sig



listOfAssociations <- lapply(rownames(BlissMat_summarized), function(y){
  
  
  geneAssociations_cor <- apply(RNAseq[commonSamples,], 2, function(x){
    a <- cor.test(x,BlissMat_summarized[y,commonSamples])
    return(c(a$estimate,a$p.value))
  })
  
  geneAssociations_cor <- t(geneAssociations_cor)
  geneAssociations_cor <- cbind(geneAssociations_cor,p.adjust(geneAssociations_cor[,2],method = "fdr"))
  colnames(geneAssociations_cor) <- c("estimate","pval","fdr")
  
  geneAssociations_cor <- geneAssociations_cor[order(geneAssociations_cor[,"fdr"]),]
  
  return(geneAssociations_cor)
})

names(listOfAssociations) <- rownames(BlissMat_summarized)



listOfAssociations_final <- lapply(listOfAssociations, function(x){
  return(data.frame(x,"Symbol"=gene_mappings[rownames(x),"Symbol"]))
})


##############
### Fig 4B ##
##############

pdf("Plots/Fig4B.pdf")

A <- cor.test(listOfAssociations$FLUVA_DP[,"estimate"],listOfAssociations$FLUVA_NFV[rownames(listOfAssociations$FLUVA_DP),"estimate"],method = "s")

df <- data.frame("Fluva_DP_PCC"=listOfAssociations$FLUVA_DP[,"estimate"], "FLUVA_NFV_PCC"=listOfAssociations$FLUVA_NFV[rownames(listOfAssociations$FLUVA_DP),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations$FLUVA_DP[,"fdr"],listOfAssociations$FLUVA_NFV[rownames(listOfAssociations$FLUVA_DP),"fdr"]),MARGIN = 1,min)
                 ,"Fluva_DP_fdr"=listOfAssociations$FLUVA_DP[,"fdr"], "FLUVA_NFV_fdr"=listOfAssociations$FLUVA_NFV[rownames(listOfAssociations$FLUVA_DP),"fdr"])

df <- df[!is.na(df$Fluva_DP_fdr) & !is.na(df$FLUVA_NFV_fdr), ]
df$fdr2=ifelse(listOfAssociations$FLUVA_DP[rownames(df),"fdr"]<0.1 & listOfAssociations$FLUVA_NFV[rownames(df),"fdr"]<0.1,T,F)


df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_mappings[rownames(df),"Symbol"]

df$d = densCols(df$Fluva_DP_PCC, df$FLUVA_NFV_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))


df_pos <- df[order(df$Fluva_DP_PCC),]
df_pos <- df_pos[df_pos$Fluva_DP_PCC > 0 & df_pos$FLUVA_NFV_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(as.numeric(x[c("Fluva_DP_PCC","FLUVA_NFV_PCC")]),na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]
df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_DP_PCC) >= 0.1 & abs(df_pos$FLUVA_NFV_PCC) >= 0.1,][1:5,]

df_pos <- df_pos[!is.na(df_pos$Fluva_DP_PCC),]

df_neg <- df[order(df$Fluva_DP_PCC),]
df_neg <- df_neg[df_neg$Fluva_DP_PCC < 0 & df_neg$FLUVA_NFV_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(as.numeric(x[c("Fluva_DP_PCC","FLUVA_NFV_PCC")]),na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]
df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_DP_PCC) >= 0.1 & abs(df_neg$FLUVA_NFV_PCC) >= 0.1,][1:5,]

df_neg <- df_neg[!is.na(df_neg$Fluva_DP_PCC),]


label <- df %>%
  summarise(
    Fluva_DP_PCC = min(Fluva_DP_PCC)+ 0.1,
    FLUVA_NFV_PCC = max(FLUVA_NFV_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )
ggplot(df) +
  geom_point(aes(Fluva_DP_PCC, FLUVA_NFV_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_DP gene associations") + ylab("FLUVA_NFV gene associations") +
  geom_text_repel(data = rbind.data.frame(df_pos,df_neg),mapping = aes(Fluva_DP_PCC,FLUVA_NFV_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_pos,df_neg),aes(Fluva_DP_PCC, FLUVA_NFV_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_DP_PCC,y=FLUVA_NFV_PCC)) 

df <- df[order(df$fdr),]


A = cor.test(listOfAssociations$FLUVA_DP[,"estimate"],listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_DP),"estimate"],method = "s")


df <- data.frame("Fluva_DP_PCC"=listOfAssociations$FLUVA_DP[,"estimate"], "FLUVA_HNK_PCC"=listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_DP),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations$FLUVA_DP[,"fdr"],listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_DP),"fdr"]),MARGIN = 1,min)
                 ,"Fluva_DP_fdr"=listOfAssociations$FLUVA_DP[,"fdr"], "FLUVA_HNK_fdr"=listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_DP),"fdr"])

df <- df[!is.na(df$fdr), ]
df$fdr2=ifelse(listOfAssociations$FLUVA_DP[rownames(df),"fdr"]<0.1 & listOfAssociations$FLUVA_HNK[rownames(df),"fdr"]<0.1,T,F)


df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_mappings[rownames(df),"Symbol"]

df$d = densCols(df$Fluva_DP_PCC, df$FLUVA_HNK_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))

df_pos <- df[order(df$Fluva_DP_PCC),]
df_pos <- df_pos[df_pos$Fluva_DP_PCC > 0 & df_pos$FLUVA_HNK_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(as.numeric(x[c("Fluva_DP_PCC","FLUVA_HNK_PCC")]),na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]

df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_DP_PCC) >= 0.1 & abs(df_pos$FLUVA_HNK_PCC) >= 0.1,][1:5,]
df_pos <- df_pos[!is.na(df_pos$Fluva_DP_PCC),]

df_neg <- df[order(df$Fluva_DP_PCC),]
df_neg <- df_neg[df_neg$Fluva_DP_PCC < 0 & df_neg$FLUVA_HNK_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(as.numeric(x[c("Fluva_DP_PCC","FLUVA_HNK_PCC")]),na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]
df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_DP_PCC) >= 0.1 & abs(df_neg$FLUVA_HNK_PCC) >= 0.1,][1:5,]

df_neg <- df_neg[!is.na(df_neg$Fluva_DP_PCC),]


label <- df %>%
  summarise(
    Fluva_DP_PCC = min(Fluva_DP_PCC) + 0.1,
    FLUVA_HNK_PCC = max(FLUVA_HNK_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )
ggplot(df) +
  geom_point(aes(Fluva_DP_PCC, FLUVA_HNK_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_DP gene associations") + ylab("FLUVA_HNK gene associations") +
  geom_text_repel(data = rbind.data.frame(df_pos,df_neg),mapping = aes(Fluva_DP_PCC,FLUVA_HNK_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_pos,df_neg),aes(Fluva_DP_PCC, FLUVA_HNK_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_DP_PCC,y=FLUVA_HNK_PCC))

df <- df[order(df$fdr),]



A <- cor.test(listOfAssociations$FLUVA_HNK[,"estimate"],listOfAssociations$FLUVA_NFV[rownames(listOfAssociations$FLUVA_HNK),"estimate"],method = "s")


df <- data.frame("Fluva_NFV_PCC"=listOfAssociations$FLUVA_NFV[,"estimate"], "FLUVA_HNK_PCC"=listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_NFV),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations$FLUVA_NFV[,"fdr"],listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_NFV),"fdr"]),MARGIN = 1,min)
                 ,"Fluva_NFV_fdr"=listOfAssociations$FLUVA_NFV[,"fdr"], "FLUVA_HNK_fdr"=listOfAssociations$FLUVA_HNK[rownames(listOfAssociations$FLUVA_NFV),"fdr"])

df <- df[!is.na(df$fdr), ]
df$fdr2=ifelse(listOfAssociations$FLUVA_NFV[rownames(df),"fdr"]<0.1 & listOfAssociations$FLUVA_HNK[rownames(df),"fdr"]<0.1,T,F)


df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_mappings[rownames(df),"Symbol"]

df$d = densCols(df$Fluva_NFV_PCC, df$FLUVA_HNK_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))

df_pos <- df[order(df$Fluva_NFV_PCC),]
df_pos <- df_pos[df_pos$Fluva_NFV_PCC > 0 & df_pos$FLUVA_HNK_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(as.numeric(x[c("Fluva_NFV_PCC","FLUVA_HNK_PCC")]),na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]
df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_NFV_PCC) >= 0.1 & abs(df_pos$FLUVA_HNK_PCC) >= 0.1,][1:5,]

df_pos <- df_pos[!is.na(df_pos$Fluva_NFV_PCC),]

df_neg <- df[order(df$Fluva_NFV_PCC),]
df_neg <- df_neg[df_neg$Fluva_NFV_PCC < 0 & df_neg$FLUVA_HNK_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(as.numeric(x[c("Fluva_NFV_PCC","FLUVA_NFV_PCC")]),na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]

df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_NFV_PCC) >= 0.1 & abs(df_neg$FLUVA_HNK_PCC) >= 0.1,][1:5,]
df_neg <- df_neg[!is.na(df_neg$Fluva_NFV_PCC),]

label <- df %>%
  summarise(
    Fluva_NFV_PCC = min(Fluva_NFV_PCC)+ 0.1,
    FLUVA_HNK_PCC = max(FLUVA_HNK_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )

ggplot(df) +
  geom_point(aes(Fluva_NFV_PCC, FLUVA_HNK_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_NFV gene associations") + ylab("FLUVA_HNK gene associations") +
  geom_text_repel(data = rbind.data.frame(df_pos,df_neg),mapping = aes(Fluva_NFV_PCC,FLUVA_HNK_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_pos,df_neg),aes(Fluva_NFV_PCC, FLUVA_HNK_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_NFV_PCC,y=FLUVA_HNK_PCC))

df <- df[order(df$fdr),]

dev.off()



# pathway PCL analysis

nPerm <- 10000
gsc1 <- loadGSC("data/pathwaysInfo/h.all.v6.2.symbols.gmt")

listOfPathways_enriched <- lapply(names(listOfAssociations), function(combo){
  
  sigAssoc_final_all <- listOfAssociations[[combo]]
  
  genesSymbols <- gene_mappings[rownames(sigAssoc_final_all),"Symbol"]
  ibx <- which(duplicated(genesSymbols))
  
  
  gene_stat_levels <- sigAssoc_final_all[-ibx,"estimate"]
  names(gene_stat_levels) <- genesSymbols[-ibx]
  
  gene_stat_levels <- gene_stat_levels[!is.na(gene_stat_levels)]
  
  gene_stat_levels <- sort(gene_stat_levels,decreasing = T)
  set.seed(3425)
  gseaRes <- piano::runGSA(geneLevelStats = gene_stat_levels,gsc = gsc1,adjMethod = "none",nPerm = nPerm,geneSetStat = "fgsea"
                           ,ncpus = 4)
  
  gseaResSummary <- piano::GSAsummaryTable(gseaRes)
  gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
  gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
  gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(nPerm+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
  gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr") # fdr correction
  gseares1 <- gseares1[,-(4:5)]
  
  gseares1 <- gseares1[order(gseares1$FDR,na.last = T),]
  gseares1 <- gseares1[order(gseares1$FDR,gseares1$pval,-abs(gseares1$`Stat (dist.dir)`),na.last = T),]
  

  gseares1Sig <- gseares1[gseares1$FDR<=0.05,]
  gseares1Sig <- gseares1Sig[order((gseares1Sig$`Stat (dist.dir)`)),]
  
 # pdf(paste(combo,"_Pathways_HALLMARK.pdf",sep = ""),height = 16,width = 20)
#  par(mai=c(1,8,1,1))
#  barplot(gseares1Sig$`Stat (dist.dir)`,horiz = T,names.arg = gseares1Sig$Name,las=2
#          ,main = paste("Pathways enriched in genes associated with response to",combo),xlab = "Enrichment score")
#  legend("topleft",legend = "FDR < 0.05",bty="n")
#  dev.off()
  
  return(gseares1)
})


names(listOfPathways_enriched) <- names(listOfAssociations)

listOfPathways_enriched_sig <- lapply(listOfPathways_enriched, function(x){
  gseares1Sig <- x[!is.na(x$FDR) & x$FDR<=0.05,]
  gseares1Sig <- gseares1Sig[order((gseares1Sig$`Stat (dist.dir)`)),]
  return(gseares1Sig)
})


unionPathways <- unique(unlist(lapply(listOfPathways_enriched_sig, function(x){return(x$Name)})))


df <- do.call(rbind,lapply(names(listOfPathways_enriched_sig), function(x){
  y <- listOfPathways_enriched_sig[[x]]
  
  otherPathways <- setdiff(unionPathways,y$Name)
  y <- y[,c("Name","Stat (dist.dir)","FDR")]
  colnames(y) <- c("Pathway","ES","FDR")
  y <- rbind(y,data.frame("Pathway"=otherPathways,"ES"=rep(NA,length(otherPathways)),"FDR"=rep(NA,length(otherPathways))))
  
  return(data.frame("Pathway"=y[,"Pathway"],"ES"=y[,"ES"],"FDR"=y[,"FDR"],"Combo"=rep(x,dim(y)[1]),stringsAsFactors = F))
}))




df_ES <- pivot_wider(df[,c("Pathway","ES","Combo")], names_from = Combo, values_from = ES)
df_FDR <- pivot_wider(df[,c("Pathway","FDR","Combo")], names_from = Combo, values_from = FDR)



rownames(gseares1Sig_Fluva_mono) <- gseares1Sig_Fluva_mono$Name

# Generate data
set.seed(12345);

dotmap.data <- cbind.data.frame("FLUVA"=gseares1Sig_Fluva_mono[unlist(lapply(df_ES$Pathway, as.character)),"Stat (dist.dir)"],df_ES[,c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")])
rownames(dotmap.data) <- df_ES$Pathway
bg.data <-  cbind.data.frame("FLUVA"=gseares1Sig_Fluva_mono[unlist(lapply(df_ES$Pathway, as.character)),"FDR"],df_FDR[,c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")])
rownames(bg.data) <- df_FDR$Pathway

dotmap.data <- dotmap.data[order(dotmap.data$FLUVA,na.last = T),]
bg.data <- bg.data[rownames(dotmap.data),]
spot.size.function <- function(x){
  0.1 + (4 * abs(x));
  }

spot.colour.function <- function(x){
  colours <- rep("white", length(x));
  colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1];
  colours[sign(x) ==  1] <- default.colours(2, palette.type = "dotmap")[2];
  return(colours);
  }

rownames(dotmap.data) <- gsub("_"," ",gsub("HALLMARK","",rownames(dotmap.data)))

ibx <- apply(dotmap.data, 1, function(x){ sum(!is.na(x))>1 & !is.na(x[1])})
ibx <- ibx[ibx]; ibx <- names(ibx)

idx <- apply(bg.data, 1, function(x){ sum(!is.na(x))>1 & !is.na(x[1])})
idx <- idx[idx]; idx <- names(idx)


##############
### Fig 4C ###
##############

pdf("Plots/Fig4C.pdf",width = 24,height = 8)
create.dotmap(x = t(dotmap.data[ibx,])
     , yaxis.cex = 1.5
     ,xaxis.cex = 1.5
     ,spot.size.function = spot.size.function
     ,spot.colour.function = spot.colour.function
     ,na.pch = 4,na.spot.size = 3
     ,key = list(space = "right",points = list(title="ES",cex = spot.size.function(seq(-.6,.6, 0.2)),col = spot.colour.function(seq(-.6,.6, 0.2)),pch = 19),text = list(lab = c("-0.6", "-0.4", "-0.2", " 0.0", "0.2","0.4", "0.6"),cex = 1.5,adj = 1.0,fontface = "bold")),
     # control spacing at top of key
     key.top = 1.5,
     # add borders to points
     pch = 21,
     pch.border.col = "white",
     # add the background
     bg.data = t(bg.data[idx,]),# add a colourkey
     colourkey = TRUE,
    # set colour scheme for background data
    colour.scheme = rev(c("white", "black")),
    # make bg colour scheme a discrete colour scheme, with breaks at these places
    at = c(0,0.0005,0.005,0.05,0.1),
    colourkey.labels = c("",0.0005,0.005,0.05,0.1),
    colourkey.labels.at = c(0,0.0005,0.005,0.05,0.1),xaxis.rot = 90,xaxis.tck = 0,yaxis.tck = 0);

dev.off()




#######################################

RPPA_UHN <- read.csv("data/CL_data/breast_rppa_UHN_NormLog2.tsv",sep = "\t",header = F,stringsAsFactors = F)


RPPA_UHN_final <- (RPPA_UHN[9:dim(RPPA_UHN)[1],11:dim(RPPA_UHN)[2]])

rownames(RPPA_UHN_final) <- RPPA_UHN[9:dim(RPPA_UHN)[1],7]

colnames(RPPA_UHN_final) <- RPPA_UHN[5,11:dim(RPPA_UHN)[2]]

RPPA_UHN_final <- as.matrix(RPPA_UHN_final)
class(RPPA_UHN_final) <- "numeric"


gene_names_RPPA <- as.character(RPPA_UHN[3,11:dim(RPPA_UHN)[2]])
names(gene_names_RPPA) <- colnames(RPPA_UHN_final)
source("R/fixCellLinesNames.R")

cells <- fixCellLinesNames(rownames(RPPA_UHN_final),"data/CL_data/cell_annotation_all.csv")


rownames(RPPA_UHN_final) <- cells


RPPA_UHN_mapping <- as.character(RPPA_UHN[3,11:dim(RPPA_UHN)[2]])
names(RPPA_UHN_mapping) <- RPPA_UHN[5,11:dim(RPPA_UHN)[2]]


commonSamples <- intersect(rownames(RPPA_UHN_final),colnames(BlissMat_summarized))


listOfAssociations_RPPA <- lapply(rownames(BlissMat_summarized), function(y){
  
  RPPA_Associations_cor <- apply(RPPA_UHN_final[commonSamples,], 2, function(x){
    a <- cor.test(x,BlissMat_summarized[y,commonSamples])
    return(c(a$estimate,a$p.value))
  })
  
  RPPA_Associations_cor <- t(RPPA_Associations_cor)
  RPPA_Associations_cor <- cbind(RPPA_Associations_cor,p.adjust(RPPA_Associations_cor[,2],method = "fdr"))
  colnames(RPPA_Associations_cor) <- c("estimate","pval","fdr")
  
  RPPA_Associations_cor <- RPPA_Associations_cor[order(RPPA_Associations_cor[,"fdr"]),]
  
})

names(listOfAssociations_RPPA) <- rownames(BlissMat_summarized)


##############
### Fig S6C ##
##############

pdf("Plots/FigS6C.pdf")
par(mai=c(2,2,2,1))
A <- cor.test(listOfAssociations_RPPA$FLUVA_DP[,"estimate"],listOfAssociations_RPPA$FLUVA_NFV[rownames(listOfAssociations_RPPA$FLUVA_DP),"estimate"],method = "s")

df <- data.frame("Fluva_DP_PCC"=listOfAssociations_RPPA$FLUVA_DP[,"estimate"], "FLUVA_NFV_PCC"=listOfAssociations_RPPA$FLUVA_NFV[rownames(listOfAssociations_RPPA$FLUVA_DP),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations_RPPA$FLUVA_DP[,"fdr"],listOfAssociations_RPPA$FLUVA_NFV[rownames(listOfAssociations_RPPA$FLUVA_DP),"fdr"]),MARGIN = 1,min))

df <- df[!is.na(df$fdr), ]
df$fdr2=ifelse(listOfAssociations_RPPA$FLUVA_DP[rownames(df),"fdr"]<0.05 & listOfAssociations_RPPA$FLUVA_NFV[rownames(df),"fdr"]<0.05,T,F)

library(scales)
library(ggrepel)
df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_names_RPPA[rownames(df)]

df$d = densCols(df$Fluva_DP_PCC, df$FLUVA_NFV_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))

df[df$Fluva_DP_fdr<0.05 & df$FLUVA_NFV_fdr<0.05,"d"] = "red"

df_pos <- df[order(df$Fluva_DP_PCC),]
df_pos <- df_pos[df_pos$Fluva_DP_PCC > 0 & df_pos$FLUVA_NFV_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(x[c("Fluva_DP_PCC","FLUVA_NFV_PCC")],na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]
df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_DP_PCC) >= 0.1 & abs(df_pos$FLUVA_NFV_PCC) >= 0.1,][1:5,]
#df_pos <- df_pos[df_pos$fdr2,]
df_pos <- df_pos[!is.na(df_pos$Fluva_DP_PCC),]

df_neg <- df[order(df$Fluva_DP_PCC),]
df_neg <- df_neg[df_neg$Fluva_DP_PCC < 0 & df_neg$FLUVA_NFV_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(x[c("Fluva_DP_PCC","FLUVA_NFV_PCC")],na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]
df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_DP_PCC) >= 0.1 & abs(df_neg$FLUVA_NFV_PCC) >= 0.1,][1:5,]
#df_neg <- df_neg[df_neg$fdr2,]
df_neg <- df_neg[!is.na(df_neg$Fluva_DP_PCC),]


label <- df %>%
  summarise(
    Fluva_DP_PCC = min(Fluva_DP_PCC)+ 0.05,
    FLUVA_NFV_PCC = max(FLUVA_NFV_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )
ggplot(df) +
  geom_point(aes(Fluva_DP_PCC, FLUVA_NFV_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(plot.margin=unit(c(1,1,1.1,1.1),"cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_DP protein associations") + ylab("FLUVA_NFV protein associations") +
  geom_text_repel(data = rbind.data.frame(df_neg,df_pos),mapping = aes(Fluva_DP_PCC,FLUVA_NFV_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_neg,df_pos),aes(Fluva_DP_PCC, FLUVA_NFV_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_DP_PCC,y=FLUVA_NFV_PCC))


A = cor.test(listOfAssociations_RPPA$FLUVA_DP[,"estimate"],listOfAssociations_RPPA$FLUVA_HNK[rownames(listOfAssociations_RPPA$FLUVA_DP),"estimate"],method = "s")

df <- data.frame("Fluva_DP_PCC"=listOfAssociations_RPPA$FLUVA_DP[,"estimate"], "FLUVA_HNK_PCC"=listOfAssociations_RPPA$FLUVA_HNK[rownames(listOfAssociations_RPPA$FLUVA_DP),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations_RPPA$FLUVA_DP[,"fdr"],listOfAssociations_RPPA$FLUVA_HNK[rownames(listOfAssociations_RPPA$FLUVA_DP),"fdr"]),MARGIN = 1,min))

df <- df[!is.na(df$fdr), ]
df$fdr2=ifelse(listOfAssociations_RPPA$FLUVA_DP[rownames(df),"fdr"]<0.05 & listOfAssociations_RPPA$FLUVA_HNK[rownames(df),"fdr"]<0.05,T,F)


df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_names_RPPA[rownames(df)]

df$d = densCols(df$Fluva_DP_PCC, df$FLUVA_HNK_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))
df[df$Fluva_DP_fdr<0.05 & df$FLUVA_HNK_fdr<0.05,"d"] = "red"

df_pos <- df[order(df$Fluva_DP_PCC),]
df_pos <- df_pos[df_pos$Fluva_DP_PCC > 0 & df_pos$FLUVA_HNK_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(x[c("Fluva_DP_PCC","FLUVA_HNK_PCC")],na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]
df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_DP_PCC) >= 0.1 & abs(df_pos$FLUVA_HNK_PCC) >= 0.1,][1:5,]
#df_pos <- df_pos[df_pos$fdr2,]
df_pos <- df_pos[!is.na(df_pos$Fluva_DP_PCC),]

df_neg <- df[order(df$Fluva_DP_PCC),]
df_neg <- df_neg[df_neg$Fluva_DP_PCC < 0 & df_neg$FLUVA_HNK_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(x[c("Fluva_DP_PCC","FLUVA_HNK_PCC")],na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]
df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_DP_PCC) >= 0.1 & abs(df_neg$FLUVA_HNK_PCC) >= 0.1,][1:5,]
#df_neg <- df_neg[df_neg$fdr2,]
df_neg <- df_neg[!is.na(df_neg$Fluva_DP_PCC),]


label <- df %>%
  summarise(
    Fluva_DP_PCC = min(Fluva_DP_PCC) + 0.05,
    FLUVA_HNK_PCC = max(FLUVA_HNK_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )
ggplot(df) +
  geom_point(aes(Fluva_DP_PCC, FLUVA_HNK_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(plot.margin=unit(c(1,1,1.1,1.1),"cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_DP protein associations") + ylab("FLUVA_HNK protein associations") +
  #geom_text_repel(data = df[!is.na(df$fdr) & df$fdr<0.05 & abs(df$Fluva_DP_PCC) >= 0.4 & abs(df$FLUVA_HNK_PCC) >= 0.4,],mapping = aes(Fluva_DP_PCC,FLUVA_HNK_PCC,label=Symbol),col="red") +
  geom_text_repel(data = rbind.data.frame(df_neg,df_pos),mapping = aes(Fluva_DP_PCC,FLUVA_HNK_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_neg,df_pos),aes(Fluva_DP_PCC, FLUVA_HNK_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_DP_PCC,y=FLUVA_HNK_PCC))





A <- cor.test(listOfAssociations_RPPA$FLUVA_HNK[,"estimate"],listOfAssociations_RPPA$FLUVA_NFV[rownames(listOfAssociations_RPPA$FLUVA_HNK),"estimate"],method = "s")

df <- data.frame("Fluva_NFV_PCC"=listOfAssociations_RPPA$FLUVA_NFV[,"estimate"], "FLUVA_HNK_PCC"=listOfAssociations_RPPA$FLUVA_HNK[rownames(listOfAssociations_RPPA$FLUVA_NFV),"estimate"]
                 ,"fdr"=apply(cbind(listOfAssociations_RPPA$FLUVA_NFV[,"fdr"],listOfAssociations_RPPA$FLUVA_HNK[rownames(listOfAssociations_RPPA$FLUVA_NFV),"fdr"]),MARGIN = 1,min))

df <- df[!is.na(df$fdr), ]
df$fdr2=ifelse(listOfAssociations_RPPA$FLUVA_NFV[rownames(df),"fdr"]<0.05 & listOfAssociations_RPPA$FLUVA_HNK[rownames(df),"fdr"]<0.05,T,F)



df$FDR = rescale(-log10(df$fdr),to = c(0,1))

df$Symbol <- gene_names_RPPA[rownames(df)]

df$d = densCols(df$Fluva_NFV_PCC, df$FLUVA_HNK_PCC, colramp = colorRampPalette(rev(c("black","darkgray"))))
df[df$Fluva_DP_fdr<0.05 & df$FLUVA_HNK_fdr<0.05,"d"] = "red"

df_pos <- df[order(df$Fluva_NFV_PCC),]
df_pos <- df_pos[df_pos$Fluva_NFV_PCC > 0 & df_pos$FLUVA_HNK_PCC > 0,]
df_pos$max <- apply(df_pos,1,function(x){max(x[c("Fluva_NFV_PCC","FLUVA_HNK_PCC")],na.rm = T)})
df_pos <- df_pos[order(as.numeric(df_pos$max),decreasing=T),]
df_pos <- df_pos[!is.na(df_pos$fdr) & df_pos$fdr<0.1 & abs(df_pos$Fluva_NFV_PCC) >= 0.1 & abs(df_pos$FLUVA_HNK_PCC) >= 0.1,][1:5,]
#df_pos <- df_pos[df_pos$fdr2,]
df_pos <- df_pos[!is.na(df_pos$Fluva_NFV_PCC),]

df_neg <- df[order(df$Fluva_NFV_PCC),]
df_neg <- df_neg[df_neg$Fluva_NFV_PCC < 0 & df_neg$FLUVA_HNK_PCC < 0,]
df_neg$max <- apply(df_neg,1,function(x){min(x[c("Fluva_NFV_PCC","FLUVA_HNK_PCC")],na.rm = T)})
df_neg <- df_neg[order(as.numeric(df_neg$max),decreasing=F),]
df_neg <- df_neg[!is.na(df_neg$fdr) & df_neg$fdr<0.1 & abs(df_neg$Fluva_NFV_PCC) >= 0.1 & abs(df_neg$FLUVA_HNK_PCC) >= 0.1,][1:5,]
#df_neg <- df_neg[df_neg$fdr2,]
df_neg <- df_neg[!is.na(df_neg$Fluva_NFV_PCC),]

label <- df %>%
  summarise(
    Fluva_NFV_PCC = min(Fluva_NFV_PCC)+ 0.05,
    FLUVA_HNK_PCC = max(FLUVA_HNK_PCC),
    label = paste("R: ",sprintf("%.2f",A$estimate))#,", P-val: ",sprintf("%.2E",A$p.value),sep = " ")
  )

ggplot(df) +
  geom_point(aes(Fluva_NFV_PCC, FLUVA_HNK_PCC, col = d), size = 1) +
  scale_color_identity() +
  theme_bw() + theme(plot.margin=unit(c(1,1,1.1,1.1),"cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12),axis.title = element_text(size = 14,face = "bold")) + xlab("FLUVA_NFV protein associations") + ylab("FLUVA_HNK protein associations") +
  geom_text_repel(data = rbind.data.frame(df_neg,df_pos),mapping = aes(Fluva_NFV_PCC,FLUVA_HNK_PCC,label=Symbol),col="red") +
  geom_point(data = rbind.data.frame(df_neg,df_pos),aes(Fluva_NFV_PCC, FLUVA_HNK_PCC), col="red") +
  geom_text(data = label,mapping = aes(label=label,x = Fluva_NFV_PCC,y=FLUVA_HNK_PCC))




dev.off()




##################################
##################################
##################################


gsc1 <- loadGSC("data/pathwaysInfo/EMT_genesets_all.gmt")

listOfPathways_enriched_EMT <- lapply(names(listOfAssociations), function(combo){
  
  sigAssoc_final_all <- listOfAssociations[[combo]]
  
  genesSymbols <- gene_mappings[rownames(sigAssoc_final_all),"Symbol"]
  ibx <- which(duplicated(genesSymbols))
  
  
  gene_stat_levels <- sigAssoc_final_all[-ibx,"estimate"]
  names(gene_stat_levels) <- genesSymbols[-ibx]
  
  gene_stat_levels <- gene_stat_levels[!is.na(gene_stat_levels)]
  
  gene_stat_levels <- sort(gene_stat_levels,decreasing = T)
  set.seed(3425)
  gseaRes <- piano::runGSA(geneLevelStats = gene_stat_levels,gsc = gsc1,adjMethod = "none",nPerm = nPerm,geneSetStat = "fgsea"
                           ,ncpus = 4)
  
  gseaResSummary <- piano::GSAsummaryTable(gseaRes)
  gseaResSummary <- piano::GSAsummaryTable(gseaRes)
  if(length(colnames(gseaResSummary))<8){
    if("p (dist.dir.up)" %in% colnames(gseaResSummary) & !("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
      gseaResSummary[,"p (dist.dir.dn)"] <- NA
    }else if(!("p (dist.dir.up)" %in% colnames(gseaResSummary)) & ("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
      gseaResSummary[,"p (dist.dir.up)"] <- NA
    }
  }
  gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
  gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
  gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(nPerm+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
  gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr") # fdr correction
  gseares1 <- gseares1[,-(4:5)]
  
  gseares1 <- gseares1[order(gseares1$FDR,na.last = T),]
  gseares1 <- gseares1[order(gseares1$FDR,gseares1$pval,-abs(gseares1$`Stat (dist.dir)`),na.last = T),]
  
  #write.table(gseares1,paste(directory,drug,"_pathways_Analysis_HALLMARK.csv"),quote = F,sep = ",",row.names = T,col.names = NA)
  
  #gseares1 <- read.csv2(file = paste(directory,drug,"_pathways_Analysis_HALLMARK.csv"),header = T,row.names = 1,sep = ",",quote = F)
  
  gseares1Sig <- gseares1[gseares1$FDR<=0.05,]
  gseares1Sig <- gseares1Sig[order((gseares1Sig$`Stat (dist.dir)`)),]
  
#  pdf(paste(combo,"_Pathways_HALLMARK.pdf",sep = ""),height = 16,width = 20)
#  par(mai=c(1,8,1,1))
#  barplot(gseares1Sig$`Stat (dist.dir)`,horiz = T,names.arg = gseares1Sig$Name,las=2
#          ,main = paste("Pathways enriched in genes associated with response to",combo),xlab = "Enrichment score")
#  legend("topleft",legend = "FDR < 0.05",bty="n")
#  dev.off()
  
  return(gseares1)
})

names(listOfPathways_enriched_EMT) <- names(listOfAssociations)


tmp <- unique(unlist(lapply(listOfPathways_enriched_EMT, function(x){
  return(x[x$FDR<0.05,"Name"])
})))



keepPathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN",
  "GO_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GO_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
  "GO_EPITHELIAL_MESENCHYMAL_CELL_SIGNALING",
  "EPITHELIAL_TO_MESENCHYMAL_TRANSITION" 
)

listOfPathways_enriched_EMT_final <- lapply(listOfPathways_enriched_EMT, function(x){
  tmp <- x[match(c(keepPathways),x[,"Name"]),]
  tmp$FDR <- p.adjust(tmp$pval)
  return(tmp)
})






df_ES <- do.call(cbind,lapply(listOfPathways_enriched_EMT_final, function(x){
  rownames(x) <- x$Name
  x[x$FDR>=0.05,"Stat (dist.dir)"] <- NA
  return(x[,"Stat (dist.dir)",drop=F])
}))
colnames(df_ES) <- names(listOfPathways_enriched_EMT_final)

df_FDR <- do.call(cbind,lapply(listOfPathways_enriched_EMT_final, function(x){
  rownames(x) <- x$Name
  x[x$FDR>=0.05,"FDR"] <- NA
  return(x[,"FDR",drop=F])
}))
colnames(df_FDR) <- names(listOfPathways_enriched_EMT_final)




genesSymbols <- gene_mappings[rownames(geneAssociations_cor_Mono),"Symbol"]
ibx <- which(duplicated(genesSymbols))


gene_stat_levels <-  -geneAssociations_cor_Mono[-ibx,"estimate"]
names(gene_stat_levels) <- genesSymbols[-ibx]

gene_stat_levels <- gene_stat_levels[!is.na(gene_stat_levels)]

set.seed(3425)



nPerm <- 10000

gseaRes <- piano::runGSA(geneLevelStats = gene_stat_levels,gsc = gsc1,adjMethod = "none",nPerm = nPerm,geneSetStat = "fgsea"
                         ,ncpus = 4)

gseaResSummary <- piano::GSAsummaryTable(gseaRes)
if(length(colnames(gseaResSummary))<8){
  if("p (dist.dir.up)" %in% colnames(gseaResSummary) & !("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
    gseaResSummary[,"p (dist.dir.dn)"] <- NA
  }else if(!("p (dist.dir.up)" %in% colnames(gseaResSummary)) & ("p (dist.dir.dn)" %in% colnames(gseaResSummary))){
    gseaResSummary[,"p (dist.dir.up)"] <- NA
  }
}
gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(nPerm+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr") # fdr correction
gseares1 <- gseares1[,-(4:5)]

gseares1 <- gseares1[order(gseares1$FDR,na.last = T),]
gseares1 <- gseares1[order(gseares1$FDR,gseares1$pval,-abs(gseares1$`Stat (dist.dir)`),na.last = T),]

gseares1Sig <- gseares1
gseares1Sig[gseares1Sig$FDR>=0.05,c("Stat (dist.dir)","FDR")] <- c(NA,NA)
gseares1Sig <- gseares1Sig[order(abs(gseares1Sig$`Stat (dist.dir)`),decreasing=T),]
gseares1Sig_Fluva_mono_EMT <- gseares1Sig




# Generate data
set.seed(12345);
rownames(gseares1Sig_Fluva_mono_EMT) <- gseares1Sig_Fluva_mono_EMT$Name

dotmap.data <- df_ES[,c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")] 
dotmap.data <- cbind.data.frame("FLUVA"=gseares1Sig_Fluva_mono_EMT[rownames(df_ES),"Stat (dist.dir)"],df_ES[,c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")])
rownames(dotmap.data) <- rownames(df_ES)
bg.data <-  df_FDR[,c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")]
bg.data <-  cbind.data.frame("FLUVA"=gseares1Sig_Fluva_mono_EMT[rownames(df_ES),"FDR"],df_FDR[rownames(df_ES),c("FLUVA_DP","FLUVA_NFV","FLUVA_HNK")])
rownames(bg.data) <- rownames(df_ES)

dotmap.data <- dotmap.data[order(dotmap.data$FLUVA_DP,na.last = T),]
bg.data <- bg.data[rownames(dotmap.data),]
spot.size.function <- function(x){
  0.1 + (4 * abs(x));
}

spot.colour.function <- function(x){
  colours <- rep("white", length(x));
  colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1];
  colours[sign(x) ==  1] <- default.colours(2, palette.type = "dotmap")[2];
  return(colours);
}

rownames(dotmap.data) <- gsub("_"," ",rownames(dotmap.data))

dotmap.data <- dotmap.data[-c(6,7),]
bg.data <- bg.data[-c(6,7),]


##############
### Fig S6D ##
##############

pdf("Plots/FigS6D.pdf",width = 7,height = 7)
create.dotmap(x = t(dotmap.data)
              , yaxis.cex = 1.5
              ,xaxis.cex = 0.6
              ,spot.size.function = spot.size.function
              ,spot.colour.function = spot.colour.function
              ,na.pch = 4,na.spot.size = 3
              ,key = list(space = "right",points = list(title="ES",cex = spot.size.function(seq(-.6,.6, 0.2)),col = spot.colour.function(seq(-.6,.6, 0.2)),pch = 19),text = list(lab = c("-0.6", "-0.4", "-0.2", " 0.0", "0.2","0.4", "0.6"),cex = 1.5,adj = 1.0,fontface = "bold")),
              # control spacing at top of key
              key.top = 1.5,
              # add borders to points
              pch = 21,
              pch.border.col = "white",
              # add the background
              bg.data = t(bg.data),# add a colourkey
              colourkey = TRUE,
              # set colour scheme for background data
              colour.scheme = rev(c("white", "black")),
              # make bg colour scheme a discrete colour scheme, with breaks at these places
              at = c(0,0.0005,0.005,0.05,0.5),
              colourkey.labels = c("",0.0005,0.005,0.05,0.1),
              colourkey.labels.at = c(0,0.0005,0.005,0.05,0.1),xaxis.rot = 90,xaxis.tck = 0,yaxis.tck = 0);

dev.off()






EMT_sets_genes <- as.matrix(read.csv(file = "data/pathwaysInfo/EMT_sets_genes.csv",header = T,sep = "\t",stringsAsFactors = F,na.strings = NA))

EMT_sets_genes[EMT_sets_genes==""] <- NA

EMT_sets_genes_union <- unique(as.character(EMT_sets_genes))
EMT_sets_genes_union <- EMT_sets_genes_union[-which(is.na(EMT_sets_genes_union))]

data <- matrix(data=NA,nrow = length(EMT_sets_genes_union)
               ,ncol = dim(EMT_sets_genes)[2]
               ,dimnames = list(EMT_sets_genes_union,colnames(EMT_sets_genes)))


for (gene in rownames(data)) {
  for (pathway in colnames(data)) {
    if(gene %in% EMT_sets_genes[,pathway])
      data[gene,pathway] <- 1
    else
      data[gene,pathway] <- 0
  }
}

data <- as.data.frame(data)

data <- data[,names(sort(apply(data, 2, sum),decreasing = T))]


##############
### Fig S7 ###
##############

pdf("Plots/FigS7.pdf",width = 12,height = 8)
upset(data = data[,c(1,2,3,4,5,8)],nsets = 6,order.by = "freq")
dev.off()






