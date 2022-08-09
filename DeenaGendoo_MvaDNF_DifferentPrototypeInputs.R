# Deena M.A. Gendoo
# Updated March 2022 
# To demonstrate the robustness of the pipeline, check what drugs will be discovered by using nelfinavir or honokiol as the input 
# Does iteratively applying the MVA-DNF pipeline to the discovered drugs converge to better drugs?
################################################################################################################################################

#Load the Run conducted for DP, NFV, and NNK
load("Data/TopHits_MvaDNF_DifferentPrototypeInputs.RData")

#Collate all stats scores for the top drugs that were tested, in addition to DP
AllDrugs<-rbind(Nelfinavir,Honokiol,CleanDPlist)
AllDrugs$RanksPVal<-as.numeric(AllDrugs$RanksPVal)
AllDrugs$RanksZScore<-as.numeric(AllDrugs$RanksZScore)

#Keep only hits that are p<=0.05 and z-score <-1.8
AllDrugs<-AllDrugs[AllDrugs$RanksPVal<=0.05,]
AllDrugs<-AllDrugs[AllDrugs$RanksZScore<=(-1.8),]

#hits that are p<=0.05 AND now have z-score cutoff
table(AllDrugs$DrugAnalyzed) 
# DIPYRIDAMOLE     HONOKIOL   NELFINAVIR 
#         19           20           12 

##### SHOW OVERLAPS BETWEEN TOP HITS USING DIFFERENT DRUGS AS INPUT ######
Lister<-AllDrugs[,1:2]
BinMat<-matrix(nrow=length(unique(Lister$DrugAnalyzed)),ncol=length(unique(Lister$MVAmatch)),data = 0)
colnames(BinMat)<-unique(Lister$MVAmatch)
rownames(BinMat)<-unique(Lister$DrugAnalyzed)

for(sample in 1:nrow(Lister))
{
  Source<-Lister$DrugAnalyzed[sample]
  Target<-Lister$MVAmatch[sample]
  
  BinMat[Source,Target]<-1
}

BinMat<-as.data.frame(BinMat)
Binny<-as.data.frame(t(BinMat))

# PLOT overlaps using Upset Plots

library(ComplexHeatmap)

pdf("./REVISION_DF_Output/MvaDNF_DrugHits_Overlaps_March2022.pdf",width = 15,height = 5)
m3 = ComplexHeatmap::make_comb_mat(Binny, mode = "intersect")
UpSet(m3,set_order = c("DIPYRIDAMOLE", "NELFINAVIR", "HONOKIOL"),
      comb_order = order(comb_size(m3),decreasing = T))
m3= m3[comb_degree(m3) > 0] #Keep only combinations where intersections>0!
ht = draw(UpSet(m3))
od = column_order(ht)
cs = comb_size(m3)
UpSet(m3, pt_size = unit(5, "mm"), lwd = 3,
      comb_col = c("red", "blue", "black","yellow","green","purple","orange")[comb_degree(m3)])
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 10))})

dev.off()


# Check what drugs are overlapping in the intersects!
m3 #List of all intersects
extract_comb(m3, "111") #Intersection of NFV, HNK, and DP
extract_comb(m3, "011") #Intersection of HNK and DP
extract_comb(m3, "101") #Intersection of NFV and DP
extract_comb(m3, "110") #Intersection of NFV and HNK
extract_comb(m3, "001") #DP-like drugs
extract_comb(m3, "010") #Honokiol-like drugs
extract_comb(m3, "100") #Nelfinavir-like drugs
