# Deena M.A. Gendoo
# March-7-2022
# Boxplots of time and memory complexity for MVA-DNF

################################################################################################################################################

# pdf("TimeMemory_Calculations.pdf",onefile = T)
Time<-read.csv("Data/TimeMemory_Time.csv",row.names = 1,header = T)
colnames(Time)<-c("6_Genes","10_Genes","100_Genes","500_Genes","900_Genes")
Time2<-Time/60
boxplot(Time,ylim=c(1200,1800),col = "white",outline = F,xlab="Gene size",ylab="Runtime (seconds)")
# Points
stripchart(Time,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)     

# boxplot(Time2,ylim=c(20,30),col = "white",outline = F,xlab="Gene size",ylab="Runtime (minutes)")
# # Points
# stripchart(Time2,              # Data
#            method = "jitter", # Random noise
#            pch = 19,          # Pch symbols
#            col = 4,           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)     
# 

MemoryTest<-read.csv("Data/TimeMemory_Memory.csv",row.names = 1,header = T)
colnames(MemoryTest)<-c("6_Genes","10_Genes","100_Genes","500_Genes","900_Genes")
boxplot(MemoryTest,ylim=c(510,560),col = "white",outline = F,
        xlab="Gene size",ylab="Memory (megabytes)")
# Points
stripchart(MemoryTest,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 2,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)     

# dev.off()
