# ============================================================
# GWAS using Sommer package
# ============================================================

# --- Setup ------------------------------------------------------------------
setwd()

# --- Load Required Packages -------------------------------------------------
# Set CRAN mirror
#options(repos = c(CRAN = "https://cloud.r-project.org/"))

#install.packages("sommer")
library(sommer)
#install.packages("data.table")
library(data.table)

# --- Data -------------------------------------------------------------------
#myY<-read.table("../../data/Cimmyt/PH_LODG.txt",sep =" ",header = T)
myY<-read.csv("../../data/Cimmyt/Qtraits.csv")
colnames(myY)[1]<-"FID"
#myGD<-read.table("../new_data/myGD_PH.hmp.txt",sep = ",",header = T)
myGD<-read.table("../new_data/myGD_TKW.hmp.txt",sep = ",",header = T)
myGD<-as.data.frame(myGD)
taxa<-myGD$taxa
rownames(myGD)<-taxa
myGD<-myGD[,-1]


relevant<-which(myY$FID %in% rownames(myGD))
myY<-myY[relevant,]

#myGM<-read.csv("../Gapit_PH/myGM.csv")
myGM<-read.csv("../new_data/myGM_TKW.csv")

myGD_matrix <- as.matrix(myGD)
# Convert 0 -> -1, 1 -> 0, and 2 -> 1
myGD_matrix <- myGD_matrix - 1
A <- A.mat(myGD_matrix) # additive relationship matrix
 
colnames(myGM)[1]<-"Locus"
colnames(myGM)[2]<-"Chrom"
# ============================================================
# 1. GWAS
# ============================================================
Sys.time()
mix1 <- GWAS(TKW~1,
             random=~vsr(FID,Gu=A),
             rcov=~units,
             data=myY, nIters=15,n.PC = 2,
             M=myGD_matrix, gTerm = "u:FID",
             verbose = FALSE)

Sys.time()


# ============================================================
# 2. Save Result
# ============================================================
#save.image("Lodg.RData")


ms <- as.data.frame(mix1$scores)
ms$Locus <- rownames(ms)
MP2 <- merge(myGM,ms,by="Locus",all.x = TRUE);
names(MP2)[names(MP2) == "Chromosome"] <- "Chrom"
names(MP2)[names(MP2) == "TKW"] <- "p.val"
MP2$p.val[MP2$p.val==Inf] <- NA

colnames(MP2)[1]<-"SNP"
colnames(MP2)[2]<-"Chrom"
#colnames(MP2)[3]<-"Position"
#colnames(MP2)[4]<-"p.val"


manhattan(MP2)

write.csv(MP2, "Sommer.TKW2.csv",row.names = F)

