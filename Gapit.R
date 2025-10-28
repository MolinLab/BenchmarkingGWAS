setwd()

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))
library(data.table)
library(readr)
# GAPIT source
#devtools::install_github("https://github.com/SFUStatgen/LDheatmap")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#############################################################################
myY<-read.csv("../../data/Cimmyt/Qtraits.csv")
myG<-read.delim("../new_data/geno_Qtraits_filtered_40_imputed_KNN_with_accuracy2.hmp.txt",head = F)

myGAPIT=GAPIT(#SNP.P3D =True
  Y=myY[,c(1,5)], 
  G=myG,
  #KI=myKI,
  #CV=myCV,
  PCA.total=2,
  model=c("GLM","MLM","MLMM","SUPER","FarmCPU","BLINK","CMLM"),
  Multiple_analysis=T)
