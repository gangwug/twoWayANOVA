#Analysis explanation
#Two-way ANOVA for each gene, having as factors the group and the time.
#The data file has expression values for each gene in different time points, for two groups.
#Afer applying ANOVA, results are then corrected for multiple comparisons (multiple genes).

##Loading packages
rm(list=ls())
library(dplyr)
library(tidyr)

##Reading the raw data
totalData<-read.csv("dataTable.csv")

##General process
#For each gene, data will be read from the original table, stored in a temporary input
#object, run through the ANOVA and the results will be stored in a final output ANOVA table.
#The process will be repeated for all genes.

##PREPARING THE TEMPORARY INPUT TABLE
groups<-c("g1","g2")
nGroups<-2
nSamples<-9  #number samples collected for each group
nReplicates<-3 #number of replicate samples for each time point in each group
#Object to store the data for a single gene in every iteration of the ANOVA
dataTr<-data.frame(group=rep(groups,each=nSamples),
                   time=rep(rep(c("t1","t2","t3"),each=nReplicates),nGroups),
                   expression=rep(NA,nSamples*nGroups))
##PREPARING THE OUTPUT TABLE
#Object to store two-way ANOVA result
nGenes<-length(totalData$CycID)
anova2way_data<-data.frame(
  CydID=totalData$CycID,
  Fgroup=rep(NA,nGenes),
  Ftime=rep(NA,nGenes),
  FtimeXgroup=rep(NA,nGenes),
  Pgroup=rep(NA,nGenes),
  Ptime=rep(NA,nGenes),
  PtimeXgroup=rep(NA,nGenes)
)

##Performing two-way ANOVA for each gene, using a for loop
for(i in seq(nGenes)){   #For each gene
  dataCols<-(2:19)   #columns in the original table that have expression data
  #Gets the data in the original table and stores in the input table dataTr.
  dataTr$expression <- unlist(totalData[i,dataCols])
  #Performs two-way ANOVA for the given gene
  d_obj <- lm(expression ~ group*time , data=dataTr)
  d_anova <- anova(d_obj)
  #Stores the ANOVA result in the output object
  anova2way_data$Fgroup[i] <- d_anova$`F value`[1]
  anova2way_data$Ftime[i] <- d_anova$`F value`[2]
  anova2way_data$FtimeXgroup[i] <- d_anova$`F value`[3]
  anova2way_data$Pgroup[i] <- d_anova$`Pr(>F)`[1]
  anova2way_data$Ptime[i] <-  d_anova$`Pr(>F)`[2]
  anova2way_data$PtimeXgroup[i] <-  d_anova$`Pr(>F)`[3]
}

##Correcting the p-values with BH (fdr) method
for(i in 5:7){
  param <- names(anova2way_data)[i]
  p <- anova2way_data[,i]
  p_fdr <- p.adjust(p, method = "fdr")
  #Storing the corrected p-value back in the ANOVA object, in a new column
  colN<-length(anova2way_data)+1
  anova2way_data[,colN]<-p_fdr
  names(anova2way_data)[colN]<-paste0(param,"FDR")
}




