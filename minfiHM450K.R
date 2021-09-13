library(data.table)
library(tidyverse)
library(minfi)
library(EnhancedVolcano)

#Download HM450K hg19 annotation file on https://zwdzwd.github.io/InfiniumAnnotation
hg19annotaion <- read.csv("hg19_HumanMethylation450_annotation.csv", header = TRUE)
hg19annotaion_CGI <- hg19annotaion[,c(1,25)]

#Download HM450K data from TCGA using manifest in Supplementary Table 5
files = list.files(path = "Data/.", pattern = "hg38.txt$", recursive = TRUE) 

# read data
for(f in files){
  print(f)
  temp = fread(file = paste("Data/",f, sep=""))
  colnames(temp)[2] = strsplit(f, "/")[[1]][2]
  if(!exists("tcga_betas")){
    tcga_betas = temp
  } else{
    tcga_betas = merge(tcga_betas, temp[,1:2], by="Composite Element REF")
  }
}

#beta value data.frame
tcga_betas_matrix <- tcga_betas[,-c(3:11)]
tcga_betas_matrix <- data.frame(tcga_betas_matrix, check.names = F)#without check.names = F, column name will change
rownames(tcga_betas_matrix) <- tcga_betas_matrix[,1];tcga_betas_matrix <- tcga_betas_matrix[,-1]
colname <- data.frame(colnames(tcga_betas_matrix))
sample_sheet = read.csv("gdc_sample_sheet.2021-03-15_filtered.txt",sep = "\t")
condition <- merge(colname, sample_sheet[,c(2,8)],by.x = "colnames.tcga_betas_matrix.", by.y = "File.Name",sort=FALSE)
condition$newname<-gsub('-','.',condition$colnames.tcga_betas_matrix.)

#minfi Primary tumor vs Normal adjacent tissue
Type <- condition$Sample.Type
tcga_betas_matrix <- data.matrix(tcga_betas_matrix)
dmp <- dmpFinder(tcga_betas_matrix, pheno = Type  , type = "continuous")
dmp$beta <- dmp$beta*(-1)
FDR <- dmp %>% filter(qval < 0.1)
DMC <- subset(FDR, beta >0.2)

##DMC annotation
Annotated_DMC <- merge(DMC, hg19annotaion_CGI,by.x=0,by.y="IlmnID")
Annotated_DMC <- subset(Annotated_DMC, UCSC_CpG_Islands_Name!= "")#DMCs on CpG islands
Annotated_DMC <- Annotated_DMC %>% separate(UCSC_CpG_Islands_Name, c("chr", "start", "end"))
Annotated_DMC <- transform(Annotated_DMC, start = as.numeric(start))
Annotated_DMC <- transform(Annotated_DMC, end = as.numeric(end))
Annotated_DMC$start <- Annotated_DMC$start+1
Annotated_DMC <-unite(Annotated_DMC, "position", c("chr","start"), sep = ":", remove = TRUE)
Annotated_DMC <-unite(Annotated_DMC, "position", c("position","end"), sep = "-", remove = TRUE)#20008
write.table(Annotated_DMC,"Annotated_DMC_COAD_q0.1_delta0.2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
Annotated_DMC_CGI <- Annotated_DMC[!duplicated(Annotated_DMC$position), ]#unique position#4630
write.table(Annotated_DMC_CGI,"Annotated_DMC_CGI_COAD_q0.1_delta0.2.txt", quote = FALSE, row.names = FALSE, sep = "\t")

##minfi DMC dmpFinder -----Primary tumor vs PBMCs
#PBMCs data can be found on GSE53045
PBMC_data <- read.table("PBMCs61_df.txt", sep = "\t", header = TRUE)
rownames(PBMC_data)<-PBMC_data[,1];PBMC_data<-PBMC_data[,-1]
tcga_betas_matrix <- data.frame(tcga_betas_matrix)
include<-subset(condition,Sample.Type=="Primary Tumor")[3]
include <- tcga_betas_matrix[,include$newname]
PBMC_df <- cbind(include,PBMC_data)
PBMC_condition <- factor(c(rep("COAD",ncol(include)), rep("PBMC",ncol(PBMC_data))))
PBMC_df <- data.matrix(PBMC_df)
PBMC_dmp <- dmpFinder(PBMC_df, pheno = PBMC_condition  , type = "continuous")
PBMC_dmp$beta <- PBMC_dmp$beta*(-1)
PBMC_FDR <- subset(PBMC_dmp, qval < 0.1)
PBMC_DMC <- subset(PBMC_FDR, beta > 0.2)
Annotated_PBMC_DMC <- merge(PBMC_DMC, hg19annotaion_CGI,by.x=0,by.y="IlmnID")
Annotated_PBMC_DMC <- subset(Annotated_PBMC_DMC, UCSC_CpG_Islands_Name!= "")#PBMC_DMCs on CpG islands
Annotated_PBMC_DMC <- Annotated_PBMC_DMC %>% separate(UCSC_CpG_Islands_Name, c("chr", "start", "end"))
Annotated_PBMC_DMC <- transform(Annotated_PBMC_DMC, start = as.numeric(start))
Annotated_PBMC_DMC <- transform(Annotated_PBMC_DMC, end = as.numeric(end))
Annotated_PBMC_DMC$start <- Annotated_PBMC_DMC$start+1
Annotated_PBMC_DMC <-unite(Annotated_PBMC_DMC, "position", c("chr","start"), sep = ":", remove = TRUE)
Annotated_PBMC_DMC <-unite(Annotated_PBMC_DMC, "position", c("position","end"), sep = "-", remove = TRUE)
write.table(Annotated_PBMC_DMC,"Annotated_DMC_COADvsPBMC_q0.1_delta0.2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
Annotated_PBMC_DMC_CGI <- Annotated_PBMC_DMC[!duplicated(Annotated_PBMC_DMC$position), ]
write.table(Annotated_PBMC_DMC_CGI,"Annotated_DMC_CGI_COADvsPBMC_q0.1_delta0.2.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Volcano plot
png(file="Volcano-COAD vs normal (DMCs).png",width=500, height=500)
EnhancedVolcano(dmp, lab = NA, x = "beta", y = "qval", xlab = "Delta beta value", ylab = expression('-Log'[10]*' q value'), title = "COAD primary tumor vs normal tissue", subtitle = "",
                pCutoff = 0.1, FCcutoff = 0.2, xlim = c(-1,1), ylim = c(0,30), pointSize = 1, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
png(file="Volcano-COAD vs PBMC (DMCs).png",width=500, height=500)
EnhancedVolcano(PBMC_dmp, lab = NA, x = "beta", y = "qval", xlab = "Delta beta value", ylab = expression('-Log'[10]*' q value'), title = "COAD primary tumor vs PBMC", subtitle = "",
                pCutoff = 0.1, FCcutoff = 0.2, xlim = c(-1,1), ylim = c(0,80), pointSize = 1, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
