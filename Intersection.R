library(data.table)
library(tidyverse)
#Intersection 
Colorectal_Tissue_DMC <- read.table("Annotated_DMC_CGI_COAD_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Colorectal_DMR <- read.table("CpGislands_DMRs_Colorectal.txt", sep = "\t", header = TRUE)
Colorectal_DMR <- subset(Colorectal_DMR, log2FoldChange >=1)
Colorectal_Intersection <- merge(Colorectal_DMR, Colorectal_Tissue_DMC,  by.x = "X", by.y = "position")
Colorectal_PBMC_DMC <- read.table("Annotated_DMC_CGI_COADvsPBMC_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Colorectal_Intersection <- merge(Colorectal_Intersection, Colorectal_PBMC_DMC,  by.x = "X", by.y = "position")
write.table(Colorectal_Intersection, "Tissue_Colorectal_Intersection.txt",row.names = FALSE,sep = "\t",quote = FALSE)
table(Colorectal_DMR$X %in% Colorectal_PBMC_DMC$position)
table(Colorectal_Tissue_DMC$position %in% Colorectal_PBMC_DMC$position)
Lung_Tissue_DMC <- read.table("Annotated_DMC_CGI_LUAD_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Lung_DMR <- read.table("CpGislands_DMRs_Lung.txt", sep = "\t", header = TRUE)
Lung_DMR <- subset(Lung_DMR, log2FoldChange >=1)
Lung_Intersection <- merge(Lung_DMR, Lung_Tissue_DMC,  by.x = "X", by.y = "position")
Lung_PBMC_DMC <- read.table("Annotated_DMC_CGI_LUADvsPBMC_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Lung_Intersection <- merge(Lung_Intersection, Lung_PBMC_DMC,  by.x = "X", by.y = "position")
write.table(Lung_Intersection, "Tissue_Lung_Intersection.txt",row.names = FALSE,sep = "\t",quote = FALSE)
table(Lung_DMR$X %in% Lung_PBMC_DMC$position)
table(Lung_Tissue_DMC$position %in% Lung_PBMC_DMC$position)
Pancreas_Tissue_DMC <- read.table("Annotated_DMC_CGI_PAAD_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Pancreas_DMR <- read.table("CpGislands_DMRs_Pancreas.txt", sep = "\t", header = TRUE)
Pancreas_DMR <- subset(Pancreas_DMR, log2FoldChange >=1)
Pancreas_Intersection <- merge(Pancreas_DMR, Pancreas_Tissue_DMC,  by.x = "X", by.y = "position")
Pancreas_PBMC_DMC <- read.table("Annotated_DMC_CGI_PAADvsPBMC_q0.1_delta0.2.txt", sep = "\t", header = TRUE)
Pancreas_Intersection <- merge(Pancreas_Intersection, Pancreas_PBMC_DMC,  by.x = "X", by.y = "position")
write.table(Pancreas_Intersection, "Tissue_Pancreas_Intersection.txt",row.names = FALSE,sep = "\t",quote = FALSE)
table(Pancreas_DMR$X %in% Pancreas_PBMC_DMC$position)
table(Pancreas_Tissue_DMC$position %in% Pancreas_PBMC_DMC$position)
#Overlap 
Merge_CL <- merge(Colorectal_Intersection, Lung_Intersection, by = "X")
Merge_CP <- merge(Colorectal_Intersection, Pancreas_Intersection, by = "X")
Merge_LP <- merge(Lung_Intersection, Pancreas_Intersection, by = "X")
Merge_CLP <- merge(Merge_CL, Pancreas_Intersection, by = "X")
#Tissue specific
Colorectal_specific <- Colorectal_Intersection[!(Colorectal_Intersection$X %in% Merge_CL$X),]
Colorectal_specific <- Colorectal_specific[!(Colorectal_specific$X %in% Merge_CP$X),]
write.table(Colorectal_specific, "Tissue_Colorectal_specific.txt",row.names = FALSE,sep = "\t",quote = FALSE)
Lung_specific <- Lung_Intersection[!(Lung_Intersection$X %in% Merge_CL$X),]
Lung_specific <- Lung_specific[!(Lung_specific$X %in% Merge_LP$X),]
write.table(Lung_specific, "Tissue_Lung_specific.txt",row.names = FALSE,sep = "\t",quote = FALSE)
Pancreas_specific <- Pancreas_Intersection[!(Pancreas_Intersection$X %in% Merge_LP$X),]
Pancreas_specific <- Pancreas_specific[!(Pancreas_specific$X %in% Merge_CP$X),]
write.table(Pancreas_specific, "Tissue_Pancreas_specific.txt",row.names = FALSE,sep = "\t",quote = FALSE)

# Veen diagram
library(eulerr)
library(Rcpp)
fit1 <- euler(c("A" = 1783, "B" = 2588, "C" = 4906,
               "A&B" = 994, "B&C" = 1994, "A&C" = 1340,
               "A&B&C" = 939))
png(file="Venn different technique.png",width=500, height=500)
plot(fit1, fills = list(fill = c("firebrick1", "cyan","lawngreen")),labels = c("cfMBD-seq\nCase vs Control","HM450K\nTumor vs Normal tissue","HM450K\nTumor tissue vs Normal blood cell"))
dev.off()

fit2 <- euler(c("L" = 939, "C" = 1486, "P" = 896,
               "C&L" = 425, "C&P" = 589, "L&P" = 410,
               "C&L&P" = 266))
png(file="Venn cancer specific DMCGIs.png",width=300, height=300)
plot(fit2, fills = list(fill = c("red", "green","cyan")),labels = c("Lung\n370","Colorectal\n738","Pancreatic\n163"))
dev.off()
