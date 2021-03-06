library(DESeq2)
library(tidyverse)
library(data.table)

#CpG_Islands & CpG_Islands_extention
CpG_Islands <- filter(df, Type == "island")
CpG_Islands <- CpG_Islands[,-c(4:6)]
CpG_Islands <-unite(CpG_Islands, "position",c("Chr","Start"),sep = ":",remove = TRUE)
CpG_Islands <-unite(CpG_Islands, "position",c("position","End"),sep = "-",remove = TRUE)
CpG_Islands <- data.frame(CpG_Islands[,-1], row.names=CpG_Islands[,1])
write.table(CpG_Islands,"CpG_Islands_Raw_counts.txt",sep = "\t")
CpG_Islands_extention <- filter(df, Type != "inter")
CpG_Islands_extention <- CpG_Islands_extention[,-c(4:6)]
CpG_Islands_extention <-unite(CpG_Islands_extention, "position",c("Chr","Start"),sep = ":",remove = TRUE)
CpG_Islands_extention <-unite(CpG_Islands_extention, "position",c("position","End"),sep = "-",remove = TRUE)
CpG_Islands_extention <- data.frame(CpG_Islands_extention[,-1], row.names=CpG_Islands_extention[,1])
write.table(CpG_Islands_extention,"CpG_Islands_extention_Raw_counts.txt",sep = "\t")

colData <- read.table("colData.txt", header = T)
#Colorectal
colData_Colorectal <- subset(colData,Type == "Non_Cancer" | Type == "Cancer_Colorectal")
CpG_Islands_filter_Colorectal <- CpG_Islands_filter %>% select(colData_Colorectal[,1])
#write.table(CpG_Islands_filter_Colorectal,"CpG_Islands_filter_Colorectal_Raw_counts.txt", col.name = NA, sep = "\t", quote = FALSE)
CpG_Islands_extention_filter_Colorectal <- CpG_Islands_extention_filter %>% select(colData_Colorectal[,1])
colData_Colorectal$Type <- as.factor(colData_Colorectal$Type)
colData_Colorectal$Batch <- as.factor(colData_Colorectal$Batch)
dds_Colorectal <- DESeqDataSetFromMatrix(CpG_Islands_filter_Colorectal,colData_Colorectal, ~ Type)
dds_Colorectal <- DESeq(dds_Colorectal)
Colorectal = results(dds_Colorectal, contrast = c("Type","Cancer_Colorectal","Non_Cancer"))
summary(Colorectal)
Colorectal_islands_DMRs <-  subset(data.frame(Colorectal), padj < 0.1)
write.table(Colorectal_islands_DMRs, "CpGislands_DMRs_Colorectal.txt",sep = "\t", col.names = NA, quote = FALSE)

#Lung
colData_Lung <- subset(colData,Type == "Non_Cancer" | Type == "Cancer_Lung")
CpG_Islands_filter_Lung <- CpG_Islands_filter %>% select(colData_Lung[,1])
#write.table(CpG_Islands_filter_Lung,"CpG_Islands_filter_Lung_Raw_counts.txt", col.name = NA, sep = "\t", quote = FALSE)
CpG_Islands_extention_filter_Lung <- CpG_Islands_extention_filter %>% select(colData_Lung[,1]) 
colData_Lung$Type <- as.factor(colData_Lung$Type)
colData_Lung$Batch <- as.factor(colData_Lung$Batch)
dds_Lung <- DESeqDataSetFromMatrix(CpG_Islands_filter_Lung,colData_Lung, ~ Type)
dds_Lung <- DESeq(dds_Lung)
Lung = results(dds_Lung, contrast = c("Type","Cancer_Lung","Non_Cancer"))
summary(Lung)
Lung_islands_DMRs <-  subset(data.frame(Lung), padj < 0.1)
write.table(Lung_islands_DMRs, "CpGislands_DMRs_Lung.txt",sep = "\t", col.names = NA, quote = FALSE)

#Pancreas
colData_Pancreas <- subset(colData,Type == "Non_Cancer" | Type == "Cancer_Pancreas")
CpG_Islands_filter_Pancreas <- CpG_Islands_filter %>% select(colData_Pancreas[,1])
#write.table(CpG_Islands_filter_Pancreas,"CpG_Islands_filter_Pancreas_Raw_counts.txt", col.name = NA, sep = "\t", quote = FALSE)
CpG_Islands_extention_filter_Pancreas <- CpG_Islands_extention_filter %>% select(colData_Pancreas[,1]) 
colData_Pancreas$Type <- as.factor(colData_Pancreas$Type)
colData_Pancreas$Batch <- as.factor(colData_Pancreas$Batch)
dds_Pancreas <- DESeqDataSetFromMatrix(CpG_Islands_filter_Pancreas,colData_Pancreas, ~ Type)
dds_Pancreas <- DESeq(dds_Pancreas)
Pancreas = results(dds_Pancreas, contrast = c("Type","Cancer_Pancreas","Non_Cancer"))
summary(Pancreas)
Pancreas_islands_DMRs <-  subset(data.frame(Pancreas), padj < 0.1)
#nrow(subset(Pancreas_islands_DMRs, log2FoldChange>=1))
write.table(Pancreas_islands_DMRs, "CpGislands_DMRs_Pancreas.txt",sep = "\t", col.names = NA, quote = FALSE)

#MA-plot
png(file="MA-Colorectal vs Non-cancer (CGI).png",width=500, height=500)
plotMA(Colorectal, ylim=c(-3,3), main = "Colorectal vs Non-cancer (CGI)");dev.off()
png(file="MA-Lung vs Non-cancer (CGI).png",width=500, height=500)
plotMA(Lung, ylim=c(-3,3), main = "Lung vs Non-cancer (CGI)");dev.off()
png(file="MA-Pancreas vs Non-cancer (CGI).png",width=500, height=500)
plotMA(Pancreas, ylim=c(-3,3), main = "Pancreas vs Non-cancer (CGI)");dev.off()
png(file="MA-Colorectal vs Non-cancer (extended CGI).png",width=500, height=500)

#volcano plot
library(EnhancedVolcano)
png(file="Volcano-Lung vs Non-cancer (CGI).png",width=500, height=500)
EnhancedVolcano(Lung, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Lung vs Non-cancer (CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
png(file="Volcano-Colorectal vs Non-cancer (CGI).png",width=500, height=500)
EnhancedVolcano(Colorectal, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Colorectal vs Non-cancer (CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
png(file="Volcano-Pancreas vs Non-cancer (CGI).png",width=500, height=500)
EnhancedVolcano(Pancreas, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Pancreas vs Non-cancer (CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()

png(file="Volcano-Lung vs Non-cancer (extended CGI).png",width=500, height=500)
EnhancedVolcano(Lung_extended, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Lung vs Non-cancer (extended CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
png(file="Volcano-Colorectal vs Non-cancer (extended CGI).png",width=500, height=500)
EnhancedVolcano(Colorectal_extended, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Colorectal vs Non-cancer (extended CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()
png(file="Volcano-Pancreas vs Non-cancer (extended CGI).png",width=500, height=500)
EnhancedVolcano(Pancreas_extended, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' q value'), title = "Pancreas vs Non-cancer (extended CGI)", subtitle = "",
                pCutoff = 0.1, FCcutoff = 1, xlim = c(-3,3), ylim = c(0,5), pointSize = 2, col = c("black","black","blue","red3"), legendLabels = c("NS", "FC", "FDR", "FDR and FC"))
dev.off()

#vsd
vsd_Colorectal <- vst(dds_Colorectal, blind=FALSE)
vsd_Colorectal_df<- data.frame(assay(vsd_Colorectal),check.names = FALSE)
vsd_Lung <- vst(dds_Lung, blind=FALSE)
vsd_Lung_df<- data.frame(assay(vsd_Lung),check.names = FALSE)
vsd_Pancreas <- vst(dds_Pancreas, blind=FALSE)
vsd_Pancreas_df<- data.frame(assay(vsd_Pancreas),check.names = FALSE)

#PCA
library("devtools")
library("pcaExplorer")

png(file="PCA12-Colorectal vs Non-cancer (CGI).png",width=800, height=400)
PCA_Colorectal12 <- pcaplot(vsd_Colorectal,intgroup = "Type",ntop = 1000, ellipse = TRUE, returnData = FALSE,
                            pcX = 1, pcY = 2, text_labels = FALSE, title = "Colorectal DMCGIs PC1 vs PC2")
PCA_Colorectal12+ylim(-10,10)+xlim(-20,32) + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        legend.title = element_blank(),legend.text = element_text(colour = "black",size = 18),legend.position = c(.99, .01),legend.justification = c("right", "bottom"),legend.background=element_blank())
dev.off()
png(file="PCA12-Lung vs Non-cancer (CGI).png",width=800, height=400)
PCA_Lung12 <- pcaplot(vsd_Lung,intgroup = "Type",ntop = 1000, ellipse = TRUE, returnData = FALSE,
                      pcX = 1, pcY = 2, text_labels = FALSE, title = "Lung DMCGIs PC1 vs PC2")
PCA_Lung12+ylim(-10,10)+xlim(-20,32) + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        legend.title = element_blank(),legend.text = element_text(colour = "black",size = 18),legend.position = c(.99, .01),legend.justification = c("right", "bottom"),legend.background=element_blank())
dev.off()
png(file="PCA12-Pancreas vs Non-cancer (CGI).png",width=800, height=400)
PCA_Pancreas12 <- pcaplot(vsd_Pancreas,intgroup = "Type",ntop = 1000, ellipse = TRUE, returnData = FALSE,
                          pcX = 1, pcY = 2, text_labels = FALSE, title = "Pancreatic DMCGIs PC1 vs PC2")
PCA_Pancreas12+ylim(-10,10)+xlim(-30,22) + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        legend.title = element_blank(),legend.text = element_text(colour = "black",size = 18),legend.position = c(.99, .01),legend.justification = c("right", "bottom"),legend.background=element_blank())
dev.off()
pca_proportion_plot = function(dst, title = NULL){
  rv = rowVars(assay(dst))
  select = order(rv, decreasing = T)[seq_len(min(1000, length(rv)))]
  pca = prcomp(t(assay(dst)[select,]))
  res = pcascree(pca, type = "pev", pc_nr = 5, paste("Proportion of variance across principal components - ", title, sep="")) + theme_bw()
  print(res)
}
png(file="PCA Proportion-Colorectal vs Non-cancer.png",width=500, height=500)
pca_proportion_plot(vsd_Colorectal, "Colorectal")+ylim(0,0.6)+ ggtitle("Proportion of PCs - Colorectal")+
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
png(file="PCA Proportion-Lung vs Non-cancer.png",width=500, height=500)
pca_proportion_plot(vsd_Lung, "Lung")+ylim(0,0.6)+ ggtitle("Proportion of PCs - Lung")+
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
png(file="PCA Proportion-Pancreas vs Non-cancer.png",width=500, height=500)
pca_proportion_plot(vsd_Pancreas, "Pancreas")+ylim(0,0.6)+ ggtitle("Proportion of PCs - Pancreatic")+
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()

#Normalized counts
normalized_Colorectal <- counts(dds_Colorectal,normalized=TRUE)
normalized_Colorectal <- normalized_Colorectal[rownames(subset(Colorectal_islands_DMRs, log2FoldChange>=0)),]
normalized_Colorectal <- as.data.frame(t(normalized_Colorectal))
normalized_Colorectal <- setDT(normalized_Colorectal, keep.rownames = TRUE)[]
normalized_Colorectal$rn <- colData_Colorectal$Type
write.table(normalized_Colorectal, "Normalized_count_matrix_Colorectal_heatmap.txt",sep = "\t", row.names = FALSE, quote = FALSE)
normalized_Lung <- counts(dds_Lung,normalized=TRUE)
normalized_Lung <- normalized_Lung[rownames(subset(Lung_islands_DMRs, log2FoldChange>=0)),]
normalized_Lung <- as.data.frame(t(normalized_Lung))
normalized_Lung <- setDT(normalized_Lung, keep.rownames = TRUE)[]
normalized_Lung$rn <- colData_Lung$Type
write.table(normalized_Lung, "Normalized_count_matrix_Lung_heatmap.txt",sep = "\t", row.names = FALSE, quote = FALSE)
normalized_Pancreas <- counts(dds_Pancreas,normalized=TRUE)
normalized_Pancreas <- normalized_Pancreas[rownames(subset(Pancreas_islands_DMRs, log2FoldChange>=0)),]
normalized_Pancreas <- as.data.frame(t(normalized_Pancreas))
normalized_Pancreas <- setDT(normalized_Pancreas, keep.rownames = TRUE)[]
normalized_Pancreas$rn <- colData_Pancreas$Type
write.table(normalized_Pancreas, "Normalized_count_matrix_Pancreas_heatmap.txt",sep = "\t", row.names = FALSE, quote = FALSE)
