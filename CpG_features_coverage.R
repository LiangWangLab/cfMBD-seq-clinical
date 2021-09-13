library(tidyverse)
library(data.table)

#Listing files with .txt extension in the `data` folder
files = list.files(path = "data/.", pattern = ".txt") 
file_names = unlist(strsplit(files, ".txt")) #getting file basename
#Parallel read all files in data directory
data <- lapply(paste("data/", files, sep=""), data.table::fread)
df <- dplyr::bind_cols(data)
df <- data.frame(df)[,c(1:4,seq(5, ncol(df), by = 5))] #autoscale
colnames(df) = c("Chr", "Start", "End", "Annotation", file_names) #rename columns by file name
#Add length and type column
df <- add_column(df, Length = df$End - df$Start + 1, .after = "End")
df <- separate(df, "Annotation", "Type", sep = ":", remove = FALSE)
#df is raw sequence read for DESeq2

#TPM normalization
df_TPM <- df
df_TPM[-(1:6)] <- df_TPM[-(1:6)]/(df_TPM[,4]/1000)[row(df_TPM[-(1:6)])]
df_TPM[,-(1:6)] = apply(df_TPM[,-(1:6)],2,function(x){x/sum(x)*1e6})
head(df_TPM)
df_TPM <- df_TPM[,-c(4:5)]
df_TPM <-unite(df_TPM, "position",c("Chr","Start"),sep = ":",remove = TRUE)
df_TPM <-unite(df_TPM, "position",c("position","End"),sep = "-",remove = TRUE)
df_TPM <- data.frame(df_TPM[,-1], row.names=df_TPM[,1])
df_TPM <- df_TPM %>% select(c(colData_Colorectal[,1],subset(colData_Lung,Type=="Cancer_Lung")[,1],subset(colData_Pancreas,Type=="Cancer_Pancreas")[,1]))
df_TPM$Type <- df$Type
write.table(df_TPM[,-ncol(df_TPM)],"CpG_annotation_coverage_TPM.txt", sep = "\t", col.names = NA, quote = FALSE)

#Sum of each feature
Island_sum = apply(subset(df_TPM, Type == "island")[,-ncol(df_TPM)], 2, sum)
Shore_sum = apply(subset(df_TPM, Type == "shore")[,-ncol(df_TPM)], 2, sum)
Shelf_sum = apply(subset(df_TPM, Type == "shelf")[,-ncol(df_TPM)], 2, sum)
Inter_sum = apply(subset(df_TPM, Type == "inter")[,-ncol(df_TPM)], 2, sum)
Extended_sum = apply(subset(df_TPM, Type != "inter")[,-ncol(df_TPM)], 2, sum)
#To data frame
All <- data.frame(Island_sum,Shore_sum,Shelf_sum,Inter_sum,Extended_sum)
All <- round(All/10000,digits = 2)
colnames(All)[1:4] <- c("Island","Shore","Shelf","Inter")
All <- merge(All,colData,by.x=0,by.y="ExperimentID")
All$Type <- gsub("Non_Cancer", "Healthy", All$Type)
All$Type <- gsub("Cancer_Colorectal", "Colorectal", All$Type)
All$Type <- gsub("Cancer_Lung", "Lung", All$Type)
All$Type <- gsub("Cancer_Pancreas", "Pancreas", All$Type)
All$Type <- factor(All$Type, levels=c("Healthy", "Colorectal", "Lung", "Pancreas"))
summary(All$Island)
summary(All$Extended_sum)
write.table(All,"CpG_annotation_coverage_summary.txt", sep = "\t", row.names = F, quote = FALSE)

###ggplot2
library(ggplot2)
png(file="CpG annotation coverage summary(violin).png",width=500, height=500)
ggplot(reshape2::melt(All[,2:5]), aes(x=variable, y=value))+ ylim(0,60) + labs(title="CpG features coverage",x="", y = "Percentage %") + 
  geom_violin(trim=FALSE, fill="firebrick") + geom_boxplot(width=0.1, fill="steelblue",outlier.colour=NA) + theme_classic() + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
png(file="CpG islands coverage(violin).png",width=500, height=500)
ggplot(All, aes(x=Type, y=Island)) + ylim(20,70) + labs(title="CpG islands coverage",x="", y = "Percentage %") + 
  geom_violin(trim=FALSE, fill="firebrick") + geom_boxplot(width=0.1, fill="steelblue",outlier.colour=NA) + theme_classic() + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
png(file="Extended CpG islands coverage(violin).png",width=500, height=500)
ggplot(All, aes(x=Type, y=Extended_sum)) + ylim(87.5,97.5) + labs(title="Extended CpG islands coverage",x="", y = "Percentage %") + 
  geom_violin(trim=FALSE, fill="firebrick") + geom_boxplot(width=0.1, fill="steelblue",outlier.colour=NA) + theme_classic() + 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
