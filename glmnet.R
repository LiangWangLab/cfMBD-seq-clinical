library(data.table)
library(tidyverse)
library(Rtsne)
library(caret)
library(glmnet)
library(ROCR)
library(ggplot2)

#Tissue_specific DMR from cfDNA
Tissue_specific_colorectal <- read.table("Tissue_Colorectal_specific.txt",header = T)
Tissue_specific_lung <- read.table("Tissue_Lung_specific.txt",header = T)
Tissue_specific_pancreas <- read.table("Tissue_Pancreas_specific.txt",header = T)
#Top cfDNA DMRs
Tissue_specific_colorectal_top<-top_n(Tissue_specific_colorectal, 100, log2FoldChange)[,1:7]
Tissue_specific_lung_top<-top_n(Tissue_specific_lung, 100, log2FoldChange)[,1:7]
Tissue_specific_pancreas_top<-top_n(Tissue_specific_pancreas, 100, log2FoldChange)[,1:7]
#450k locus annotation
hg19_450k_annotaion <- read.csv("hg19_HumanMethylation450_annotation.csv", header = TRUE)
hg19_450k_CGI <- hg19_450k_annotaion[,c(1,25)]
hg19_450k_CGI <- subset(hg19_450k_CGI, UCSC_CpG_Islands_Name!= "")#CpG sites on CpG islands
hg19_450k_CGI <- hg19_450k_CGI %>% separate(UCSC_CpG_Islands_Name, c("chr", "start", "end"))
hg19_450k_CGI <- transform(hg19_450k_CGI, start = as.numeric(start))
hg19_450k_CGI <- transform(hg19_450k_CGI, end = as.numeric(end))
hg19_450k_CGI$start <- hg19_450k_CGI$start+1
hg19_450k_CGI <-unite(hg19_450k_CGI, "position", c("chr","start"), sep = ":", remove = TRUE)
hg19_450k_CGI <-unite(hg19_450k_CGI, "position", c("position","end"), sep = "-", remove = TRUE)
#tansform 450k data to CGI data
All_450k_data <- read.table("TCGA_450k_757_beta_value.txt", sep = "\t", header = TRUE)
rownames(All_450k_data)<-All_450k_data[,1];All_450k_data<-All_450k_data[,-1]
All_450k_data_CGI <- All_450k_data[(row.names(All_450k_data) %in% hg19_450k_CGI$IlmnID),]
All_450k_data_CGI <- merge(All_450k_data_CGI, hg19_450k_CGI, by.x = 0, by.y = "IlmnID")
All_CGI_data <- aggregate(All_450k_data_CGI[,2:(ncol(All_450k_data_CGI)-1)], by=list(All_450k_data_CGI$position), FUN=mean, na.rm=TRUE, na.action=NULL)

####Predictive modeling
##Colorectal
Top_colorectal <- merge(Tissue_specific_colorectal_top, All_CGI_data,  by.x = "X", by.y = "Group.1")
Top_colorectal <- Top_colorectal[,-(2:7)]
rownames(Top_colorectal)<-Top_colorectal[,1];Top_colorectal<-Top_colorectal[,-1]
x_colorectal<-t(Top_colorectal)
y_colorectal = factor(ifelse(grepl("COAD", row.names(x_colorectal)), "Colorectal", "NONcolorectal"))
rm(coef_colorectal_df)
coef_colorectal_list <- list()
AUC_colorectal_list <- double()
for(i in 1:200){
set.seed(i)
splitSample_colorectal <- createDataPartition(y_colorectal, p = 0.8, list = FALSE)
training_x_colorectal <- x_colorectal[splitSample_colorectal,]
training_y_colorectal <- y_colorectal[splitSample_colorectal]
validation_x_colorectal <- x_colorectal[-splitSample_colorectal,]
validation_y_colorectal <- y_colorectal[-splitSample_colorectal]
#glmnet
colorectal_cv.lasso <- cv.glmnet(training_x_colorectal, training_y_colorectal, alpha = 1, family = "binomial")
best_colorectal_lambda <- colorectal_cv.lasso$lambda.min
colorectal_model <- glmnet(training_x_colorectal, training_y_colorectal, alpha = 1, family = "binomial", lambda = best_colorectal_lambda)
coef_colorectal <- coef(colorectal_model)
prob_colorectal <- predict(colorectal_model,newx = validation_x_colorectal)
pred_colorectal <- prediction(prob_colorectal, validation_y_colorectal)
AUC_colorectal<- attr(performance(pred_colorectal, "auc"), "y.values")[[1]]
coef_colorectal_list[[i]] = coef_colorectal
AUC_colorectal_list[[i]] = AUC_colorectal
if(!exists("coef_colorectal_df")){
  coef_colorectal_df = as.data.frame(as.matrix(coef_colorectal_list[[i]]))
} else {
  coef_colorectal_df = cbind(coef_colorectal_df, as.data.frame(as.matrix(coef_colorectal_list[[i]])))
}
}
summary(AUC_colorectal_list)
coef_colorectal_df_filter <- coef_colorectal_df[apply(coef_colorectal_df, 1, function(x) !any(x==0)),]
coef_colorectal_df_filter <- rownames(coef_colorectal_df_filter[-1,])
##Lung
Top_lung <- merge(Tissue_specific_lung_top, All_CGI_data,  by.x = "X", by.y = "Group.1")
Top_lung <- Top_lung[,-(2:7)]
rownames(Top_lung)<-Top_lung[,1];Top_lung<-Top_lung[,-1]
x_lung<-t(Top_lung)
y_lung = factor(ifelse(grepl("LUAD", row.names(x_lung)), "Lung", "NONlung"))

rm(coef_lung_df)
coef_lung_list <- list()
AUC_lung_list <- double()
for(i in 1:200){
  set.seed(i)
  splitSample_lung <- createDataPartition(y_lung, p = 0.8, list = FALSE)
  training_x_lung <- x_lung[splitSample_lung,]
  training_y_lung <- y_lung[splitSample_lung]
  validation_x_lung <- x_lung[-splitSample_lung,]
  validation_y_lung <- y_lung[-splitSample_lung]
#glmnet
  lung_cv.lasso <- cv.glmnet(training_x_lung, training_y_lung, alpha = 1, family = "binomial")
  best_lung_lambda <- lung_cv.lasso$lambda.min
  lung_model <- glmnet(training_x_lung, training_y_lung, alpha = 1, family = "binomial", lambda = best_lung_lambda)
  coef_lung <- coef(lung_model)
  prob_lung <- predict(lung_model,newx = validation_x_lung)
  pred_lung <- prediction(prob_lung, validation_y_lung)
  AUC_lung<- attr(performance(pred_lung, "auc"), "y.values")[[1]]
  coef_lung_list[[i]] = coef_lung
  AUC_lung_list[[i]] = AUC_lung
  if(!exists("coef_lung_df")){
    coef_lung_df = as.data.frame(as.matrix(coef_lung_list[[i]]))
  } else {
    coef_lung_df = cbind(coef_lung_df, as.data.frame(as.matrix(coef_lung_list[[i]])))
  }
}
summary(AUC_lung_list)
coef_lung_df_filter <- coef_lung_df[apply(coef_lung_df, 1, function(x) !any(x==0)),]
coef_lung_df_filter <- rownames(coef_lung_df_filter[-1,])      
##Pancreas
Top_pancreas <- merge(Tissue_specific_pancreas_top, All_CGI_data,  by.x = "X", by.y = "Group.1")
Top_pancreas <- Top_pancreas[,-(2:7)]
rownames(Top_pancreas)<-Top_pancreas[,1];Top_pancreas<-Top_pancreas[,-1]
x_pancreas<-t(Top_pancreas)
y_pancreas = factor(ifelse(grepl("PAAD", row.names(x_pancreas)), "Pancreas", "NONpancreas"))
rm(coef_pancreas_df)
coef_pancreas_list <- list()
AUC_pancreas_list <- double()
for(i in 1:200){
  set.seed(i)
  splitSample_pancreas <- createDataPartition(y_pancreas, p = 0.8, list = FALSE)
  training_x_pancreas <- x_pancreas[splitSample_pancreas,]
  training_y_pancreas <- y_pancreas[splitSample_pancreas]
  validation_x_pancreas <- x_pancreas[-splitSample_pancreas,]
  validation_y_pancreas <- y_pancreas[-splitSample_pancreas]
#glmnet
  pancreas_cv.lasso <- cv.glmnet(training_x_pancreas, training_y_pancreas, alpha = 1, family = "binomial")
  best_pancreas_lambda <- pancreas_cv.lasso$lambda.min
  pancreas_model <- glmnet(training_x_pancreas, training_y_pancreas, alpha = 1, family = "binomial", lambda = best_pancreas_lambda)
  coef_pancreas <- coef(pancreas_model)
  prob_pancreas <- predict(pancreas_model,newx = validation_x_pancreas)
  pred_pancreas <- prediction(prob_pancreas, validation_y_pancreas)
  AUC_pancreas<- attr(performance(pred_pancreas, "auc"), "y.values")[[1]]
  coef_pancreas_list[[i]] = coef_pancreas
  AUC_pancreas_list[[i]] = AUC_pancreas
  if(!exists("coef_pancreas_df")){
    coef_pancreas_df = as.data.frame(as.matrix(coef_pancreas_list[[i]]))
  } else {
    coef_pancreas_df = cbind(coef_pancreas_df, as.data.frame(as.matrix(coef_pancreas_list[[i]])))
  }
}
summary(AUC_pancreas_list)
coef_pancreas_df_filter <- coef_pancreas_df[apply(coef_pancreas_df, 1, function(x) !any(x==0)),]
coef_pancreas_df_filter <- rownames(coef_pancreas_df_filter[-1,])
##AUC
AUC <- data.frame(AUC=c(AUC_colorectal_list,AUC_lung_list,AUC_pancreas_list),
                  Type=c(rep("COAD",200),rep("LUAD",200),rep("PAAD",200)))
png(file="AUC.png",width=500, height=500)
ggplot(AUC, aes(x=Type, y= AUC)) + labs(title="HM450K testing set",x="", y = "AUC")+
  geom_boxplot(width=0.1, fill="steelblue",outlier.colour=NA) + geom_jitter(shape=16, position=position_jitter(0.1),colour = "red")+ ylim(0.95,1) + theme_classic()+ 
  theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2),legend.position = "none")
dev.off()
                                                  
####t-sne plot
coef<-c(coef_colorectal_df_filter,coef_lung_df_filter,coef_pancreas_df_filter)                                                  
##for HM450k data    
PBMC_data <- read.table("PBMCs61_df.txt", sep = "\t", header = TRUE)
rownames(PBMC_data)<-PBMC_data[,1];PBMC_data<-PBMC_data[,-1]
PBMC_data_CGI <- PBMC_data[(row.names(PBMC_data) %in% hg19_450k_CGI$IlmnID),]
PBMC_data_CGI <- merge(PBMC_data_CGI, hg19_450k_CGI, by.x = 0, by.y = "IlmnID")      
PBMC_data_CGI <- aggregate(PBMC_data_CGI[,2:(ncol(PBMC_data_CGI)-1)], by=list(PBMC_data_CGI$position), FUN=mean, na.rm=TRUE, na.action=NULL)
rownames(PBMC_data_CGI)<-PBMC_data_CGI[,1];PBMC_data_CGI<-PBMC_data_CGI[,-1] 
All_CGI_df<-All_CGI_data[,-1];rownames(All_CGI_df)<-All_CGI_data[,1]
All_CGI_df<-cbind(All_CGI_df,PBMC_data_CGI)
All_CGI_df<-All_CGI_df[coef,]
All_CGI_df<-as.data.frame(t(All_CGI_df))
condition_450k <- c(substr(rownames(All_CGI_df[1:757,]),13,16), substr(rownames(All_CGI_df[758:818,]),1,3))
condition_450k <-gsub("GSM","PBMC",condition_450k)
All_CGI_df$Type <-condition_450k
colors_450k <- c("green","red","yellow","blue")
names(colors_450k)<- sort(unique(All_CGI_df$Type))
tsne_450k <- Rtsne(All_CGI_df, perplexity=10, check_duplicates = FALSE)
png(file="tsne-450k.png",width=500, height=500)
plot(tsne_450k$Y, col=colors_450k[All_CGI_df$Type], xlab="tSNE 1", ylab="tSNE 2", pch =20, main = "HM450K data set", cex=2, cex.lab = 1.5, cex.axis = 1.7, cex.main = 2)
legend("topright", legend = unique(All_CGI_df$Type), fill = colors_450k[unique(All_CGI_df$Type)],cex=1.5)
dev.off()                                                 
##for cfMBD data   
cfMBD_coldata <- read.table("cfMBD_coldata.txt",header = T)
cfMBD_coldata$Type <- replace(cfMBD_coldata$Type, cfMBD_coldata$Type=="Colorectal", "Cancer-Colorectal")
cfMBD_coldata$Type <- replace(cfMBD_coldata$Type, cfMBD_coldata$Type=="Lung", "Cancer-Lung")
cfMBD_coldata$Type <- replace(cfMBD_coldata$Type, cfMBD_coldata$Type=="Pancreas", "Cancer-Pancreas")
cfMBD_coldata$Type <- replace(cfMBD_coldata$Type, cfMBD_coldata$Type=="Healthy", "Non-Cancer")
cfMBD_data <- read.table("cfMBD_vsd_df.txt",header = T)
cfMBD_data_classifier <- cfMBD_data[coef,]
cfMBD_data_classifier <- as.data.frame(t(cfMBD_data_classifier))
cfMBD_data_classifier$Type <-cfMBD_coldata$Type
colors_cfMBD <- c("green","red","yellow","blue")
names(colors_cfMBD)<- sort(unique(cfMBD_data_classifier$Type))
tsne_cfMBD <- Rtsne(cfMBD_data_classifier, perplexity=16, check_duplicates = FALSE)
png(file="tsne-cfMBD.png",width=500, height=500)
plot(tsne_cfMBD$Y, col=colors_cfMBD[cfMBD_data_classifier$Type], xlab="tSNE 1", ylab="tSNE 2", pch =20, main = "cfMBD-seq data set", cex=2.5, cex.lab = 1.5, cex.axis = 1.7, cex.main = 2)
legend("topright", legend = sort(unique(cfMBD_data_classifier$Type)), fill = colors_cfMBD[sort(unique(cfMBD_data_classifier$Type))],cex=1.5)
dev.off()                                                  
