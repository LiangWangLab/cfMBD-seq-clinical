library(ramwas)

#run once
library(BSgenome.Hsapiens.UCSC.hg19)
cpgset = getCpGsetCG(BSgenome.Hsapiens.UCSC.hg19)
saveRDS(file = "UCSC.hg19_cpgset.rds", object = cpgset)
#Put all bam files in a 'bams' folder
files = list.files(path = "bams/", pattern = ".bam$") 
files = gsub(".bam","",files)
writeLines(files,"bam_list.txt")
#Parameters
param = ramwasParameters(
  filebamlist = "bam_list.txt",
  filecpgset = "UCSC.hg19_cpgset.rds",
  cputhreads = 10,
  scoretag = "MAPQ",
  minscore = 4,
  minfragmentsize = 100,
  maxfragmentsize = 600,
  minavgcpgcoverage = 0.3,
  minnonzerosamples = 0.3,
  filecovariates = "covariates.txt",
  modelcovariates = "Batch",
  modeloutcome = "Group",
  modelPCs = 0
)
##Quality control
ramwas1scanBams(param)
ramwas2collectqc(param)
#Open /qc/summary_bams/Summary_QC.txt to see summary QC.

#figure
summary <- read.csv("qc/summary_bams/Summary_QC.txt",sep = "\t")
coldata <- read.csv("coldata.txt",sep = "\t")#metadata of your samples
summary$Type <- coldata$GroupName
png(file="total reads.png",width=500, height=500)
boxplot(as.numeric(gsub("," ,"", summary$Total.reads))/1000000~summary$Type, main = "Total sequenced reads", xlab = "",ylab = "Millions",col = "steelblue",border = "black", cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
dev.off()

png(file="High quality reads.png",width=500, height=500)
boxplot(as.numeric(gsub("," ,"", summary$Reads.used.for.coverage))/1000000~summary$Type, main = "High quality reads", xlab = "",ylab = "Millions",col = "steelblue",border = "black", cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
dev.off()

png(file="noise.png",width=500, height=500)
boxplot(summary$Non.Cpg.CpG.coverage.ratio~summary$Type, main = "Noise", xlab = "",ylab = "non-CpG/CpG coverage ratio",col = "steelblue",border = "black", cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
dev.off()

png(file="noCpG.png",width=500, height=500)
boxplot(as.numeric(substr(summary$Non.CpG.reads....,1,4))~summary$Type, main = "Non-CpG reads percentage", xlab = "",ylab = "Percentage %",col = "steelblue",border = "black", cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
dev.off()

png(file="peak.png",width=500, height=500)
boxplot(summary$Peak.SQRT^2~summary$Type, main = "CpG density at peak", xlab = "",ylab = "CpG density",col = "steelblue",border = "black", cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
dev.off()
