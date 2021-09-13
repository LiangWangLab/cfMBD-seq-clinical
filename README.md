# cfMBD-seq-clinical sample
Author: Jinyong Huang
# Introduction
The code within the LiangWangLab/cfMBD-seq-clinical repository is used for quick quality controls and differential methylation analysis of cfMBD-seq data.  
Bash programming and pre-installation of software are required. Access to a high-performance computing cluster is also recommended. 
# Requirements
Computer running a Linux system (≥ 8 GB RAM) Cluster computing is highly recommended when working with FASTQ/BAM files  
R (Version 4.1.1 or greater) 
# Linux Modules: 
fastp (Version 0.20.1);  
Bowtie2 (Version 2.4.2);   
Reference genome-human hg19;   
samtools (Version 1.11);   
bedtools (Version 2.28.0);   
Integrative Genomics Viewer.
# R Packages: 
BSgenome.hsapiens.UCSC.hg19;   
RaMWAS (Version 1.12.0);   
MEDIPS (Version 1.40.0);  
DESeq2 (Version 1.32.0);  
EnhancedVolcano (Version 1.10.0);  
pcaExplorer (Version 2.18.0);  
minfi (Version 1.36.0);  
annotatr (Version 1.16.0);  
Caret (Version 6.0-88);  
glmnet (Version 4.1-2);  
Rtsne (Version 0.15).
# Procedure
1. Bash to process raw data FASTQ files to BAM files.  
```./fastq_to_bam.sh```  
2. Quality control using RaMWAS to generate summary QC.  
```Rscript RaMWAS.R```    
3. Generate hg19 CpG features annotation reference.  
```Rscript CpG_annotations_reference.R```  
4. Call CpG annotations coverage using bedtools and normalize the reads.   
```./Call_CpG_coverage.sh```  
```Rscript CpG_features_coverage.R```  
5. Differential methylation analysis
```Rscript DESeq2.R```  
