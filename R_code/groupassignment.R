# Question: What is the difference in survivability between old, mid, and young patients 
# with low vs. high gene expression of TP53, PIK3CA, MUC16?

# Installation
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
library(devtools)
library(robustbase)
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("arsenal"))install.packages(c("arsenal"))
if(!requireNamespace("survival"))install.packages(c("survival"))
if(!requireNamespace("survminer"))install.packages(c("survminer"))

# Loading packages
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(arsenal)
library(survival)
library(survminer)

# Barcodes
barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
                     "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  barcode = barcodes_rnaseq)
# GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", barcode=barcodes_clinic)
# GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
shorter_barcodes <- substr(colnames(htseq_counts), 1, 12)
colnames(htseq_counts) <- shorter_barcodes

# Accessing counts data for specific genes
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ TP53_mask ]
TP53_counts <- htseq_counts[TP53_ENSG_name, ]

# Reorder RNAseq barcodes to match clinical barcodes
TP53_transpose = t(TP53_counts)
htseq_counts = htseq_counts[match(colnames(htseq_counts),rownames(clinic))]

# Adding age to clinical data
age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

