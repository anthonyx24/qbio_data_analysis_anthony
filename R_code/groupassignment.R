# Question: What is the difference in survivability between old, mid, and young patients 
# with low vs. high gene expression of TP53, PIK3CA, MUC16?

# Installation
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
# if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
# library(devtools)
# library(robustbase)
# if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
# if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
# install.packages(c("arsenal"), dependencies=TRUE)
# install.packages(c("survival"), dependencies=TRUE)
# install.packages(c("survminer"), dependencies=TRUE)

# Loading packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

# Barcodes
# barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
#                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
#                    "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
# barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
#                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
# barcode = barcodes_rnaseq
# GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type = "xml")
# barcode=barcodes_clinic
# GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient

# Matching clinical data sample order to RNAseq sample order
row_order <- match(colnames(htseq_counts), clinic$bcr_patient_barcode)
clinic_ordered  <- clinic[row_order, ]

# Get rid of nonmatching samples in clinical and htseq
matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))
clinic_matched <- clinic_ordered[matching,]

# Adding age to clinical data
age_clinical = clinic_matched$age_at_initial_pathologic_diagnosis
clinic_matched$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

# Accessing counts data for TP53, categorize expression, and add to clinical data
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ TP53_mask ]
TP53_counts <- htseq_counts[TP53_ENSG_name, clinic_matched$bcr_patient_barcode]
TP53_quartiles <- quantile(TP53_counts)
TP53_expression_level <- ifelse(TP53_counts > TP53_quartiles[4], "High", ifelse(TP53_counts < TP53_quartiles[2], "Low", "Mid"))
clinic_matched$TP53_expression = TP53_expression_level

# Accessing counts data for PIK3CA, categorize expression, and add to clinical data
PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ PIK3CA_mask ]
PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, clinic_matched$bcr_patient_barcode]
PIK3CA_quartiles <- quantile(PIK3CA_counts)
PIK3CA_expression_level <- ifelse(PIK3CA_counts > PIK3CA_quartiles[4], "High", ifelse(PIK3CA_counts < PIK3CA_quartiles[2], "Low", "Mid"))
clinic_matched$PIK3CA_expression = PIK3CA_expression_level

# Accessing counts data for MUC16, categorize expression, and add to clinical data
MUC16_mask <- rowData(sum_exp)$external_gene_name == "MUC16"
MUC16_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ MUC16_mask ]
MUC16_counts <- htseq_counts[MUC16_ENSG_name, clinic_matched$bcr_patient_barcode]
MUC16_quartiles <- quantile(MUC16_counts)
MUC16_expression_level <- ifelse(MUC16_counts > MUC16_quartiles[4], "High", ifelse(MUC16_counts < MUC16_quartiles[2], "Low", "Mid"))
clinic_matched$MUC16_expression = MUC16_expression_level

# Splitting clinic_matched based on age
clinic_old <- subset(clinic_matched, clinic_matched$age_category=="Old")
clinic_mid <- subset(clinic_matched, clinic_matched$age_category=="Mid")
clinic_young <- subset(clinic_matched, clinic_matched$age_category=="Young")

# Making KM plots for TP53
TCGAanalyze_survival(clinic_old, "TP53_expression", main="Kaplan-Meier Survival Curves for Old Patients with Varying TP53 Expression", filename = "TP53old.pdf")
TCGAanalyze_survival(clinic_mid, "TP53_expression", main="Kaplan-Meier Survival Curves for Middle-Aged Patients with Varying TP53 Expression", filename = "TP53mid.pdf")
TCGAanalyze_survival(clinic_young, "TP53_expression", main="Kaplan-Meier Survival Curves for Young Patients with Varying TP53 Expression", filename = "TP53young.pdf")

# Making KM plots for PIK3CA
TCGAanalyze_survival(clinic_old, "PIK3CA_expression", main="Kaplan-Meier Survival Curves for Old Patients with Varying PIK3CA Expression", filename = "PIK3CAold.pdf")
TCGAanalyze_survival(clinic_mid, "PIK3CA_expression", main="Kaplan-Meier Survival Curves for Middle-Aged Patients with Varying PIK3CA Expression", filename = "PIK3CAmid.pdf")
TCGAanalyze_survival(clinic_young, "PIK3CA_expression", main="Kaplan-Meier Survival Curves for Young Patients with Varying PIK3CA Expression", filename = "PIK3CAyoung.pdf")

# Making KM plots for MUC16
TCGAanalyze_survival(clinic_old, "MUC16_expression", main="Kaplan-Meier Survival Curves for Old Patients with Varying MUC16 Expression", filename = "MUC16old.pdf")
TCGAanalyze_survival(clinic_mid, "MUC16_expression", main="Kaplan-Meier Survival Curves for Middle-Aged Patients with Varying MUC16 Expression", filename = "MUC16mid.pdf")
TCGAanalyze_survival(clinic_young, "MUC16_expression", main="Kaplan-Meier Survival Curves for Young Patients with Varying MUC16 Expression", filename = "MUC16young.pdf")
