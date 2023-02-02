# Preprocessing Script for the GSE64930
# Raw data files could be obtained from the corresponding repository at GSE64930 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64930
# Note 1: The authors used a last updated version of GSE64930 obtained directly from the MPIP (please contact MPIP to obtain raw files)
# Note 2: The equal sign = was used as an assignment operator as authors don't buy the idea of using <- for typing/productivity reasons

Working_directory = "..." # Replace with an appropriate directory
setwd(Working_directory)

################### Importing packages ###################
library(stringr)
library(stringi)
library(data.table)
library(GEOquery)
library(Biobase)
library(utils)


################### Defining functions ###################
# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 14, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}


################### Importing data ###################
# Getting suppl. files for GSE64930
getGEOSuppFiles("GSE64930")
files = list.files("GSE64930")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE64930", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE64930_raw"))

files = list.files("GSE64930")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE64930", "/", files)
lapply(files, gunzip)

files = list.files("GSE64930_raw")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE64930_raw", "/", files)
lapply(files, gunzip)

# Using the last updated version of the data (obtained from MPIP based on a request)
GSE64930_data_processed = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/Private_data/GSE64930_eData.txt")
GSE64930_annotation = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/GSE64930_raw/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt",
                                  skip = "[Probes]")
all(GSE64930_data_processed$ID %in% GSE64930_annotation$Probe_Id)
rownames(GSE64930_annotation) = GSE64930_annotation$Probe_Id
GSE64930_annotation = GSE64930_annotation[GSE64930_data_processed$ID,]
rownames(GSE64930_data_processed) = GSE64930_data_processed$ID
GSE64930_data_processed$ID = NULL
GSE64930_phenotypes =  smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/Private_data/GSE64930_pData.csv")
all(GSE64930_phenotypes$RNA_ID %in% colnames(GSE64930_data_processed)) # Everything is presented

# Curating phenotypes
GSE64930_phenotypes_baseline = GSE64930_phenotypes[GSE64930_phenotypes$Dex == 0,]
GSE64930_phenotypes_baseline$Sex = ifelse(GSE64930_phenotypes_baseline$Sex == 1, "M", "F")
GSE64930_phenotypes_baseline$Status = ifelse(GSE64930_phenotypes_baseline$Status == 1, "case", "control")
GSE64930_phenotypes_baseline$Status = factor(GSE64930_phenotypes_baseline$Status , levels = c("case", "control"))
GSE64930_phenotypes_baseline$Sex = factor(GSE64930_phenotypes_baseline$Sex)
GSE64930_data_processed_baseline = GSE64930_data_processed[,GSE64930_phenotypes_baseline$RNA_ID]
all(colnames(GSE64930_data_processed_baseline) == GSE64930_phenotypes_baseline$RNA_ID) # The order is matching

# Writing the output
write.csv(GSE64930_annotation, "GSE64930_annotation.csv")
write.csv(GSE64930_phenotypes_baseline, "GSE64930_phenotypes_baseline.csv")
write.csv(GSE64930_data_processed_baseline, "GSE64930_data_processed_baseline.csv")
