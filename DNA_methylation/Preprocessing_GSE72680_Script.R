# Preprocessing Script for the GSE72680
# Raw data files could be obtained from corresponding repository at GSE72680 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680
# Note 1: To execute the script several raw files used have to be downloaded from GEO
# Note 2: Computer-specific file paths are shown as "..."
# Note 3: The equal sign = was used as an assignment for typing/productivity reasons
# Note 4: In many cases loops were deliberately used instead of apply functions to enable better control of the variables

Working_directory = "..." # Replace with an appropriate directory
setwd(Working_directory)

# Importing packages (not all of them are required)
library(stringr)
library(dplyr)
library(openxlsx)
library(grid)
library(gdata)
library(magrittr)
library(data.table)
library(XML)
library(RCurl)
library(stringi)
library(httr)
library(FlowSorted.Blood.450k)
library(minfi)
library(wateRmelon)
library(lumi)
library(lmtest)
library(sva)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(GEOquery)
library(meffil)

# Setting options
options(stringsAsFactors = FALSE)


################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to quickly read data
smart_fread = function(x){
  x = as.data.frame(fread(x, nThread = 14))
  rownames(x) = x$V1
  x$V1 = NULL
  return(x)
}

################### Importing data ###################

# Obtaining data from GEO
gse = getGEO("GSE72680")
gse = gse[[1]]
experiment = gse@experimentData
GSE72680_phenotypes = gse@phenoData@data
# Phenotype data appeared to be incorrect since columns were skewed, continuing only with Beta values

# Working only with BETA values
Raw_betas = smart_fread(".../GSE72680_beta_values.txt") # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680
pvals = Raw_betas[,seq(from = 2, to = ncol(Raw_betas), by = 2)]
Raw_betas = Raw_betas[,seq(from = 1, to = ncol(Raw_betas)-1, by = 2)]

# Betas were not normalized but with corrected background
# Performing quantile normalization as in other cohorts
betQN = betaqn(Raw_betas)
betQN[betQN == 1] = 1 - min(betQN)

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # rows are not fully matching
betQN = betQN[rownames(betQN) %in% rownames(annot),]
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # now all rows are matching

# Beta-Mixture Quantile (BMIQ) Normalization to correct for probe type bias
betQN.BMIQ = BMIQ(betQN, annot$Type)
betas = as.data.frame(betQN.BMIQ$nbeta)

################### Filtering probes ###################

# First of all, filter the samples (we have to keep only the samples that are of good quality)
cols_to_select_samples = colSums(pvals<=0.00005)>=0.75*nrow(pvals) #Select samples with more than 75% of probes with detection P-value less than 0.00005

# Then, filter the probes (we have to keep only the probes that are of good quality)
lines_to_select_pvals = rownames(pvals)[rowSums(pvals<=0.01)>=0.75*ncol(pvals)] #Select probes with more than 75% of samples with detection P-value less than 0.01

# Remove sites on X and Y chromosomes
indexes_to_select_X = which(annot[,'chr'] != 'chrX')
indexes_to_select_Y = which(annot[,'chr'] != 'chrY')
indexes_to_select_XY = intersect(indexes_to_select_X, indexes_to_select_Y)
lines_to_select_XY = rownames(annot[indexes_to_select_XY, ])

# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch=rownames(annot)[which(substr(rownames(annot), 1, 2) == 'cg')]

# Remove sites with SNPs within probe with maf >5% 
indexes_to_select_SNP = which(is.na(annot[,'Probe_maf']) == T | annot[,'Probe_maf'] < 0.05)
lines_to_select_SNP = rownames(annot[indexes_to_select_SNP, ])

# Remove sites with SNPs within the SBE and CpG
indexes_to_select_SNP_SBE_CPG = which(is.na(annot[,'CpG_rs']) == T & is.na(annot[,'SBE_rs']) == T)
lines_to_select_SNP_SBE_CPG = rownames(annot[indexes_to_select_SNP_SBE_CPG, ])

# Inspecting removed probes
test = annot[indexes_to_select_SNP_SBE_CPG,]
table(test$CpG_rs) # questionable CpGs have been removed
table(test$SBE_rs) # questionable probes with SBE SNP have been removed
rm(test)

# Intersect all these probes that were selected, to only obtain the probes that pass all these filtering steps
lines_to_select = intersect(lines_to_select_XY,lines_to_select_pvals)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP)
lines_to_select = intersect(lines_to_select,lines_to_select_ch)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP_SBE_CPG)

# Get the beta-value matrix without probes that didn't pass the quality control and samples that didn't
betas2=betas[lines_to_select,cols_to_select_samples]

# Additional filtering for Illumina 450K
Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_MILES_BENTON$V1,]

# Remove sites with missing beta values
betas2=na.omit(betas2)

# Transform beta to M values
M.val=beta2m(betas2)

# Note 5: Batch effect for GSE72680 is corrected in the final model during analysis
# Note 6: Cell-type composition correction is performed in the final model during the analysis based on calculated proportions from GSE phenotpyes data. 
# -> The RGset object with all CpGs to correct for cell types is not available and minfi package cannot be used.

# Saving output
write.table(betas2, "beta_corrected_GSE72680.txt",sep="\t")
write.table(M.val, "Mval_corrected_GSE72680.txt",sep="\t")



