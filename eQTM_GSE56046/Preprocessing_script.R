Working_directory = "/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/eQTM_GSE56046" # Replace with an appropriate path
setwd(Working_directory)

# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)


################### Package import ###################
# Importing that might be used (Note: Not all of them may be required)
library(fun)
library(stringr)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(grid)
library(gdata)
library(RColorBrewer)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(magrittr)
library(igraph)
library(visNetwork)
library(data.table)
library(XML)
library(rvest)
library(RCurl)
library(HGNChelper)
library(stringi)
library(httr)
library(lubridate)
library(rjson)
library(rtracklayer)
library(rstudioapi)
library(tidyr)
library(Gviz)
library(limma)
library(FactoMineR)
library(ggthemes)
library(igraph)
library(RSelenium)
library(lumi)
library(outliers)
library(svglite)
library(scatterplot3d)
library(sva)
library(jsonlite)
library(ggrepel)
library(parallel)
library(bacon)
library(gridExtra)
library(ggplotify)
library(HGNChelper)
library(jetset)
library(GEOquery)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(chromoMap)
library(RIdeogram)
library(ggVennDiagram)
library(seqinr)
library(Biostrings)
library(dbparser)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to convert list to data frame
list_to_df = function(data_list){
  if (length(data_list) > 1){
    data_list = do.call(rbind, data_list)
  } else {
    data_list = data_list[[1]]
  }
  return(data_list)
}

# A function to replace multiple patterns by multiple replacements in a string
multiple_stri_replacer = function(string, pattern_vector, replacement_vector){
  
  # Pattern_vector and replacement_vector should have the same length
  for (i in 1:length(pattern_vector)){
    string = stri_replace_all_fixed(str = string, pattern = pattern_vector[i], replacement = replacement_vector[i])
  }
  return(string)
}

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 10, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

# A function to get value table from GEO series file (sometimes it is different from the standard import)
geo_data_table_gsm_extract = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "") #useful!
  Data_full = parseGEO(Data_GEO)
  file.remove(Data_GEO)
  GSMS = Data_full@gsms
  Full_GSMS = lapply(GSMS, function(x) x@dataTable@table)
  Full_Columns = lapply(GSMS, function(x) x@dataTable@columns)
  Full_Columns = Full_Columns[[1]]
  Full_GSMS = do.call(cbind, Full_GSMS)
  Output = list()
  Output[[1]] = Full_Columns
  Output[[2]] = Full_GSMS
  return(Output)
}

################### Reading raw file ###################
raw_expression = smart_fread("expression_GSE56045_non_normalized.txt")
raw_expression = raw_expression[stri_detect_fixed(raw_expression$ID_REF, "ILMN"),]
rownames(raw_expression) = raw_expression$ID_REF
raw_expression$ID_REF = NULL

raw_expression_vals = raw_expression[,stri_detect_fixed(colnames(raw_expression), pattern = "intensity")]
raw_expression_P = raw_expression[,stri_detect_fixed(colnames(raw_expression), pattern = "detectionPval")]

# -> Do we have a reason not to use already normalized data?
# -> Using normalized values as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280798/

################### Getting phenotypes ###################
GSE56045 = getGEO("GSE56045")
GSE56045 = GSE56045[[1]]
GSE56045_pheno = GSE56045@phenoData@data
colnames(GSE56045_pheno) = stri_replace_all_fixed(colnames(GSE56045_pheno), pattern = ":ch1", replacement = "")
table(GSE56045_pheno$source_name_ch1) # all monocytes

################### Getting normalized expression (as in the initial publication) ###################
normalized_expression = geo_data_table_gsm_extract("GSE56045")
normalized_expression_vals = normalized_expression[[2]]

IDs = normalized_expression_vals[, stri_detect_fixed(colnames(normalized_expression_vals), "ID_REF")]

# Row-wise check is all IDs match
check = apply(IDs, 1, function(x) length(unique(x)) == 1)
all(check) # TRUE

# Select eset
normalized_eset = normalized_expression_vals[, stri_detect_fixed(colnames(normalized_expression_vals), "VALUE")]
rownames(normalized_eset) = IDs$GSM1352002.ID_REF
colnames(normalized_eset) = stri_replace_all_fixed(colnames(normalized_eset), pattern = ".VALUE", replacement = "")
all(colnames(normalized_eset)  %in% GSE56045_pheno$geo_accession) # TRUE

# Getting Detection P
normalized_expression_P =  normalized_expression_vals[, stri_detect_fixed(colnames(normalized_expression_vals), "detectionPval")]
rownames(normalized_expression_P) = IDs$GSM1352002.ID_REF

# Filtering expression based on detection (similar to GSE64930) to obtain roughly similar number of probes
filter_rows = rownames(normalized_expression_P)[rowSums(normalized_expression_P <= 0.05) >= 0.5 * ncol(normalized_expression_P)]

normalized_eset_filtered = normalized_eset[filter_rows,]

################### Methylation data ###################
# Methylation data will be used as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280798/

# GSE56046 (methylation)
GSE56046 = getGEO("GSE56046")
GSE56046 = GSE56046[[1]]
GSE56046_pheno = GSE56046@phenoData@data
colnames(GSE56046_pheno) = stri_replace_all_fixed(colnames(GSE56046_pheno), pattern = ":ch1", replacement = "")
table(GSE56046_pheno$source_name_ch1) # all monocytes

GSE56046_pheno$geo_accession == GSE56045_pheno$geo_accession # not equal order
rownames(GSE56046_pheno) == GSE56046_pheno$geo_accession # TRUE
rownames(GSE56045_pheno) == GSE56045_pheno$geo_accession # TRUE

# We will use the order as in GSE56046_pheno
# IDs do not match
GSE56045_pheno$geo_accession %in% GSE56046_pheno$geo_accession

################### Obtain participant IDs ###################
GSE56045_pheno$participant_ID = stri_replace_all_fixed(GSE56045_pheno$title, pattern = "_peripheral_CD14", replacement = "")
GSE56046_pheno$participant_ID = stri_replace_all_fixed(GSE56046_pheno$title, pattern = "_peripheral_CD14 [methylation]", replacement = "")

all(GSE56045_pheno$participant_ID == GSE56046_pheno$participant_ID) # All IDs match!
all(GSE56045_pheno$age == GSE56046_pheno$age) # Ages match (as check)

################### Inspecting methylation ###################
normalized_methylation = smart_fread("GSE56046_methylome_normalized.txt")
rownames(normalized_methylation) = normalized_methylation$ID_REF
normalized_methylation$ID_REF = NULL

normalized_methylation_P =  normalized_methylation[, stri_detect_fixed(colnames(normalized_methylation), "detectionPval")]
normalized_methylation_mval =  normalized_methylation[, stri_detect_fixed(colnames(normalized_methylation), "Mvalue")]

# Filtering methylation similar to main methyl cohorts
annot = minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot = annot[rownames(normalized_methylation_mval),]

# First of all, filter the samples (we have to keep only the samples that are of good quality)
cols_to_select_samples = colSums(normalized_methylation_P<=0.00005)>=0.75*nrow(normalized_methylation_P) #Select samples with more than 75% of probes with detection P-value less than 0.00005

# Then, filter the probes (we have to keep only the probes that are of good quality)
lines_to_select_pvals = rownames(normalized_methylation_P)[rowSums(normalized_methylation_P<=0.01)>=0.75*ncol(normalized_methylation_P)] #Select probes with more than 75% of samples with detection P-value less than 0.01

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

# Intersect all these probes that were selected, to only obtain the probes that pass all these filtering steps
lines_to_select = intersect(lines_to_select_XY,lines_to_select_pvals)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP)
lines_to_select = intersect(lines_to_select,lines_to_select_ch)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP_SBE_CPG)

# Get values that passed filtering
normalized_methylation_mval_filtered = normalized_methylation_mval[lines_to_select,cols_to_select_samples]

# Additional filtering for Illumina 450K
Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
normalized_methylation_mval_filtered = normalized_methylation_mval_filtered[rownames(normalized_methylation_mval_filtered) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
normalized_methylation_mval_filtered = normalized_methylation_mval_filtered[rownames(normalized_methylation_mval_filtered) %!in% Cross_reactive_probes_MILES_BENTON$V1,]

# Rename columns
colnames(normalized_methylation_mval_filtered) = stri_replace_all_fixed(colnames(normalized_methylation_mval_filtered), 
                                                                        pattern = ".Mvalue", 
                                                                        replacement = "")

# check
all(colnames(normalized_methylation_mval_filtered) == GSE56046_pheno$participant_ID) # All IDs match!
all(GSE56045_pheno$participant_ID == GSE56046_pheno$participant_ID) # All IDs match!

all(colnames(normalized_eset_filtered) == GSE56045_pheno$geo_accession) # All IDs match!

# Replace column names in expression set
colnames(normalized_eset_filtered) = GSE56045_pheno$participant_ID

all(colnames(normalized_methylation_mval_filtered) == GSE56046_pheno$participant_ID) # All IDs match!
all(colnames(normalized_eset_filtered) == GSE56046_pheno$participant_ID) # All IDs match!
# Everything is matching

################### Saving data  ###################
fwrite(normalized_eset_filtered, "normalized_eset_filtered.csv", sep = ",", row.names = TRUE)
fwrite(normalized_methylation_mval_filtered, "normalized_methylation_mval_filtered.csv", sep = ",", row.names = TRUE)
fwrite(GSE56046_pheno, "GSE56046_pheno.csv", sep = ",", row.names = TRUE)


