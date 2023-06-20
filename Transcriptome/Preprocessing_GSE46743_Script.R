# Preprocessing Script for the GSE46743
# Raw data files could be obtained from the corresponding repository at GSE46743 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46743
# Note 1: To execute the script, the user has to download files from repository
# Note 2: The equal sign = was used as an assignment operator for typing/productivity reasons

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

# A function to process a variable from GSE dataset (in the form of "Sex: Male" or "Age: 20", etc)
preprocess_GSE_var = function(x){
  x = sapply(x, function(z){
    z = stri_split_fixed(z, pattern = ":")
    z = unlist(z)
    z = z[2]
    z = str_trim(z)
    return(z)
  })
  return(x)
}

# A function to process GSE variable names
preprocess_GSE_var_names = function(x){
  x = sapply(x, function(z){
    z = stri_split_fixed(z, pattern = ":")
    z = unlist(z)
    z = z[1]
    z = str_trim(z)
    return(z)
  })
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

# A function to get phenotype table from GEO series file (sometimes it is different from the standard import)
geo_data_table_gsm_extract_phenotypes = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO, remove = TRUE, overwrite = TRUE)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "")
  Data_full = parseGEO(Data_GEO)
  GSMS = Data_full@gsms
  file.remove(Data_GEO)
  Full_GSMS = lapply(GSMS, function(x) x@header)
  Colnames_GSMS = names(Full_GSMS[[1]])
  Full_GSMS_DFs = lapply(Full_GSMS, function(x){
    x = lapply(x, function(j){
      if (length(j) > 1){
        j = paste0(j, collapse = "___")
      }
      return(j)
    })
    x = data.frame(x)
    char_vector = x$characteristics_ch1
    char_vector = unlist(stri_split_fixed(char_vector, pattern =  "___"))
    vars_GSM = preprocess_GSE_var(char_vector)
    names_GSM = preprocess_GSE_var_names(char_vector)
    vars_GSM = as.data.frame(vars_GSM)
    vars_GSM = t(vars_GSM)
    vars_GSM = as.data.frame(vars_GSM)
    colnames(vars_GSM) = names_GSM
    x$characteristics_ch1 = NULL
    x$data_row_count = NULL
    x = x[1,]
    x = cbind(x, vars_GSM)
    rownames(x) = NULL
    return(x)
  })
  Data_pheno_GSMS = as.data.frame(do.call(rbind, Full_GSMS_DFs))
  return(Data_pheno_GSMS)
}


################### Importing data ###################

# Getting data table and phenotype table
GSE46743_phenotypes_table = geo_data_table_gsm_extract_phenotypes(GEO_ID = "GSE46743")
GSE46743_data_table = geo_data_table_gsm_extract(GEO_ID = "GSE46743")
GSE46743_data_table = GSE46743_data_table[[2]]
GSE46743_data_table_IDs = GSE46743_data_table[,seq(from = 1, to = ncol(GSE46743_data_table), by = 2)]
all(apply(GSE46743_data_table_IDs, 1, function(x) length(unique(x)) == 1)) # All values are matching
GSE46743_data_table = GSE46743_data_table[,seq(from = 2, to = ncol(GSE46743_data_table), by = 2)]
rownames(GSE46743_data_table) = GSE46743_data_table_IDs$GSM1137699.ID_REF
colnames(GSE46743_data_table) = stri_replace_all_fixed(colnames(GSE46743_data_table), pattern = ".VALUE", replacement = "")

# Getting suppl. files for GSE46743
getGEOSuppFiles("GSE46743")
files = list.files("GSE46743")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE46743", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE46743_raw"))

files = list.files("GSE46743_raw")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE46743_raw", "/", files)
lapply(files, gunzip)
GSE46743_annotation = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE46743_MDD_control_dexamethasone_exposure/GSE46743_raw/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx",
                                  skip = "[Probes]")
all(rownames(GSE46743_data_table) %in% GSE46743_annotation$Probe_Id)
rownames(GSE46743_annotation) = GSE46743_annotation$Probe_Id
GSE46743_annotation = GSE46743_annotation[rownames(GSE46743_data_table),]

# Fixing phenotypes
GSE46743_phenotypes_table = GSE46743_phenotypes_table[match(colnames(GSE46743_data_table), GSE46743_phenotypes_table$geo_accession),]
colnames(GSE46743_phenotypes_table)[30] = "Depr_status"
all(GSE46743_phenotypes_table$geo_accession == colnames(GSE46743_data_table)) # Everything is matching

# Getting baseline data
GSE46743_phenotypes_table_baseline = GSE46743_phenotypes_table[GSE46743_phenotypes_table$stimulus == "baseline",]
GSE46743_data_processed_baseline = GSE46743_data_table[, colnames(GSE46743_data_table) %in% GSE46743_phenotypes_table_baseline$geo_accession]
all(GSE46743_phenotypes_table_baseline$geo_accession == colnames(GSE46743_data_processed_baseline)) # Everything is matching

# Fixing phenotype variables
GSE46743_phenotypes_table_baseline$Depr_status = factor(GSE46743_phenotypes_table_baseline$Depr_status, levels = c("case", "control"))
GSE46743_phenotypes_table_baseline$age = as.numeric(GSE46743_phenotypes_table_baseline$age)
GSE46743_phenotypes_table_baseline$bmi = as.numeric(GSE46743_phenotypes_table_baseline$bmi)

# Writing the output
write.csv(GSE46743_phenotypes_table_baseline, "GSE46743_phenotypes_table_baseline.csv")
write.csv(GSE46743_data_processed_baseline, "GSE46743_data_processed_baseline.csv")
write.csv(GSE46743_annotation, "GSE46743_annotation.csv")







