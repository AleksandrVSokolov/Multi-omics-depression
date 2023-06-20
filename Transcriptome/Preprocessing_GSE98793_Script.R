# Preprocessing Script for the GSE98793
# Raw data files could be obtained from corresponding repository at GSE98793 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98793
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
library(affy)
library(hgu133plus2.db)
library(sva)
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


################### Importing data ###################

# Getting supplementary files
getGEOSuppFiles("GSE98793")
files = list.files("GSE98793")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE98793", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE98793_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE98793_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE98793_CEL", "/", files)
lapply(files, gunzip)
# All files are fixed

# Getting GEO repository and inspecting data
GSE98793 = getGEO("GSE98793")
GSE98793 = GSE98793[[1]]
GSE98793_phenotypes = pData(GSE98793)
colnames(GSE98793_phenotypes) = stri_replace_all_fixed(colnames(GSE98793_phenotypes), pattern = ":ch1", "")
colnames(GSE98793_phenotypes) = stri_replace_all_fixed(colnames(GSE98793_phenotypes), pattern = " ", "_")
table(GSE98793_phenotypes$subject_group, GSE98793_phenotypes$anxiety)
GSE98793_probes = featureData(GSE98793)@data
GSE98793_expression = assayData(GSE98793)[["exprs"]]

# Getting tabular data with values
GSE98793_data_tables = geo_data_table_gsm_extract("GSE98793")
GSE98793_data_tables = GSE98793_data_tables[[2]]

# Getting correct CEL files
GSE98793_phenotypes$Corrected_CEL = sapply(GSE98793_phenotypes$supplementary_file, function(x){
  x = unlist(stri_split_fixed(x, pattern = "/"))
  x = x[length(x)]
  x = str_trim(x)
  x = stri_replace_all_fixed(x, pattern = ".gz", replacement = "")
  return(x)
})

# Check 
affy::list.celfiles("GSE98793_CEL") %in% GSE98793_phenotypes$Corrected_CEL

# Reading CEL files
Cel_files = affy::list.celfiles("GSE98793_CEL")
Cel_files = paste0("GSE98793_CEL", "/", Cel_files)
rawdata_GSE98793 = affy::ReadAffy(filenames = Cel_files)

# QC commands
affy::hist(rawdata_GSE98793) # Variance because of batching, need to use Combat from sva
affy::boxplot(rawdata_GSE98793) # To preview boxplots
head(intensity(rawdata_GSE53987)) # To preview intensities
image(rawdata_GSE53987, transfo=log) # To view images

# Correcting background, Normalizing, Summarizing
Eset_GSE98793  = affy::expresso(afbatch = rawdata_GSE98793, bgcorrect.method = "rma",
                                normalize.method = "quantiles",
                                pmcorrect.method = "pmonly",
                                summary.method = "medianpolish")
Eset_GSE98793_expression = affy::exprs(Eset_GSE98793)
affy::boxplot(Eset_GSE98793)
affy::hist(Eset_GSE98793)

# Writing normalized data
write.csv(Eset_GSE98793_expression, file = "Eset_GSE98793_expression.csv")


################### Adjusting for batch effects ###################
# Check and preparing phenotypes (we will use complex model adjustment with covariates)
all(colnames(Eset_GSE98793_expression) == GSE98793_phenotypes$Corrected_CEL) # The order is matching
GSE98793_phenotypes$batch = factor(GSE98793_phenotypes$batch)
GSE98793_phenotypes$age = as.numeric(GSE98793_phenotypes$age)
GSE98793_phenotypes$anxiety = factor(GSE98793_phenotypes$anxiety, levels = c("yes", "no"))
GSE98793_phenotypes$gender = factor(GSE98793_phenotypes$gender)
GSE98793_phenotypes$subject_group = ifelse(GSE98793_phenotypes$subject_group == "CASE; major depressive disorder (MDD) patient", "MDD", "Control")
GSE98793_phenotypes$subject_group =  factor(GSE98793_phenotypes$subject_group, levels = c("MDD", "Control"))

# Batch correction (no covariates)
Eset_GSE98793_expression_combat_simple = ComBat(Eset_GSE98793_expression, batch = GSE98793_phenotypes$batch)
write.csv(Eset_GSE98793_expression_combat_simple, file = "Eset_GSE98793_expression_combat_simple.csv")

# Batch correction (with covariates)
Batch.model.matrix = model.matrix(~ subject_group + anxiety + gender + age, data = GSE98793_phenotypes)
Eset_GSE98793_expression_combat_model = ComBat(Eset_GSE98793_expression, batch = GSE98793_phenotypes$batch, mod = Batch.model.matrix)
write.csv(Eset_GSE98793_expression_combat_model, file = "Eset_GSE98793_expression_combat_model.csv")

# All order is matching
all(colnames(Eset_GSE98793_expression_combat_model) == GSE98793_phenotypes$Corrected_CEL)

# Writing the files
write.csv(GSE98793_phenotypes, file = "GSE98793_phenotypes.csv")
write.csv(GSE98793_probes, file = "GSE98793_probes.csv")

