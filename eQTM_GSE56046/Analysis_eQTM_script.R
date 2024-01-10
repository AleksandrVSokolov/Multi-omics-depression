Working_directory = "/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/eQTM_GSE56046" # Replace with an appropriate path
setwd(Working_directory)

# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)


################### Package import ###################
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
  x = as.data.frame(fread(x, nThread = 10, header = TRUE, ...))
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


################### Importing data ###################
GSE56046_pheno = smart_fread("GSE56046_pheno.csv")
normalized_methylation_mval_filtered = smart_fread("normalized_methylation_mval_filtered.csv")
normalized_eset_filtered = smart_fread("normalized_eset_filtered.csv")

CpG_anno = minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
CpG_anno = as.data.frame(CpG_anno)
CpG_anno = CpG_anno[rownames(normalized_methylation_mval_filtered),]
Illumina_HT_12_V3_V4_probes_Alignments_PASSED = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/Transcriptome/Illumina_HT_12_V3_V4_probes_Alignments_PASSED.csv")
Illumina_HT_12_V3_V4_probes_Alignments_PASSED = Illumina_HT_12_V3_V4_probes_Alignments_PASSED[Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID %in% rownames(normalized_eset_filtered),]

normalized_eset_filtered = normalized_eset_filtered[rownames(normalized_eset_filtered) %in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
################### Processing pheno data ###################
GSE56046_pheno$age # numeric
GSE56046_pheno$bcell # numeric
GSE56046_pheno$ChIP = factor(GSE56046_pheno$ChIP) # factor
GSE56046_pheno$neutro # numeric
GSE56046_pheno$nkcell # numeric
GSE56046_pheno$tcell # numeric
GSE56046_pheno$racegendersite # factor
GSE56046_pheno$racegendersite = factor(GSE56046_pheno$racegendersite)

################### Analysis ###################
short_CpG_anno = CpG_anno[,c( "Name", "chr", "pos")]

CpGs = CpG_anno$Name

dir.create("results")
CpG_chunks = split(CpGs, ceiling(seq_along(CpGs) / 10000))

for (i in 1:length(CpG_chunks)){
  
  writeLines(paste0("Chunk: ", i))
  print(Sys.time())
  
  selected_CpGs = CpG_chunks[[i]]
  
  cis_eQTM_analysis = mclapply(selected_CpGs, function(x){
    
    select_anno_CpG = short_CpG_anno[x,]
    select_anno_expr = Illumina_HT_12_V3_V4_probes_Alignments_PASSED[Illumina_HT_12_V3_V4_probes_Alignments_PASSED$chrom == select_anno_CpG$chr,]
    select_anno_expr = select_anno_expr[,1:6]
    
    # Selecting CpG sector
    position_CpG = select_anno_CpG$pos
    
    if (position_CpG > 1000000){
      min_area = position_CpG-1000000
      max_area = position_CpG+1000000
    } else {
      min_area = 0
      max_area = position_CpG+1000000
    }
    
    # Selecting related genes
    selection_idx = mapply(function(x,y){
      
      if (x >= min_area & x<= max_area){
        return(TRUE)
      }
      
      if (y >= min_area & y<= max_area){
        return(TRUE)
      }
      
      return(FALSE)
    }, select_anno_expr$txStart, select_anno_expr$txEnd)
    candidate_expression = select_anno_expr[selection_idx, ]
    
    if (nrow(candidate_expression) < 1){
      
      return(NA)
      
    }
    
    candidate_expression_values = normalized_eset_filtered[candidate_expression$ID,]
    
    df = cbind(t(normalized_methylation_mval_filtered[x,]), t(candidate_expression_values), GSE56046_pheno)
    
    model_formulas = paste0(" ~ ", x , " + age + bcell + ChIP + neutro + nkcell + tcell + racegendersite")
    model_formulas = paste0(rownames(candidate_expression_values), model_formulas)
    
    # running models
    models = lapply(model_formulas, function(z){
      
      lin_model = lm(formula = as.formula(z), data = df)
      
      coefs = summary(lin_model)$coefficients
      coefs = coefs[2,]
      coefs = as.data.frame(t(coefs))
      
      return(coefs)
    })
    models = list_to_df(models)
    
    models = cbind(models, candidate_expression)
    models$CpG = select_anno_CpG$Name
    models$CpG_pos = select_anno_CpG$pos
    
    return(models)
    
  }, mc.cores = 10)
  cis_eQTM_analysis = cis_eQTM_analysis[sapply(cis_eQTM_analysis, is.data.frame)]
  cis_eQTM_analysis = list_to_df(cis_eQTM_analysis)
  
  file_path = paste0("results/", i, "_chunk.csv")
  fwrite(cis_eQTM_analysis, file_path, sep = ",", row.names = TRUE)
  gc()
  
}

files = list.files("results")
files = paste0("results/", files)
all_eQTMs = lapply(files, smart_fread)
all_eQTMs = list_to_df(all_eQTMs)
rownames(all_eQTMs) = NULL
fwrite(all_eQTMs, "all_eQTMs.csv", sep = ",", row.names = FALSE)
