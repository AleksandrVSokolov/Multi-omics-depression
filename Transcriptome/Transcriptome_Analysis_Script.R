################# This script shows analysis for Transcriptome data ################# 
# Note 1: The equal sign = was used as an assignment operator as authors don't buy the idea of using <- for typing/productivity reasons
# Note 2: In many cases loops were deliberately used instead of apply functions to enable better control of the variables (even though loops in R are slow and computationally inefficient)
# Note 3: Some variables in the loops contain prefixes to enable easy cleanup of the environment once the loop is executed
# Note 4: The whole running environment gets heavy once the data is imported (It is not advised to run a script on a machine with less than 32 GB of RAM)
# Note 5: Some of the raw files (primarily an updated data for GSE64930) could not be deposited publicly since was obtained with a request
# -> To get an updated data for GSE64930, please contact MPIP or a corresponding author for further instructions

Working_directory = "~/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/Transcriptome" # Replace with an appropriate path
setwd(Working_directory)

# Setting options
getOption("scipen") # default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)

################### Package import ###################
# Importing thet might be used (Note: Not all of them may be required)
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
library(UniprotR)
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
library(UniprotR)
library(RIdeogram)
library(ggVennDiagram)
library(seqinr)
library(Biostrings)


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
  
  # pattern_vector and replacement_vector should have the same length
  for (i in 1:length(pattern_vector)){
    string = stri_replace_all_fixed(str = string, pattern = pattern_vector[i], replacement = replacement_vector[i])
  }
  return(string)
}

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 14, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

# A function to detect at least one pattern in a string
multiple_stri_detector = function(string, pattern_vector){
  output_list = list()
  for (i in 1:length(pattern_vector)){
    output_list[[i]] = stri_detect_fixed(str = string, pattern = pattern_vector[i])
  }
  output_list = do.call(rbind, output_list)
  apply(output_list, 2, any)
}

# A function to expand a data frame where several columns contain condensed cells with a specified separator
multiple_expander = function(df, cols_to_expand, pattern){
  #
  orig_colnames = colnames(df)
  df_modif = df[, cols_to_expand, drop = FALSE]
  df_const = df[,-cols_to_expand, drop = FALSE]
  orig_colnames_modif = colnames(df_modif)
  
  # running expansion
  df_list = list()
  for (i in 1:nrow(df_const)){
    print(i)
    curr_df_const = df_const[i,, drop = FALSE]
    curr_df_modif = df_modif[i,, drop = FALSE]
    curr_df_modif = apply(curr_df_modif, 2, function(x) unlist(stri_split_fixed(x, pattern = pattern)))
    
    if (length(cols_to_expand) > 1){
      curr_df_modif = do.call(cbind, curr_df_modif)
    } else {
      curr_df_modif = as.character(curr_df_modif)
    }
    
    if (is.matrix(curr_df_modif)){
      for (b in 1:nrow(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    } else {
      for (b in 1:length(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    }
    
    curr_df_combined = cbind(curr_df_const, curr_df_modif)
    
    if (length(cols_to_expand) == 1){
      colnames(curr_df_combined)[ncol(curr_df_combined)] = orig_colnames[cols_to_expand]
    }
    df_list[[i]] = curr_df_combined
  }
  
  if (length(df_list) > 1){
    df_list = do.call(rbind, df_list)
  } else {
    df_list = df_list[[1]]
  }
  
  df_list = df_list[, orig_colnames]
  return(df_list)
}

# A function to check gene symbols
# Importing NIH dataset
Homo_Sapiens_Gene_info_NIH = smart_fread("/home/aleksandr/Desktop/WORK/UCSC_ID_MAP/Homo_sapiens.gene_info") # https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/ (replace with an appropriate path)
Homo_Sapiens_Gene_info_NIH_expanded = multiple_expander(df = Homo_Sapiens_Gene_info_NIH, cols_to_expand = 5, pattern = "|")
check_gene_symbol_NIH = function(PRF_gene_symbols, PRF_ref_NIH_expanded, PRF_replace_NA_with_old = FALSE){
  PRF_gene_symbols_check = lapply(PRF_gene_symbols, function(x){
    if (x %in% PRF_ref_NIH_expanded$Symbol_from_nomenclature_authority){
      Curr_gene = x
      Approved = TRUE
      Suggested.Symbol = x
    } else if (x %in% PRF_ref_NIH_expanded$Symbol){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Symbol == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
        if (Suggested.Symbol == x){
          Approved = TRUE
        }
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else if (x %in% PRF_ref_NIH_expanded$Synonyms){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Synonyms == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else {
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = NA
      if (PRF_replace_NA_with_old){
        Suggested.Symbol = x
      }
    }
    Dataset = data.frame(x = Curr_gene, Approved = Approved, Suggested.Symbol = Suggested.Symbol)
  })
  PRF_gene_symbols_check = list_to_df(PRF_gene_symbols_check)
  return(PRF_gene_symbols_check)
}

# This function performs differential expression analysis for Affymetrix-based arrays
# It automatically removes all probes that contain more that a single gene in the annotation
# The function variables are labelled with the prefix PREFIX_ (to easily identify all variables if debugging/cleaning the environment)
# The LogFC threshold is used for Volcano plot
test_genes_affymetrix_limma_generelized = function(PREFIX_gene_names,
                                                   PREFIX_remove_non_specific = TRUE,
                                                   PREFIX_Pheno_df,
                                                   PREFIX_Expression_matrix,
                                                   PREFIX_Probes_df,
                                                   PREFIX_contrast_col_number, 
                                                   PREFIX_participants_col_number, 
                                                   PREFIX_contrast_vector, 
                                                   PREFIX_model.covariates,
                                                   PREFIX_plots_folder = NULL,
                                                   PREFIX_plots_str = NULL,
                                                   PREFIX_logFC_threshold = 0.1){
  # All parameters have prefix so it is easier to remove after making a function
  # gene names: a list of gene names to test or NULL
  # Pheno_df: a data frame with phenotypes (must be ordered as Expression_matrix)
  # Expression_matrix: a matrix with log2 expression values
  # Probes_df: a data frame with probes
  # contrast_col_number: a column, specifying the main factor (contrast for analysis)
  # participants_col_number: a column number that matches IDs for participants in Pheno_df and Expression_matrix
  # contrast_vector: a vector of contrasts to compare
  # model.covariates: covariates of the model (besides the main contrast)
  # plots_folder folder for volcano plots
  
  if (is.null(PREFIX_gene_names)){
    writeLines("Testing all available probes")
    writeLines("Filtering probes")
    PREFIX_Probes_df = PREFIX_Probes_df[!is.na(PREFIX_Probes_df$ID),]
    if (PREFIX_remove_non_specific){
      writeLines("Removing non-specific probes")
      PREFIX_Probes_df = PREFIX_Probes_df[!stri_detect_fixed(PREFIX_Probes_df$`Gene Symbol`, pattern = " /// "),]
    }
    writeLines("Removing non-mapped probes")
    PREFIX_Probes_df = PREFIX_Probes_df[PREFIX_Probes_df$`Gene Symbol` != "",]
    PREFIX_Probes_df = PREFIX_Probes_df[!is.na(PREFIX_Probes_df$`Gene Symbol`),]
    PREFIX_Probes_test = PREFIX_Probes_df
    PREFIX_Probes_test$Fixed_ID = paste0("ID_",PREFIX_Probes_test$ID)
  } else {
    
    # Checking gene names in probes
    PREFIX_probes_genes = lapply(PREFIX_Probes_df$`Gene Symbol`, function(x){
      PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = " /// "))
      PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
      return(PREFIX_Current_genes)
    })
    PREFIX_probes_genes = unlist(PREFIX_probes_genes)
    PREFIX_probes_genes = unique(PREFIX_probes_genes)
    if (all(PREFIX_gene_names %in% PREFIX_probes_genes)){
      writeLines("All genes are presented")
    } else {
      PREFIX_Not_present_genes = PREFIX_gene_names[PREFIX_gene_names %!in% PREFIX_probes_genes]
      PREFIX_Not_present_genes = paste0(PREFIX_Not_present_genes, collapse = "; ")
      PREFIX_Not_present_genes = paste0("Several genes are not presented: ", PREFIX_Not_present_genes)
      writeLines(PREFIX_Not_present_genes)
      PREFIX_gene_names = PREFIX_gene_names[PREFIX_gene_names %in% PREFIX_probes_genes]
      if (length(PREFIX_gene_names) < 1){
        writeLines("No genes to test, returning NA")
        return(NA)
      }
    }
    
    # Extracting and ordering probes
    PREFIX_gene_names = sort(PREFIX_gene_names)
    PREFIX_probes_genes = lapply(PREFIX_Probes_df$`Gene Symbol`, function(x){
      PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = " /// "))
      PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
      return(PREFIX_Current_genes)
    })
    PREFIX_Probes_test = lapply(PREFIX_gene_names, function(x){
      PREFIX_Indeces = sapply(PREFIX_probes_genes, function(j){
        if (x %in% j){
          return(TRUE)
        } else {
          return(FALSE)
        }
      })
      PREFIX_gene_probes = PREFIX_Probes_df[PREFIX_Indeces,]
      return(PREFIX_gene_probes)
    })
    
    # Probes df small
    PREFIX_Probes_test = do.call(rbind, PREFIX_Probes_test)
    PREFIX_Probes_test = dplyr::distinct(PREFIX_Probes_test, ID, .keep_all = TRUE)
    PREFIX_Probes_test$Fixed_ID = paste0("ID_",PREFIX_Probes_test$ID)
  }
  
  # Preparing the data
  colnames(PREFIX_Pheno_df)[PREFIX_contrast_col_number] = "Case_Control"
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Case_Control %in% PREFIX_contrast_vector, ]
  
  # Limma reorders the factors
  PREFIX_Pheno_df$Case_Control = factor(PREFIX_Pheno_df$Case_Control, levels = c(PREFIX_contrast_vector[2], PREFIX_contrast_vector[1]))
  colnames(PREFIX_Pheno_df)[PREFIX_participants_col_number] = "Participant_ID_FUN"
  
  # Getting expression values
  PREFIX_Tested_expression = PREFIX_Expression_matrix[PREFIX_Probes_test$ID,]
  PREFIX_Tested_expression = PREFIX_Tested_expression[,colnames(PREFIX_Tested_expression) %in% PREFIX_Pheno_df$Participant_ID_FUN]
  
  # Running the analysis
  PREFIX_Model.formula.string = paste0(PREFIX_model.covariates, collapse = " + ")
  PREFIX_Model.formula.string = paste0("~ ", "Case_Control + ", PREFIX_Model.formula.string)
  PREFIX_Design.matrix = model.matrix(as.formula(PREFIX_Model.formula.string), data = PREFIX_Pheno_df)
  PREFIX_fit = lmFit(PREFIX_Tested_expression, PREFIX_Design.matrix)
  PREFIX_fitE = eBayes(PREFIX_fit)
  PREFIX_Top_table = limma::topTable(fit = PREFIX_fitE, coef = 2, adjust.method = "fdr", number = Inf)
  
  # Adding Variables
  PREFIX_Probes_test = PREFIX_Probes_test[rownames(PREFIX_Top_table),]
  PREFIX_Top_table$ID = rownames(PREFIX_Top_table)
  PREFIX_Top_table$Model.formula.string = PREFIX_Model.formula.string
  PREFIX_Top_table$Contrast = paste0(PREFIX_contrast_vector, collapse = "/")
  PREFIX_Top_table$Sequence_Type = PREFIX_Probes_test$`Sequence Type`
  PREFIX_Top_table$Gene = PREFIX_Probes_test$`Gene Symbol`
  
  # Checking gene symbols
  PREFIX_tt_genes = PREFIX_Top_table$Gene
  PREFIX_tt_genes = unlist(stri_split_fixed(PREFIX_tt_genes, pattern = " /// "))
  PREFIX_tt_genes = unique(PREFIX_tt_genes)
  PREFIX_tt_genes = PREFIX_tt_genes[PREFIX_tt_genes != ""]
  PREFIX_tt_genes = PREFIX_tt_genes[!is.na(PREFIX_tt_genes)]
  PREFIX_probes_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = PREFIX_tt_genes, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                    PRF_replace_NA_with_old = TRUE)
  PREFIX_probes_genes_check_dict = PREFIX_probes_genes_check$Suggested.Symbol
  names(PREFIX_probes_genes_check_dict) = PREFIX_probes_genes_check$x
  PREFIX_Top_table$Updated_gene_names = sapply(PREFIX_Top_table$Gene, function(x){
    if (is.na(x)){return(NA)}
    if (x == ""){return(NA)}
    PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = " /// "))
    PREFIX_Current_genes = unique(PREFIX_Current_genes)
    PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
    PREFIX_New_names = PREFIX_probes_genes_check_dict[PREFIX_Current_genes]
    PREFIX_New_names = paste0(PREFIX_New_names, collapse = ";")
    return(PREFIX_New_names)
  })
  PREFIX_Top_table = PREFIX_Top_table[,c("Contrast", "Model.formula.string", "Gene", "Updated_gene_names", "ID", "Sequence_Type", "logFC", "AveExpr", "t" ,"P.Value" , "adj.P.Val","B")]
  
  # Adjusting p-values with bacon
  PREFIX_Degr.freedom = PREFIX_fitE$df.total[1]
  PREFIX_BACON_stats = bacon(PREFIX_Top_table$t)
  PREFIX_Inflation.Bacon = inflation(PREFIX_BACON_stats)
  PREFIX_Bias.Bacon  = bias(PREFIX_BACON_stats)
  PREFIX_T.Bacon = bacon::tstat(PREFIX_BACON_stats)
  
  # Getting p-values for T.bacon
  PREFIX_P.val.Bacon = sapply(PREFIX_T.Bacon, function(x){
    if (x < 0){
      stats = pt(q = x, df = PREFIX_Degr.freedom, lower.tail = TRUE)*2
    } 
    if (x > 0){
      stats = pt(q = x, df = PREFIX_Degr.freedom, lower.tail = FALSE)*2
    }
    return(stats)
  })
  PREFIX_P.val.Bacon = unlist(PREFIX_P.val.Bacon)
  PREFIX_Adj.P.val.Bacon = p.adjust(p = PREFIX_P.val.Bacon, method = "fdr")
  PREFIX_Top_table = cbind(PREFIX_Top_table,
                           PREFIX_Inflation.Bacon,
                           PREFIX_Bias.Bacon, 
                           PREFIX_T.Bacon,
                           PREFIX_P.val.Bacon,
                           PREFIX_Adj.P.val.Bacon)
  colnames(PREFIX_Top_table) = stri_replace_all_fixed(colnames(PREFIX_Top_table), pattern = "PREFIX_", replacement = "")
  PREFIX_Top_table_output = PREFIX_Top_table

  # Making the volcano plot
  # Adding highlight for transcript based on the fold change
  PREFIX_Top_table$is.highlight = sapply(PREFIX_Top_table$logFC, function(x){
    if (x >  PREFIX_logFC_threshold){
      x = "Up"
    } else if (x < - PREFIX_logFC_threshold){
      x = "Down"
    } else {
      x = "None"
    }
  })
  
  # Modification of Calculated_stats
  Pval_treshold = PREFIX_Top_table[PREFIX_Top_table$adj.P.Val < 0.05,]
  if (nrow(Pval_treshold) < 1){
    writeLines("No Significant P-values")
    PREFIX_Top_table$is_annotate = "no"
    Pval_treshold = NULL
  } else {
    Pval_treshold = arrange(Pval_treshold, -P.Value)
    Pval_treshold = Pval_treshold$P.Value[1]
    Pval_treshold = -log10(Pval_treshold)
    PREFIX_Top_table = mutate(PREFIX_Top_table, is_annotate = ifelse(-log10(P.Value)>Pval_treshold & is.highlight != "None", "yes", "no"))
  }
  
  # Make the plot
  plot = ggplot(PREFIX_Top_table, aes(x=logFC, y=-log10(P.Value))) +
    
    # Show all points
    geom_point(aes(color= factor(is.highlight, levels = c("None", "Down", "Up"))), alpha=0.6, size=4) +
    scale_color_manual(values = c("grey", "skyblue", "red")) 
  
  # Add pval line FDR
  if (!is.null(Pval_treshold)){
    plot = plot + geom_hline(yintercept=Pval_treshold, linetype="dashed", 
                             color = "red", size=0.5)
  }
  
  # Add pval line 0.05
  plot = plot + geom_hline(yintercept= -log10(0.05), linetype="dashed", 
                           color = "blue", size=0.5)
  # Add logFC lines
  if (min(PREFIX_Top_table$logFC) < -PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= -PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  if (max(PREFIX_Top_table$logFC) > PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  # Add label using ggrepel to avoid overlapping
  if (any(PREFIX_Top_table$is_annotate == "yes")){
    plot = plot + 
      geom_label_repel(data = subset(PREFIX_Top_table, is_annotate == "yes"), aes(label=ID), size=4, force = 10, 
                       max.overlaps = 50)
  } 
  
  plot = plot +
    labs(x = "log2 fold change", y = "-log10 p-value") +
    # Customizing the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
      panel.background = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  dir.create(PREFIX_plots_folder)
  
  file_name = paste0(PREFIX_plots_folder, "/", PREFIX_plots_str, "Volcano", "_",
                     paste0(PREFIX_contrast_vector, collapse = "_"), ".png")
  ggsave(file = file_name, plot = plot, width=2560, height=1440, units = "px", scale = 2)
  
  return(PREFIX_Top_table_output)
}

# This function performs differential expression analysis for Illumina-based arrays
# It automatically removes all probes that contain more that a single gene in the annotation
# The function variables are labelled with the prefix PREFIX_ (to easily identify all variables if debugging/cleaning the environment)
# The LogFC threshold is used for Volcano plot
test_genes_Illumina_limma_generelized = function(PREFIX_gene_names,
                                                 PREFIX_remove_non_specific = TRUE,
                                                 PREFIX_Pheno_df,
                                                 PREFIX_Expression_matrix,
                                                 PREFIX_Probes_df,
                                                 PREFIX_contrast_col_number, 
                                                 PREFIX_participants_col_number, 
                                                 PREFIX_contrast_vector, 
                                                 PREFIX_model.covariates,
                                                 PREFIX_plots_folder = NULL,
                                                 PREFIX_plots_str = NULL,
                                                 PREFIX_Remove_NA_predictors = TRUE,
                                                 PREFIX_logFC_threshold = 0.1){
  # All parameters have prefix so it is easier to remove after making a function
  # gene names: a list of gene names to test or NULL
  # Pheno_df: a data frame with phenotypes (must be ordered as Expression_matrix)
  # Expression_matrix: a matrix with log2 expression values
  # Probes_df: a data frame with probes
  # contrast_col_number: a column, specifying the main factor (contrast for analysis)
  # participants_col_number: a column number that matches IDs for participants in Pheno_df and Expression_matrix
  # contrast_vector: a vector of contrasts to compare
  # model.covariates: covariates of the model (besides the main contrast)
  # plots_folder folder for volcano plots
  
  if (is.null(PREFIX_gene_names)){
    writeLines("Testing all available probes")
    writeLines("Filtering probes")
    PREFIX_Probes_df = PREFIX_Probes_df[!is.na(PREFIX_Probes_df$Probe_Id),]
    if (PREFIX_remove_non_specific){
      writeLines("Removing non-specific probes")
      PREFIX_Probes_df = PREFIX_Probes_df[!stri_detect_fixed(PREFIX_Probes_df$Symbol, pattern = ", "),]
      PREFIX_Probes_df = PREFIX_Probes_df[!stri_detect_fixed(PREFIX_Probes_df$Symbol, pattern = "; "),]
      PREFIX_Probes_df = PREFIX_Probes_df[!stri_detect_fixed(PREFIX_Probes_df$Symbol, pattern = ";"),]
    }
    writeLines("Removing non-mapped probes")
    PREFIX_Probes_df = PREFIX_Probes_df[PREFIX_Probes_df$Symbol != "",]
    PREFIX_Probes_df = PREFIX_Probes_df[!is.na(PREFIX_Probes_df$Symbol),]
    PREFIX_Probes_df_symbol_check = check_gene_symbol_NIH(PRF_gene_symbols = PREFIX_Probes_df$Symbol, 
                                                          PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded, 
                                                          PRF_replace_NA_with_old = TRUE)
    PREFIX_Probes_df$Symbol_upd = PREFIX_Probes_df_symbol_check$Suggested.Symbol
    PREFIX_Probes_df$All_symb = paste0(PREFIX_Probes_df$Symbol,"; ", PREFIX_Probes_df$Synonyms,"; ", PREFIX_Probes_df$Obsolete_Probe_Id)
    PREFIX_Probes_test = PREFIX_Probes_df
  } else {
    
    #Checking gene names in probes
    PREFIX_probes_genes = lapply(PREFIX_Probes_df$All_symb, function(x){
      PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = "; "))
      PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
      return(PREFIX_Current_genes)
    })
    PREFIX_probes_genes = unlist(PREFIX_probes_genes)
    PREFIX_probes_genes = unique(PREFIX_probes_genes)
    if (all(PREFIX_gene_names %in% PREFIX_probes_genes)){
      writeLines("All genes are presented")
    } else {
      PREFIX_Not_present_genes = PREFIX_gene_names[PREFIX_gene_names %!in% PREFIX_probes_genes]
      PREFIX_Not_present_genes = paste0(PREFIX_Not_present_genes, collapse = "; ")
      PREFIX_Not_present_genes = paste0("Several genes are not presented: ", PREFIX_Not_present_genes)
      writeLines(PREFIX_Not_present_genes)
      PREFIX_gene_names = PREFIX_gene_names[PREFIX_gene_names %in% PREFIX_probes_genes]
      if (length(PREFIX_gene_names) < 1){
        writeLines("No genes to test, returning NA")
        return(NA)
      }
    }
    
    # Extracting and ordering probes
    PREFIX_gene_names = sort(PREFIX_gene_names)
    PREFIX_probes_genes = lapply(PREFIX_Probes_df$All_symb, function(x){
      PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = "; "))
      PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
      return(PREFIX_Current_genes)
    })
    PREFIX_Probes_test = lapply(PREFIX_gene_names, function(x){
      PREFIX_Indeces = sapply(PREFIX_probes_genes, function(j){
        if (x %in% j){
          return(TRUE)
        } else {
          return(FALSE)
        }
      })
      PREFIX_gene_probes = PREFIX_Probes_df[PREFIX_Indeces,]
      return(PREFIX_gene_probes)
    })
    PREFIX_Probes_test = do.call(rbind, PREFIX_Probes_test)
    PREFIX_Probes_test = dplyr::distinct(PREFIX_Probes_test, Probe_Id, .keep_all = TRUE)
  }
  
  # Preparing the data
  colnames(PREFIX_Pheno_df)[PREFIX_contrast_col_number] = "Case_Control"
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Case_Control %in% PREFIX_contrast_vector, ]
  
  # Limma reorders the factor levels
  PREFIX_Pheno_df$Case_Control = factor(PREFIX_Pheno_df$Case_Control, levels = c(PREFIX_contrast_vector[2], PREFIX_contrast_vector[1]))
  colnames(PREFIX_Pheno_df)[PREFIX_participants_col_number] = "Participant_ID_FUN"
  
  # Getting expression values
  PREFIX_Tested_expression = PREFIX_Expression_matrix[PREFIX_Probes_test$Probe_Id,]
  PREFIX_Tested_expression = PREFIX_Tested_expression[,colnames(PREFIX_Tested_expression) %in% PREFIX_Pheno_df$Participant_ID_FUN]
  
  # Preparing pheno df before model (removing samples with missing values if needed)
  PREFIX_Pheno_df = PREFIX_Pheno_df[,c("Participant_ID_FUN", "Case_Control", PREFIX_model.covariates)]
  if (PREFIX_Remove_NA_predictors){
    writeLines("Removing NA predictors")
    PREFIX_Participants_missing_data = apply(PREFIX_Pheno_df, 1, function(x){
      x = as.character(x)
      if (any(is.na(x))){
        return(TRUE)
      } else{
        return(FALSE)
      }
    })
    PREFIX_Participants_missing_data = PREFIX_Pheno_df$Participant_ID_FUN[PREFIX_Participants_missing_data]
  }
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Participant_ID_FUN %!in% PREFIX_Participants_missing_data,]
  PREFIX_Tested_expression = PREFIX_Tested_expression[,colnames(PREFIX_Tested_expression) %in% PREFIX_Pheno_df$Participant_ID_FUN]
  PREFIX_Tested_expression = PREFIX_Tested_expression[,PREFIX_Pheno_df$Participant_ID_FUN]
  print(table(PREFIX_Pheno_df$Case_Control))
  
  # Running the analysis
  PREFIX_Model.formula.string = paste0(PREFIX_model.covariates, collapse = " + ")
  PREFIX_Model.formula.string = paste0("~ ", "Case_Control + ", PREFIX_Model.formula.string)
  PREFIX_Design.matrix = model.matrix(as.formula(PREFIX_Model.formula.string), data = PREFIX_Pheno_df)
  PREFIX_fit = lmFit(PREFIX_Tested_expression, PREFIX_Design.matrix)
  PREFIX_fitE = eBayes(PREFIX_fit)
  PREFIX_Top_table = limma::topTable(fit = PREFIX_fitE, coef = 2, adjust.method = "fdr", number = Inf)
  
  # Adding Variables
  PREFIX_Probes_test = PREFIX_Probes_test[rownames(PREFIX_Top_table),]
  PREFIX_Top_table$ID = rownames(PREFIX_Top_table)
  PREFIX_Top_table$Model.formula.string = PREFIX_Model.formula.string
  PREFIX_Top_table$Contrast = paste0(PREFIX_contrast_vector, collapse = "/")
  PREFIX_Top_table$Probe_Type = PREFIX_Probes_test$Probe_Type
  PREFIX_Top_table$Gene = PREFIX_Probes_test$Symbol_upd
  PREFIX_Top_table = PREFIX_Top_table[,c("Contrast", "Model.formula.string", "Gene", "ID", "Probe_Type", "logFC", "AveExpr", "t" ,"P.Value" , "adj.P.Val","B")]
  
  # Adjusting p-values with bacon
  PREFIX_Degr.freedom = PREFIX_fitE$df.total[1]
  PREFIX_BACON_stats = bacon(PREFIX_Top_table$t)
  PREFIX_Inflation.Bacon = inflation(PREFIX_BACON_stats)
  PREFIX_Bias.Bacon  = bias(PREFIX_BACON_stats)
  PREFIX_T.Bacon = bacon::tstat(PREFIX_BACON_stats)
  
  # Getting p-values for T.bacon
  PREFIX_P.val.Bacon = sapply(PREFIX_T.Bacon, function(x){
    if (x < 0){
      stats = pt(q = x, df = PREFIX_Degr.freedom, lower.tail = TRUE)*2
    } 
    if (x > 0){
      stats = pt(q = x, df = PREFIX_Degr.freedom, lower.tail = FALSE)*2
    }
    return(stats)
  })
  PREFIX_P.val.Bacon = unlist(PREFIX_P.val.Bacon)
  PREFIX_Adj.P.val.Bacon = p.adjust(p = PREFIX_P.val.Bacon, method = "fdr")
  PREFIX_Top_table = cbind(PREFIX_Top_table,
                           PREFIX_Inflation.Bacon,
                           PREFIX_Bias.Bacon, 
                           PREFIX_T.Bacon,
                           PREFIX_P.val.Bacon,
                           PREFIX_Adj.P.val.Bacon)
  colnames(PREFIX_Top_table) = stri_replace_all_fixed(colnames(PREFIX_Top_table), pattern = "PREFIX_", replacement = "")
  PREFIX_Top_table_output = PREFIX_Top_table
  
  # Making volcano plot
  # Adding the highlight for transcript based on the fold change
  PREFIX_Top_table$is.highlight = sapply(PREFIX_Top_table$logFC, function(x){
    if (x > PREFIX_logFC_threshold){
      x = "Up"
    } else if (x < -PREFIX_logFC_threshold){
      x = "Down"
    } else {
      x = "None"
    }
  })
  # Modification of Calculated_stats
  Pval_treshold = PREFIX_Top_table[PREFIX_Top_table$adj.P.Val < 0.05,]
  if (nrow(Pval_treshold) < 1){
    writeLines("No Significant P-values")
    PREFIX_Top_table$is_annotate = "no"
    Pval_treshold = NULL
  } else {
    Pval_treshold = arrange(Pval_treshold, -P.Value)
    Pval_treshold = Pval_treshold$P.Value[1]
    Pval_treshold = -log10(Pval_treshold)
    PREFIX_Top_table = mutate(PREFIX_Top_table, is_annotate = ifelse(-log10(P.Value)>Pval_treshold & is.highlight != "None", "yes", "no"))
  }
  
  # Make the plot
  plot = ggplot(PREFIX_Top_table, aes(x=logFC, y=-log10(P.Value))) +
    
    # Show all points
    geom_point(aes(color= factor(is.highlight, levels = c("None", "Down", "Up"))), alpha=0.6, size=4) +
    scale_color_manual(values = c("grey", "skyblue", "red")) 
  
  # Add pval line
  if (!is.null(Pval_treshold)){
    plot = plot + geom_hline(yintercept=Pval_treshold, linetype="dashed", 
                             color = "red", size=0.5)
  }
  
  # Add pval line 0.05
  plot = plot + geom_hline(yintercept= -log10(0.05), linetype="dashed", 
                           color = "blue", size=0.5)
  # Add logFC lines
  if (min(PREFIX_Top_table$logFC) < -PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= -PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  if (max(PREFIX_Top_table$logFC) > PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  # Add label using ggrepel to avoid overlapping
  if (any(PREFIX_Top_table$is_annotate == "yes")){
    plot = plot + 
      geom_label_repel(data=subset(PREFIX_Top_table, is_annotate=="yes"), aes(label=ID), size=4, force = 10, 
                       max.overlaps = 50)
  } 
  
  plot = plot +
    labs(x = "log2 fold change", y = "-log10 p-value") +
    
    # Customizing the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
      panel.background = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  dir.create(PREFIX_plots_folder)
  
  file_name = paste0(PREFIX_plots_folder, "/", PREFIX_plots_str, "Volcano", "_",
                     paste0(PREFIX_contrast_vector, collapse = "_"), ".png")
  ggsave(file = file_name, plot = plot, width=2560, height=1440, units = "px", scale = 2)
  
  return(PREFIX_Top_table_output)
}

# A function to make Venn diagrams from named lists (uses ggplot2 and ggVennDiagram)
make_Venn_digram_list = function(named_list, plot_full_path = NULL, ...){
  
  # Customizable Venn diagram
  venn = Venn(named_list)
  data = process_data(venn)
  data@region$full_lable = sapply(data@region$count, function(x){
    Number = x
    Percent = x/sum(data@region$count)
    Percent = Percent*100
    Percent = round(Percent, digits = 1)
    Label_full = paste0(Number, "\n","(", Percent, "%)")
    return(Label_full)
  })
  
  # to see available shapes plot_shapes()
  plot = ggplot() +
    
    # 1. region count layer
    geom_sf(aes(fill = id), data = venn_region(data), ...) +
    
    # 2. set edge layer
    geom_sf(color="black", size = 0.5, data = venn_setedge(data), show.legend = FALSE, ...) +
    
    # 3. set label layer
    geom_sf_text(aes(label = name), data = venn_setlabel(data), ...) +
    
    # 4. region label layer
    geom_sf_label(aes(label = full_lable), data = venn_region(data), ...) +
    scale_fill_brewer(...) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(1,1,1,1, "cm")
    )
  print(plot)
  if (is.character(plot_full_path)){
    pdf(plot_full_path, width = 10, height = 10)
    print(plot)
    dev.off()
  }
}

# A function to get table with information on chromosomes
chromosome_table_getter = function(){
  URL = "https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37"
  HTML = read_html(URL)
  Tables = html_table(HTML)
  Tables = Tables[[1]]
  Tables$`Total length (bp)` = stri_replace_all_fixed(Tables$`Total length (bp)`, pattern = ',', replacement = '')
  Tables$`Total length (bp)` = as.numeric(Tables$`Total length (bp)`)
  write.csv(Tables, "Chromosome_coord_table.csv")
  return(Tables)
}

# function to perform enrichment for GO terms and KEGG and create figures that may be useful
# all results are saved in the specified folder
# genes and universe are accepted as vectors of Entrez IDs
# some of the ploting is performed within tryCatch since these plots may not be rendered depending on the enrichment results
run_enrichment_GO_KEGG_gene_set = function(genes, universe, categories_to_show = 30, folder, plot_name_pref){
  
  # genes should be submitted as Entrez IDs
  if (length(genes) < 10) {
    Minsize = length(genes)
  } else {
    Minsize = 10
  }
  
  # IDs should be characters
  genes = as.character(genes)
  universe = as.character(universe)
  
  # Running GO analysis
  GO_String = c("BP", "MF", "CC")
  GO_Enrichment = lapply(GO_String, function(x){
    result = clusterProfiler::enrichGO(gene = genes, OrgDb = 'org.Hs.eg.db', ont = x, universe = universe, minGSSize = Minsize)
    return(result)
  })
  KEGG_Enrichment = clusterProfiler::enrichKEGG(gene = genes,  organism = 'hsa' , universe = universe, minGSSize = Minsize)
  
  # Saving GO results
  dir.create(folder)
  combined_enrich = list()
  for (i in 1:length(GO_String)){
    filename = paste0(plot_name_pref, "_GO_", GO_String[i], ".xlsx")
    filename =  paste0(folder, "/", filename)
    Curr_GO = GO_Enrichment[[i]]
    Curr_GO_readable = setReadable(Curr_GO, 'org.Hs.eg.db', 'ENTREZID')
    Curr_GO_result = Curr_GO@result
    Curr_GO_result$geneSymbol = Curr_GO_readable@result$geneID
    combined_enrich[[i]] = Curr_GO_result
    openxlsx::write.xlsx(x = Curr_GO_result, file = filename, overwrite = TRUE)
  }
  
  # Saving KEGG results
  if (!is.null(KEGG_Enrichment)){
    KEGG_Enrichment_readable = setReadable(KEGG_Enrichment, 'org.Hs.eg.db', 'ENTREZID')
    KEGG_Enrichment_result = KEGG_Enrichment@result
    KEGG_Enrichment_result$geneSymbol = KEGG_Enrichment_readable@result$geneID
    combined_enrich[[4]] = KEGG_Enrichment_result
    openxlsx::write.xlsx(x = KEGG_Enrichment_result, file = paste0(folder, "/", plot_name_pref, "_KEGG.xlsx"), overwrite = TRUE)
  } else {
    combined_enrich[[4]] = NA
  }
  # Saving full results
  openxlsx::write.xlsx(x = combined_enrich, file = paste0(folder, "/", plot_name_pref, "_enrichment_combined_sheet.xlsx"), overwrite = TRUE)
  
  # Making plots
  GO_Enrichment_readable = lapply(GO_Enrichment, function(x) setReadable(x, 'org.Hs.eg.db', 'ENTREZID'))
  
  if (!is.null(KEGG_Enrichment)){
    GO_Enrichment_readable[[4]] = KEGG_Enrichment_readable
  } else {
    GO_Enrichment_readable[[4]] = NA
  }
  
  GO_String = c("BP", "MF", "CC", "KEGG")
  gc()
  
  for (i in 1:length(GO_String)){
    Curr_GO = GO_Enrichment_readable[[i]]
    
    if (is.na(Curr_GO)){
      writeLines(paste0(filename, " has NO results and analysis was not conducted"))
      next
    }
    
    filename = paste0(plot_name_pref, "_GO_", GO_String[i])
    filename =  paste0(folder, "/", filename)
    
    if (GO_String[i] == "KEGG"){
      # Figures for KEGG
      dimensions = c(19*nrow(summary(Curr_GO))/categories_to_show, 11*nrow(summary(Curr_GO))/categories_to_show)
      
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      if (nrow(summary(Curr_GO)) < 1){
        writeLines(paste0(filename, " has NO significant results"))
        next
      }
      
      # Barplot
      pdf(paste0(filename, "_bar.pdf")
          , width = 10, height = 10)
      print(barplot(Curr_GO, showCategory = categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_bar_full.pdf")
          , width = 10, height = 10*nrow(summary(Curr_GO))/categories_to_show)
      print(barplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Dotplot
      pdf(paste0(filename, "_dot.pdf")
          , width = 10, height = 10)
      print(dotplot(Curr_GO, showCategory= categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_dot_full.pdf")
          , width = 10, height = dimensions[2])
      print(dotplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Gene-Concept Network
      pdf(paste0(filename, "_gene_conc.pdf")
          , width = 19, height = 11)
      print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                     showCategory = categories_to_show, cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
      dev.off()
      
      tryCatch({
        pdf(paste0(filename, "_gene_conc_full.pdf")
            , width = dimensions[1], height = dimensions[2])
        print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                       showCategory = nrow(summary(Curr_GO)), cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
        dev.off()
      }, error = function(e) {writeLines(paste0("Full gene conc plot is not available for ", filename))})
      
    } else {
      # Figures for GO
      dimensions = c(19*nrow(summary(Curr_GO))/categories_to_show, 11*nrow(summary(Curr_GO))/categories_to_show)
      
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      if (nrow(summary(Curr_GO)) < 1){
        writeLines(paste0(filename, " has NO significant results"))
        next
      }
      
      # Barplot
      pdf(paste0(filename, "_bar.pdf")
          , width = 10, height = 10)
      print(barplot(Curr_GO, showCategory = categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_bar_full.pdf")
          , width = 10, height = dimensions[2])
      print(barplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Dotplot
      pdf(paste0(filename, "_dot.pdf")
          , width = 10, height = 10)
      print(dotplot(Curr_GO, showCategory= categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_dot_full.pdf")
          , width = 10, height = dimensions[2])
      print(dotplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Gene-Concept Network
      pdf(paste0(filename, "_gene_conc.pdf")
          , width = 19, height = 11)
      print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                     showCategory = categories_to_show, cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
      dev.off()
      
      tryCatch({
        pdf(paste0(filename, "_gene_conc_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                       showCategory = nrow(summary(Curr_GO)), cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
        dev.off()
      }, error = function(e) {writeLines(paste0("Full gene conc plot is not available for ", filename))})
      
      #GO induced graph
      tryCatch({
        pdf(paste0(filename, "_go_graph.pdf"), 
            width = 19, height = 11)
        print(goplot(Curr_GO, showCategory = categories_to_show))
        dev.off()
      }, error = function(e) {writeLines(paste0("Go graph is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_go_graph_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(goplot(Curr_GO, showCategory = nrow(summary(Curr_GO))))
        dev.off()
      }, error = function(e) {writeLines(paste0("Go graph (big) is not available for ", filename))})
      
      #Tree plot
      Curr_GO_pairwise = enrichplot::pairwise_termsim(Curr_GO)
      dimensions = c(19*nrow(Curr_GO_pairwise@termsim)/categories_to_show, 11*nrow(Curr_GO_pairwise@termsim)/categories_to_show)
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      tryCatch({
        pdf(paste0(filename, "_tree.pdf"), 
            width = 19, height = 11)
        print(treeplot(Curr_GO_pairwise, showCategory = categories_to_show, nCluster = 10))
        dev.off()
      }, error = function(e) {writeLines(paste0("Tree plot is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_tree_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(treeplot(Curr_GO_pairwise, showCategory = nrow(Curr_GO_pairwise@termsim), nCluster = 10))
        dev.off()
      }, error = function(e) {writeLines(paste0("Tree plot (big) is not available for ", filename))})
      
      
      #Enrichment map
      tryCatch({
        pdf(paste0(filename, "_enr_map.pdf"), 
            width = 19, height = 11)
        print(emapplot(Curr_GO_pairwise, showCategory = categories_to_show, layout = "fr", cex_label_category = 0.9, 
                       cex_line = 0.8, repel = TRUE))
        dev.off()
      }, error = function(e) {writeLines(paste0("Enrich. map is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_enr_map_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(emapplot(Curr_GO_pairwise, showCategory = nrow(Curr_GO_pairwise@termsim), layout = "fr", cex_label_category = 0.9, 
                       cex_line = 0.8, repel = TRUE))
        dev.off()
      }, error = function(e) {writeLines(paste0("Enrich. map is not available for ", filename))})
      
    }
  }
  
  # preparing outputs to load into the environment
  GO_Enrichment[[4]] = KEGG_Enrichment
  
  if (is.null(KEGG_Enrichment)){
    GO_Enrichment[[4]] = NA
  }
  
  Output = list()
  Output[[1]] = GO_Enrichment
  Output[[2]] = GO_Enrichment_readable
  Output[[3]] = combined_enrich
  names(Output) = c("Enrichment", "Mapped Enrichment", "Combined Results")
  
  for (i in 1:length(Output)){
    names(Output[[i]]) = GO_String
  }
  
  return(Output)
}

# A function to create demographic characteristics for a dataset
characterize_dataset_generelized_two_subgroups = function(dataset,
                                                          study_char,
                                                          contrast_col_number,
                                                          contrast_vector,
                                                          participants_col_number,
                                                          model_covariates_col_vector,
                                                          columns_to_characterise_vector,
                                                          Remove_NA_predictors = TRUE,
                                                          drop_P = TRUE,
                                                          simplif_P = 3){
  
  # Running Raw Statistics
  total_number = nrow(dataset)
  total_number_message = paste0("Initial dataset includes ", total_number, " participants")
  writeLines(total_number_message)
  first_row_output = data.frame(Name = total_number_message,Category.1 = "", Category.2 = "", P.value = NA)
  
  # Removing missing data
  colnames(dataset)[participants_col_number] = "Participant_ID_FUN"
  indeces = c(participants_col_number, contrast_col_number, model_covariates_col_vector)
  indeces = unique(indeces)
  Model_df = dataset[,indeces]
  if (Remove_NA_predictors){
    
    writeLines("Removing participants with missing data")
    Participants_missing_data = apply(Model_df, 1, function(x){
      x = as.character(x)
      
      if (any(is.na(x))){
        return(TRUE)
      } else{
        return(FALSE)
      }
      
    })
    Participants_missing_data = Model_df$Participant_ID_FUN[Participants_missing_data]
    
    if (length(Participants_missing_data) < 1){
      excluded_number = 0
      
    } else {
      
      excluded_number = length(Participants_missing_data)
    }
    
    Participants_excluded_df = dataset[dataset$Participant_ID_FUN %in% Participants_missing_data, ]
    dataset = dataset[dataset$Participant_ID_FUN %!in% Participants_missing_data, ]
    excluded_message = paste0("Participants with missing data excluded: ", excluded_number,"\n",
                              "Resulting number of participants: ", nrow(dataset))
    second_row_output = data.frame(Name = excluded_message,Category.1 = "", Category.2 = "", P.value = NA)
    writeLines(excluded_message)
    
  } else {
    
    writeLines("Exclusion of participants was not performed")
    second_row_output = data.frame(Name = "Exclusion of participants was not performed", Category.1 = "", Category.2 = "", P.value = NA)
    Participants_excluded_df = NA
  }
  
  main_indeces = c(contrast_col_number, columns_to_characterise_vector)
  main_indeces = unique(main_indeces)
  Main_dataset = list()
  
  
  # Making categories
  Initial_names = colnames(dataset)[main_indeces]
  colnames(dataset)[contrast_col_number] = "variable_to_split"
  levels_var = dataset[,"variable_to_split"]
  levels_var = levels(levels_var)
  names(levels_var) = NULL
  Category_1 = levels_var[1]
  Category_2 = levels_var[2]
  
  for (i in 1:length(main_indeces)){
    curr_variable = dataset[, main_indeces[i]]
    if (is.factor(curr_variable)){
      Characterised_var = characterize_categorical_variable(df = dataset,
                                                            Category_1 = Category_1, 
                                                            Category_2 = Category_2, 
                                                            Variable_name = colnames(dataset)[main_indeces[i]], 
                                                            keep_missing = TRUE)
    } else if (is.character(curr_variable)){
      stop("All columns should be either a Factor or Numeric")
    } else if (is.numeric(curr_variable)){
      Characterised_var = characterize_numeric_variable(df = dataset, 
                                                        Category_1 = Category_1, 
                                                        Category_2 = Category_2, 
                                                        Variable_name = colnames(dataset)[main_indeces[i]], 
                                                        Mention_NAs = TRUE)
    }
    Current_DF = data.frame(Name = Initial_names[i], 
                            Category.1 = Characterised_var[[1]], 
                            Category.2 = Characterised_var[[2]], 
                            P.value = Characterised_var[[3]])
    Main_dataset[[i]] = Current_DF
  }
  Main_dataset = list_to_df(data_list = Main_dataset)
  Main_dataset$P.value = round(Main_dataset$P.value, digits = simplif_P)
  
  # Compiling full dataframe
  header = data.frame(Name = study_char, Category.1 = "", Category.2 = "", P.value = NA)
  Full_table = rbind(header, first_row_output, second_row_output, Main_dataset)
  
  if(drop_P){
    Full_table$P.value = NULL
  }
  
  output = list()
  output$Table = Full_table
  output$Excluded = Participants_excluded_df
  return(output)
}
# Several helper functions for characterization (used internally in characterize_dataset_generelized_two_subgroups)
characterize_categorical_variable = function(df, Category_1, Category_2, Variable_name, P.val.valid = TRUE, keep_missing = FALSE){
  Total_variable = df[, Variable_name]
  
  if (!is.ordered(Variable_name)){
    df[, Variable_name] = ordered(Total_variable, levels = sort(unique(Total_variable)))
    Total_variable = df[, Variable_name]
  }
  
  Contingency_table = table(df[,Variable_name], df[,"variable_to_split"])
  
  if (keep_missing){
    Contingency_table = table(df[,Variable_name], df[,"variable_to_split"], exclude = NULL)
  }
  
  Contingency_table_percents_1 = (Contingency_table[,Category_1]/(sum(Contingency_table[,Category_1])) * 100) %>% round(., digits = 1)
  Contingency_table_percents_2 = (Contingency_table[,Category_2]/(sum(Contingency_table[,Category_2])) * 100) %>% round(., digits = 1)
  Report_1 = paste0(names(Contingency_table_percents_1), ": ")
  Report_1 = paste0(Report_1, Contingency_table[,Category_1], " (",Contingency_table_percents_1, "%)", collapse = "\n")
  Report_2 = paste0(names(Contingency_table_percents_2), ": ")
  Report_2 = paste0(Report_2, Contingency_table[,Category_2], " (",Contingency_table_percents_2, "%)", collapse = "\n")
  
  Contingency_table_pval = table(df[,Variable_name], df[,"variable_to_split"])
  if (P.val.valid){
    tryCatch({
      Var.pval <<- chisq.test(Contingency_table_pval)
      Var.pval <<- Var.pval$p.value}, warning = function(w){
        message = paste0(Variable_name, " has small counts in some groups, using P.val from Fisher's Exact Test")
        writeLines(message)
        Var.pval <<- fisher.test(Contingency_table_pval)
        Var.pval <<- Var.pval$p.value
      })
  } else {
    Var.pval = "Not valid"
  }
  
  Output_list = list()
  Output_list[[1]] = Report_1
  Output_list[[2]] = Report_2
  Output_list[[3]] = Var.pval
  return(Output_list)
}

characterize_numeric_variable = function(df, Category_1, Category_2, Variable_name, Mention_NAs = FALSE){
  Total_variable = df[, Variable_name]
  Variable_1 = df[df$variable_to_split == Category_1, Variable_name]
  Variable_1_report = describe_vector_numeric(Variable_1, Mention_NAs = Mention_NAs)
  Variable_2 = df[df$variable_to_split == Category_2, Variable_name]
  Variable_2_report = describe_vector_numeric(Variable_2, Mention_NAs = Mention_NAs)
  Shapiro_normality_check = shapiro.test(Total_variable)
  Shapiro_normality_check = Shapiro_normality_check$p.value
  if (Shapiro_normality_check < 0.05){
    message = paste0(Variable_name, " is not normally distributed. Using P-val from Mann Whitney U Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Total_pval = wilcox.test(formula = as.formula(formula_test), data = df)
    Total_pval = Total_pval$p.value
  } else {
    message = paste0(Variable_name, " is normally distributed. Using P-val from T Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Var_test_check = var.test(formula	= as.formula(formula_test), data = df)
    Var_test_check = Var_test_check$p.value 
    if (Var_test_check >= 0.05){
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = TRUE)
    } else  {
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = FALSE)
    }
    Total_pval = Total_pval$p.value
  }
  Output_list = list()
  Output_list[[1]] = Variable_1_report
  Output_list[[2]] = Variable_2_report
  Output_list[[3]] = Total_pval
  return(Output_list)
}

describe_vector_numeric = function(x, Mention_NAs, show_as_percent = FALSE, Describe_min_max = TRUE){
  Mean_var = mean(x, na.rm = TRUE) %>% round(., digits = 2)
  SD_var = sd(x, na.rm = TRUE) %>% round(., digits = 2)
  Var_report = paste0(Mean_var, " ± ", SD_var)
  if (Mention_NAs){
    Missing_count = length(x[is.na(x)])
    Missing_val_percent = (Missing_count/length(x) * 100) %>% round(., digits = 2)
    if (Missing_count > 0){
      Var_report = paste0(Var_report, "\nMissing val: ", Missing_count, " (", Missing_val_percent, "%)")
    }
  }
  if (show_as_percent){
    Mean_var = Mean_var * 100
    SD_var = SD_var * 100
    Var_report = paste0(Mean_var, " ± ", SD_var, "%")
  }
  if (Describe_min_max){
    Min_var = min(x, na.rm = TRUE) %>% round(., digits = 2)
    Max_var = max(x, na.rm = TRUE) %>% round(., digits = 2)
    if (show_as_percent){
      Min_var = Min_var * 100
      Max_var = Max_var * 100
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, "%, Max: ",Max_var, "%")
    } else {
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, ", Max: ",Max_var)
    }
  }
  return(Var_report)
}


################### Importing and preparing GSE98793 data ###################

# Reading files
Eset_GSE98793_expresion_full = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE98793_MDD_blood/Eset_GSE98793_expression_combat_model.csv") # Replace with an appropriate path
GSE98793_phenotypes = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE98793_MDD_blood/GSE98793_phenotypes.csv") # Replace with an appropriate path
GSE98793_probes = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE98793_MDD_blood/GSE98793_probes.csv") # Replace with an appropriate path

# Correcting phenotypes
GSE98793_phenotypes$batch = factor(GSE98793_phenotypes$batch)
GSE98793_phenotypes$age = as.numeric(GSE98793_phenotypes$age)
GSE98793_phenotypes$anxiety = factor(GSE98793_phenotypes$anxiety, levels = c("yes", "no"))
GSE98793_phenotypes$gender = factor(GSE98793_phenotypes$gender)
GSE98793_phenotypes$subject_group =  factor(GSE98793_phenotypes$subject_group, levels = c("MDD", "Control"))
str(GSE98793_phenotypes)


################### Importing and preparing GSE46743 data ###################

# Reading files
GSE46743_data_processed_baseline = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE46743_MDD_control_dexamethasone_exposure/GSE46743_data_processed_baseline.csv") # Replace with an appropriate path
GSE46743_phenotypes_table_baseline = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE46743_MDD_control_dexamethasone_exposure/GSE46743_phenotypes_table_baseline.csv") # Replace with an appropriate path
GSE46743_annotation = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE46743_MDD_control_dexamethasone_exposure/GSE46743_annotation.csv") # Replace with an appropriate path

# Correcting phenotypes
GSE46743_phenotypes_table_baseline$Depr_status = factor(GSE46743_phenotypes_table_baseline$Depr_status, levels = c("case", "control"))
str(GSE46743_phenotypes_table_baseline)


################### Importing and preparing GSE64930 data ###################

# Reading files
GSE64930_data_processed_baseline = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/GSE64930_data_processed_baseline.csv") # Replace with an appropriate path
GSE64930_phenotypes_baseline = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/GSE64930_phenotypes_baseline.csv") # Replace with an appropriate path
GSE64930_annotation = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/GSE64930_MDD_control_dexamethasone_exposure/GSE64930_annotation.csv") # Replace with an appropriate path
table(GSE64930_phenotypes_baseline$Sex, GSE64930_phenotypes_baseline$Status)

# Correcting phenotypes
GSE64930_phenotypes_baseline$Status = factor(GSE64930_phenotypes_baseline$Status , levels = c("case", "control"))
GSE64930_phenotypes_baseline$Sex = factor(GSE64930_phenotypes_baseline$Sex)
str(GSE64930_phenotypes_baseline)


################### Importing and preparing UCSC tracks ###################

# KnownGene track
UCSC_Known_Gene = smart_fread("/home/aleksandr/Desktop/WORK/KNOWN_GENE_TRACK_UCSC/Known_gene_track.txt") # Replace with another path if needed
Mapped_UCSC_Known_Gene = smart_fread("/home/aleksandr/Desktop/WORK/KNOWN_GENE_TRACK_UCSC/Known_to_Symbol_kgxref.txt") # Replace with another path if needed
UCSC_Known_Gene_full = inner_join(UCSC_Known_Gene, Mapped_UCSC_Known_Gene, by = c("#name" = "#kgID"))

# NCBI RefSeq track
UCSC_NCBI_RefSeq = smart_fread("/home/aleksandr/Desktop/WORK/KNOWN_GENE_TRACK_UCSC/ncbiRefSeq.txt") # Replace with another path if needed

# GENCODE lifted annotations from V38lift37 (Ensembl 104) track
UCSC_wgEncodeGencodeCompV38lift37 = smart_fread("/home/aleksandr/Downloads/wgEncodeGencodeCompV38lift37.txt") # Replace with another path if needed
UCSC_wgEncodeGencodeCompV38lift37_pseudo = smart_fread("/home/aleksandr/Downloads/wgEncodeGencodeCompV38lift37_pseudo.txt") # Replace with another path if needed
UCSC_wgEncodeGencodeCompV38lift37 = rbind(UCSC_wgEncodeGencodeCompV38lift37, UCSC_wgEncodeGencodeCompV38lift37_pseudo)


################### Filtering and mapping U133 Plus 2.0 Array probes ###################

# Note: The associated product files could be downloaded from https://www.thermofisher.com/order/catalog/product/900466

# Getting Probe sequences
U133_Plus_2_Probe_Seq = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/Expression_Cohorts_Analysis/Affymextrix_hgu_133_plus_files/HG-U133_Plus_2.probe_tab") # Replace with an appropriate path
U133_Plus_2_FASTA_targets = seqinr::read.fasta("/home/aleksandr/Desktop/WORK/Expression Cohorts/Expression_Cohorts_Analysis/Affymextrix_hgu_133_plus_files/HG-U133_Plus_2.target") # Replace with an appropriate path
U133_Plus_2_FASTA_Consensus_seq = seqinr::read.fasta("/home/aleksandr/Desktop/WORK/Expression Cohorts/Expression_Cohorts_Analysis/Affymextrix_hgu_133_plus_files/HG-U133_Plus_2.consensus") # Replace with an appropriate path
U133_Plus_2_FASTA_Exemplar_seq = seqinr::read.fasta("/home/aleksandr/Desktop/WORK/Expression Cohorts/Expression_Cohorts_Analysis/Affymextrix_hgu_133_plus_files/HG-U133_Plus_2_exemplar") # Replace with an appropriate path
U133_Plus_2_FASTA_targets_names = names(U133_Plus_2_FASTA_targets)
U133_Plus_2_FASTA_targets_names = sapply(U133_Plus_2_FASTA_targets_names, function(x){
  x = unlist(stri_split_fixed(x, pattern = ":"))
  x = x[3]
  x = stri_replace_all_fixed(x, pattern = ";", replacement = "")
  x = str_trim(x)
  return(x)
})

# DNA Sequence Alignment Scoring Matrix
Scoring_EDNAFULL_NUC4.4 = smart_fread("/home/aleksandr/Desktop/WORK/Expression Cohorts/Expression_Cohorts_Analysis/Affymextrix_hgu_133_plus_files/NUC.4.4.txt") # Replace with an appropriate path
Scoring_EDNAFULL_NUC4.4 = as.matrix(Scoring_EDNAFULL_NUC4.4)

# Preparing annotation for U133 Plus 2.0 Array probes
U133_Plus_2_Probes_Annotation = GSE98793_probes
U133_Plus_2_Probes_Annotation = U133_Plus_2_Probes_Annotation[U133_Plus_2_Probes_Annotation$`Gene Symbol` != "",]
U133_Plus_2_Probes_Annotation = U133_Plus_2_Probes_Annotation[!is.na(U133_Plus_2_Probes_Annotation$`Gene Symbol` != ""),]
U133_Plus_2_Probes_Annotation = U133_Plus_2_Probes_Annotation[!stri_detect_fixed(U133_Plus_2_Probes_Annotation$`Gene Symbol`, pattern = " /// "),]
# 42700 Updated probes

# Updating gene names
U133_Plus_2_Probes_Annotation_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = U133_Plus_2_Probes_Annotation$`Gene Symbol`, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                                  PRF_replace_NA_with_old = TRUE)
probes_genes_check_dict = U133_Plus_2_Probes_Annotation_genes_check$Suggested.Symbol
names(probes_genes_check_dict) = U133_Plus_2_Probes_Annotation_genes_check$x
U133_Plus_2_Probes_Annotation$Updated_gene_names = sapply(U133_Plus_2_Probes_Annotation$`Gene Symbol`, function(x){
  if (is.na(x)){return(NA)}
  if (x == ""){return(NA)}
  Current_genes = unlist(stri_split_fixed(x, pattern = " /// "))
  Current_genes = unique(Current_genes)
  Current_genes = str_trim(Current_genes)
  New_names = probes_genes_check_dict[Current_genes]
  New_names = paste0(New_names, collapse = ";")
  return(New_names)
})

# Mapping loop
dir.create("Mapping_Affymetrix_Probes_Alignments") # A folder to store alignments

# Note: The loop is using UCSC API and downloading sequences -> constant internet connection is required
# Note: The loop takes several hours to execute
# Note: The loop variables are labelled with the prefix "PRF_" to enable easy cleanup
i = 1
while (i <= nrow(U133_Plus_2_Probes_Annotation)){
  print(i)
  PRF_Current_gene = U133_Plus_2_Probes_Annotation$Updated_gene_names[i]
  PRF_Current_UCSC_Known_Gene_df = UCSC_Known_Gene_full[UCSC_Known_Gene_full$geneSymbol == PRF_Current_gene,]
  PRF_Ref_DS = "knownGene"
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    PRF_Current_UCSC_Known_Gene_df = UCSC_NCBI_RefSeq[UCSC_NCBI_RefSeq$name2 == PRF_Current_gene,]
    PRF_Ref_DS = "ncbiRefSeq"
  }
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    PRF_Current_UCSC_Known_Gene_df = UCSC_wgEncodeGencodeCompV38lift37[UCSC_wgEncodeGencodeCompV38lift37$name2 == PRF_Current_gene,]
    PRF_Ref_DS = "wgEncodeGencodeCompV38lift37"
  }
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    PRF_Ref_DS = "NOT IDENTIFIED"
    Output_Data = data.frame(
      ID = U133_Plus_2_Probes_Annotation$ID[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = NA,
      txStart = NA,
      txEnd = NA,
      PRF_URL = NA,
      PRF_DNA_init = NA,
      PRF_DNA_bio_str = NA,
      PRF_Strand = NA,
      PRF_Target_seq = NA,
      PRF_DNA_transcript = NA,
      PRF_DNA_R_C_str = NA,
      PRF_DNA_transcript_bio_R_C_str = NA,
      PRF_Align_Start_Gene = NA,
      PRF_Align_Start_Gene_genomic = NA,
      PRF_Align_Start_Target = NA,
      PRF_Align_End_Target = NA,
      PRF_Align_End_Target_Expected = NA,
      PRF_Align_End_Target_genomic = NA,
      PRF_Alignment_score = NA
    )
    write.csv(Output_Data, file = paste0("Mapping_Affymetrix_Probes_Alignments/",i, "__",U133_Plus_2_Probes_Annotation$ID[i], ".csv"))
    i = i + 1
    next
  }
  
  # Selecting Best Gene representation
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[PRF_Current_UCSC_Known_Gene_df$exonCount == max(PRF_Current_UCSC_Known_Gene_df$exonCount),]
  PRF_Current_UCSC_Known_Gene_df_lengths = mapply(function(x,y){
    Len = x:y
    Len = length(Len)
    return(Len)
  }, PRF_Current_UCSC_Known_Gene_df$cdsStart, PRF_Current_UCSC_Known_Gene_df$cdsEnd)
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[which.max(PRF_Current_UCSC_Known_Gene_df_lengths),]
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[1,]
  
  # Getting Sequence
  PRF_URL = paste0("https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=",
                   PRF_Current_UCSC_Known_Gene_df$chrom,
                   ";start=",
                   PRF_Current_UCSC_Known_Gene_df$txStart,
                   ";end=",
                   PRF_Current_UCSC_Known_Gene_df$txEnd)
  PRF_EROR_FLAG = FALSE
  PRF_EROR_Counter = 0
  
  tryCatch({
    PRF_curr_DNA = jsonlite::fromJSON(PRF_URL)
    PRF_curr_DNA = PRF_curr_DNA$dna
  }
  , error = function(e) {PRF_EROR_FLAG <<-TRUE})
  
  while(PRF_EROR_FLAG){
    PRF_EROR_Counter = PRF_EROR_Counter + 1
    warning(paste0('Error detected  during downloading file from   ', PRF_URL, '  ', Sys.time(), ' Trying to reopen link  #', PRF_EROR_Counter), immediate. = TRUE)
    if (PRF_EROR_Counter > 20){
      PRF_curr_DNA = NA
      break
    }
    if (PRF_EROR_Counter > 15) {
      Sys.sleep(sample(1200:3600, 1))
      closeAllConnections()
    }
    Sys.sleep(5*PRF_EROR_Counter)
    tryCatch({
      PRF_curr_DNA = jsonlite::fromJSON(PRF_URL)
      PRF_curr_DNA = PRF_curr_DNA$dna
      PRF_EROR_FLAG = FALSE}
      , error = function(e) {PRF_EROR_FLAG <<-TRUE})
  }
  
  if (is.na(PRF_curr_DNA)|PRF_curr_DNA == ""|PRF_curr_DNA == " "){
    Output_Data = data.frame(
      ID = U133_Plus_2_Probes_Annotation$ID[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = PRF_Current_UCSC_Known_Gene_df$chrom,
      txStart = PRF_Current_UCSC_Known_Gene_df$txStart,
      txEnd = PRF_Current_UCSC_Known_Gene_df$txEnd,
      PRF_URL = PRF_URL,
      PRF_DNA_init = NA,
      PRF_DNA_bio_str = NA,
      PRF_Strand = NA,
      PRF_Target_seq = NA,
      PRF_DNA_transcript = NA,
      PRF_DNA_R_C_str = NA,
      PRF_DNA_transcript_bio_R_C_str = NA,
      PRF_Align_Start_Gene = NA,
      PRF_Align_Start_Gene_genomic = NA,
      PRF_Align_Start_Target = NA,
      PRF_Align_End_Target = NA,
      PRF_Align_End_Target_Expected = NA,
      PRF_Align_End_Target_genomic = NA,
      PRF_Alignment_score = NA
    )
    write.csv(Output_Data, file = paste0("Mapping_Affymetrix_Probes_Alignments/",i, "__",U133_Plus_2_Probes_Annotation$ID[i], ".csv"))
    i = i + 1
  } else {
    PRF_curr_DNA_init = PRF_curr_DNA
    PRF_Strand = PRF_Current_UCSC_Known_Gene_df$strand
    PRF_curr_Target_seq = U133_Plus_2_Probes_Annotation$ID[i]
    PRF_curr_Target_seq = U133_Plus_2_FASTA_targets[U133_Plus_2_FASTA_targets_names == PRF_curr_Target_seq]
    PRF_curr_Target_seq = PRF_curr_Target_seq[[1]]
    PRF_curr_Target_seq = toupper(c2s(PRF_curr_Target_seq))
    PRF_curr_Target_seq_bio = DNAString(PRF_curr_Target_seq)
    
    # Getting Full Transcripts and Spliced Transcripts
    PRF_Exon_Starts = unlist(stri_split_fixed(PRF_Current_UCSC_Known_Gene_df$exonStarts, pattern = ","))
    PRF_Exon_Starts = PRF_Exon_Starts[PRF_Exon_Starts != ""]
    PRF_Exon_Ends = unlist(stri_split_fixed(PRF_Current_UCSC_Known_Gene_df$exonEnds, pattern = ","))
    PRF_Exon_Ends = PRF_Exon_Ends[PRF_Exon_Ends != ""]
    
    # Getting Spliced Sequence
    PRF_Exon_Coords = mapply(function(x,y){
      Vector = x:y
      return(Vector)
    }, PRF_Exon_Starts, PRF_Exon_Ends)
    PRF_Exon_Coords = unlist(PRF_Exon_Coords) # it is sorted
    names(PRF_Exon_Coords) = NULL
    
    # Getting coordinates  for a sequence
    PRF_curr_DNA_init_vector = unlist(strsplit(PRF_curr_DNA_init, split = "", fixed = TRUE))
    PRF_Gene_Coords = 1:length(PRF_curr_DNA_init_vector)
    PRF_Gene_Coords_inverse = length(PRF_curr_DNA_init_vector):1
    PRF_Gene_Coords_Genomic = PRF_Gene_Coords + (PRF_Current_UCSC_Known_Gene_df$txStart)
    
    # Getting spliced coordinates
    PRF_Splice_Coords_Ticker = 0
    PRF_Splice_Coords_Vector = vector()
    
    for (i_spl in 1:length(PRF_Gene_Coords_Genomic)){
      PRF_Curr_coord = PRF_Gene_Coords_Genomic[i_spl]
      if (PRF_Curr_coord %in% PRF_Exon_Coords){
        PRF_Spl_coord = PRF_Splice_Coords_Ticker[length(PRF_Splice_Coords_Ticker)] + 1
        PRF_Splice_Coords_Vector[i_spl] = PRF_Spl_coord
        PRF_Splice_Coords_Ticker = c(PRF_Splice_Coords_Ticker, PRF_Spl_coord)
      } else {
        PRF_Splice_Coords_Vector[i_spl] = 0
      }
    }
    
    PRF_Coord_matrix = cbind(PRF_Gene_Coords, PRF_Gene_Coords_inverse, PRF_Gene_Coords_Genomic, PRF_Splice_Coords_Vector)
    PRF_Coord_matrix = as.data.frame(PRF_Coord_matrix)
    PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse = NA
    PRF_TR_len = length(which(PRF_Coord_matrix$PRF_Splice_Coords_Vector != 0))
    PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse = sapply(PRF_Coord_matrix$PRF_Splice_Coords_Vector, function(x){
      if(x == 0){
        return(0)
      }
      x = abs(x-PRF_TR_len) + 1
      return(x)
    })
    PRF_Coord_matrix$DNA_init = PRF_curr_DNA_init_vector
    PRF_Coord_matrix$DNA_UP = toupper(PRF_curr_DNA_init_vector)
    PRF_curr_DNA_transcript = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector != 0, "DNA_UP"]
    PRF_curr_DNA_transcript = paste0(PRF_curr_DNA_transcript, collapse = "")
    if (PRF_Strand == "-"){
      PRF_curr_DNA_bio_str = toupper(PRF_curr_DNA_init)
      PRF_curr_DNA_bio = DNAString(PRF_curr_DNA_bio_str)
      PRF_curr_DNA_R_C = reverseComplement(PRF_curr_DNA_bio)
      PRF_curr_DNA_transcript_bio = DNAString(PRF_curr_DNA_transcript)
      PRF_curr_DNA_transcript_bio_R_C = reverseComplement(PRF_curr_DNA_transcript_bio)
      PRF_Alignment = pairwiseAlignment(PRF_curr_Target_seq_bio, PRF_curr_DNA_transcript_bio_R_C, substitutionMatrix = Scoring_EDNAFULL_NUC4.4, type = "overlap", gapOpening = 10, gapExtension = 0.001)
      
      # Returning to characters
      PRF_curr_DNA_R_C_str = as.character(PRF_curr_DNA_R_C)
      PRF_curr_DNA_transcript_bio_R_C_str = as.character(PRF_curr_DNA_transcript_bio_R_C)
      
      # Calculating appropriate position
      PRF_Align_Start_Gene = PRF_Alignment@subject@range@start
      PRF_Align_Start_Gene_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == PRF_Align_Start_Gene,"PRF_Gene_Coords_Genomic"]
      PRF_Align_Start_Target = PRF_Alignment@pattern@range@start
      PRF_Align_End_Target = PRF_Alignment@pattern@range@width
      PRF_Align_End_Target_Expected = nchar(PRF_curr_Target_seq)
      PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == (PRF_Alignment@pattern@range@width + PRF_Align_Start_Gene - 1), "PRF_Gene_Coords_Genomic"]
      
      if (length(PRF_Align_End_Target_genomic) < 1){
        PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == max(PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse), "PRF_Gene_Coords_Genomic"]
      }
      
      # Getting Alignment Scores
      PRF_Alignment_score = PRF_Alignment@score
    } else {
      PRF_curr_DNA_bio_str = toupper(PRF_curr_DNA_init)
      PRF_curr_DNA_bio = DNAString(PRF_curr_DNA_bio_str)
      PRF_curr_DNA_R_C = NA
      PRF_curr_DNA_transcript_bio = DNAString(PRF_curr_DNA_transcript)
      PRF_curr_DNA_transcript_bio_R_C = NA
      PRF_Alignment = pairwiseAlignment(PRF_curr_Target_seq_bio, PRF_curr_DNA_transcript_bio, substitutionMatrix = Scoring_EDNAFULL_NUC4.4, type = "overlap", gapOpening = 10, gapExtension = 0.001)
      
      # Returning to characters
      PRF_curr_DNA_R_C_str = NA
      PRF_curr_DNA_transcript_bio_R_C_str = NA
      
      # Calculating appropriate position
      PRF_Align_Start_Gene = PRF_Alignment@subject@range@start
      PRF_Align_Start_Gene_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == PRF_Align_Start_Gene,"PRF_Gene_Coords_Genomic"]
      PRF_Align_Start_Target = PRF_Alignment@pattern@range@start
      PRF_Align_End_Target = PRF_Alignment@pattern@range@width
      PRF_Align_End_Target_Expected = nchar(PRF_curr_Target_seq)
      PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == (PRF_Alignment@pattern@range@width + PRF_Align_Start_Gene - 1), "PRF_Gene_Coords_Genomic"]
      if (length(PRF_Align_End_Target_genomic) < 1){
        PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == max(PRF_Coord_matrix$PRF_Splice_Coords_Vector), "PRF_Gene_Coords_Genomic"]
      }
      
      # Getting Alignment Scores
      PRF_Alignment_score = PRF_Alignment@score
    }
    Output_Data = data.frame(
      ID = U133_Plus_2_Probes_Annotation$ID[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = PRF_Current_UCSC_Known_Gene_df$chrom,
      txStart = PRF_Current_UCSC_Known_Gene_df$txStart,
      txEnd = PRF_Current_UCSC_Known_Gene_df$txEnd,
      PRF_URL = PRF_URL,
      PRF_DNA_init = PRF_curr_DNA_init,
      PRF_DNA_bio_str = PRF_curr_DNA_bio_str,
      PRF_Strand = PRF_Strand,
      PRF_Target_seq = PRF_curr_Target_seq,
      PRF_DNA_transcript = PRF_curr_DNA_transcript,
      PRF_DNA_R_C_str = PRF_curr_DNA_R_C_str,
      PRF_DNA_transcript_bio_R_C_str = PRF_curr_DNA_transcript_bio_R_C_str,
      PRF_Align_Start_Gene = PRF_Align_Start_Gene,
      PRF_Align_Start_Gene_genomic = PRF_Align_Start_Gene_genomic,
      PRF_Align_Start_Target = PRF_Align_Start_Target,
      PRF_Align_End_Target = PRF_Align_End_Target,
      PRF_Align_End_Target_Expected = PRF_Align_End_Target_Expected,
      PRF_Align_End_Target_genomic = PRF_Align_End_Target_genomic,
      PRF_Alignment_score = PRF_Alignment_score
    )
    write.csv(Output_Data, file = paste0("Mapping_Affymetrix_Probes_Alignments/",i, "__",U133_Plus_2_Probes_Annotation$ID[i], ".csv"))
    i = i + 1
  }
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

#Reading through files and creating the dataset
PRF_files_affymetrix_probes = list.files("Mapping_Affymetrix_Probes_Alignments")
PRF_files_affymetrix_probes = paste0("Mapping_Affymetrix_Probes_Alignments/", PRF_files_affymetrix_probes)
U133_Plus_2_Probes_Alignments = lapply(PRF_files_affymetrix_probes, function(x){
  PRF_data = smart_fread(x)
  PRF_data = PRF_data[,colnames(PRF_data) %!in% c("PRF_DNA_init", "PRF_DNA_bio_str", "PRF_DNA_R_C_str")]
  return(PRF_data)
})
U133_Plus_2_Probes_Alignments = list_to_df(U133_Plus_2_Probes_Alignments)
U133_Plus_2_Probes_Alignments = U133_Plus_2_Probes_Alignments[match(U133_Plus_2_Probes_Annotation$ID, U133_Plus_2_Probes_Alignments$ID),]
table(U133_Plus_2_Probes_Alignments$PRF_Strand, useNA = "always")
U133_Plus_2_Probes_Alignments_NA = U133_Plus_2_Probes_Alignments[is.na(U133_Plus_2_Probes_Alignments$PRF_Strand),]
U133_Plus_2_Probes_Annotation_Unalingned_1 = U133_Plus_2_Probes_Annotation[U133_Plus_2_Probes_Annotation$ID %in% U133_Plus_2_Probes_Alignments_NA$ID,]
# A part of these probes could be excluded since cannot be adequately mapped (Part of the probes is mitochondrial and pseudogenes)

# Match check
all(U133_Plus_2_Probes_Alignments$ID == U133_Plus_2_Probes_Annotation$ID) # All are matching

# Additional stats and modifications of the dataset
U133_Plus_2_Probes_Alignments$Align_Coord_Start = NA
U133_Plus_2_Probes_Alignments$Align_Coord_End = NA
U133_Plus_2_Probes_Alignments$Sequence_matching = NA

# Full match is 5 based on the NUC4.4 scoring matrix
for (i in 1:nrow(U133_Plus_2_Probes_Alignments)){
  if (is.na(U133_Plus_2_Probes_Alignments$PRF_Strand[i])){
    U133_Plus_2_Probes_Alignments$Align_Coord_Start[i] = NA
    U133_Plus_2_Probes_Alignments$Align_Coord_End[i] = NA
  } else if (U133_Plus_2_Probes_Alignments$PRF_Strand[i] == "-"){
    U133_Plus_2_Probes_Alignments$Align_Coord_Start[i] = U133_Plus_2_Probes_Alignments$PRF_Align_End_Target_genomic[i]
    U133_Plus_2_Probes_Alignments$Align_Coord_End[i] = U133_Plus_2_Probes_Alignments$PRF_Align_Start_Gene_genomic[i]
  } else {
    U133_Plus_2_Probes_Alignments$Align_Coord_Start[i] = U133_Plus_2_Probes_Alignments$PRF_Align_Start_Gene_genomic[i]
    U133_Plus_2_Probes_Alignments$Align_Coord_End[i] = U133_Plus_2_Probes_Alignments$PRF_Align_End_Target_genomic[i]
  }
  if (!is.na(U133_Plus_2_Probes_Alignments$PRF_Alignment_score[i])){
    U133_Plus_2_Probes_Alignments$Sequence_matching[i] = U133_Plus_2_Probes_Alignments$PRF_Alignment_score[i]/(U133_Plus_2_Probes_Alignments$PRF_Align_End_Target[i]*5)
  }
}

hist(U133_Plus_2_Probes_Alignments$Sequence_matching,
     main = "Sequence matching frequency",
     xlab = "Sequence matching")

png("Sequence_matching Affymetrix U133 Plus 2.png", width = 1000, height = 1000)
hist(U133_Plus_2_Probes_Alignments$Sequence_matching,
     main = "Sequence matching frequency",
     xlab = "Sequence matching") # Keep only scores above 0.8
dev.off()

# Creating and writing outputs
U133_Plus_2_Probes_Alignments_PASSED = U133_Plus_2_Probes_Alignments[!is.na(U133_Plus_2_Probes_Alignments$PRF_Strand), ]
U133_Plus_2_Probes_Alignments_PASSED = U133_Plus_2_Probes_Alignments_PASSED[U133_Plus_2_Probes_Alignments_PASSED$Sequence_matching >= 0.8,]
U133_Plus_2_Probes_Alignments_EXCLUDED = U133_Plus_2_Probes_Alignments[U133_Plus_2_Probes_Alignments$ID %!in% U133_Plus_2_Probes_Alignments_PASSED$ID,]
U133_Plus_2_Probes_Annotation_PASSED = U133_Plus_2_Probes_Annotation[U133_Plus_2_Probes_Annotation$ID %in% U133_Plus_2_Probes_Alignments_PASSED$ID,]
U133_Plus_2_Probes_Annotation_EXCLUDED = U133_Plus_2_Probes_Annotation[U133_Plus_2_Probes_Annotation$ID %in% U133_Plus_2_Probes_Alignments_EXCLUDED$ID,]
write.csv(x = U133_Plus_2_Probes_Alignments_PASSED, file = "U133_Plus_2_Probes_Alignments_PASSED.csv")
write.csv(x = U133_Plus_2_Probes_Alignments_EXCLUDED, file = "U133_Plus_2_Probes_Alignments_EXCLUDED.csv")
write.csv(x = U133_Plus_2_Probes_Annotation_PASSED, file = "U133_Plus_2_Probes_Annotation_PASSED.csv")
write.csv(x = U133_Plus_2_Probes_Annotation_EXCLUDED, file = "U133_Plus_2_Probes_Annotation_EXCLUDED.csv")


################### Filtering and mapping Illumina HumanHT-12 v3 & Illumina HumanHT-12 v4 probes ###################

# Creating_Combined annotation file
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined = rbind(GSE64930_annotation, GSE46743_annotation)
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined = distinct(Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined,  Probe_Id, .keep_all = TRUE) #14613 unique probes
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined[Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined$Source != "ILMN_Controls",]
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr[Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Symbol != "",]
table(Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Source) # RefSeq; 13879 probes in total
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_CONTROLS = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined[Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined$Probe_Id 
                                                                                                           %!in% Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id,]
write.csv(x = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr, file = "Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr.csv")
write.csv(x = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_CONTROLS, file = "Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_CONTROLS.csv")

# Updating gene names
Illumina_HT_12_V3_V4_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Symbol, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                         PRF_replace_NA_with_old = TRUE)
probes_genes_check_dict = Illumina_HT_12_V3_V4_genes_check$Suggested.Symbol
names(probes_genes_check_dict) = Illumina_HT_12_V3_V4_genes_check$x
Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Updated_gene_names = sapply(Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Symbol, function(x){
  if (is.na(x)){return(NA)}
  if (x == ""){return(NA)}
  Current_genes = unlist(stri_split_fixed(x, pattern = " /// "))
  Current_genes = unique(Current_genes)
  Current_genes = str_trim(Current_genes)
  New_names = probes_genes_check_dict[Current_genes]
  New_names = paste0(New_names, collapse = ";")
  return(New_names)
})

# Mapping probes
dir.create("Mapping_Illumina_HT_12_V3_V4") # A folder to store alignments

# Note: The loop is using UCSC API and downloading sequences -> constant internet connection is required
# Note: The loop takes several hours to execute
# Note: The loop variables are labelled with the prefix "PRF_" to enable easy cleanup
i = 1
while (i <= nrow(Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr)){
  print(i)
  PRF_Current_gene = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Updated_gene_names[i]
  PRF_Current_UCSC_Known_Gene_df = UCSC_Known_Gene_full[UCSC_Known_Gene_full$geneSymbol == PRF_Current_gene,]
  PRF_Ref_DS = "knownGene"
  
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    PRF_Current_UCSC_Known_Gene_df = UCSC_NCBI_RefSeq[UCSC_NCBI_RefSeq$name2 == PRF_Current_gene,]
    PRF_Ref_DS = "ncbiRefSeq"
  }
  
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    PRF_Current_UCSC_Known_Gene_df = UCSC_wgEncodeGencodeCompV38lift37[UCSC_wgEncodeGencodeCompV38lift37$name2 == PRF_Current_gene,]
    PRF_Ref_DS = "wgEncodeGencodeCompV38lift37"
  }
  
  if (nrow(PRF_Current_UCSC_Known_Gene_df) < 1){
    
    # A gene was not mapped in any of the datasets
    PRF_Ref_DS = "NOT IDENTIFIED"
    Output_Data = data.frame(
      ID = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = NA,
      txStart = NA,
      txEnd = NA,
      PRF_URL = NA,
      PRF_DNA_init = NA,
      PRF_DNA_bio_str = NA,
      PRF_Strand = NA,
      PRF_Probe_seq = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Sequence[i],
      PRF_Target_seq = NA,
      PRF_DNA_transcript = NA,
      PRF_DNA_R_C_str = NA,
      PRF_DNA_transcript_bio_R_C_str = NA,
      PRF_Align_Start_Gene = NA,
      PRF_Align_Start_Gene_genomic = NA,
      PRF_Align_Start_Target = NA,
      PRF_Align_End_Target = NA,
      PRF_Align_End_Target_Expected = NA,
      PRF_Align_End_Target_genomic = NA,
      PRF_Alignment_score = NA
    )
    write.csv(Output_Data, file = paste0("Mapping_Illumina_HT_12_V3_V4/",i, "__",Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i], ".csv"))
    i = i + 1
    next
  }
  
  # Selecting the best gene representation
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[PRF_Current_UCSC_Known_Gene_df$exonCount == max(PRF_Current_UCSC_Known_Gene_df$exonCount),]
  PRF_Current_UCSC_Known_Gene_df_lengths = mapply(function(x,y){
    Len = x:y
    Len = length(Len)
    return(Len)
  }, PRF_Current_UCSC_Known_Gene_df$cdsStart, PRF_Current_UCSC_Known_Gene_df$cdsEnd)
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[which.max(PRF_Current_UCSC_Known_Gene_df_lengths),]
  PRF_Current_UCSC_Known_Gene_df = PRF_Current_UCSC_Known_Gene_df[1,]
  
  # Getting the sequence
  PRF_URL = paste0("https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=",
                   PRF_Current_UCSC_Known_Gene_df$chrom,
                   ";start=",
                   PRF_Current_UCSC_Known_Gene_df$txStart,
                   ";end=",
                   PRF_Current_UCSC_Known_Gene_df$txEnd)
  PRF_EROR_FLAG = FALSE
  PRF_EROR_Counter = 0
  
  tryCatch({
    PRF_curr_DNA = jsonlite::fromJSON(PRF_URL)
    PRF_curr_DNA = PRF_curr_DNA$dna
  }
  , error = function(e) {PRF_EROR_FLAG <<-TRUE})
  
  while(PRF_EROR_FLAG){
    PRF_EROR_Counter = PRF_EROR_Counter + 1
    warning(paste0('Error detected  during downloading file from   ', PRF_URL, '  ', Sys.time(), ' Trying to reopen link  #', PRF_EROR_Counter), immediate. = TRUE)
    if (PRF_EROR_Counter > 20){
      PRF_curr_DNA = NA
      break
    }
    if (PRF_EROR_Counter > 15) {
      Sys.sleep(sample(1200:3600, 1))
      closeAllConnections()
    }
    Sys.sleep(5*PRF_EROR_Counter)
    tryCatch({
      PRF_curr_DNA = jsonlite::fromJSON(PRF_URL)
      PRF_curr_DNA = PRF_curr_DNA$dna
      PRF_EROR_FLAG = FALSE}
      , error = function(e) {PRF_EROR_FLAG <<-TRUE})
  }
  
  if (is.na(PRF_curr_DNA)|PRF_curr_DNA == ""|PRF_curr_DNA == " "){
    Output_Data = data.frame(
      ID = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = PRF_Current_UCSC_Known_Gene_df$chrom,
      txStart = PRF_Current_UCSC_Known_Gene_df$txStart,
      txEnd = PRF_Current_UCSC_Known_Gene_df$txEnd,
      PRF_URL = PRF_URL,
      PRF_DNA_init = NA,
      PRF_DNA_bio_str = NA,
      PRF_Strand = NA,
      PRF_Probe_seq = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Sequence[i],
      PRF_Target_seq = NA,
      PRF_DNA_transcript = NA,
      PRF_DNA_R_C_str = NA,
      PRF_DNA_transcript_bio_R_C_str = NA,
      PRF_Align_Start_Gene = NA,
      PRF_Align_Start_Gene_genomic = NA,
      PRF_Align_Start_Target = NA,
      PRF_Align_End_Target = NA,
      PRF_Align_End_Target_Expected = NA,
      PRF_Align_End_Target_genomic = NA,
      PRF_Alignment_score = NA
    )
    write.csv(Output_Data, file = paste0("Mapping_Illumina_HT_12_V3_V4/",i, "__",Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i], ".csv"))
    i = i + 1
  } else {
    PRF_curr_DNA_init = PRF_curr_DNA
    PRF_Strand = PRF_Current_UCSC_Known_Gene_df$strand
    PRF_Probe_seq = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Sequence[i]
    PRF_Probe_seq_bio = DNAString(PRF_Probe_seq)
    PRF_curr_Target_seq_bio = reverseComplement(DNAString(PRF_Probe_seq))
    PRF_curr_Target_seq = as.character(PRF_curr_Target_seq_bio)
    
    # Getting Transcripts and Spliced transcripts
    PRF_Exon_Starts = unlist(stri_split_fixed(PRF_Current_UCSC_Known_Gene_df$exonStarts, pattern = ","))
    PRF_Exon_Starts = PRF_Exon_Starts[PRF_Exon_Starts != ""]
    PRF_Exon_Ends = unlist(stri_split_fixed(PRF_Current_UCSC_Known_Gene_df$exonEnds, pattern = ","))
    PRF_Exon_Ends = PRF_Exon_Ends[PRF_Exon_Ends != ""]
    
    # Getting Spliced Sequence
    PRF_Exon_Coords = mapply(function(x,y){
      Vector = x:y
      return(Vector)
    }, PRF_Exon_Starts, PRF_Exon_Ends)
    PRF_Exon_Coords = unlist(PRF_Exon_Coords) # it is sorted
    names(PRF_Exon_Coords) = NULL
    
    # Getting coordinates for sequence
    PRF_curr_DNA_init_vector = unlist(strsplit(PRF_curr_DNA_init, split = "", fixed = TRUE))
    PRF_Gene_Coords = 1:length(PRF_curr_DNA_init_vector)
    PRF_Gene_Coords_inverse = length(PRF_curr_DNA_init_vector):1
    PRF_Gene_Coords_Genomic = PRF_Gene_Coords + (PRF_Current_UCSC_Known_Gene_df$txStart)
    
    # Getting spliced coordinates
    PRF_Splice_Coords_Ticker = 0
    PRF_Splice_Coords_Vector = vector()
    for (i_spl in 1:length(PRF_Gene_Coords_Genomic)){
      PRF_Curr_coord = PRF_Gene_Coords_Genomic[i_spl]
      if (PRF_Curr_coord %in% PRF_Exon_Coords){
        PRF_Spl_coord = PRF_Splice_Coords_Ticker[length(PRF_Splice_Coords_Ticker)] + 1
        PRF_Splice_Coords_Vector[i_spl] = PRF_Spl_coord
        PRF_Splice_Coords_Ticker = c(PRF_Splice_Coords_Ticker, PRF_Spl_coord)
      } else {
        PRF_Splice_Coords_Vector[i_spl] = 0
      }
    }
    PRF_Coord_matrix = cbind(PRF_Gene_Coords, PRF_Gene_Coords_inverse, PRF_Gene_Coords_Genomic, PRF_Splice_Coords_Vector)
    PRF_Coord_matrix = as.data.frame(PRF_Coord_matrix)
    PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse = NA
    PRF_TR_len = length(which(PRF_Coord_matrix$PRF_Splice_Coords_Vector != 0))
    PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse = sapply(PRF_Coord_matrix$PRF_Splice_Coords_Vector, function(x){
      if(x == 0){
        return(0)
      }
      x = abs(x-PRF_TR_len) + 1
      return(x)
    })
    PRF_Coord_matrix$DNA_init = PRF_curr_DNA_init_vector
    PRF_Coord_matrix$DNA_UP = toupper(PRF_curr_DNA_init_vector)
    PRF_curr_DNA_transcript = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector != 0, "DNA_UP"]
    PRF_curr_DNA_transcript = paste0(PRF_curr_DNA_transcript, collapse = "")
    
    if (PRF_Strand == "-"){
      PRF_curr_DNA_bio_str = toupper(PRF_curr_DNA_init)
      PRF_curr_DNA_bio = DNAString(PRF_curr_DNA_bio_str)
      PRF_curr_DNA_R_C = reverseComplement(PRF_curr_DNA_bio)
      PRF_curr_DNA_transcript_bio = DNAString(PRF_curr_DNA_transcript)
      PRF_curr_DNA_transcript_bio_R_C = reverseComplement(PRF_curr_DNA_transcript_bio)
      PRF_Alignment = pairwiseAlignment(PRF_Probe_seq_bio, PRF_curr_DNA_transcript_bio_R_C, substitutionMatrix = Scoring_EDNAFULL_NUC4.4, type = "overlap", gapOpening = 10, gapExtension = 0.001)
      # Illumina reports probe seq as a target in reality
      
      # Returning to characters
      PRF_curr_DNA_R_C_str = as.character(PRF_curr_DNA_R_C)
      PRF_curr_DNA_transcript_bio_R_C_str = as.character(PRF_curr_DNA_transcript_bio_R_C)
      
      # Calculating the appropriate position
      PRF_Align_Start_Gene = PRF_Alignment@subject@range@start
      PRF_Align_Start_Gene_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == PRF_Align_Start_Gene,"PRF_Gene_Coords_Genomic"]
      PRF_Align_Start_Target = PRF_Alignment@pattern@range@start
      PRF_Align_End_Target = PRF_Alignment@pattern@range@width
      PRF_Align_End_Target_Expected = nchar(PRF_curr_Target_seq)
      PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == (PRF_Alignment@pattern@range@width + PRF_Align_Start_Gene - 1), "PRF_Gene_Coords_Genomic"]
      
      if (length(PRF_Align_End_Target_genomic) < 1){
        PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse == max(PRF_Coord_matrix$PRF_Splice_Coords_Vector_Inverse), "PRF_Gene_Coords_Genomic"]
      }
      
      # Getting Alignment Scores
      PRF_Alignment_score = PRF_Alignment@score
    } else {
      PRF_curr_DNA_bio_str = toupper(PRF_curr_DNA_init)
      PRF_curr_DNA_bio = DNAString(PRF_curr_DNA_bio_str)
      PRF_curr_DNA_R_C = NA
      PRF_curr_DNA_transcript_bio = DNAString(PRF_curr_DNA_transcript)
      PRF_curr_DNA_transcript_bio_R_C = NA
      PRF_Alignment = pairwiseAlignment(PRF_Probe_seq_bio, PRF_curr_DNA_transcript_bio, substitutionMatrix = Scoring_EDNAFULL_NUC4.4, type = "overlap", gapOpening = 10, gapExtension = 0.001)
      # Probe seq is in fact the target sequence, the logic is maintained
      
      # Returning to characters
      PRF_curr_DNA_R_C_str = NA
      PRF_curr_DNA_transcript_bio_R_C_str = NA
      
      # Calculating the appropriate position
      PRF_Align_Start_Gene = PRF_Alignment@subject@range@start
      PRF_Align_Start_Gene_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == PRF_Align_Start_Gene,"PRF_Gene_Coords_Genomic"]
      PRF_Align_Start_Target = PRF_Alignment@pattern@range@start
      PRF_Align_End_Target = PRF_Alignment@pattern@range@width
      PRF_Align_End_Target_Expected = nchar(PRF_curr_Target_seq)
      PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == (PRF_Alignment@pattern@range@width + PRF_Align_Start_Gene - 1), "PRF_Gene_Coords_Genomic"]
      
      if (length(PRF_Align_End_Target_genomic) < 1){
        PRF_Align_End_Target_genomic = PRF_Coord_matrix[PRF_Coord_matrix$PRF_Splice_Coords_Vector == max(PRF_Coord_matrix$PRF_Splice_Coords_Vector), "PRF_Gene_Coords_Genomic"]
      }
      
      # Getting Alignment Scores
      PRF_Alignment_score = PRF_Alignment@score
    }
    Output_Data = data.frame(
      ID = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i],
      PRF_Gene_symbol = PRF_Current_gene,
      PRF_Ref_DS = PRF_Ref_DS,
      chrom = PRF_Current_UCSC_Known_Gene_df$chrom,
      txStart = PRF_Current_UCSC_Known_Gene_df$txStart,
      txEnd = PRF_Current_UCSC_Known_Gene_df$txEnd,
      PRF_URL = PRF_URL,
      PRF_DNA_init = PRF_curr_DNA_init,
      PRF_DNA_bio_str = PRF_curr_DNA_bio_str,
      PRF_Strand = PRF_Strand,
      PRF_Probe_seq = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Sequence[i],
      PRF_Target_seq = PRF_curr_Target_seq,
      PRF_DNA_transcript = PRF_curr_DNA_transcript,
      PRF_DNA_R_C_str = PRF_curr_DNA_R_C_str,
      PRF_DNA_transcript_bio_R_C_str = PRF_curr_DNA_transcript_bio_R_C_str,
      PRF_Align_Start_Gene = PRF_Align_Start_Gene,
      PRF_Align_Start_Gene_genomic = PRF_Align_Start_Gene_genomic,
      PRF_Align_Start_Target = PRF_Align_Start_Target,
      PRF_Align_End_Target = PRF_Align_End_Target,
      PRF_Align_End_Target_Expected = PRF_Align_End_Target_Expected,
      PRF_Align_End_Target_genomic = PRF_Align_End_Target_genomic,
      PRF_Alignment_score = PRF_Alignment_score
    )
    write.csv(Output_Data, file = paste0("Mapping_Illumina_HT_12_V3_V4/",i, "__",Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id[i], ".csv"))
    i = i + 1
  }
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Reading files
PRF_files_Illumina_HT_12_V3_V4_probes = list.files("Mapping_Illumina_HT_12_V3_V4")
PRF_files_Illumina_HT_12_V3_V4_probes = paste0("Mapping_Illumina_HT_12_V3_V4/", PRF_files_Illumina_HT_12_V3_V4_probes)
Illumina_HT_12_V3_V4_probes_Alignments = lapply(PRF_files_Illumina_HT_12_V3_V4_probes, function(x){
  PRF_data = smart_fread(x)
  PRF_data = PRF_data[,colnames(PRF_data) %!in% c("PRF_DNA_init", "PRF_DNA_bio_str", "PRF_DNA_R_C_str")]
  return(PRF_data)
})
Illumina_HT_12_V3_V4_probes_Alignments = list_to_df(Illumina_HT_12_V3_V4_probes_Alignments)

# Additional stats and modifications of the dataset
Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_Start = NA
Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_End = NA
Illumina_HT_12_V3_V4_probes_Alignments$Sequence_matching = NA

# Full match is 5 based on the NUC4.4 scoring matrix
for (i in 1:nrow(Illumina_HT_12_V3_V4_probes_Alignments)){
  if (is.na(Illumina_HT_12_V3_V4_probes_Alignments$PRF_Strand[i])){
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_Start[i] = NA
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_End[i] = NA
  } else if (Illumina_HT_12_V3_V4_probes_Alignments$PRF_Strand[i] == "-"){
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_Start[i] = Illumina_HT_12_V3_V4_probes_Alignments$PRF_Align_End_Target_genomic[i]
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_End[i] = Illumina_HT_12_V3_V4_probes_Alignments$PRF_Align_Start_Gene_genomic[i]
  } else {
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_Start[i] = Illumina_HT_12_V3_V4_probes_Alignments$PRF_Align_Start_Gene_genomic[i]
    Illumina_HT_12_V3_V4_probes_Alignments$Align_Coord_End[i] = Illumina_HT_12_V3_V4_probes_Alignments$PRF_Align_End_Target_genomic[i]
  }
  if (!is.na(Illumina_HT_12_V3_V4_probes_Alignments$PRF_Alignment_score[i])){
    Illumina_HT_12_V3_V4_probes_Alignments$Sequence_matching[i] = Illumina_HT_12_V3_V4_probes_Alignments$PRF_Alignment_score[i]/(Illumina_HT_12_V3_V4_probes_Alignments$PRF_Align_End_Target[i]*5)
  }
}

Illumina_HT_12_V3_V4_probes_Alignments = Illumina_HT_12_V3_V4_probes_Alignments[match(Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id, Illumina_HT_12_V3_V4_probes_Alignments$ID),]

hist(Illumina_HT_12_V3_V4_probes_Alignments$Sequence_matching,
     main = "Sequence matching frequency",
     xlab = "Sequence matching")

png("Sequence_matching Illumina HT-12 V3 V4.png", width = 1000, height = 1000)
hist(Illumina_HT_12_V3_V4_probes_Alignments$Sequence_matching,
     main = "Sequence matching frequency",
     xlab = "Sequence matching") # Keep only scores above 0.8
dev.off()

# Creating and writing outputs
Illumina_HT_12_V3_V4_probes_Alignments_PASSED = Illumina_HT_12_V3_V4_probes_Alignments[!is.na(Illumina_HT_12_V3_V4_probes_Alignments$PRF_Strand), ]
Illumina_HT_12_V3_V4_probes_Alignments_PASSED = Illumina_HT_12_V3_V4_probes_Alignments_PASSED[Illumina_HT_12_V3_V4_probes_Alignments_PASSED$Sequence_matching >= 0.8,] #11705 probes passed out of 13879
Illumina_HT_12_V3_V4_Annotation_PASSED = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr[Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id %in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
Illumina_HT_12_V3_V4_probes_Alignments_EXCLUDED = Illumina_HT_12_V3_V4_probes_Alignments[Illumina_HT_12_V3_V4_probes_Alignments$ID %!in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
Illumina_HT_12_V3_V4_Annotation_EXCLUDED = Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr[Illumina_HT_12_GSE64930_GSE46743_V3_V4_combined_no_ctr$Probe_Id %in% Illumina_HT_12_V3_V4_probes_Alignments_EXCLUDED$ID,]
write.csv(x = Illumina_HT_12_V3_V4_probes_Alignments_PASSED, file = "Illumina_HT_12_V3_V4_probes_Alignments_PASSED.csv")
write.csv(x = Illumina_HT_12_V3_V4_probes_Alignments_EXCLUDED, file = "Illumina_HT_12_V3_V4_probes_Alignments_EXCLUDED.csv")
write.csv(x = Illumina_HT_12_V3_V4_Annotation_PASSED, file = "Illumina_HT_12_V3_V4_Annotation_PASSED.csv")
write.csv(x = Illumina_HT_12_V3_V4_Annotation_EXCLUDED, file = "Illumina_HT_12_V3_V4_Annotation_EXCLUDED.csv")


################### Transcriptome profiling analysis ###################

# Preparing data after filtering probes

# Curating affymetrix GSE98793 
Eset_GSE98793_expresion_full = Eset_GSE98793_expresion_full[rownames(Eset_GSE98793_expresion_full) %in% U133_Plus_2_Probes_Alignments_PASSED$ID,]
GSE98793_probes = GSE98793_probes[GSE98793_probes$ID %in% U133_Plus_2_Probes_Alignments_PASSED$ID, ]
all(GSE98793_probes$ID == rownames(Eset_GSE98793_expresion_full)) # Not a full match
all(GSE98793_probes$ID %in% rownames(Eset_GSE98793_expresion_full))
all(rownames(Eset_GSE98793_expresion_full)  %in% GSE98793_probes$ID)
GSE98793_probes = GSE98793_probes[match(rownames(Eset_GSE98793_expresion_full), GSE98793_probes$ID),]
all(GSE98793_probes$ID == rownames(Eset_GSE98793_expresion_full)) # Full match

# Curating GSE46743
GSE46743_data_processed_baseline = GSE46743_data_processed_baseline[rownames(GSE46743_data_processed_baseline) %in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
GSE46743_annotation = GSE46743_annotation[GSE46743_annotation$Probe_Id %in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
all(GSE46743_annotation$Probe_Id == rownames(GSE46743_data_processed_baseline))# Full match

# Curating GSE64930
GSE64930_data_processed_baseline = GSE64930_data_processed_baseline[rownames(GSE64930_data_processed_baseline) %in%  Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
GSE64930_annotation = GSE64930_annotation[GSE64930_annotation$Probe_Id %in% Illumina_HT_12_V3_V4_probes_Alignments_PASSED$ID,]
all(GSE64930_annotation$Probe_Id == rownames(GSE64930_data_processed_baseline)) #Full match

# Differential expression analysis in GSE98793
Analysis_limma_affymetrix_depr_broad_GSE98793 = test_genes_affymetrix_limma_generelized(PREFIX_gene_names = NULL,
                                                                                        PREFIX_remove_non_specific = TRUE,
                                                                                        PREFIX_Pheno_df  = GSE98793_phenotypes,
                                                                                        PREFIX_Expression_matrix  = Eset_GSE98793_expresion_full,
                                                                                        PREFIX_Probes_df = GSE98793_probes,
                                                                                        PREFIX_contrast_col_number  =  41,
                                                                                        PREFIX_participants_col_number  = 43,
                                                                                        PREFIX_contrast_vector  =  c("MDD", "Control"),
                                                                                        PREFIX_model.covariates  = c("anxiety", "gender", "age", "batch"),
                                                                                        PREFIX_plots_folder = "Broad_depr_GSE98793_plots",
                                                                                        PREFIX_plots_str = "Depression_genes_",
                                                                                        PREFIX_logFC_threshold = 0.1)
Analysis_limma_affymetrix_depr_broad_GSE98793$Contrast = paste0("GSE98793: ",Analysis_limma_affymetrix_depr_broad_GSE98793$Contrast)
Analysis_limma_affymetrix_depr_broad_GSE98793_signif = Analysis_limma_affymetrix_depr_broad_GSE98793[Analysis_limma_affymetrix_depr_broad_GSE98793$P.Value < 0.05,]
Analysis_limma_affymetrix_depr_broad_GSE98793_signif = Analysis_limma_affymetrix_depr_broad_GSE98793_signif[abs(Analysis_limma_affymetrix_depr_broad_GSE98793_signif$logFC) >= 0.1,]
write.csv(Analysis_limma_affymetrix_depr_broad_GSE98793, "Analysis_limma_affymetrix_depr_broad_GSE98793.csv")
write.csv(Analysis_limma_affymetrix_depr_broad_GSE98793_signif, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif.csv") # 3708 probes

# Significant probes at FDR
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR = Analysis_limma_affymetrix_depr_broad_GSE98793[Analysis_limma_affymetrix_depr_broad_GSE98793$adj.P.Val < 0.05,]
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR = Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR[abs(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR$logFC) >= 0.1,] #70 probes
write.csv(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR.csv")

# Getting matching direction probes for GSE98793 (Probes are considered matching if more than 80% had a concordant direction)
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs = list()

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:length(unique(Analysis_limma_affymetrix_depr_broad_GSE98793_signif$Updated_gene_names))){
  PRF_current_gene = unique(Analysis_limma_affymetrix_depr_broad_GSE98793_signif$Updated_gene_names)[i]
  PRF_curr_df = Analysis_limma_affymetrix_depr_broad_GSE98793_signif[Analysis_limma_affymetrix_depr_broad_GSE98793_signif$Updated_gene_names == PRF_current_gene,]
  PRF_curr_LFs = PRF_curr_df$logFC
  if (length(PRF_curr_LFs) == 1){
    Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[[i]] = PRF_curr_df
  } else {
    PRF_Idx_more_than_0 = length(which(PRF_curr_LFs > 0))/length(PRF_curr_LFs)
    PRF_Idx_more_less_than_0 = length(which(PRF_curr_LFs < 0))/length(PRF_curr_LFs)
    PRF_idx_Total = c(PRF_Idx_more_than_0, PRF_Idx_more_less_than_0)
    if (any(PRF_idx_Total >= 0.8)){
      
      # Determining incorrect direction
      PRF_false = which(PRF_idx_Total < 0.8)
      
      if (PRF_false == 1){
        # False gene LF is > 0
        # All valid probes should have a negative logFC
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC < 0,]
      }
      
      if (PRF_false == 2){
        # False gene is < 0
        # All valid probes should have a positive logFC
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC > 0,]
      }
      
      Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[[i]] = PRF_curr_df
    } else {
      Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[[i]] = NA
    }
  }
}
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs = Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[sapply(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs, is.data.frame)]
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs = list_to_df(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs)
write.csv(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs.csv") #3702 probes

# Differential expression analysis in GSE46743
Analysis_limma_illumina_depr_broad_GSE46743 = test_genes_Illumina_limma_generelized(PREFIX_gene_names = NULL,
                                                                                    PREFIX_remove_non_specific = TRUE,
                                                                                    PREFIX_Pheno_df  = GSE46743_phenotypes_table_baseline,
                                                                                    PREFIX_Expression_matrix  = GSE46743_data_processed_baseline,
                                                                                    PREFIX_Probes_df = GSE46743_annotation,
                                                                                    PREFIX_contrast_col_number  = 30,
                                                                                    PREFIX_participants_col_number  = 11,
                                                                                    PREFIX_contrast_vector  =  c("case", "control"),
                                                                                    PREFIX_model.covariates  = c("age", "bmi"),
                                                                                    PREFIX_plots_folder = "Broad_depr_GSE46743_plots",
                                                                                    PREFIX_plots_str = "Broad_depr_GSE46743_",
                                                                                    PREFIX_Remove_NA_predictors = TRUE,
                                                                                    PREFIX_logFC_threshold = 0.1)
Analysis_limma_illumina_depr_broad_GSE46743$Contrast = paste0("GSE46743: ",Analysis_limma_illumina_depr_broad_GSE46743$Contrast)
Analysis_limma_illumina_depr_broad_GSE46743_signif = Analysis_limma_illumina_depr_broad_GSE46743[Analysis_limma_illumina_depr_broad_GSE46743$P.Value < 0.05,]
Analysis_limma_illumina_depr_broad_GSE46743_signif = Analysis_limma_illumina_depr_broad_GSE46743_signif[abs(Analysis_limma_illumina_depr_broad_GSE46743_signif$logFC) >= 0.1,]
write.csv(Analysis_limma_illumina_depr_broad_GSE46743, "Analysis_limma_illumina_depr_broad_GSE46743.csv") # 11352 analyses
write.csv(Analysis_limma_illumina_depr_broad_GSE46743_signif, "Analysis_limma_illumina_depr_broad_GSE46743_signif.csv") # 814 nominally significant probes

# Checking the uniquness of genes
Analysis_GSE46743_genes = as.data.frame(table(Analysis_limma_illumina_depr_broad_GSE46743_signif$Gene)) # Some genes have 3 probes

# Getting matching direction probes for GSE46743 (Probes are considered matching if more than 80% had a concordant direction)
Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs = list()

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:length(unique(Analysis_limma_illumina_depr_broad_GSE46743_signif$Gene))){
  PRF_current_gene = unique(Analysis_limma_illumina_depr_broad_GSE46743_signif$Gene)[i]
  PRF_curr_df = Analysis_limma_illumina_depr_broad_GSE46743_signif[Analysis_limma_illumina_depr_broad_GSE46743_signif$Gene == PRF_current_gene,]
  PRF_curr_LFs = PRF_curr_df$logFC
  if (length(PRF_curr_LFs) == 1){
    Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[[i]] = PRF_curr_df
  } else {
    PRF_Idx_more_than_0 = length(which(PRF_curr_LFs > 0))/length(PRF_curr_LFs)
    PRF_Idx_more_less_than_0 = length(which(PRF_curr_LFs < 0))/length(PRF_curr_LFs)
    PRF_idx_Total = c(PRF_Idx_more_than_0, PRF_Idx_more_less_than_0)
    if (any(PRF_idx_Total >= 0.8)){
      # Determining the false direction
      PRF_false = which(PRF_idx_Total < 0.8)
      
      if (PRF_false == 1){
        # False gene is > 0
        # All genes should be less than 0
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC < 0,]
      }
      
      if (PRF_false == 2){
        #False gene is < 0
        # All genes should be more than 0
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC > 0,]
      }
      
      Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[[i]] = PRF_curr_df
    } else {
      Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[[i]] = NA
    }
  }
}
Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs = Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[sapply(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs, is.data.frame)]
Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs = list_to_df(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs)
write.csv(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs, "Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs.csv") # 814 passed, ALL are matching

# Differential expression analysis in GSE64930
Analysis_limma_illumina_depr_broad_GSE64930 = test_genes_Illumina_limma_generelized(PREFIX_gene_names = NULL,
                                                                                    PREFIX_remove_non_specific = TRUE,
                                                                                    PREFIX_Pheno_df  = GSE64930_phenotypes_baseline,
                                                                                    PREFIX_Expression_matrix  = GSE64930_data_processed_baseline,
                                                                                    PREFIX_Probes_df = GSE64930_annotation,
                                                                                    PREFIX_contrast_col_number  = 5,
                                                                                    PREFIX_participants_col_number  = 1,
                                                                                    PREFIX_contrast_vector  =  c("case", "control"),
                                                                                    PREFIX_model.covariates  = c("Sex", "Age", "BMI", "RIN", "SV1", "SV2", "SV3"),
                                                                                    PREFIX_plots_folder = "Broad_depr_GSE64930_plots",
                                                                                    PREFIX_plots_str = "Broad_depr_GSE64930_",
                                                                                    PREFIX_logFC_threshold = 0.1)
Analysis_limma_illumina_depr_broad_GSE64930$Contrast = paste0("GSE64930: ",Analysis_limma_illumina_depr_broad_GSE64930$Contrast)
Analysis_limma_illumina_depr_broad_GSE64930_signif = Analysis_limma_illumina_depr_broad_GSE64930[Analysis_limma_illumina_depr_broad_GSE64930$P.Value < 0.05,]
Analysis_limma_illumina_depr_broad_GSE64930_signif = Analysis_limma_illumina_depr_broad_GSE64930_signif[abs(Analysis_limma_illumina_depr_broad_GSE64930_signif$logFC) >= 0.1,]
write.csv(Analysis_limma_illumina_depr_broad_GSE64930, "Analysis_limma_illumina_depr_broad_GSE64930.csv") # 10846 analyses
write.csv(Analysis_limma_illumina_depr_broad_GSE64930_signif, "Analysis_limma_illumina_depr_broad_GSE64930_signif.csv") # 148 passed

# Getting matching direction probes for GSE46743 (Probes are considered matching if more than 80% had a concordant direction)
Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs = list()

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:length(unique(Analysis_limma_illumina_depr_broad_GSE64930_signif$Gene))){
  PRF_current_gene = unique(Analysis_limma_illumina_depr_broad_GSE64930_signif$Gene)[i]
  PRF_curr_df = Analysis_limma_illumina_depr_broad_GSE64930_signif[Analysis_limma_illumina_depr_broad_GSE64930_signif$Gene == PRF_current_gene,]
  PRF_curr_LFs = PRF_curr_df$logFC
  if (length(PRF_curr_LFs) == 1){
    Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[[i]] = PRF_curr_df
  } else {
    PRF_Idx_more_than_0 = length(which(PRF_curr_LFs > 0))/length(PRF_curr_LFs)
    PRF_Idx_more_less_than_0 = length(which(PRF_curr_LFs < 0))/length(PRF_curr_LFs)
    PRF_idx_Total = c(PRF_Idx_more_than_0, PRF_Idx_more_less_than_0)
    if (any(PRF_idx_Total >= 0.8)){
      # Determining the false direction
      PRF_false = which(PRF_idx_Total < 0.8)
      
      if (PRF_false == 1){
        # False direction is > 0
        # All genes should be less than 0
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC < 0,]
      }
      
      if (PRF_false == 2){
        # False direction is < 0
        # All genes should be more than 0
        PRF_curr_df = PRF_curr_df[PRF_curr_df$logFC > 0,]
      }
      
      Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[[i]] = PRF_curr_df
    } else {
      Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[[i]] = NA
    }
  }
}
Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs = Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[sapply(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs, is.data.frame)]
Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs = list_to_df(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs)
write.csv(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs, "Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs.csv") #148 passed, ALL passed
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])


################### Making Venn diagrams ###################

# Covered genes
Venn_covered_genes = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930$Gene)
)
make_Venn_digram_list(named_list = Venn_covered_genes, palette = 2, plot_full_path = "Venn_covered_genes_expression_depression.pdf")

# Covered genes significant
Venn_covered_genes_signif = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs$Gene)
)
make_Venn_digram_list(named_list = Venn_covered_genes_signif, palette = 3, plot_full_path = "Venn_expression_depression_sign_cohort_matching.pdf")

# Covered genes significant (present in all cohorts)
Venn_covered_genes = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930$Gene)
)
Transcriptome_intersect_all = Reduce(intersect, list(Venn_covered_genes[[1]], Venn_covered_genes[[2]], Venn_covered_genes[[3]]))
Venn_covered_genes_signif_intersect = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs$Gene)
)
Venn_covered_genes_signif_intersect = lapply(Venn_covered_genes_signif_intersect, function(x){
  x = x[x %in% Transcriptome_intersect_all]
  return(x)
})
make_Venn_digram_list(named_list = Venn_covered_genes_signif_intersect, palette = 4, plot_full_path = "Venn_expression_depression_sign_cohort_matching_intersect.pdf")

# only 40.7% of genes are common for all studies
Expression_Genes_intersect_all = Reduce(intersect, list(Venn_covered_genes[[1]], Venn_covered_genes[[2]], Venn_covered_genes[[3]]))

# Making combined expression dataset
Colnames_to_select_GSE98793 = c("Contrast",
                                "Model.formula.string",
                                "Updated_gene_names",
                                "ID",
                                "logFC",
                                "AveExpr",
                                "t",
                                "P.Value",
                                "adj.P.Val",
                                "B",
                                "Inflation.Bacon",
                                "Bias.Bacon",
                                "T.Bacon",
                                "P.val.Bacon",
                                "Adj.P.val.Bacon")

Colnames_to_select_GSE46743_GSE64930 = c("Contrast",
                                         "Model.formula.string",
                                         "Gene",
                                         "ID",
                                         "logFC",
                                         "AveExpr",
                                         "t",
                                         "P.Value",
                                         "adj.P.Val",
                                         "B",
                                         "Inflation.Bacon",
                                         "Bias.Bacon",
                                         "T.Bacon",
                                         "P.val.Bacon",
                                         "Adj.P.val.Bacon")

Combined_expression_df = list(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[,Colnames_to_select_GSE98793],
                              Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930],
                              Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930])
Combined_expression_df = lapply(Combined_expression_df, function(x){
  colnames(x) = Colnames_to_select_GSE46743_GSE64930
  return(x)
})
Combined_expression_df = do.call(rbind, Combined_expression_df)
Combined_expression_df = Combined_expression_df[Combined_expression_df$Gene %in% Expression_Genes_intersect_all,]

# Getting significant probe IDs
PRF_One_gene_study_df = list()
PRF_study_levels = unique(Combined_expression_df$Contrast)
for (index2 in 1:length(PRF_study_levels)){
  PRF_curr_study = Combined_expression_df[Combined_expression_df$Contrast == PRF_study_levels[index2],]
  # Simplify the data to the first probe (they all are matching in direction)
  PRF_curr_study = distinct(PRF_curr_study, Gene, .keep_all = TRUE)
  PRF_One_gene_study_df[[index2]] = PRF_curr_study
}
PRF_One_gene_study_df = list_to_df(PRF_One_gene_study_df)

# Getting genes with matching direction
PRF_matching_dir_index = sapply(unique(PRF_One_gene_study_df$Gene), function(PRF_x){
  PRF_LFs = PRF_One_gene_study_df[PRF_One_gene_study_df$Gene == PRF_x, "logFC"]
  PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
  if (length(PRF_LFs) < 2){
    return(FALSE)
  }
  PRF_probes_greater_0 = which(PRF_LFs > 0)
  PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
  PRF_probes_smaller_0 = which(PRF_LFs < 0)
  PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
  PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
  if (any(PRF_combined == 1)){
    return(TRUE)
  }
  return(FALSE)
})
Genes_matching_dir = unique(PRF_One_gene_study_df$Gene)[PRF_matching_dir_index]
Genes_matching_dir = Genes_matching_dir[Genes_matching_dir %in% Transcriptome_intersect_all] #230

# Venn diagram
Venn_covered_genes_signif_intersect = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs$Gene)
)
Venn_covered_genes_signif_intersect_full_match = lapply(Venn_covered_genes_signif_intersect, function(x){
  x = x[x %in% Transcriptome_intersect_all]
  x = x[x %in% Genes_matching_dir]
  return(x)
})
make_Venn_digram_list(named_list = Venn_covered_genes_signif_intersect_full_match, palette = 5, plot_full_path = "Venn_expression_depression_sign_full_matching_intersect.pdf")
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])


################### Making Step-based plots ###################

# Significant probes total vs P-value
# Note: Temporary variables are labelled with the prefix PRF_ to enable easy cleanup after execution
Dataset_Step_P_val = list()
Combined_expression_df = list(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[,Colnames_to_select_GSE98793],
                              Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930],
                              Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930])
Combined_expression_df = lapply(Combined_expression_df, function(x){
  colnames(x) = Colnames_to_select_GSE46743_GSE64930
  return(x)
})
Combined_expression_df = do.call(rbind, Combined_expression_df)
Combined_expression_df = Combined_expression_df[Combined_expression_df$Gene %in% Expression_Genes_intersect_all,]
# 3157 probes are significant and presented in the overlap

# Preparing steps in P
Combined_expression_df$P.Value.log = -log10(Combined_expression_df$P.Value)
PRF_Max_P_val_log = max(Combined_expression_df$P.Value.log)
PRF_Max_P_val_log = round(PRF_Max_P_val_log, digits = 1)
PRF_Min_P_val_log = min(Combined_expression_df$P.Value.log)
PRF_Min_P_val_log = round(PRF_Min_P_val_log, digits = 1)
PRF_P_val_seq = seq(from = PRF_Min_P_val_log, to = PRF_Max_P_val_log, by = 0.05)

# loop
for (i in 1:length(PRF_P_val_seq)){
  PRF_current_P_val_log = PRF_P_val_seq[i]
  PRF_current_df = Combined_expression_df[Combined_expression_df$P.Value.log >= PRF_current_P_val_log,]
  PRF_current_df_GSE98793 = PRF_current_df[PRF_current_df$Contrast == "GSE98793: MDD/Control",]
  PRF_current_df_GSE46743 = PRF_current_df[PRF_current_df$Contrast == "GSE46743: case/control",]
  PRF_current_df_GSE64930 = PRF_current_df[PRF_current_df$Contrast == "GSE64930: case/control",]
  
  # Counting significant Genes
  PRF_sign_GSE98793 = length(unique(PRF_current_df_GSE98793$Gene))
  PRF_sign_GSE46743 = length(unique(PRF_current_df_GSE46743$Gene))
  PRF_sign_GSE64930 = length(unique(PRF_current_df_GSE64930$Gene))
  PRF_sign_sum = sum(PRF_sign_GSE98793, PRF_sign_GSE46743, PRF_sign_GSE64930)
  PRF_sign_unqiue = length(unique(PRF_current_df$Gene))
  
  # Getting overlapping genes
  PRF_sign_GSE98793_genes = unique(PRF_current_df_GSE98793$Gene)
  PRF_sign_GSE46743_genes = unique(PRF_current_df_GSE46743$Gene)
  PRF_sign_GSE64930_genes = unique(PRF_current_df_GSE64930$Gene)
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes))
  
  PRF_intersect_GSE98793_GSE46743 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes)
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE98793_GSE64930 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE46743_GSE64930 = intersect(PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_GSE98793_GSE46743_count = length(PRF_intersect_GSE98793_GSE46743)
  PRF_intersect_GSE98793_GSE64930_count = length(PRF_intersect_GSE98793_GSE64930)
  PRF_intersect_GSE46743_GSE64930_count = length(PRF_intersect_GSE46743_GSE64930)
  
  # Percent total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_all = round((PRF_intersect_GSE98793_GSE46743_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_all = round((PRF_intersect_GSE98793_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_all = round((PRF_intersect_GSE46743_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  
  # Percent intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_GSE98793_GSE46743, PRF_intersect_GSE98793_GSE64930, PRF_intersect_GSE46743_GSE64930)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_intersect = round((PRF_intersect_GSE98793_GSE46743_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_intersect = round((PRF_intersect_GSE98793_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_intersect = round((PRF_intersect_GSE46743_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_output_df = data.frame(
    PRF_sign_GSE98793 = PRF_sign_GSE98793,
    PRF_sign_GSE46743 = PRF_sign_GSE46743,
    PRF_sign_GSE64930 = PRF_sign_GSE64930,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_GSE98793_GSE46743_count = PRF_intersect_GSE98793_GSE46743_count,
    PRF_intersect_GSE98793_GSE64930_count = PRF_intersect_GSE98793_GSE64930_count,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_GSE98793_GSE46743_count_percent_all = PRF_intersect_GSE98793_GSE46743_count_percent_all,
    PRF_intersect_GSE98793_GSE64930_count_percent_all = PRF_intersect_GSE98793_GSE64930_count_percent_all,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_GSE98793_GSE46743_count_percent_intersect = PRF_intersect_GSE98793_GSE46743_count_percent_intersect,
    PRF_intersect_GSE98793_GSE64930_count_percent_intersect = PRF_intersect_GSE98793_GSE64930_count_percent_intersect,
    PRF_intersect_GSE46743_GSE64930_count_percent_intersect = PRF_intersect_GSE46743_GSE64930_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(P.Val.10.log = PRF_current_P_val_log, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_P_val[[i]] = PRF_output_df
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Dataset_Step_P_val = list_to_df(Dataset_Step_P_val)
Dataset_Step_P_val_2 = Dataset_Step_P_val[stri_detect_fixed(Dataset_Step_P_val$Score_name, pattern = "_sign_"),]
Dataset_Step_P_val_2$Score_name = factor(Dataset_Step_P_val_2$Score_name, levels = unique(Dataset_Step_P_val_2$Score_name), 
                                         labels = c("GSE98793",
                                                    "GSE46743",
                                                    "GSE64930",
                                                    "Total",
                                                    "Total unique"))
plot = ggplot(data = Dataset_Step_P_val_2, aes(x = P.Val.10.log, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 50, expand = c(0, 0), limits = c(0,3500)) +
  labs(y = "Signif. proteins count", x = "-log10 p-value", col = "Set", shape = "Set", linetype = "Set") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", 
             color = "red", size=0.5) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17), # Customizing the plot title
    axis.title = element_text(face = "bold", size = 14), # Customizing the axes titles
    legend.title = element_text(face = "bold", size = 14), # Customizing the legend title
    axis.text = element_text(size = 13, colour = "black"), # Customizing the axes text
    legend.text = element_text(size = 13, colour = "black"), # Customizing the legend text
    panel.background = element_blank(), # Removing ugly gray background
    axis.line = element_line(size = 0.5), # Adding axis lines
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
    strip.background = element_blank(),  # Removes facet background
    strip.text = element_text(face = "bold", size = 14, colour = "blue") # Customizing facet text
  )
plot
pdf("Total_Prot_expr_vs_P_val.pdf", width = 10, height = 10)
plot
dev.off()

# Significant probes total vs logFC
# Note: Temporary variables are labelled with the prefix PRF_ to enable easy cleanup after execution
Dataset_Step_Lf = list()

# Preparing logFC steps
PRF_Max_LF = max(abs(Combined_expression_df$logFC))
PRF_Max_LF = round(PRF_Max_LF, digits = 1)
PRF_Max_LF_seq = seq(from = 0.1, to = PRF_Max_LF, by = 0.02)

# loop
for (i in 1:length(PRF_Max_LF_seq)){
  PRF_current_LF = PRF_Max_LF_seq[i]
  PRF_current_df = Combined_expression_df[abs(Combined_expression_df$logFC) >= PRF_current_LF,]
  PRF_current_df_GSE98793 = PRF_current_df[PRF_current_df$Contrast == "GSE98793: MDD/Control",]
  PRF_current_df_GSE46743 = PRF_current_df[PRF_current_df$Contrast == "GSE46743: case/control",]
  PRF_current_df_GSE64930 = PRF_current_df[PRF_current_df$Contrast == "GSE64930: case/control",]
  
  # Counting significant genes
  PRF_sign_GSE98793 = length(unique(PRF_current_df_GSE98793$Gene))
  PRF_sign_GSE46743 = length(unique(PRF_current_df_GSE46743$Gene))
  PRF_sign_GSE64930 = length(unique(PRF_current_df_GSE64930$Gene))
  PRF_sign_sum = sum(PRF_sign_GSE98793, PRF_sign_GSE46743, PRF_sign_GSE64930)
  PRF_sign_unqiue = length(unique(PRF_current_df$Gene))
  
  # Getting overlapping gens
  PRF_sign_GSE98793_genes = unique(PRF_current_df_GSE98793$Gene)
  PRF_sign_GSE46743_genes = unique(PRF_current_df_GSE46743$Gene)
  PRF_sign_GSE64930_genes = unique(PRF_current_df_GSE64930$Gene)
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes))
  
  PRF_intersect_GSE98793_GSE46743 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes)
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE98793_GSE64930 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE46743_GSE64930 = intersect(PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930 %!in% PRF_intersect_all]
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_GSE98793_GSE46743_count = length(PRF_intersect_GSE98793_GSE46743)
  PRF_intersect_GSE98793_GSE64930_count = length(PRF_intersect_GSE98793_GSE64930)
  PRF_intersect_GSE46743_GSE64930_count = length(PRF_intersect_GSE46743_GSE64930)
  
  # Percent_total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_all = round((PRF_intersect_GSE98793_GSE46743_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_all = round((PRF_intersect_GSE98793_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_all = round((PRF_intersect_GSE46743_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  
  #Percent_intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_GSE98793_GSE46743, PRF_intersect_GSE98793_GSE64930, PRF_intersect_GSE46743_GSE64930)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_intersect = round((PRF_intersect_GSE98793_GSE46743_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_intersect = round((PRF_intersect_GSE98793_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_intersect = round((PRF_intersect_GSE46743_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_output_df = data.frame(
    PRF_sign_GSE98793 = PRF_sign_GSE98793,
    PRF_sign_GSE46743 = PRF_sign_GSE46743,
    PRF_sign_GSE64930 = PRF_sign_GSE64930,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_GSE98793_GSE46743_count = PRF_intersect_GSE98793_GSE46743_count,
    PRF_intersect_GSE98793_GSE64930_count = PRF_intersect_GSE98793_GSE64930_count,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_GSE98793_GSE46743_count_percent_all = PRF_intersect_GSE98793_GSE46743_count_percent_all,
    PRF_intersect_GSE98793_GSE64930_count_percent_all = PRF_intersect_GSE98793_GSE64930_count_percent_all,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_GSE98793_GSE46743_count_percent_intersect = PRF_intersect_GSE98793_GSE46743_count_percent_intersect,
    PRF_intersect_GSE98793_GSE64930_count_percent_intersect = PRF_intersect_GSE98793_GSE64930_count_percent_intersect,
    PRF_intersect_GSE46743_GSE64930_count_percent_intersect = PRF_intersect_GSE46743_GSE64930_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(LF = PRF_current_LF, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_Lf[[i]] = PRF_output_df
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Dataset_Step_Lf = list_to_df(Dataset_Step_Lf)
Dataset_Step_Lf_2 = Dataset_Step_Lf[stri_detect_fixed(Dataset_Step_Lf$Score_name, pattern = "_sign_"),]
Dataset_Step_Lf_2$Score_name = factor(Dataset_Step_Lf_2$Score_name, levels = unique(Dataset_Step_Lf_2$Score_name), 
                                      labels = c("GSE98793",
                                                 "GSE46743",
                                                 "GSE64930",
                                                 "Total",
                                                 "Total unique"))
plot = ggplot(data = Dataset_Step_Lf_2, aes(x = LF, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 50, expand = c(0, 0), limits = c(0,3500)) +
  labs(y = "Signif. proteins count", x = "|Log2 fold change|", col = "Set", shape = "Set", linetype = "Set") +
  geom_vline(xintercept = 0.1, linetype="dashed", 
             color = "red", size=0.5) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17), # Customizing the plot title
    axis.title = element_text(face = "bold", size = 14), # Customizing the axes titles
    legend.title = element_text(face = "bold", size = 14), # Customizing the legend title
    axis.text = element_text(size = 13, colour = "black"), # Customizing the axes text
    legend.text = element_text(size = 13, colour = "black"), # Customizing the legend text
    panel.background = element_blank(), # Removing ugly gray background
    axis.line = element_line(size = 0.5), # Adding axis lines
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
    strip.background = element_blank(),  # Removes facet background
    strip.text = element_text(face = "bold", size = 14, colour = "blue") # Customizing facet text
  )
plot
pdf( "Total_Prot_expr_vs_LF.pdf", width = 10, height = 10)
plot
dev.off()

# Significant probes (similar direction) vs P-value
Dataset_Step_P_val = list()

# Making steps for P
Combined_expression_df$P.Value.log = -log10(Combined_expression_df$P.Value)
PRF_Max_P_val_log = max(Combined_expression_df$P.Value.log)
PRF_Max_P_val_log = round(PRF_Max_P_val_log, digits = 1)
PRF_Min_P_val_log = min(Combined_expression_df$P.Value.log)
PRF_Min_P_val_log = round(PRF_Min_P_val_log, digits = 1)
PRF_P_val_seq = seq(from = PRF_Min_P_val_log, to = PRF_Max_P_val_log, by = 0.05)

# loop
for (i in 1:length(PRF_P_val_seq)){
  PRF_current_P_val_log = PRF_P_val_seq[i]
  PRF_current_df = Combined_expression_df[Combined_expression_df$P.Value.log >= PRF_current_P_val_log,]
  PRF_current_df_GSE98793 = PRF_current_df[PRF_current_df$Contrast == "GSE98793: MDD/Control",]
  PRF_current_df_GSE46743 = PRF_current_df[PRF_current_df$Contrast == "GSE46743: case/control",]
  PRF_current_df_GSE64930 = PRF_current_df[PRF_current_df$Contrast == "GSE64930: case/control",]
  
  # Counting significant Genes
  PRF_sign_GSE98793 = length(unique(PRF_current_df_GSE98793$Gene))
  PRF_sign_GSE46743 = length(unique(PRF_current_df_GSE46743$Gene))
  PRF_sign_GSE64930 = length(unique(PRF_current_df_GSE64930$Gene))
  PRF_sign_sum = sum(PRF_sign_GSE98793, PRF_sign_GSE46743, PRF_sign_GSE64930)
  PRF_sign_unqiue = length(unique(PRF_current_df$Gene))
  # Getting overlapping Genes
  PRF_sign_GSE98793_genes = unique(PRF_current_df_GSE98793$Gene)
  PRF_sign_GSE46743_genes = unique(PRF_current_df_GSE46743$Gene)
  PRF_sign_GSE64930_genes = unique(PRF_current_df_GSE64930$Gene)
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes))
  PRF_intersect_all_index = sapply(PRF_intersect_all, function(PRF_x){
    PRF_LFs = PRF_current_df[PRF_current_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_all_index) < 1){
    PRF_intersect_all_index = FALSE
  }
  PRF_intersect_all = PRF_intersect_all[PRF_intersect_all_index]
  
  PRF_intersect_GSE98793_GSE46743 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes)
  PRF_intersect_GSE98793_GSE46743_index = sapply(PRF_intersect_GSE98793_GSE46743, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE98793, PRF_current_df_GSE46743)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE98793_GSE46743_index) < 1){
    PRF_intersect_GSE98793_GSE46743_index = FALSE
  }
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743_index]
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE98793_GSE64930 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE98793_GSE64930_index = sapply(PRF_intersect_GSE98793_GSE64930, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE98793, PRF_current_df_GSE64930)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE98793_GSE64930_index) < 1){
    PRF_intersect_GSE98793_GSE64930_index = FALSE
  }
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930_index]
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE46743_GSE64930 = intersect(PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE46743_GSE64930_index = sapply(PRF_intersect_GSE46743_GSE64930, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE46743, PRF_current_df_GSE64930)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE46743_GSE64930_index) < 1){
    PRF_intersect_GSE46743_GSE64930_index = FALSE
  }
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930_index]
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_GSE98793_GSE46743_count = length(PRF_intersect_GSE98793_GSE46743)
  PRF_intersect_GSE98793_GSE64930_count = length(PRF_intersect_GSE98793_GSE64930)
  PRF_intersect_GSE46743_GSE64930_count = length(PRF_intersect_GSE46743_GSE64930)
  
  # Percent total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_all = round((PRF_intersect_GSE98793_GSE46743_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_all = round((PRF_intersect_GSE98793_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_all = round((PRF_intersect_GSE46743_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  # Percent intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_GSE98793_GSE46743, PRF_intersect_GSE98793_GSE64930, PRF_intersect_GSE46743_GSE64930)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_intersect = round((PRF_intersect_GSE98793_GSE46743_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_intersect = round((PRF_intersect_GSE98793_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_intersect = round((PRF_intersect_GSE46743_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_output_df = data.frame(
    PRF_sign_GSE98793 = PRF_sign_GSE98793,
    PRF_sign_GSE46743 = PRF_sign_GSE46743,
    PRF_sign_GSE64930 = PRF_sign_GSE64930,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_GSE98793_GSE46743_count = PRF_intersect_GSE98793_GSE46743_count,
    PRF_intersect_GSE98793_GSE64930_count = PRF_intersect_GSE98793_GSE64930_count,
    PRF_intersect_GSE46743_GSE64930_count = PRF_intersect_GSE46743_GSE64930_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_GSE98793_GSE46743_count_percent_all = PRF_intersect_GSE98793_GSE46743_count_percent_all,
    PRF_intersect_GSE98793_GSE64930_count_percent_all = PRF_intersect_GSE98793_GSE64930_count_percent_all,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_GSE98793_GSE46743_count_percent_intersect = PRF_intersect_GSE98793_GSE46743_count_percent_intersect,
    PRF_intersect_GSE98793_GSE64930_count_percent_intersect = PRF_intersect_GSE98793_GSE64930_count_percent_intersect,
    PRF_intersect_GSE46743_GSE64930_count_percent_intersect = PRF_intersect_GSE46743_GSE64930_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(P.Val.10.log = PRF_current_P_val_log, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_P_val[[i]] = PRF_output_df
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Dataset_Step_P_val = list_to_df(Dataset_Step_P_val)
Dataset_Step_P_val_2 = Dataset_Step_P_val[stri_detect_fixed(Dataset_Step_P_val$Score_name, pattern = "_count_percent_all"),]
Dataset_Step_P_val_2$Score_name = factor(Dataset_Step_P_val_2$Score_name, levels = unique(Dataset_Step_P_val_2$Score_name),
                                         labels = c("Intresect all cohorts",
                                                    "GSE98793 & GSE46743",
                                                    "GSE98793 & GSE64930",
                                                    "GSE46743 & GSE64930",
                                                    "Intersect in any"))
plot = ggplot(data = Dataset_Step_P_val_2, aes(x = P.Val.10.log, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) +
  scale_color_brewer(palette="Set2") +
  scale_x_continuous(n.breaks = 10, limits = c(NA,4.5)) +
  scale_y_continuous(n.breaks = 50, limits = c(NA,11), expand = c(0, 0)) +
  labs(y = "% of overlapping proteins (similar dir.)", x = "-log10 p-value", col = "Overlap type", shape = "Overlap type", linetype = "Overlap type") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed",
             color = "red", size=0.5) +
  theme(
    legend.position = c(0.8, 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17), # Customizing the plot title
    axis.title = element_text(face = "bold", size = 14), # Customizing the axes titles
    legend.title = element_text(face = "bold", size = 14), # Customizing the legend title
    axis.text = element_text(size = 13, colour = "black"), # Customizing the axes text
    legend.text = element_text(size = 13, colour = "black"), # Customizing the legend text
    panel.background = element_blank(), # Removing ugly gray background
    axis.line = element_line(size = 0.5), # Adding axis lines
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
    strip.background = element_blank(),  # Removes facet background
    strip.text = element_text(face = "bold", size = 14, colour = "blue") # Customizing facet text
  )
plot
pdf("Overlap_Gene_Expr_vs_P_val_similar_dir.pdf", width = 10, height = 10)
plot
dev.off()

# Significant probes (similar direction) vs logFC
Dataset_Step_Lf = list()

# Making steps in logFC
PRF_Max_LF = max(abs(Combined_expression_df$logFC))
PRF_Max_LF = round(PRF_Max_LF, digits = 1)
PRF_Max_LF_seq = seq(from = 0.1, to = PRF_Max_LF, by = 0.02)

# loop
for (i in 1:length(PRF_Max_LF_seq)){
  PRF_current_LF = PRF_Max_LF_seq[i]
  PRF_current_df = Combined_expression_df[abs(Combined_expression_df$logFC) >= PRF_current_LF,]
  PRF_current_df_GSE98793 = PRF_current_df[PRF_current_df$Contrast == "GSE98793: MDD/Control",]
  PRF_current_df_GSE46743 = PRF_current_df[PRF_current_df$Contrast == "GSE46743: case/control",]
  PRF_current_df_GSE64930 = PRF_current_df[PRF_current_df$Contrast == "GSE64930: case/control",]
  
  # Counting significant genes
  PRF_sign_GSE98793 = length(unique(PRF_current_df_GSE98793$Gene))
  PRF_sign_GSE46743 = length(unique(PRF_current_df_GSE46743$Gene))
  PRF_sign_GSE64930 = length(unique(PRF_current_df_GSE64930$Gene))
  PRF_sign_sum = sum(PRF_sign_GSE98793, PRF_sign_GSE46743, PRF_sign_GSE64930)
  PRF_sign_unqiue = length(unique(PRF_current_df$Gene))
  
  # Getting overlapping genes
  PRF_sign_GSE98793_genes = unique(PRF_current_df_GSE98793$Gene)
  PRF_sign_GSE46743_genes = unique(PRF_current_df_GSE46743$Gene)
  PRF_sign_GSE64930_genes = unique(PRF_current_df_GSE64930$Gene)
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes))
  PRF_intersect_all_index = sapply(PRF_intersect_all, function(PRF_x){
    PRF_LFs = PRF_current_df[PRF_current_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_all_index) < 1){
    PRF_intersect_all_index = FALSE
  }
  PRF_intersect_all = PRF_intersect_all[PRF_intersect_all_index]
  
  PRF_intersect_GSE98793_GSE46743 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE46743_genes)
  PRF_intersect_GSE98793_GSE46743_index = sapply(PRF_intersect_GSE98793_GSE46743, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE98793, PRF_current_df_GSE46743)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE98793_GSE46743_index) < 1){
    PRF_intersect_GSE98793_GSE46743_index = FALSE
  }
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743_index]
  PRF_intersect_GSE98793_GSE46743 = PRF_intersect_GSE98793_GSE46743[PRF_intersect_GSE98793_GSE46743 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE98793_GSE64930 = intersect(PRF_sign_GSE98793_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE98793_GSE64930_index = sapply(PRF_intersect_GSE98793_GSE64930, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE98793, PRF_current_df_GSE64930)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE98793_GSE64930_index) < 1){
    PRF_intersect_GSE98793_GSE64930_index = FALSE
  }
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930_index]
  PRF_intersect_GSE98793_GSE64930 = PRF_intersect_GSE98793_GSE64930[PRF_intersect_GSE98793_GSE64930 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE46743_GSE64930 = intersect(PRF_sign_GSE46743_genes, PRF_sign_GSE64930_genes)
  PRF_intersect_GSE46743_GSE64930_index = sapply(PRF_intersect_GSE46743_GSE64930, function(PRF_x){
    PRF_curr_small_df = rbind(PRF_current_df_GSE46743, PRF_current_df_GSE64930)
    PRF_LFs = PRF_curr_small_df[PRF_curr_small_df$Gene == PRF_x, "logFC"]
    PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
    PRF_probes_greater_0 = which(PRF_LFs > 0)
    PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
    PRF_probes_smaller_0 = which(PRF_LFs < 0)
    PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
    PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
    if (any(PRF_combined == 1)){
      return(TRUE)
    }
    return(FALSE)
  })
  if (length(PRF_intersect_GSE46743_GSE64930_index) < 1){
    PRF_intersect_GSE46743_GSE64930_index = FALSE
  }
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930_index]
  PRF_intersect_GSE46743_GSE64930 = PRF_intersect_GSE46743_GSE64930[PRF_intersect_GSE46743_GSE64930 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_GSE98793_GSE46743_count = length(PRF_intersect_GSE98793_GSE46743)
  PRF_intersect_GSE98793_GSE64930_count = length(PRF_intersect_GSE98793_GSE64930)
  PRF_intersect_GSE46743_GSE64930_count = length(PRF_intersect_GSE46743_GSE64930)
  
  # Percent total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_all = round((PRF_intersect_GSE98793_GSE46743_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_all = round((PRF_intersect_GSE98793_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_all = round((PRF_intersect_GSE46743_GSE64930_count/PRF_sign_unqiue)*100, digits = 2)
  # Percent intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_GSE98793_GSE46743, PRF_intersect_GSE98793_GSE64930, PRF_intersect_GSE46743_GSE64930)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE46743_count_percent_intersect = round((PRF_intersect_GSE98793_GSE46743_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE98793_GSE64930_count_percent_intersect = round((PRF_intersect_GSE98793_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE46743_GSE64930_count_percent_intersect = round((PRF_intersect_GSE46743_GSE64930_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_output_df = data.frame(
    PRF_sign_GSE98793 = PRF_sign_GSE98793,
    PRF_sign_GSE46743 = PRF_sign_GSE46743,
    PRF_sign_GSE64930 = PRF_sign_GSE64930,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_GSE98793_GSE46743_count = PRF_intersect_GSE98793_GSE46743_count,
    PRF_intersect_GSE98793_GSE64930_count = PRF_intersect_GSE98793_GSE64930_count,
    PRF_intersect_GSE46743_GSE64930_count = PRF_intersect_GSE46743_GSE64930_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_GSE98793_GSE46743_count_percent_all = PRF_intersect_GSE98793_GSE46743_count_percent_all,
    PRF_intersect_GSE98793_GSE64930_count_percent_all = PRF_intersect_GSE98793_GSE64930_count_percent_all,
    PRF_intersect_GSE46743_GSE64930_count_percent_all = PRF_intersect_GSE46743_GSE64930_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_GSE98793_GSE46743_count_percent_intersect = PRF_intersect_GSE98793_GSE46743_count_percent_intersect,
    PRF_intersect_GSE98793_GSE64930_count_percent_intersect = PRF_intersect_GSE98793_GSE64930_count_percent_intersect,
    PRF_intersect_GSE46743_GSE64930_count_percent_intersect = PRF_intersect_GSE46743_GSE64930_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(LF = PRF_current_LF, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_Lf[[i]] = PRF_output_df
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Dataset_Step_Lf = list_to_df(Dataset_Step_Lf)
Dataset_Step_Lf_2 = Dataset_Step_Lf[stri_detect_fixed(Dataset_Step_Lf$Score_name, pattern = "_count_percent_all"),]
Dataset_Step_Lf_2$Score_name = factor(Dataset_Step_Lf_2$Score_name, levels = unique(Dataset_Step_Lf_2$Score_name),
                                      labels = c("Intresect all cohorts",
                                                 "GSE98793 & GSE46743",
                                                 "GSE98793 & GSE64930",
                                                 "GSE46743 & GSE64930",
                                                 "Intersect in any"))

plot = ggplot(data = Dataset_Step_Lf_2, aes(x = LF, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Dark2") +
  scale_y_continuous(n.breaks = 50, limits = c(NA,11), expand = c(0, 0)) +
  scale_x_continuous(n.breaks = 10, limits = c(NA,0.5)) +
  geom_vline(xintercept = 0.1, linetype="dashed",
             color = "red", size=0.5) +
  labs(y = "% of overlapping genes (similar dir.)", x = "|Log2 fold change|", col = "Overlap type", shape = "Overlap type", linetype = "Overlap type") +
  theme(
    legend.position = c(0.8, 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17), # Customizing the plot title
    axis.title = element_text(face = "bold", size = 14), # Customizing the axes titles
    legend.title = element_text(face = "bold", size = 14), # Customizing the legend title
    axis.text = element_text(size = 13, colour = "black"), # Customizing the axes text
    legend.text = element_text(size = 13, colour = "black"), # Customizing the legend text
    panel.background = element_blank(), # Removing ugly gray background
    axis.line = element_line(size = 0.5), # Adding axis lines
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
    strip.background = element_blank(),  # Removes facet background
    strip.text = element_text(face = "bold", size = 14, colour = "blue") # Customizing facet text
  )
plot
pdf("Overlap_Gene_Expr_vs_LF_similar_dir.pdf", width = 10, height = 10)
plot
dev.off()


################### Making Chromosome maps ###################

# Getting cytoBand track from UCSC
Chromosome_coord_table = chromosome_table_getter()
Full_cytoband_track = list()
mySession = browserSession("UCSC")
genome(mySession) = "hg19"
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_CytoBand_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_Cytoband_track = getTable(mySession, table = "cytoBand", range =  PRF_CytoBand_grange)
  Full_cytoband_track[[i]] = PRF_Cytoband_track
}
Full_cytoband_track = do.call(rbind, Full_cytoband_track)
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Preparing Chromosome table for maps
Chromosome_coord_table_prepared = Chromosome_coord_table
Chromosome_coord_table_prepared$Chromosome = paste0("chr", Chromosome_coord_table_prepared$Chromosome)
Chromosome_coord_table_prepared = Chromosome_coord_table_prepared[,1:2]
Chromosome_coord_table_prepared$start = 1
Centromeres = Full_cytoband_track[Full_cytoband_track$gieStain == "acen",]
Centromere_Starts = Centromeres$chromStart[seq(from= 1, by = 2, to = length(Centromeres$chromStart))]
Centromere_Ends = Centromeres$chromEnd[seq(from= 2, by = 2, to = length(Centromeres$chromEnd))]
Centromeres = Centromeres$chromEnd
Centromeres = Centromeres[seq(from= 1, by = 2, to = length(Centromeres))]
Chromosome_coord_table_prepared = cbind(Chromosome_coord_table_prepared, Centromeres)
colnames(Chromosome_coord_table_prepared) = c("Chromosome", "End", "Start", "Centromere")
Chromosome_coord_table_prepared = Chromosome_coord_table_prepared[c("Chromosome", "Start", "End", "Centromere")]

# Making chromosome dataframe
Chromosome_map_Rideogram = Chromosome_coord_table_prepared
Chromosome_map_Rideogram$Start = 0
Chromosome_map_Rideogram$Centromere = NULL
Chromosome_map_Rideogram$Centromere_St = Centromere_Starts
Chromosome_map_Rideogram$Centromere_End = Centromere_Ends
colnames(Chromosome_map_Rideogram) = c("Chr","Start","End","CE_start","CE_end")

# Preparing combined gene reference track
Names_Ref_1 = c("#name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd","exonCount","exonStarts","exonEnds", "geneSymbol")
Ref_1 = UCSC_Known_Gene_full[,Names_Ref_1]
Names_Ref_2_3 = c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd","exonCount","exonStarts","exonEnds", "name2")
Ref_2 = UCSC_NCBI_RefSeq[,Names_Ref_2_3]
Ref_3 = UCSC_wgEncodeGencodeCompV38lift37[, Names_Ref_2_3]
Ref_4 = UCSC_wgEncodeGencodeCompV38lift37_pseudo[, Names_Ref_2_3]
colnames(Ref_1) = colnames(Ref_2)
Ref_Gene_DS = rbind(Ref_1, Ref_2, Ref_3, Ref_4)
Ref_Gene_DS_list = list()
for (i in 1:nrow(Chromosome_map_Rideogram)){
  PRF_Curr_df = Ref_Gene_DS[Ref_Gene_DS$chrom == Chromosome_map_Rideogram$Chr[i],]
  Ref_Gene_DS_list[[i]] = PRF_Curr_df
}
Ref_Gene_DS_list = list_to_df(Ref_Gene_DS_list)
write.csv(Ref_Gene_DS_list, "Ref_Gene_DS_list.csv")

# Preparing Annotation datasets
Ref_annot_Affymetrix_U133 = U133_Plus_2_Probes_Alignments_PASSED[,c(-9)]
Ref_annot_Illumina_HT_12 = Illumina_HT_12_V3_V4_probes_Alignments_PASSED[,c(-9, -10)]
Ref_annot_Affymetrix_U133 = Ref_annot_Affymetrix_U133[,colnames(Ref_annot_Illumina_HT_12)]
Ref_annot_expression = rbind(Ref_annot_Affymetrix_U133, Ref_annot_Illumina_HT_12)

# Combined expression dataset
Combined_expression_df = list(Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs[,Colnames_to_select_GSE98793],
                              Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930],
                              Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs[,Colnames_to_select_GSE46743_GSE64930])
Combined_expression_df = lapply(Combined_expression_df, function(x){
  colnames(x) = Colnames_to_select_GSE46743_GSE64930
  return(x)
})
Combined_expression_df = do.call(rbind, Combined_expression_df) # 4664 objects in total (not limiting to full intersect, otherwise 3157)

# loop
Chrom_df_Transcriptome_map = list()
for (i in 1:nrow(Chromosome_map_Rideogram)){
  PRF_curr_chrom = Chromosome_map_Rideogram$Chr[i]
  writeLines(PRF_curr_chrom)
  PRF_Curr_chrom_df_Ref_genes = Ref_Gene_DS_list[Ref_Gene_DS_list$chrom == PRF_curr_chrom,]
  PRF_curr_chrom_map = Chromosome_map_Rideogram[i,]
  PRF_curr_chrom_expr_probes = Ref_annot_expression[Ref_annot_expression$chrom == PRF_curr_chrom,]
  
  # Getting the largest trnascript
  PRF_Curr_chrom_df_Ref_genes_2 = list()
  PRF_Curr_Genes_Chromosome = unique(PRF_Curr_chrom_df_Ref_genes$name2)
  for (i_2 in 1:length(PRF_Curr_Genes_Chromosome)){
    PRF_gene_df = PRF_Curr_chrom_df_Ref_genes[PRF_Curr_chrom_df_Ref_genes$name2 == PRF_Curr_Genes_Chromosome[i_2],]
    PRF_gene_df = PRF_gene_df[PRF_gene_df$exonCount == max(PRF_gene_df$exonCount),]
    PRF_gene_df_lengths = mapply(function(x,y){
      PRF_vect = x:y
      PRF_vect = length(PRF_vect)
      return(PRF_vect)
    }, PRF_gene_df$cdsStart, PRF_gene_df$cdsEnd)
    PRF_gene_df_length_max = which.max(PRF_gene_df_lengths)
    PRF_gene_df_length_max = PRF_gene_df_length_max[1]
    PRF_gene_df = PRF_gene_df[PRF_gene_df_length_max,]
    PRF_Curr_chrom_df_Ref_genes_2[[i_2]] = PRF_gene_df
  }
  PRF_Curr_chrom_df_Ref_genes_2 = list_to_df(PRF_Curr_chrom_df_Ref_genes_2)
  
  # Current chromosome gene sequences
  PRF_Curr_chrom_gene_seq = mapply(function(x,y){
    PRF_vct = x:y
    return(PRF_vct)
  }, PRF_Curr_chrom_df_Ref_genes_2$txStart, PRF_Curr_chrom_df_Ref_genes_2$txEnd)
  
  # Getting vectors of transcript sequences 
  PRF_Curr_chrom_transcript_seq = mapply(function(x,y){
    PRF_vct = x:y
    return(PRF_vct)
  }, PRF_curr_chrom_expr_probes$Align_Coord_Start, PRF_curr_chrom_expr_probes$Align_Coord_End)
  
  # Making chromosome vector
  PRF_curr_chrom_vector = seq(from = PRF_curr_chrom_map$Start, by = 1000000, to = PRF_curr_chrom_map$End)
  if (PRF_curr_chrom_vector[length(PRF_curr_chrom_vector)] < PRF_curr_chrom_map$End){
    PRF_curr_chrom_vector = c(PRF_curr_chrom_vector, PRF_curr_chrom_map$End)
  }
  PRF_Curr_Chrom_Table = list()
  
  # Running through intervals (inner chromosomal loop)
  for (index in 1:(length(PRF_curr_chrom_vector)-1)){
    
    if (index == 1){
      
      # The first row
      PRF_Start = PRF_curr_chrom_vector[1]
      PRF_End = PRF_curr_chrom_vector[2]
    } else {
      
      # The second and further rows
      PRF_Start = PRF_curr_chrom_vector[index] + 1
      PRF_End = PRF_curr_chrom_vector[index + 1]
    }
    RPF_Interval = PRF_Start:PRF_End
    PRF_gene_overlap_idx = sapply(PRF_Curr_chrom_gene_seq, function(x){
      if (x[1] > RPF_Interval[length(RPF_Interval)]){
        return(FALSE)
      }
      if (x[length(x)] < RPF_Interval[1]){
        return(FALSE)
      }
      if (any(x %in% RPF_Interval)){
        return(TRUE)
      }
      return(FALSE)
    })
    PRF_identified_ref_gene = PRF_Curr_chrom_df_Ref_genes_2[PRF_gene_overlap_idx,]
    
    # General stats
    PRF_identified_genes = paste0(PRF_identified_ref_gene$name2, collapse = ";")
    PRF_identified_genes_count = length(PRF_identified_ref_gene$name2)
    PRF_transcript_overlap_idx = sapply(PRF_Curr_chrom_transcript_seq, function(x){
      if (x[1] > RPF_Interval[length(RPF_Interval)]){
        return(FALSE)
      }
      if (x[length(x)] < RPF_Interval[1]){
        return(FALSE)
      }
      if (any(x %in% RPF_Interval)){
        return(TRUE)
      }
      return(FALSE)
    })
    PRF_identified_transcripts = PRF_curr_chrom_expr_probes[PRF_transcript_overlap_idx,]
    PRF_identified_transcripts_names = unique(PRF_identified_transcripts$ID)
    PRF_identified_transcripts_names = paste0(PRF_identified_transcripts_names, collapse = ";")
    PRF_identified_transcripts_count = nrow(PRF_identified_transcripts)
    PRF_identified_transcripts_genes_unique = unique(PRF_identified_transcripts$PRF_Gene_symbol)
    PRF_identified_transcripts_genes_unique_count = length(PRF_identified_transcripts_genes_unique)
    PRF_identified_transcripts_genes_unique = paste0(PRF_identified_transcripts_genes_unique, collapse = ";")
    
    # Specific stats
    PRF_signif_tarnscripts = Combined_expression_df[Combined_expression_df$ID %in% PRF_identified_transcripts$ID,]
    if (nrow(PRF_signif_tarnscripts) < 1){
      PRF_signif_tarnscripts_names = NA
      PRF_signif_tarnscripts_count_unique = 0
      PRF_signif_tarnscripts_genes_unique = NA
      PRF_signif_tarnscripts_genes_unique_count = 0
      PRF_Cross_signif_genes = NA
      PRF_Cross_signif_genes_count = 0
      PRF_Cross_signif_genes_similar_dir = NA
      PRF_Cross_signif_genes_similar_dir_count = 0
    } else {
      PRF_signif_tarnscripts_names = unique(PRF_signif_tarnscripts$ID)
      PRF_signif_tarnscripts_names = paste0(PRF_signif_tarnscripts_names, collapse = ";")
      PRF_signif_tarnscripts_count_unique = length(unique(PRF_signif_tarnscripts$ID))
      PRF_signif_tarnscripts_genes_unique = unique(PRF_signif_tarnscripts$Gene)
      PRF_signif_tarnscripts_genes_unique_count = length(PRF_signif_tarnscripts_genes_unique)
      PRF_signif_tarnscripts_genes_unique = paste0(PRF_signif_tarnscripts_genes_unique, collapse = ";")
      
      # Identifying cross-valid genes
      PRF_One_gene_study_df = list()
      PRF_study_levels = unique(PRF_signif_tarnscripts$Contrast)
      for (index2 in 1:length(PRF_study_levels)){
        PRF_curr_study = PRF_signif_tarnscripts[PRF_signif_tarnscripts$Contrast == PRF_study_levels[index2],]
        PRF_curr_study = distinct(PRF_curr_study, Gene, .keep_all = TRUE)
        PRF_One_gene_study_df[[index2]] = PRF_curr_study
      }
      PRF_One_gene_study_df = list_to_df(PRF_One_gene_study_df)
      PRF_genes_table = as.data.frame(table(PRF_One_gene_study_df$Gene))
      PRF_Cross_signif_genes = PRF_genes_table[PRF_genes_table$Freq > 1, "Var1"]
      if (length(PRF_Cross_signif_genes) < 1){
        PRF_Cross_signif_genes = NA
        PRF_Cross_signif_genes_count = 0
        PRF_Cross_signif_genes_similar_dir = NA
        PRF_Cross_signif_genes_similar_dir_count = 0
      } else {
        PRF_Cross_signif_genes_count = length(PRF_Cross_signif_genes)
        PRF_Cross_signif_df = PRF_One_gene_study_df[PRF_One_gene_study_df$Gene %in% PRF_Cross_signif_genes,]
        PRF_Cross_signif_genes = paste0(PRF_Cross_signif_genes, collapse = ";")
        PRF_matching_dir_index = sapply(unique(PRF_Cross_signif_df$Gene), function(PRF_x){
          PRF_LFs = PRF_Cross_signif_df[PRF_Cross_signif_df$Gene == PRF_x, "logFC"]
          PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
          PRF_probes_greater_0 = which(PRF_LFs > 0)
          PRF_probes_greater_0 = length(PRF_probes_greater_0)/length(PRF_LFs)
          PRF_probes_smaller_0 = which(PRF_LFs < 0)
          PRF_probes_smaller_0 = length(PRF_probes_smaller_0)/length(PRF_LFs)
          PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
          if (any(PRF_combined == 1)){
            return(TRUE)
          }
          return(FALSE)
        })
        if (length(PRF_matching_dir_index) < 1){
          PRF_matching_dir_index = FALSE
        }
        PRF_Cross_signif_genes_similar_dir = unique(PRF_Cross_signif_df$Gene)[PRF_matching_dir_index]
        PRF_Cross_signif_genes_similar_dir_count = length(PRF_Cross_signif_genes_similar_dir)
        PRF_Cross_signif_genes_similar_dir = paste0(PRF_Cross_signif_genes_similar_dir, collapse = ";")
      }
    }
    PRF_Current_identified_df = data.frame(
      Chrom = PRF_curr_chrom,Start = PRF_Start, End = PRF_End,
      PRF_identified_genes = PRF_identified_genes,
      PRF_identified_genes_count = PRF_identified_genes_count,
      PRF_identified_transcripts_names = PRF_identified_transcripts_names,
      PRF_identified_transcripts_count = PRF_identified_transcripts_count,
      PRF_identified_transcripts_genes_unique = PRF_identified_transcripts_genes_unique,
      PRF_identified_transcripts_genes_unique_count = PRF_identified_transcripts_genes_unique_count,
      PRF_signif_tarnscripts_names = PRF_signif_tarnscripts_names,
      PRF_signif_tarnscripts_count_unique = PRF_signif_tarnscripts_count_unique,
      PRF_signif_tarnscripts_genes_unique = PRF_signif_tarnscripts_genes_unique,
      PRF_signif_tarnscripts_genes_unique_count = PRF_signif_tarnscripts_genes_unique_count,
      PRF_Cross_signif_genes = PRF_Cross_signif_genes,
      PRF_Cross_signif_genes_count = PRF_Cross_signif_genes_count,
      PRF_Cross_signif_genes_similar_dir = PRF_Cross_signif_genes_similar_dir,
      PRF_Cross_signif_genes_similar_dir_count = PRF_Cross_signif_genes_similar_dir_count
    )
    
    PRF_Curr_Chrom_Table[[index]] = PRF_Current_identified_df
  }
  
  PRF_Curr_Chrom_Table = list_to_df(PRF_Curr_Chrom_Table)
  Chrom_df_Transcriptome_map[[i]] = PRF_Curr_Chrom_Table
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Chrom_df_Transcriptome_map = do.call(rbind, Chrom_df_Transcriptome_map)

# Calculation ratios (normalization to gene counts on arrays)
Chrom_df_Transcriptome_map$Signif_genes_Ratio = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Chrom_df_Transcriptome_map$PRF_signif_tarnscripts_genes_unique_count, Chrom_df_Transcriptome_map$PRF_identified_transcripts_genes_unique_count)
Chrom_df_Transcriptome_map$Signif_genes_cross_valid_Ratio = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Chrom_df_Transcriptome_map$PRF_Cross_signif_genes_count, Chrom_df_Transcriptome_map$PRF_identified_transcripts_genes_unique_count)

# Checking the data
table(Chrom_df_Transcriptome_map$PRF_Cross_signif_genes_similar_dir_count)
sum(Chrom_df_Transcriptome_map$PRF_Cross_signif_genes_similar_dir_count) # 235 genes since this map takes all probes
Gene_names = Chrom_df_Transcriptome_map$PRF_Cross_signif_genes_similar_dir
Gene_names = Gene_names[!is.na(Gene_names)]
Gene_names = Gene_names[Gene_names != ""]
Gene_names = unlist(stri_split_fixed(Gene_names, pattern = ";"))
# Check
any(Gene_names == "NA") # No genes == "NA"
any(duplicated(Gene_names)) # No duplicated...
Gene_names[Gene_names %!in% Genes_matching_dir] # "HSPA6"  "KCNIP3" "LYPD2"  "RPS17"  "FKBP1A" these genes were detected with matching directions in 2 out of 3 cohorts but were 
# -> not presented in all 3 cohorts (thus they can't be adequately analyzed in Enrichment analysis)
Gene_names = Gene_names[Gene_names %in% Expression_Genes_intersect_all] # Now it is 230 genes!
write.csv(Chrom_df_Transcriptome_map, "Chrom_df_Transcriptome_map.csv")
write(Genes_matching_dir, "Genes_matching_dir.txt")

################### Enrichment analysis ###################

# Preparing Entrez dataset
# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/, https://www.ncbi.nlm.nih.gov/gene/
ENTREZ_genes_Homo_Sapiens = smart_fread("/home/aleksandr/Desktop/WORK/GWAS_CAT_CHARACT/data") # Replace with another path if needed
ENTREZ_genes_Homo_Sapiens = ENTREZ_genes_Homo_Sapiens[ENTREZ_genes_Homo_Sapiens$`#tax_id` == "9606",]
ENTREZ_genes_Homo_Sapiens$Ensembl = sapply(ENTREZ_genes_Homo_Sapiens$dbXrefs, function(x){
  x = unlist(stri_split_fixed(x, pattern = "|"))
  x = x[stri_detect_fixed(str = x, pattern = "Ensembl")]
  x = stri_replace_all_fixed(str = x, pattern = "Ensembl:", replacement = "")
  names(x) = NULL
  
  if (length(x) > 1){
    x = paste0(x, collapse = ";")
  }
  
  if (length(x) > 0){
    return(x)
  } else {
    return(NA)
  }
  
})
ENTREZ_genes_Homo_Sapiens$chromosome = paste0("chr", ENTREZ_genes_Homo_Sapiens$chromosome)

# Preparing a Gene Universe
Venn_covered_genes = list(
  "GSE98793" = unique(Analysis_limma_affymetrix_depr_broad_GSE98793$Updated_gene_names),
  "GSE46743" = unique(Analysis_limma_illumina_depr_broad_GSE46743$Gene),
  "GSE64930" = unique(Analysis_limma_illumina_depr_broad_GSE64930$Gene)
)
Transcriptome_intersect_all = Reduce(intersect, list(Venn_covered_genes[[1]], Venn_covered_genes[[2]], Venn_covered_genes[[3]]))
write(Transcriptome_intersect_all, "Transcriptome_intersect_all.txt")

# Mapping to Entrez
ENTREZ_Transcriptome_universe = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Transcriptome_intersect_all){
    return(TRUE)
  }
  return(FALSE)
}),]

ENTREZ_Transcriptome_universe = ENTREZ_Transcriptome_universe[ENTREZ_Transcriptome_universe$chromosome != "chr-",]
Non_mapped_genes = Transcriptome_intersect_all[Transcriptome_intersect_all %!in% ENTREZ_Transcriptome_universe$Symbol]

# Mapping to Entrez through synonyms
ENTREZ_Transcriptome_universe_idx = vector()
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(Ref_annot_expression[Ref_annot_expression$PRF_Gene_symbol %in% Splitted_synonyms, "chrom"])){
      ENTREZ_Transcriptome_universe_idx[i] = TRUE
    } else {
      ENTREZ_Transcriptome_universe_idx[i] = FALSE
    }
  } else {
    ENTREZ_Transcriptome_universe_idx[i] = FALSE
  }
}
ENTREZ_Transcriptome_universe_2 = ENTREZ_genes_Homo_Sapiens[ENTREZ_Transcriptome_universe_idx,] # No additional genes mapped
ENTREZ_Transcriptome_universe_ids = unique(ENTREZ_Transcriptome_universe$GeneID) # 7941 genes were successfully mapped in total
write.csv(ENTREZ_Transcriptome_universe, "ENTREZ_Transcriptome_universe_overlap.csv")

# Mapping genes with matching direction to Entrez
ENTREZ_matching_dir = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Genes_matching_dir){
    return(TRUE)
  }
  return(FALSE)
}),] # 230 genes
ENTREZ_matching_dir = ENTREZ_matching_dir[ENTREZ_matching_dir$chromosome != "chr-",] # 230 genes
write.csv(ENTREZ_matching_dir, "ENTREZ_matching_dir.csv")
Non_mapped_genes = Genes_matching_dir[Genes_matching_dir %!in% ENTREZ_matching_dir$Symbol] # All genes were mapped
ENTREZ_matching_dir_ids = unique(ENTREZ_matching_dir$GeneID)

ENRICHMENT_Transcriptome_Overlap_matching_dir = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_matching_dir_ids,
                                                                                universe = ENTREZ_Transcriptome_universe_ids,
                                                                                categories_to_show = 30,
                                                                                folder = "Depression_Broad_Enrichment_Transcriptome_Overlap_matching_dir",
                                                                                plot_name_pref = "Depression_Broad_Enrichment_Transcriptome_Overlap_matching_dir")

# Transcriptome overlap is strict and requires all cohorts with significant results to overlap
# Methylation overlap is not very strict since DNA methylation is generally known to be very dynamic process and could be also changed as a negative feedback mechanism ->
# -> directions may be the opposite


################### Characterization of used cohorts ###################

# Characterization of GSE98793
GSE98793_CHARACT_DF = GSE98793_phenotypes
table(GSE98793_phenotypes$subject_group)
colnames(GSE98793_CHARACT_DF)
colnames(GSE98793_CHARACT_DF) = multiple_stri_replacer(colnames(GSE98793_CHARACT_DF), 
                                                       pattern_vector = c("gender",
                                                                          "age",
                                                                          "anxiety",
                                                                          "subject_group"), 
                                                       replacement_vector = c("Gender",
                                                                              "Age",
                                                                              "Anxiety",
                                                                              "Diagnosis"))

# Exploring the variables
GSE98793_CHARACT_DF$Anxiety # 38
GSE98793_CHARACT_DF$Age # 37
GSE98793_CHARACT_DF$Gender # 40
GSE98793_CHARACT_DF$Diagnosis # 41
GSE98793_CHARACT_DF$Gender = factor(GSE98793_CHARACT_DF$Gender, levels = c("F", "M"), labels = c("Female", "Male"))
GSE98793_CHARACT_DF$Anxiety = factor(GSE98793_CHARACT_DF$Anxiety, levels = c("yes", "no"), labels = c("Yes", "No"))
GSE98793_CHARACT_DF$Diagnosis = factor(GSE98793_CHARACT_DF$Diagnosis, levels = c("MDD", "Control"), labels = c("MDD", "Healthy"))
GSE98793_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = GSE98793_CHARACT_DF,
                                                                     study_char = "GSE98793",
                                                                     contrast_col_number = 41,
                                                                     contrast_vector = c("MDD", "Healthy"),
                                                                     participants_col_number = 43,
                                                                     model_covariates_col_vector = c(40,37,38),
                                                                     columns_to_characterise_vector = c(41,40,37,38),
                                                                     Remove_NA_predictors = TRUE,
                                                                     drop_P = TRUE,
                                                                     simplif_P = 3)
openxlsx::write.xlsx(GSE98793_CHARACT_DF[["Table"]], file = "GSE98793_demographics.xlsx", overwrite = TRUE)

# Characterization of GSE46743
GSE46743_CHARACT_DF = GSE46743_phenotypes_table_baseline
colnames(GSE46743_CHARACT_DF)
colnames(GSE46743_CHARACT_DF) = multiple_stri_replacer(colnames(GSE46743_CHARACT_DF), 
                                                       pattern_vector = c("age",
                                                                          "bmi",
                                                                          "Depr_status"), 
                                                       replacement_vector = c("Age",
                                                                              "BMI",
                                                                              "Diagnosis"))
# Exploring the variables
GSE46743_CHARACT_DF$Age # 31
GSE46743_CHARACT_DF$BMI # 32
GSE46743_CHARACT_DF$Diagnosis # 30
GSE46743_CHARACT_DF$Diagnosis = factor(GSE46743_CHARACT_DF$Diagnosis, levels = c("case", "control"), labels = c("MDD", "Healthy"))
GSE46743_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = GSE46743_CHARACT_DF,
                                                                     study_char = "GSE46743",
                                                                     contrast_col_number = 30,
                                                                     contrast_vector = c("MDD", "Healthy"),
                                                                     participants_col_number = 11,
                                                                     model_covariates_col_vector = c(31,32),
                                                                     columns_to_characterise_vector = c(30,31,32),
                                                                     Remove_NA_predictors = TRUE,
                                                                     drop_P = TRUE,
                                                                     simplif_P = 3)
openxlsx::write.xlsx(GSE46743_CHARACT_DF[["Table"]], file = "GSE46743_demographics.xlsx", overwrite = TRUE)

# Characterization of GSE64930
GSE64930_CHARACT_DF = GSE64930_phenotypes_baseline
colnames(GSE64930_CHARACT_DF)

# Exploring the variables
GSE64930_CHARACT_DF$Age # 6
GSE64930_CHARACT_DF$Sex # 4
GSE64930_CHARACT_DF$Status # 5
GSE64930_CHARACT_DF$Status = factor(GSE64930_CHARACT_DF$Status, levels = c("case", "control"), labels = c("MDD", "Healthy"))
GSE64930_CHARACT_DF$BMI # 7
GSE64930_CHARACT_DF$RIN # 8
GSE64930_CHARACT_DF$HAMD # 9
GSE64930_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = GSE64930_CHARACT_DF,
                                                                     study_char = "GSE64930",
                                                                     contrast_col_number = 5,
                                                                     contrast_vector = c("MDD", "Healthy"),
                                                                     participants_col_number = 1,
                                                                     model_covariates_col_vector = c(4,6,7,8,10, 11, 12),
                                                                     columns_to_characterise_vector = c(5,4,6,7,8,9),
                                                                     Remove_NA_predictors = TRUE,
                                                                     drop_P = TRUE,
                                                                     simplif_P = 3)
openxlsx::write.xlsx(GSE64930_CHARACT_DF[["Table"]], file = "GSE64930_demographics.xlsx", overwrite = TRUE)









