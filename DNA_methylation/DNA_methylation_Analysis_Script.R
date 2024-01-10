################# This script shows analysis for DNA methylation data ################# 
# Note 1: The equal sign = was used as an assignment operator for typing/productivity reasons
# Note 2: In many cases loops were deliberately used instead of apply functions to enable better control of the variables
# Note 3: Some variables in the loops contain prefixes to enable easy cleanup of the environment once the loop is executed
# Note 4: The whole running environment gets heavy once the data is imported (It is not advised to run a script on a machine with less than 32 GB of RAM)
# Note 5: Some of the raw files (primarily phenotype data (PSY and GSE125105) and DNA methylation data (PSY)) could not be deposited publicly (for ethical reasons)
# -> Raw data files could be provided based on authorized reasonable request to a corresponding author if approved by an ethical review board of Uppsala University
# Note 6: To get additional phenotypic data for GSE125105, please contact MPIP or a corresponding author for further instructions

Working_directory = "~/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/DNA_methylation" # Replace with an appropriate path
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
library(grid)
library(ggplotify)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(chromoMap)
library(RIdeogram)
library(ggVennDiagram)


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
  
  # Running expansion
  df_list = list()
  for (i in 1:nrow(df_const)){
    print(i)
    curr_df_const = df_const[i,, drop = FALSE]
    curr_df_modif = df_modif[i,, drop = FALSE]
    curr_df_modif = apply(curr_df_modif, 2, function(x) unlist(stri_split_fixed(x, pattern = pattern)), simplify = FALSE)
    
    if (length(cols_to_expand) > 1){
      curr_df_modif = do.call(cbind, curr_df_modif)
    } else {
      curr_df_modif = unlist(curr_df_modif)
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

# This function perform differential DNA methylation analysis using limma package. Additionally, it updates all gene symbols, 
# -> performs correction for batch with Combat and removes participants with missing data
# The function variables are labelled with the prefix PREFIX_ (to easily identify all variables if debugging/cleaning the environment)
# This function depends on check_gene_symbol_NIH, make_manhattan_plot_toptable_generalized, and calculate_cumulative_coordinates
test_CpGs_limma_generalized = function(PREFIX_CpGs,
                                       PREFIX_Mvals,
                                       PREFIX_Pheno_df,
                                       PREFIX_contrast_col_number,
                                       PREFIX_participants_col_number,
                                       PREFIX_contrast_vector,
                                       PREFIX_Study,
                                       PREFIX_model.covariates = NULL,
                                       PREFIX_annotation,
                                       PREFIX_plots_folder = NULL,
                                       PREFIX_plots_str = NULL,
                                       PREFIX_use_Combat = FALSE,
                                       PREFIX_Batch_col_number = NULL,
                                       PREFIX_Remove_NA_predictors = FALSE,
                                       PREFIX_caclulate_cumul_positions = FALSE,
                                       PREFIX_Manhattan_plot_name = FALSE,
                                       PREFIX_Log_FC_threshold = 1){
  
  if (all(PREFIX_CpGs %in% rownames(PREFIX_Mvals))){
    writeLines("All CpGs are presented")
  } else {
    PREFIX_Not_present_CpGs = PREFIX_CpGs[PREFIX_CpGs %!in% rownames(PREFIX_Mvals)]
    PREFIX_Not_present_CpGs = paste0(PREFIX_Not_present_CpGs, collapse = "; ")
    PREFIX_Not_present_CpGs = paste0("Several CpGs are not presented: ", PREFIX_Not_present_CpGs)
    writeLines(PREFIX_Not_present_CpGs)
    PREFIX_CpGs = PREFIX_CpGs[PREFIX_CpGs %in% rownames(PREFIX_Mvals)]
    if (length(PREFIX_CpGs) < 1){
      writeLines("No CpGs to test, returning NA")
      return(NA)
    }
  }
  
  if (PREFIX_use_Combat){
    writeLines("Adjusting for Batches with combat")
    colnames(PREFIX_Pheno_df)[PREFIX_Batch_col_number] = "Batch"
    if (!is.factor(PREFIX_Pheno_df$Batch)){
      PREFIX_Pheno_df$Batch = as.factor(PREFIX_Pheno_df$Batch)
      writeLines("Creating Batch factor")
    }
    PREFIX_Mvals = ComBat(PREFIX_Mvals, batch = PREFIX_Pheno_df$Batch)
    PREFIX_Mvals = as.data.frame(PREFIX_Mvals)
  }
  
  # Preparing phenotypes
  colnames(PREFIX_Pheno_df)[PREFIX_contrast_col_number] = "Case_Control"
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Case_Control %in% PREFIX_contrast_vector, ]
  PREFIX_Pheno_df$Case_Control = factor(PREFIX_Pheno_df$Case_Control, levels = c(PREFIX_contrast_vector[2], PREFIX_contrast_vector[1]))
  colnames(PREFIX_Pheno_df)[PREFIX_participants_col_number] = "Participant_ID_FUN"
  
  if (!is.null(PREFIX_model.covariates)){
    PREFIX_full_model_components = c("Case_Control", PREFIX_model.covariates)
    PREFIX_Pheno_df = PREFIX_Pheno_df[,c(PREFIX_full_model_components, "Participant_ID_FUN")]
  } else {
    PREFIX_full_model_components = "Case_Control"
    PREFIX_Pheno_df = PREFIX_Pheno_df[,c(PREFIX_full_model_components, "Participant_ID_FUN")]
  }
  
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
  
  # Correcting pheno and M-values
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Participant_ID_FUN %!in% PREFIX_Participants_missing_data,]
  PREFIX_Mvals = PREFIX_Mvals[,colnames(PREFIX_Mvals) %in% PREFIX_Pheno_df$Participant_ID_FUN]
  PREFIX_Mvals = PREFIX_Mvals[PREFIX_CpGs,]
  print(table(PREFIX_Pheno_df$Case_Control))
  print(nrow(PREFIX_Pheno_df))
  
  # Removing bad predictors (with a single category)
  PREFIX_Test_levels = apply(PREFIX_Pheno_df, 2, function(x) table(x))
  PREFIX_Test_levels = sapply(PREFIX_Test_levels, function(x) length(names(x)))
  PREFIX_Predictor_to_remove = names(PREFIX_Test_levels)[PREFIX_Test_levels < 2]
  PREFIX_full_model_components = PREFIX_full_model_components[PREFIX_full_model_components %!in% PREFIX_Predictor_to_remove]
  PREFIX_Pheno_df = PREFIX_Pheno_df[,colnames(PREFIX_Pheno_df) %!in% PREFIX_Predictor_to_remove]
  
  # Running the analysis using limma
  PREFIX_Model.formula.string = paste0(PREFIX_full_model_components, collapse = " + ")
  PREFIX_Model.formula.string = paste0("~ ",PREFIX_Model.formula.string)
  PREFIX_Design.matrix = model.matrix(as.formula(PREFIX_Model.formula.string), data = PREFIX_Pheno_df)
  PREFIX_fit = lmFit(PREFIX_Mvals, PREFIX_Design.matrix)
  PREFIX_fitE = eBayes(PREFIX_fit)
  PREFIX_Top_table = limma::topTable(fit = PREFIX_fitE, coef = 2, adjust.method = "fdr", number = Inf)
  
  # Annotation
  PREFIX_Probes_test = PREFIX_annotation[rownames(PREFIX_Top_table),]
  PREFIX_Top_table$CpG = rownames(PREFIX_Top_table)
  PREFIX_Top_table$Model.formula.string = PREFIX_Model.formula.string
  PREFIX_Top_table$Contrast = paste0(PREFIX_Study, PREFIX_contrast_vector[1], "/", PREFIX_contrast_vector[2])
  PREFIX_Top_table$Relation_to_Island = PREFIX_Probes_test$Relation_to_Island
  PREFIX_Top_table$Relation_to_Gene = PREFIX_Probes_test$UCSC_RefGene_Group
  PREFIX_Top_table$Relation_to_Gene = sapply(PREFIX_Top_table$Relation_to_Gene, function(x){
    x = unlist(stri_split_fixed(x, pattern = ";"))
    x = unique(x)
    x = paste0(x, collapse = ";")
    return(x)
  })
  PREFIX_Top_table$Regulatory_feature = PREFIX_Probes_test$Regulatory_Feature_Group
  PREFIX_Top_table$Regulatory_feature = sapply(PREFIX_Top_table$Regulatory_feature, function(x){
    x = unlist(stri_split_fixed(x, pattern = ";"))
    x = unique(x)
    x = paste0(x, collapse = ";")
    return(x)
  })
  PREFIX_Top_table$Gene = PREFIX_Probes_test$UCSC_RefGene_Name
  
  # Updated gene names
  PREFIX_probes_genes = lapply(PREFIX_Probes_test$UCSC_RefGene_Name, function(x){
    PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = ";"))
    PREFIX_Current_genes = unique(PREFIX_Current_genes)
    PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
    return(PREFIX_Current_genes)
  })
  PREFIX_probes_genes = unlist(PREFIX_probes_genes)
  PREFIX_probes_genes = unique(PREFIX_probes_genes)
  
  # Checking gene symbols
  PREFIX_probes_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = PREFIX_probes_genes, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                    PRF_replace_NA_with_old = TRUE)
  PREFIX_probes_genes_check_dict = PREFIX_probes_genes_check$Suggested.Symbol
  names(PREFIX_probes_genes_check_dict) = PREFIX_probes_genes_check$x
  PREFIX_Probes_test$Updated_gene_names = sapply(PREFIX_Probes_test$UCSC_RefGene_Name, function(x){
    PREFIX_Current_genes = unlist(stri_split_fixed(x, pattern = ";"))
    PREFIX_Current_genes = unique(PREFIX_Current_genes)
    PREFIX_Current_genes = str_trim(PREFIX_Current_genes)
    PREFIX_New_names = PREFIX_probes_genes_check_dict[PREFIX_Current_genes]
    PREFIX_New_names = paste0(PREFIX_New_names, collapse = ";")
    return(PREFIX_New_names)
  })
  PREFIX_Top_table$Upd_gene_name = PREFIX_Probes_test$Updated_gene_names
  PREFIX_Top_table$Gene = sapply(PREFIX_Top_table$Gene, function(x){
    x = unlist(stri_split_fixed(x, pattern = ";"))
    x = unique(x)
    x = paste0(x, collapse = ";")
    return(x)
  })
  PREFIX_Top_table$Chromosome = PREFIX_Probes_test$chr
  PREFIX_Top_table$Position_Hg19 = PREFIX_Probes_test$pos
  PREFIX_Top_table$Strand = PREFIX_Probes_test$strand
  gc()
  
  if (PREFIX_caclulate_cumul_positions){
    writeLines(paste0("Getting CpG_cumul_pos at  ", Sys.time()))
    CpG_cumul_pos = mapply(function(x,y){
      output = calculate_cumulative_coordinates(chrom = x, pos = y)
      return(output)
    }, PREFIX_Top_table$Chromosome, PREFIX_Top_table$Position_Hg19)
  } else {
    CpG_cumul_pos = 0
  }
  
  PREFIX_Top_table$CpG_cumul_pos = CpG_cumul_pos
  PREFIX_Top_table = PREFIX_Top_table[,c("Contrast", "Model.formula.string", "Gene", "Upd_gene_name",
                                         "CpG", "Chromosome","Position_Hg19", "Strand", 
                                         "CpG_cumul_pos", "Relation_to_Island", "Relation_to_Gene", "Regulatory_feature",
                                         "logFC", "AveExpr", "t" ,"P.Value" , "adj.P.Val","B")]
  
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
  
  # Plotting
  if (is.character(PREFIX_plots_str)){
    # Making volcano plot
    # Adding highlight for transcript based on the fold change
    PREFIX_Log_FC_threshold = abs(PREFIX_Log_FC_threshold)
    PREFIX_Top_table$is.highlight = sapply(PREFIX_Top_table$logFC, function(x){
      if (x > PREFIX_Log_FC_threshold){
        x = "Up"
      } else if (x < - PREFIX_Log_FC_threshold){
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
      PREFIX_Top_table = mutate(PREFIX_Top_table, is_annotate=ifelse(-log10(P.Value)>Pval_treshold & is.highlight != "None", "yes", "no"))
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
    if (min(PREFIX_Top_table$logFC) < -PREFIX_Log_FC_threshold){
      plot = plot + geom_vline(xintercept= -PREFIX_Log_FC_threshold, linetype="dashed", 
                               color = "grey", size=0.5)
    }
    
    if (max(PREFIX_Top_table$logFC) > PREFIX_Log_FC_threshold){
      plot = plot + geom_vline(xintercept= PREFIX_Log_FC_threshold, linetype="dashed", 
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
      # Custom the theme:
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
    
  }
  gc()
  if (is.character(PREFIX_Manhattan_plot_name)){
    writeLines("Making Manhattan Plot")
    make_manhattan_plot_toptable_generalized(Toptable_df = PREFIX_Top_table_output, calculate_cumul_pos = FALSE, plot_name = PREFIX_Manhattan_plot_name, plot_folder = PREFIX_plots_folder)
  }
  gc()
  return(PREFIX_Top_table_output)
}

# This function creates the Manhattan plot and saves it in the specified folder (used internally in the function above)
make_manhattan_plot_toptable_generalized = function(Toptable_df, calculate_cumul_pos = TRUE, plot_name, plot_folder = FALSE, target_genes = FALSE){
  Toptable_df = arrange(Toptable_df, CpG_cumul_pos)
  
  if (calculate_cumul_pos){
    writeLines("Calculating CpGs cumul. pos.")
    Toptable_df$CpG_cumul_pos = mapply(function(x,y){
      output = calculate_cumulative_coordinates(chrom = x, pos = y)
      return(output)
    }, Toptable_df$Chromosome, Toptable_df$Position_Hg19)
    Toptable_df = arrange(Toptable_df, CpG_cumul_pos)
  }
  
  Toptable_df$Chromosome = stri_replace_all_fixed(Toptable_df$Chromosome, pattern = "chr", replacement = "")
  Toptable_df$Chromosome = as.numeric(Toptable_df$Chromosome)
  
  axisdf = list()
  for (i in 1:length(unique(Toptable_df$Chromosome))){
    Chrom = unique(Toptable_df$Chromosome)[i]
    Current_df = Toptable_df[Toptable_df$Chromosome == Chrom,]
    Current_df$CpG_cumul_pos = as.double(Current_df$CpG_cumul_pos)
    center = (max(Current_df$CpG_cumul_pos) + min(Current_df$CpG_cumul_pos))/2
    data_small = data.frame(Chromosome = Chrom, center = center)
    axisdf[[i]] = data_small
  }
  axisdf = list_to_df(axisdf)
  
  # Adding highlight for CpGs based on the gene
  if (is.character(target_genes)){
    writeLines("Adding highlight")
    Toptable_df = mutate(Toptable_df, is.highlight = ifelse(multiple_stri_detector(CpG_gene, pattern_vector = target_genes), 
                                                            "yes", "no"))
  }
  
  # Modification of CpGs DF
  Pval_treshold = Toptable_df[Toptable_df$adj.P.Val < 0.05,]
  Pval_treshold = arrange(Pval_treshold, -P.Value)
  Pval_treshold = Pval_treshold$P.Value[1]
  Pval_treshold = -log10(Pval_treshold)
  Toptable_df = mutate(Toptable_df, is_annotate= ifelse(-log10(Toptable_df$P.Value)>Pval_treshold, "yes", "no"))
  
  # Make the plot
  plot = ggplot(Toptable_df, aes(x=CpG_cumul_pos, y=-log10(P.Value))) +
    
    # Show all points
    geom_point( aes(color=as.factor(Chromosome)), alpha=0.6, size=4) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chromosome, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 1)) +     # Remove space between plot area and x axis
    
    # Add p-val line
    geom_hline(yintercept=Pval_treshold, linetype="dashed", 
               color = "red", size=0.5)
  
  if (is.character(target_genes)){
    # Add highlighted points
    plot = plot + geom_point(data=filter(Toptable_df, is.highlight =="yes"), color="orange", size = 6, alpha= 0.8)
  }
  
  plot = plot + 
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data = subset(Toptable_df, is_annotate=="yes"), aes(label = CpG), size = 8, force = 10, 
                      max.overlaps = 50) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      axis.title.x = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold")
    ) +
    labs(x = "Chromosome")
  
  # Saving the plot
  writeLines("Saving plot")
  if (is.character(plot_folder)){
    plot_name = paste0(plot_folder, "/", plot_name)
  }
  png(plot_name, width = 4096, height = 2160, units = "px")
  print(plot)
  dev.off()
  
}

# This function calculates cumulated genomic coordinates for CpGs (used internally in the function test_CpGs_limma_generalized)
Chromosome_coord_table = smart_fread("/home/aleksandr/Desktop/WORK/depression analysis/manhattan_SCR/Chromosome_coord_table.csv") # helper dataset for the function (replace with an appropriate path)
calculate_cumulative_coordinates = function(Coordinates = NULL, chrom = NULL, pos = NULL, pos_full = FALSE){
  if (!is.null(chrom) & !is.null(pos)){
    pos = as.double(pos)
    Curr_chrom = stri_replace_all_fixed(chrom, pattern = 'chr', replacement = '')
    if (Curr_chrom != "1"){
      Index = which(Chromosome_coord_table$Chromosome == Curr_chrom)
      Curr_table = Chromosome_coord_table[1:(Index-1), ]
      Combined_length = sum(Curr_table$`Total length (bp)`)
      pos = pos + Combined_length
    }
    if (pos_full){
      pos = paste0(chrom, ":", pos)
    }
    return(pos)
  }
  
  if (!is.null(Coordinates)){
    Coordinates = unlist(stri_split_fixed(Coordinates, pattern = c(':')))
    Coordinates = unlist(stri_split_fixed(Coordinates, pattern = c('-')))
    Curr_chrom = Coordinates[1]
    Curr_chrom = stri_replace_all_fixed(Curr_chrom, pattern = 'chr', replacement = '')
    if (Curr_chrom == "1"){
      Start = Coordinates[2]
      End = Coordinates[3]
    } else {
      Start = as.numeric(Coordinates[2])
      End = as.numeric(Coordinates[3])
      Index = which(Chromosome_coord_table$Chromosome == Curr_chrom)
      Curr_table = Chromosome_coord_table[1:(Index-1), ]
      Combined_length = sum(Curr_table$`Total length (bp)`)
      Start = Start + Combined_length
      End = End + Combined_length
    }
    coord_string = paste0(Start,"-", End)
    return(coord_string)
  }
  
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

# A function to filter Illumina annotation files to remove bad probes (same as in data preprocessing)
filter_annotation_illumina = function(annot_df, mode){
  
  # Standard SNP filtering
  annot_df = annot_df[stri_detect_fixed(rownames(annot_df), pattern = "cg"),]
  annot_df = annot_df[is.na(annot_df[,'Probe_maf']) | annot_df[,'Probe_maf']<0.05,]
  annot_df = annot_df[is.na(annot_df[,'CpG_rs']) & is.na(annot_df[,'SBE_rs']),]
  
  # Filtering based on the array type
  if (mode == "450K"){
    Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
    annot_df = annot_df[rownames(annot_df) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
    Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
    annot_df = annot_df[rownames(annot_df) %!in% Cross_reactive_probes_MILES_BENTON$V1,]
  }
  
  if (mode == "EPIC"){
    Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
    annot_df = annot_df[rownames(annot_df) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
    Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
    annot_df = annot_df[rownames(annot_df) %!in% Cross_reactive_probes_MILES_BENTON$V1,]
    
    # EPIC-specific filtering
    Files_EPIC_bad_probes = list.files("/home/aleksandr/Desktop/WORK/Files with pverlapping probes EPIC") #https://github.com/sirselim/illumina450k_filtering; #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1#Sec22
    Files_EPIC_bad_probes = Files_EPIC_bad_probes[stri_detect_fixed(Files_EPIC_bad_probes, pattern = ".csv")]
    Files_EPIC_bad_probes = paste0("/home/aleksandr/Desktop/WORK/Files with pverlapping probes EPIC/", Files_EPIC_bad_probes)
    Files_EPIC_bad_probes = lapply(Files_EPIC_bad_probes, function(x) as.data.frame(fread(x)))
    Bad_probes_Epic = lapply(Files_EPIC_bad_probes, function(x) x[,1])
    Bad_probes_Epic = unlist(Bad_probes_Epic)
    Bad_probes_Epic = unique(Bad_probes_Epic)
    annot_df = annot_df[rownames(annot_df) %!in% Bad_probes_Epic,]
  }
  
  # Adding updated gene names
  # Checking symbols
  probes_genes = lapply(annot_df$UCSC_RefGene_Name, function(x){
    Current_genes = unlist(stri_split_fixed(x, pattern = ";"))
    Current_genes = unique(Current_genes)
    Current_genes = str_trim(Current_genes)
    return(Current_genes)
  })
  probes_genes = unlist(probes_genes)
  probes_genes = unique(probes_genes)
  writeLines("Getting Updated Symbols")
  probes_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = probes_genes, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                    PRF_replace_NA_with_old = TRUE)
  probes_genes_check_dict = probes_genes_check$Suggested.Symbol
  names(probes_genes_check_dict) = probes_genes_check$x
  writeLines("Adding Updated Symbols")
  annot_df$Updated_gene_names = sapply(annot_df$UCSC_RefGene_Name, function(x){
    Current_genes = unlist(stri_split_fixed(x, pattern = ";"))
    Current_genes = unique(Current_genes)
    Current_genes = str_trim(Current_genes)
    New_names = probes_genes_check_dict[Current_genes]
    New_names = paste0(New_names, collapse = ";")
    return(New_names)
  })
  return(annot_df)
}

# Function to perform enrichment for GO terms and KEGG and create figures that may be useful
# All results are saved in the specified folder
# Genes and universe are accepted as vectors of Entrez IDs
# Some of the plotting is performed within tryCatch since these plots may not be rendered depending on the enrichment results
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
      
      # GO induced graph
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
      
      # Tree plot
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
      
      
      # Enrichment map
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
  
  # Preparing outputs to load into the environment
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


################### Importing and preparing PSY data ###################
# Note: data files for the PSY could be provided based on authorized reasonable request to a corresponding author 
#  -> if approved by an ethical review board of Uppsala University
M_values = read.table("/home/aleksandr/Desktop/WORK/Preprocessing_SCR_Combinded/M_values_adj_221_PSY_screen_comb_prepr.txt", sep = "\t") # Replace with an appropriate path
Beta_values = read.table("/home/aleksandr/Desktop/WORK/Preprocessing_SCR_Combinded/Beta_values_adj_221_PSY_screen_comb_prepr.txt", sep = "\t") # Replace with an appropriate path
Sample_IDs_SCR = read.csv("/home/aleksandr/Desktop/WORK/Preprocessing_SCR_Combinded/SAMPLE_ID_FOR_PSY_screen_comb_prepr.csv") # Replace with an appropriate path
Colnames_BM_values = colnames(M_values)

# Fixing column names
Colnames_BM_values = stri_replace_all_fixed(Colnames_BM_values ,pattern = "X", replacement = "")
New_Col_names = sapply(Colnames_BM_values, function(x){
  index = which(Sample_IDs_SCR$barcode == x)
  index = index[1]
  ID = Sample_IDs_SCR$Participant[index]
  return(ID)
})
any(duplicated(New_Col_names)) # no duplicated columns
colnames(M_values) = New_Col_names
colnames(Beta_values) = New_Col_names

# Preparing phenotypes data
Sample_IDs_SCR = Sample_IDs_SCR[Sample_IDs_SCR$barcode %in% Colnames_BM_values,]
Phenotypes_SCR = read.csv("/home/aleksandr/Desktop/WORK/depression analysis/Screening_pheno_small.csv", encoding = "UTF-8") # Replace with an appropriate path
Phenotypes_SCR = inner_join(Sample_IDs_SCR[Sample_IDs_SCR$barcode %in% Colnames_BM_values,], Phenotypes_SCR, by = c("Participant" = "Code"))
all(Phenotypes_SCR$Participant == New_Col_names) # order is matching
rownames(Phenotypes_SCR) = New_Col_names

# Inspecting data
str(Phenotypes_SCR)
Phenotypes_SCR = Phenotypes_SCR[colnames(M_values),]
all(Phenotypes_SCR$Participant == colnames(M_values)) # order is matching

# Setting variables
Phenotypes_SCR$Study = factor(Phenotypes_SCR$Study)
Phenotypes_SCR$Gender = factor(Phenotypes_SCR$Gender)
Phenotypes_SCR$Depr_risk = ifelse(Phenotypes_SCR$DAWBA_DEPBAND > 3, "High depr. risk", "Low depr. risk")
Phenotypes_SCR$Depr_risk = factor(Phenotypes_SCR$Depr_risk, levels = c("High depr. risk", "Low depr. risk"))
Phenotypes_SCR = Phenotypes_SCR[Phenotypes_SCR$Participant %in% colnames(M_values),]
all(Phenotypes_SCR$Participant == colnames(M_values)) # order is matching

# Additional phenotypes
Genotyped_participants_PSY = smart_fread("/home/aleksandr/Desktop/WORK/depression analysis/depression_SNPs_153_all_data_model.csv") # Replace with an appropriate path
Genotyped_participants_PSY = Genotyped_participants_PSY$ID
Phenotypes_SCR$Is_genotyped = sapply(Phenotypes_SCR$Participant, function(x){
  if (x %in% Genotyped_participants_PSY){
    return("YES")
  } else {
    return("NO")
  }
})
Phenotypes_SCR$Is_genotyped = factor(Phenotypes_SCR$Is_genotyped, levels = c("YES", "NO"))


################### Importing and preparing GSE72680 data ###################

# Reading data
E_GEOD_72680_beta_val = smart_fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/beta_corrected_GSE72680.txt") # Replace with an appropriate path
E_GEOD_72680_Mval = smart_fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/Mval_corrected_GSE72680.txt") # Replace with an appropriate path 
E_GEOD_72680_Phenotypes = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/GSE_version/data", skip = 30, nThread = 14)) # Replace with an appropriate path
# The data file was downloaded from GEO (series matrix file)
# Structure is complex and should be reformatted

# Preparing scaffold data dataframe with phenotypes
E_GEOD_72680_Phenotypes_Scaffold = E_GEOD_72680_Phenotypes[-(8:36),]
E_GEOD_72680_Phenotypes_Scaffold = t(E_GEOD_72680_Phenotypes_Scaffold)
E_GEOD_72680_Phenotypes_GEO_ID = rownames(E_GEOD_72680_Phenotypes_Scaffold)[-1]
colnames(E_GEOD_72680_Phenotypes_Scaffold) = E_GEOD_72680_Phenotypes_Scaffold[1,]
E_GEOD_72680_Phenotypes_Scaffold = E_GEOD_72680_Phenotypes_Scaffold[-1,]
E_GEOD_72680_Phenotypes_Scaffold = cbind(E_GEOD_72680_Phenotypes_GEO_ID, E_GEOD_72680_Phenotypes_Scaffold)
E_GEOD_72680_Phenotypes_Scaffold = as.data.frame(E_GEOD_72680_Phenotypes_Scaffold)

# Getting characteristics
Characteristics_matrix = E_GEOD_72680_Phenotypes[8:36, ]
Characteristics_matrix = apply(Characteristics_matrix[,-1], 2, function(x){
  All_traits = x
  All_traits = All_traits[!is.na(All_traits)]
  All_traits = All_traits[All_traits != ""]
  All_traits = sapply(All_traits, function(z){
    z = unlist(stri_split_fixed(z, pattern = ": "))
    z = z[1]
    return(z)})
  return(All_traits)
})
Characteristics_matrix = unlist(Characteristics_matrix)
names(Characteristics_matrix) = NULL
Characteristics_matrix = unique(Characteristics_matrix)
TMP_DF_E_GEOD_72680 = lapply(Characteristics_matrix, function(x){
  x = rep(NA, times = 422)
  return(x)
}) 
names(TMP_DF_E_GEOD_72680) = Characteristics_matrix
TMP_DF_E_GEOD_72680 = do.call(cbind, TMP_DF_E_GEOD_72680)
TMP_DF_E_GEOD_72680 = as.data.frame(TMP_DF_E_GEOD_72680)
colnames(TMP_DF_E_GEOD_72680) = stri_replace_all_fixed(Characteristics_matrix, pattern = " ", replacement = "_")
TMP_DF_E_GEOD_72680 = cbind(E_GEOD_72680_Phenotypes_GEO_ID, TMP_DF_E_GEOD_72680)

# Looping through participants and obtaining associated traits
i = 1 
while(i <= nrow(TMP_DF_E_GEOD_72680)){
  Current_participant = TMP_DF_E_GEOD_72680$E_GEOD_72680_Phenotypes_GEO_ID[i]
  Subset_data = E_GEOD_72680_Phenotypes[8:36, Current_participant]
  for (trait in Characteristics_matrix){
    if (any(stri_detect_fixed(str = Subset_data, pattern = trait))){
      Current_char = Subset_data[stri_detect_fixed(str = Subset_data, pattern = trait)]
      Value_to_return = unlist(stri_split_fixed(Current_char, pattern = ": "))
      Value_to_return = Value_to_return[2]
      TMP_DF_E_GEOD_72680[i,stri_replace_all_fixed(trait, pattern = " ", replacement = "_")] = Value_to_return
    }
  }
  rm(list = c("Current_participant", "Subset_data", "trait", "Current_char", "Value_to_return"))
  print(i)
  i = i + 1
}
E_GEOD_72680_Phenotypes = cbind(TMP_DF_E_GEOD_72680, E_GEOD_72680_Phenotypes_Scaffold[,-1])

# Removing temporary files
rm(list = c("E_GEOD_72680_Phenotypes_Scaffold", "TMP_DF_E_GEOD_72680", 
            "Characteristics_matrix", "E_GEOD_72680_Phenotypes_GEO_ID"))
# Matching order of beta-values, M-values, and phenotypes
E_GEOD_72680_beta_val = E_GEOD_72680_beta_val[,E_GEOD_72680_Phenotypes$`!Sample_description.1`]
E_GEOD_72680_Mval = E_GEOD_72680_Mval[, E_GEOD_72680_Phenotypes$`!Sample_description.1`]

# Check
colnames(E_GEOD_72680_beta_val) == E_GEOD_72680_Phenotypes$`!Sample_description.1`
colnames(E_GEOD_72680_Mval) == E_GEOD_72680_Phenotypes$`!Sample_description.1`

# Inspecting the data and preparing variables
E_GEOD_72680_Phenotypes$Batch = sapply(E_GEOD_72680_Phenotypes$`!Sample_description.1`, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
E_GEOD_72680_Phenotypes$Batch = factor(E_GEOD_72680_Phenotypes$Batch)
str(E_GEOD_72680_Phenotypes)
E_GEOD_72680_Phenotypes$age = as.numeric(E_GEOD_72680_Phenotypes$age)
E_GEOD_72680_Phenotypes$Sex = as.factor(E_GEOD_72680_Phenotypes$Sex)

# Fixing child abuse
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` = ifelse(E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` == "--", NA,
                                                                                       E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`)
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` = factor(E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`, 
                                                                                       levels = c("Yes", "No"))
# Fixing treatment_for_depression
E_GEOD_72680_Phenotypes$treatment_for_depression
E_GEOD_72680_Phenotypes$treatment_for_depression = ifelse(E_GEOD_72680_Phenotypes$treatment_for_depression == "--", NA,
                                                          E_GEOD_72680_Phenotypes$treatment_for_depression)
E_GEOD_72680_Phenotypes$treatment_for_depression = factor(E_GEOD_72680_Phenotypes$treatment_for_depression, levels = c("Yes", "No"))

# Fixing bipolar disorder
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder == "--", NA,
                                                                E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder)
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder, levels = c("Yes", "No"))

# Fixing posttraumatic_stress_disorder 
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder == "--", NA,
                                                                             E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder)
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder, levels = c("Yes", "No"))

# Fixing anxiety
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder == "--", NA,
                                                                E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder)
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder, levels = c("Yes", "No"))

# Numerical scores
E_GEOD_72680_Phenotypes[,9:21] = apply(E_GEOD_72680_Phenotypes[,9:21],2, function(x){
  x = str_trim(x)
  x = as.numeric(x)
  return(x)
})
E_GEOD_72680_Phenotypes$`race/ethnicity` = as.factor(E_GEOD_72680_Phenotypes$`race/ethnicity`)
E_GEOD_72680_Phenotypes[,24:30] = apply(E_GEOD_72680_Phenotypes[,24:30],2, function(x){
  x = str_trim(x)
  x = as.numeric(x)
  return(x)
})

# BDI score categorization
# https://www.ismanet.org/doctoryourspirit/pdfs/Beck-Depression-Inventory-BDI.pdf
# https://www.pnas.org/content/pnas/suppl/2019/05/20/1816847116.DCSupplemental/pnas.1816847116.sapp.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905604/

# Inspecting
E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score
summary(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score)

# Making categories based on thresholds
# Standard
E_GEOD_72680_Phenotypes$Beck_depression_binary = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=19){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary)

# All categories
E_GEOD_72680_Phenotypes$Beck_depression_full = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x %in% 0:10){
    return("Normal")
  }
  if (x %in% 11:16){
    return("Mild.mood.dist")
  }
  if (x %in% 17:20){
    return("Borderline.clin.depr")
  }
  if (x %in% 21:30){
    return("Mod.depr")
  }
  if (x %in% 31:40){
    return("Sev.depr")
  }
  if (x > 40){
    return("Extr.depr")
  }
})
table(E_GEOD_72680_Phenotypes$Beck_depression_full)
E_GEOD_72680_Phenotypes$Beck_depression_full = ordered(E_GEOD_72680_Phenotypes$Beck_depression_full, levels = c("Normal",
                                                                                                                "Mild.mood.dist",
                                                                                                                "Borderline.clin.depr",
                                                                                                                "Mod.depr",
                                                                                                                "Sev.depr",
                                                                                                                "Extr.depr"))
# Strict threshold (without borderline depression)
E_GEOD_72680_Phenotypes$Beck_depression_binary_strict = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=21){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_strict = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_strict, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_strict)

# Threshold for severe depression
E_GEOD_72680_Phenotypes$Beck_depression_binary_severe = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=31){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_severe = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_severe, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_severe)

# Threshold for extreme depression
E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=41){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme)

# Making composite depression categorization based on BDI category (strict) and different logic of interpreting missing values
E_GEOD_72680_Phenotypes$Composite_depression = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & !is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  } else if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & !is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_NA = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = NA
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]<21){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Normal"
    }
  } else if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Normal"
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_Strict = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) | is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = NA
  } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = "Depressed"
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = "Normal"
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_NA_full = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Depressed"
  } else if (!is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]) & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21) {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Depressed"
  } else if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i] < 21) {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] =  NA
  } else if (is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]) & E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No"){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] =  NA
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Normal"
  }
}

table(E_GEOD_72680_Phenotypes$Composite_depression, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_NA, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_Strict, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_NA_full, exclude = NULL)

# Preparing factors for composite scores
E_GEOD_72680_Phenotypes$Composite_depression = factor(E_GEOD_72680_Phenotypes$Composite_depression, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_NA = factor(E_GEOD_72680_Phenotypes$Composite_depression_NA, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_Strict = factor(E_GEOD_72680_Phenotypes$Composite_depression_Strict, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_NA_full = factor(E_GEOD_72680_Phenotypes$Composite_depression_NA_full, levels = c("Depressed", "Normal"))

# Fixing columns of the pheno data
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "!", replacement = "")
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "-", replacement = "_")
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "/", replacement = "_")

# Check
colnames(E_GEOD_72680_beta_val) == E_GEOD_72680_Phenotypes$Sample_description.1
colnames(E_GEOD_72680_Mval) == E_GEOD_72680_Phenotypes$Sample_description.1
summary(E_GEOD_72680_Phenotypes[,10:15]) # No missing values among the cell types


################### Importing and preparing GSE125105 data ###################

# Reading data and performin curation of series matrix file
GSE125105_pheno_init = readLines("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/data_GSE125105") # the data file was downloaded from GEO (series matrix file)
GSE125105_pheno_init_small = substr(x = GSE125105_pheno_init, 1, 200)
GSE125105_pheno_selected = GSE125105_pheno_init[c(27:46,48,49,50,53,55,56,57,58,59,60,61,62,63,67)]
GSE125105_pheno_selected = stri_split_fixed(GSE125105_pheno_selected, pattern = "\t")
Names_list_GSE125105_pheno_selected = lapply(GSE125105_pheno_selected, function(x){
  Name = x[[1]]
  if (Name == "!Sample_characteristics_ch1"){
    Name_2 = x[[2]]
    Name_2 = unlist(stri_split_fixed(Name_2, pattern = ":"))
    Name_2 = Name_2[[1]]
    Name = Name_2
  }
  Name = stri_replace_all_fixed(Name, pattern = '"', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '!', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '-', replacement = "_")
  return(Name)
})
Names_list_GSE125105_pheno_selected = unlist(Names_list_GSE125105_pheno_selected)
GSE125105_pheno_selected = lapply(GSE125105_pheno_selected, function(x){
  x = x[-1]
  x = sapply(x, function(z){
    z = unlist(stri_split_fixed(z, pattern = ":"))
    if (length(z)>1){
      z = z[2]
    }
    return(z)
  })
  x = stri_replace_all_fixed(x, pattern = '"', replacement = "")
  return(x)
})
GSE125105_pheno_selected = do.call(cbind,GSE125105_pheno_selected)
GSE125105_pheno_selected = as.data.frame(GSE125105_pheno_selected)
Names_list_GSE125105_pheno_selected[duplicated(Names_list_GSE125105_pheno_selected)] = paste0(Names_list_GSE125105_pheno_selected[duplicated(Names_list_GSE125105_pheno_selected)], "_2")
colnames(GSE125105_pheno_selected) = Names_list_GSE125105_pheno_selected

# Reading targets file
GSE125105_target = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/Targets_GSE125105.csv")) # Replace with an appropriate path
GSE125105_pheno_full = inner_join(GSE125105_target, GSE125105_pheno_selected, by = c("Participant" = "ID_REF"))

# Reading M- and beta-values
GSE125105_M_val_corrected = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/M_values_adj_combined_GSE125105.txt", nThread = 14)) # Replace with an appropriate path
rownames(GSE125105_M_val_corrected) = GSE125105_M_val_corrected$V1
GSE125105_M_val_corrected$V1 = NULL

GSE125105_Beta_corrected = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/adjusted_blood_beta_combined_GSE125105.txt", nThread = 14)) # Replace with an appropriate path
rownames(GSE125105_Beta_corrected) = GSE125105_Beta_corrected$V1
GSE125105_Beta_corrected$V1 = NULL

# Attaching calculated blood data (no really needed since values have been adjusted with regression)
GSE125105_blood = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/wbcset.txt", nThread = 14))
GSE125105_pheno_full = inner_join(GSE125105_pheno_full, GSE125105_blood, by = c("Full_file_name" = "V1"))

# Check
colnames(GSE125105_M_val_corrected) == GSE125105_pheno_full$Full_file_name
colnames(GSE125105_Beta_corrected) == GSE125105_pheno_full$Full_file_name
rm(list = c("GSE125105_pheno_init","GSE125105_pheno_init_small","GSE125105_pheno_selected", "Names_list_GSE125105_pheno_selected", "GSE125105_target", "GSE125105_blood"))

# Phenotype curation
GSE125105_pheno_full$diagnosis = str_trim(GSE125105_pheno_full$diagnosis)
GSE125105_pheno_full$diagnosis = ifelse(GSE125105_pheno_full$diagnosis == "control", "Non-depressed", "Depression")
GSE125105_pheno_full$diagnosis = factor(GSE125105_pheno_full$diagnosis, levels = c("Depression", "Non-depressed"))

GSE125105_pheno_full$age = str_trim(GSE125105_pheno_full$age)
GSE125105_pheno_full$age = as.numeric(GSE125105_pheno_full$age)

GSE125105_pheno_full$Sex = str_trim(GSE125105_pheno_full$Sex)
GSE125105_pheno_full$Sex = factor(GSE125105_pheno_full$Sex)

GSE125105_pheno_full[,stri_detect_fixed(colnames(GSE125105_pheno_full), pattern = "cellcount")] = apply(GSE125105_pheno_full[,stri_detect_fixed(colnames(GSE125105_pheno_full), pattern = "cellcount")],
                                                                                                        2, function(x){
                                                                                                          x = str_trim(x)
                                                                                                          x = as.numeric(x)
                                                                                                          return(x)
                                                                                                        })
# Getting extra data from subset GSE128235 (This cohort is a subset of GSE125105 and is deposited separately)
# The cohort GSE128235 has additional phenotypical data for participants
GSE128235_pheno_init = readLines("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/GSE128235_data") # the data file was downloaded from GEO (series matrix file) # Replace with an appropriate path
GSE128235_pheno_init_small = substr(x = GSE128235_pheno_init, 1, 200)
GSE128235_pheno_selected = GSE128235_pheno_init[c(26:49,69, 71)]
GSE128235_pheno_selected = stri_split_fixed(GSE128235_pheno_selected, pattern = "\t")
Names_list_GSE128235_pheno_selected = lapply(GSE128235_pheno_selected, function(x){
  Name = x[[1]]
  if (Name == "!Sample_characteristics_ch1"){
    Name_2 = x[[2]]
    Name_2 = unlist(stri_split_fixed(Name_2, pattern = ":"))
    Name_2 = Name_2[[1]]
    Name = Name_2
  }
  Name = stri_replace_all_fixed(Name, pattern = '"', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '!', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '-', replacement = "_")
  return(Name)
})
Names_list_GSE128235_pheno_selected = unlist(Names_list_GSE128235_pheno_selected)
GSE128235_pheno_selected = lapply(GSE128235_pheno_selected, function(x){
  x = x[-1]
  x = sapply(x, function(z){
    z = unlist(stri_split_fixed(z, pattern = ":"))
    if (length(z)>1){
      z = z[2]
    }
    z = str_trim(z)
    return(z)
  })
  x = stri_replace_all_fixed(x, pattern = '"', replacement = "")
  return(x)
})
GSE128235_pheno_selected = do.call(cbind,GSE128235_pheno_selected)
GSE128235_pheno_selected = as.data.frame(GSE128235_pheno_selected)
Names_list_GSE128235_pheno_selected[duplicated(Names_list_GSE128235_pheno_selected)] = paste0(Names_list_GSE128235_pheno_selected[duplicated(Names_list_GSE128235_pheno_selected)], "_2")
colnames(GSE128235_pheno_selected) = Names_list_GSE128235_pheno_selected
rm(list = c("GSE128235_pheno_init", "GSE128235_pheno_init_small", "Names_list_GSE128235_pheno_selected"))

# Now we need to identify the correct column to match
GSE125105_pheno_full$Participant
GSE128235_pheno_selected$Sample_relation %in% GSE125105_pheno_full$Participant
GSE128235_pheno_selected$Sample_title %in% GSE125105_pheno_full$Sample_title # Sample numbers match, but the way it is written is not

# Adding extra phenotypes
GSE125105_pheno_full$education = sapply(GSE125105_pheno_full$Participant, function(x){
  Data = GSE128235_pheno_selected[GSE128235_pheno_selected$Sample_relation == x, "education level"]
  Data = unlist(Data)
  if (length(Data) < 1){
    Data = NA
  }
  return(Data)
})
GSE125105_pheno_full$education = as.numeric(GSE125105_pheno_full$education)

GSE125105_pheno_full$myocard_inf = sapply(GSE125105_pheno_full$Participant, function(x){
  Data = GSE128235_pheno_selected[GSE128235_pheno_selected$Sample_relation == x, "history of myocardial infarction"]
  Data = unlist(Data)
  if (length(Data) < 1){
    Data = NA
  }
  return(Data)
})
GSE125105_pheno_full$myocard_inf = ifelse(GSE125105_pheno_full$myocard_inf == "NA", NA, GSE125105_pheno_full$myocard_inf)
GSE125105_pheno_full$myocard_inf = ifelse(GSE125105_pheno_full$myocard_inf == "1", "YES", "NO")
GSE125105_pheno_full$myocard_inf = factor(GSE125105_pheno_full$myocard_inf, levels = c("YES", "NO"))

GSE125105_pheno_full$SNP_PC_1 = sapply(GSE125105_pheno_full$Participant, function(x){
  Data = GSE128235_pheno_selected[GSE128235_pheno_selected$Sample_relation == x, "snp based principle component 1"]
  Data = unlist(Data)
  if (length(Data) < 1){
    Data = NA
  }
  return(Data)
})
GSE125105_pheno_full$SNP_PC_1  = as.numeric(GSE125105_pheno_full$SNP_PC_1)

GSE125105_pheno_full$SNP_PC_2 = sapply(GSE125105_pheno_full$Participant, function(x){
  Data = GSE128235_pheno_selected[GSE128235_pheno_selected$Sample_relation == x, "snp based principle component 2"]
  Data = unlist(Data)
  if (length(Data) < 1){
    Data = NA
  }
  return(Data)
})
GSE125105_pheno_full$SNP_PC_2 = as.numeric(GSE125105_pheno_full$SNP_PC_2)

# BMI phenotype data
GSE128235_bmi = smart_fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/GSE128235_BMI.txt") # This file has been provided by MPIP with a request # Replace with an appropriate path
GSE128235_bmi$ID_corrected = stri_replace_all_fixed(GSE128235_bmi$id, pattern = "sample", replacement = "")
GSE128235_pheno_selected$ID_corrected = stri_replace_all_fixed(GSE128235_pheno_selected$Sample_title, pattern = "genomic DNA from sample ", replacement = "")
GSE128235_pheno_selected = inner_join(GSE128235_pheno_selected, GSE128235_bmi, by = "ID_corrected")
GSE128235_pheno_selected$BMI_fixed = NA
# Average or available values were used for BMI
for (i in 1:nrow(GSE128235_pheno_selected)){
  BMIs = c(GSE128235_pheno_selected$bmi_00[i], GSE128235_pheno_selected$bmi_01[i])
  if (all(!is.na(BMIs))){
    GSE128235_pheno_selected$BMI_fixed[i] = mean(BMIs)
  } else if (all(is.na(BMIs))){
    GSE128235_pheno_selected$BMI_fixed[i] = NA
  } else {
    BMIs = BMIs[!is.na(BMIs)]
    GSE128235_pheno_selected$BMI_fixed[i] = BMIs
  }
  rm(BMIs)
}
GSE125105_pheno_full$BMI_fixed = sapply(GSE125105_pheno_full$Participant, function(x){
  Data = GSE128235_pheno_selected[GSE128235_pheno_selected$Sample_relation == x, "BMI_fixed"]
  Data = unlist(Data)
  if (length(Data) < 1){
    Data = NA
  }
  return(Data)
})


################### Differential DNA methylation analysis ###################

# Reading Illumina 450K annotation file
annotation = smart_fread("/home/aleksandr/Desktop/WORK/depression analysis//Annot_Illumina_450_K.csv") # Replace with another path if needed

# 1) PSY Methylation 197 normal, 24 depressed
Methyl_depression_broad_PSY_SCR = test_CpGs_limma_generalized(PREFIX_CpGs = rownames(M_values),
                                                              PREFIX_Mvals = M_values,
                                                              PREFIX_Pheno_df = Phenotypes_SCR,
                                                              PREFIX_contrast_col_number = 15,
                                                              PREFIX_participants_col_number = 6,
                                                              PREFIX_contrast_vector = c("High depr. risk", "Low depr. risk"),
                                                              PREFIX_Study = "PSY-SCR-Dawba:",
                                                              PREFIX_annotation = annotation,
                                                              PREFIX_model.covariates = c("Gender", "Age", "BMI", "Study"),
                                                              PREFIX_plots_folder = "Methyl_depression_broad",
                                                              PREFIX_plots_str = "PSY_SCR_Dawba_depr_",
                                                              PREFIX_use_Combat = FALSE,
                                                              PREFIX_Batch_col_number = NULL,
                                                              PREFIX_Remove_NA_predictors = TRUE,
                                                              PREFIX_caclulate_cumul_positions = TRUE,
                                                              PREFIX_Manhattan_plot_name = "PSY_SCR_Dawba_depr_manhattan.png", 
                                                              PREFIX_Log_FC_threshold = 0.2)
write.csv(Methyl_depression_broad_PSY_SCR, "Methyl_depression_broad_PSY_SCR.csv")
# 380756 CpGs for analysis; 0 Significant after FDR

Methyl_depression_broad_PSY_SCR_significant = Methyl_depression_broad_PSY_SCR[Methyl_depression_broad_PSY_SCR$P.Value < 0.05,]
write.csv(Methyl_depression_broad_PSY_SCR_significant, "Methyl_depression_broad_PSY_SCR_significant.csv")
# 13251 probes are nominally significant

Methyl_depression_broad_PSY_SCR_significant_LF = Methyl_depression_broad_PSY_SCR_significant[abs(Methyl_depression_broad_PSY_SCR_significant$logFC) > 0.2,]
write.csv(Methyl_depression_broad_PSY_SCR_significant_LF, "Methyl_depression_broad_PSY_SCR_significant_LF.csv")
# 10213 probes show difference in methylation 15% or higher and nominally significant

# 2) GSE125105 MPIP 167 non-depressed 324 depression
Methyl_depression_broad_GSE125105 = test_CpGs_limma_generalized(PREFIX_CpGs = rownames(GSE125105_M_val_corrected),
                                                                PREFIX_Mvals = GSE125105_M_val_corrected,
                                                                PREFIX_Pheno_df = GSE125105_pheno_full,
                                                                PREFIX_contrast_col_number = 15,
                                                                PREFIX_participants_col_number = 1,
                                                                PREFIX_contrast_vector = c("Depression", "Non-depressed"),
                                                                PREFIX_Study = "GSE125105 Depression:",
                                                                PREFIX_annotation = annotation,
                                                                PREFIX_model.covariates = c("Sex","age", "BMI_fixed", "SNP_PC_1", 
                                                                                            "SNP_PC_2"),
                                                                PREFIX_plots_folder = "Methyl_depression_broad",
                                                                PREFIX_plots_str = "GSE125105_depression_",
                                                                PREFIX_use_Combat = FALSE,
                                                                PREFIX_Batch_col_number = NULL,
                                                                PREFIX_Remove_NA_predictors = TRUE,
                                                                PREFIX_caclulate_cumul_positions = TRUE,
                                                                PREFIX_Manhattan_plot_name = "GSE125105_Depression_manhattan.png", 
                                                                PREFIX_Log_FC_threshold = 0.2)
write.csv(Methyl_depression_broad_GSE125105, "Methyl_depression_broad_GSE125105.csv")
# 378256 probes are available for analysis

Methyl_depression_broad_GSE125105_significant = Methyl_depression_broad_GSE125105[Methyl_depression_broad_GSE125105$P.Value < 0.05, ]
write.csv(Methyl_depression_broad_GSE125105_significant, "Methyl_depression_broad_GSE125105_significant.csv")
# 24767 probes are nominally significant

Methyl_depression_broad_GSE125105_significant_LF = Methyl_depression_broad_GSE125105_significant[abs(Methyl_depression_broad_GSE125105_significant$logFC) > 0.2,]
write.csv(Methyl_depression_broad_GSE125105_significant_LF, "Methyl_depression_broad_GSE125105_significant_LF.csv")
# 448 probes show difference in methylation 15% or higher and nominally significant

# 3) GSE72680, GRADY, Composite_depression_NA_full
# Analyzed Normal 161 Depressed 194 Total: 355
Methyl_depression_broad_E_GEOD_72680 = test_CpGs_limma_generalized(PREFIX_CpGs = rownames(E_GEOD_72680_Mval),
                                                                   PREFIX_Mvals = E_GEOD_72680_Mval,
                                                                   PREFIX_Pheno_df = E_GEOD_72680_Phenotypes,
                                                                   PREFIX_contrast_col_number = 69,
                                                                   PREFIX_participants_col_number = 48,
                                                                   PREFIX_contrast_vector = c("Depressed", "Normal"),
                                                                   PREFIX_Study = "E-GEOD-72680 Depression Composite:",
                                                                   PREFIX_annotation = annotation,
                                                                   PREFIX_model.covariates = c("Sex","age", "race_ethnicity", "body_mass_index", 
                                                                                               "cd8_t_cells" ,"cd4_t_cells", "natural_killer_cells", "b_cells", "monocytes", "granulocytes"), # cell-type correction
                                                                   PREFIX_plots_folder = "Methyl_depression_broad",
                                                                   PREFIX_plots_str = "E_GEOD_72680_Depr_treatment_",
                                                                   PREFIX_use_Combat = TRUE, # The adjustment for bead batch
                                                                   PREFIX_Batch_col_number = 60,
                                                                   PREFIX_Remove_NA_predictors = TRUE,
                                                                   PREFIX_caclulate_cumul_positions = TRUE,
                                                                   PREFIX_Manhattan_plot_name = "E_GEOD_72680_manhattan_",
                                                                   PREFIX_Log_FC_threshold = 0.2)
# 285900 probes were available for analysis
write.csv(Methyl_depression_broad_E_GEOD_72680, "Methyl_depression_broad_E_GEOD_72680.csv")
# One CpG passed through multiple comparisons cg20263853

Methyl_depression_broad_E_GEOD_72680_significant = Methyl_depression_broad_E_GEOD_72680[Methyl_depression_broad_E_GEOD_72680$P.Value < 0.05, ]
write.csv(Methyl_depression_broad_E_GEOD_72680_significant, "Methyl_depression_broad_E_GEOD_72680_significant.csv")
# 16887 probes are nominally significant

Methyl_depression_broad_E_GEOD_72680_significant_LF = Methyl_depression_broad_E_GEOD_72680_significant[abs(Methyl_depression_broad_E_GEOD_72680_significant$logFC) > 0.2,]
write.csv(Methyl_depression_broad_E_GEOD_72680_significant_LF, "Methyl_depression_broad_E_GEOD_72680_significant_LF.csv")
# 4344 probes show difference in methylation 15% or higher and nominally significant


################### Making Venn diagrams ###################
# Covered CpGs
Methyl_depression_broad_covered_CpGs = list(
  "PSY-SCR" = rownames(M_values),
  "GSE125105" = rownames(GSE125105_M_val_corrected),
  "E-GEOD-72680" = rownames(E_GEOD_72680_Mval)
)
make_Venn_digram_list(named_list = Methyl_depression_broad_covered_CpGs, palette = 1, plot_full_path = "Venn_CpGs_Tested.pdf")

# Significant CpGs
Methyl_depression_broad_significant_CpGs = list(
  "PSY-SCR Significant" = Methyl_depression_broad_PSY_SCR_significant$CpG,
  "GSE125105 Significant" = Methyl_depression_broad_GSE125105_significant$CpG,
  "E-GEOD-72680 Significant" = Methyl_depression_broad_E_GEOD_72680_significant$CpG
)
make_Venn_digram_list(named_list = Methyl_depression_broad_significant_CpGs, palette = 2, plot_full_path = "Venn_CpGs_Significant.pdf")

# Overlap Significant CpGs (covered in all cohorts)
Covered_CpGs_all = Reduce(intersect, Methyl_depression_broad_covered_CpGs)
Methyl_depression_broad_significant_CpGs = list(
  "PSY-SCR Significant" = Methyl_depression_broad_PSY_SCR_significant$CpG[Methyl_depression_broad_PSY_SCR_significant$CpG %in% Covered_CpGs_all],
  "GSE125105 Significant" = Methyl_depression_broad_GSE125105_significant$CpG[Methyl_depression_broad_GSE125105_significant$CpG %in% Covered_CpGs_all],
  "E-GEOD-72680 Significant" = Methyl_depression_broad_E_GEOD_72680_significant$CpG[Methyl_depression_broad_E_GEOD_72680_significant$CpG %in% Covered_CpGs_all]
)
make_Venn_digram_list(named_list = Methyl_depression_broad_significant_CpGs, palette = 7, plot_full_path = "Venn_CpGs_Significant_Covered_All.pdf")


################### Making Step-based plots ###################

# Significant probes vs P-value
# Note: Temporary variables are labelled with the prefix PRF_ to enable easy cleanup after execution
Dataset_Step_P_val = list()
PRF_combined_methyl = rbind(Methyl_depression_broad_PSY_SCR_significant, Methyl_depression_broad_GSE125105_significant, Methyl_depression_broad_E_GEOD_72680_significant)
PRF_combined_methyl = PRF_combined_methyl[PRF_combined_methyl$CpG %in% Covered_CpGs_all,]
PRF_combined_methyl$P.Value.log = -log10(PRF_combined_methyl$P.Value)
PRF_Max_P_val_log = max(PRF_combined_methyl$P.Value.log)
PRF_Max_P_val_log = round(PRF_Max_P_val_log, digits = 1)
PRF_Min_P_val_log = min(PRF_combined_methyl$P.Value.log)
PRF_Min_P_val_log = round(PRF_Min_P_val_log, digits = 1)
PRF_P_val_seq = seq(from = PRF_Min_P_val_log, to = PRF_Max_P_val_log, by = 0.05)

for (i in 1:length(PRF_P_val_seq)){
  PRF_current_P_val_log = PRF_P_val_seq[i]
  PRF_current_df = PRF_combined_methyl[PRF_combined_methyl$P.Value.log >= PRF_current_P_val_log,]
  PRF_current_df_PSY_SCR = PRF_current_df[PRF_current_df$Contrast == "PSY-SCR-Dawba:High depr. risk/Low depr. risk",]
  PRF_current_df_GSE125105 = PRF_current_df[PRF_current_df$Contrast == "GSE125105 Depression:Depression/Non-depressed",]
  PRF_current_df_GSE72680 = PRF_current_df[PRF_current_df$Contrast == "E-GEOD-72680 Depression Composite:Depressed/Normal",]
  
  # Counting significant CpGs
  PRF_sign_PSY = nrow(PRF_current_df_PSY_SCR)
  PRF_sign_GSE125105 = nrow(PRF_current_df_GSE125105)
  PRF_sign_GSE72680 = nrow(PRF_current_df_GSE72680)
  PRF_sign_sum = nrow(PRF_current_df)
  PRF_sign_unqiue = length(unique(PRF_current_df$CpG))
  
  # Getting overlapping CpGs
  PRF_sign_PSY_CpGs = PRF_current_df_PSY_SCR$CpG
  PRF_sign_GSE125105_CpGs = PRF_current_df_GSE125105$CpG
  PRF_sign_GSE72680_CpGs = PRF_current_df_GSE72680$CpG
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs))
  
  PRF_intersect_PSY_GSE125105 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs)
  PRF_intersect_PSY_GSE125105 = PRF_intersect_PSY_GSE125105[PRF_intersect_PSY_GSE125105 %!in% PRF_intersect_all]
  
  PRF_intersect_PSY_GSE72680 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_PSY_GSE72680 = PRF_intersect_PSY_GSE72680[PRF_intersect_PSY_GSE72680 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE125105_GSE72680 = intersect(PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_GSE125105_GSE72680 = PRF_intersect_GSE125105_GSE72680[PRF_intersect_GSE125105_GSE72680 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_PSY_GSE125105_count = length(PRF_intersect_PSY_GSE125105)
  PRF_intersect_PSY_GSE72680_count = length(PRF_intersect_PSY_GSE72680)
  PRF_intersect_GSE125105_GSE72680_count = length(PRF_intersect_GSE125105_GSE72680)
  
  # Percent_total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_all = round((PRF_intersect_PSY_GSE125105_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_all = round((PRF_intersect_PSY_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_all = round((PRF_intersect_GSE125105_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  
  #Percent_intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_PSY_GSE125105, PRF_intersect_PSY_GSE72680, PRF_intersect_GSE125105_GSE72680)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_intersect = round((PRF_intersect_PSY_GSE125105_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_intersect = round((PRF_intersect_PSY_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_intersect = round((PRF_intersect_GSE125105_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  
  PRF_output_df = data.frame(
    PRF_sign_PSY = PRF_sign_PSY,
    PRF_sign_GSE125105 = PRF_sign_GSE125105,
    PRF_sign_GSE72680 = PRF_sign_GSE72680,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_PSY_GSE125105_count = PRF_intersect_PSY_GSE125105_count,
    PRF_intersect_PSY_GSE72680_count = PRF_intersect_PSY_GSE72680_count,
    PRF_intersect_GSE125105_GSE72680_count = PRF_intersect_GSE125105_GSE72680_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_PSY_GSE125105_count_percent_all = PRF_intersect_PSY_GSE125105_count_percent_all,
    PRF_intersect_PSY_GSE72680_count_percent_all = PRF_intersect_PSY_GSE72680_count_percent_all,
    PRF_intersect_GSE125105_GSE72680_count_percent_all = PRF_intersect_GSE125105_GSE72680_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_PSY_GSE125105_count_percent_intersect = PRF_intersect_PSY_GSE125105_count_percent_intersect,
    PRF_intersect_PSY_GSE72680_count_percent_intersect = PRF_intersect_PSY_GSE72680_count_percent_intersect,
    PRF_intersect_GSE125105_GSE72680_count_percent_intersect = PRF_intersect_GSE125105_GSE72680_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(P.Val.10.log = PRF_current_P_val_log, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_P_val[[i]] = PRF_output_df
}

rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])
Dataset_Step_P_val = list_to_df(Dataset_Step_P_val)
Dataset_Step_P_val_2 = Dataset_Step_P_val[stri_detect_fixed(Dataset_Step_P_val$Score_name, pattern = "_sign_"),]
Dataset_Step_P_val_2$Score_name = factor(Dataset_Step_P_val_2$Score_name, levels = unique(Dataset_Step_P_val_2$Score_name), 
                                         labels = c("PSY-SCR",
                                                    "GSE125105",
                                                    "GSE72680",
                                                    "Total",
                                                    "Total unique"))
plot = ggplot(data = Dataset_Step_P_val_2, aes(x = P.Val.10.log, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 50, expand = c(0, 0), limits = c(0,47000)) +
  labs(y = "Probe count", x = "-log10 p-value", col = "Set", shape = "Set", linetype = "Set") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", 
             color = "red", size=0.5) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 13, colour = "black"),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
plot
pdf("Total_vs_P_val.pdf", width = 10, height = 10)
plot
dev.off()

# Significant probes vs LF
# Note: Temporary variables are labelled with the prefix PRF_ to enable easy cleanup after execution
Dataset_Step_Lf = list()
PRF_combined_methyl = rbind(Methyl_depression_broad_PSY_SCR_significant, Methyl_depression_broad_GSE125105_significant, Methyl_depression_broad_E_GEOD_72680_significant)
PRF_combined_methyl = PRF_combined_methyl[PRF_combined_methyl$CpG %in% Covered_CpGs_all,]
PRF_Max_LF = max(abs(PRF_combined_methyl$logFC))
PRF_Max_LF = round(PRF_Max_LF, digits = 1)
PRF_Max_LF_seq = seq(from = 0, to = PRF_Max_LF, by = 0.02)

for (i in 1:length(PRF_Max_LF_seq)){
  PRF_current_LF = PRF_Max_LF_seq[i]
  PRF_current_df = PRF_combined_methyl[abs(PRF_combined_methyl$logFC) >= PRF_current_LF,]
  PRF_current_df_PSY_SCR = PRF_current_df[PRF_current_df$Contrast == "PSY-SCR-Dawba:High depr. risk/Low depr. risk",]
  PRF_current_df_GSE125105 = PRF_current_df[PRF_current_df$Contrast == "GSE125105 Depression:Depression/Non-depressed",]
  PRF_current_df_GSE72680 = PRF_current_df[PRF_current_df$Contrast == "E-GEOD-72680 Depression Composite:Depressed/Normal",]
  
  # Counting signif CpGs
  PRF_sign_PSY = nrow(PRF_current_df_PSY_SCR)
  PRF_sign_GSE125105 = nrow(PRF_current_df_GSE125105)
  PRF_sign_GSE72680 = nrow(PRF_current_df_GSE72680)
  PRF_sign_sum = nrow(PRF_current_df)
  PRF_sign_unqiue = length(unique(PRF_current_df$CpG))
  
  # Getting overlapping CpGs
  PRF_sign_PSY_CpGs = PRF_current_df_PSY_SCR$CpG
  PRF_sign_GSE125105_CpGs = PRF_current_df_GSE125105$CpG
  PRF_sign_GSE72680_CpGs = PRF_current_df_GSE72680$CpG
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs))
  
  PRF_intersect_PSY_GSE125105 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs)
  PRF_intersect_PSY_GSE125105 = PRF_intersect_PSY_GSE125105[PRF_intersect_PSY_GSE125105 %!in% PRF_intersect_all]
  
  PRF_intersect_PSY_GSE72680 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_PSY_GSE72680 = PRF_intersect_PSY_GSE72680[PRF_intersect_PSY_GSE72680 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE125105_GSE72680 = intersect(PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_GSE125105_GSE72680 = PRF_intersect_GSE125105_GSE72680[PRF_intersect_GSE125105_GSE72680 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_PSY_GSE125105_count = length(PRF_intersect_PSY_GSE125105)
  PRF_intersect_PSY_GSE72680_count = length(PRF_intersect_PSY_GSE72680)
  PRF_intersect_GSE125105_GSE72680_count = length(PRF_intersect_GSE125105_GSE72680)
  
  # Percent_total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_all = round((PRF_intersect_PSY_GSE125105_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_all = round((PRF_intersect_PSY_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_all = round((PRF_intersect_GSE125105_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  
  #Percent_intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_PSY_GSE125105, PRF_intersect_PSY_GSE72680, PRF_intersect_GSE125105_GSE72680)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_intersect = round((PRF_intersect_PSY_GSE125105_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_intersect = round((PRF_intersect_PSY_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_intersect = round((PRF_intersect_GSE125105_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  
  PRF_output_df = data.frame(
    PRF_sign_PSY = PRF_sign_PSY,
    PRF_sign_GSE125105 = PRF_sign_GSE125105,
    PRF_sign_GSE72680 = PRF_sign_GSE72680,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_PSY_GSE125105_count = PRF_intersect_PSY_GSE125105_count,
    PRF_intersect_PSY_GSE72680_count = PRF_intersect_PSY_GSE72680_count,
    PRF_intersect_GSE125105_GSE72680_count = PRF_intersect_GSE125105_GSE72680_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_PSY_GSE125105_count_percent_all = PRF_intersect_PSY_GSE125105_count_percent_all,
    PRF_intersect_PSY_GSE72680_count_percent_all = PRF_intersect_PSY_GSE72680_count_percent_all,
    PRF_intersect_GSE125105_GSE72680_count_percent_all = PRF_intersect_GSE125105_GSE72680_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_PSY_GSE125105_count_percent_intersect = PRF_intersect_PSY_GSE125105_count_percent_intersect,
    PRF_intersect_PSY_GSE72680_count_percent_intersect = PRF_intersect_PSY_GSE72680_count_percent_intersect,
    PRF_intersect_GSE125105_GSE72680_count_percent_intersect = PRF_intersect_GSE125105_GSE72680_count_percent_intersect
  )
  PRF_output_df = as.data.frame(t(PRF_output_df))
  PRF_output_df = data.frame(LF = PRF_current_LF, Score_name = rownames(PRF_output_df), Score = PRF_output_df$V1)
  Dataset_Step_Lf[[i]] = PRF_output_df
}

rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])
Dataset_Step_Lf = list_to_df(Dataset_Step_Lf)
Dataset_Step_Lf_2 = Dataset_Step_Lf[stri_detect_fixed(Dataset_Step_Lf$Score_name, pattern = "_sign_"),]
Dataset_Step_Lf_2$Score_name = factor(Dataset_Step_Lf_2$Score_name, levels = unique(Dataset_Step_Lf_2$Score_name), 
                                      labels = c("PSY-SCR",
                                                 "GSE125105",
                                                 "GSE72680",
                                                 "Total",
                                                 "Total unique"))
plot = ggplot(data = Dataset_Step_Lf_2, aes(x = LF, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Paired") +
  scale_y_continuous(n.breaks = 50, expand = c(0, 0), limits = c(0,47000)) +
  scale_x_continuous(n.breaks = 15) +
  labs(y = "Probe count", x = "|Log2 fold change|", col = "Set", shape = "Set", linetype = "Set") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 13, colour = "black"),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
plot
pdf("Total_vs_LF.pdf", width = 10, height = 10)
plot
dev.off()

# LogFC Steps
Dataset_Step_Lf = list()
PRF_combined_methyl = rbind(Methyl_depression_broad_PSY_SCR_significant, Methyl_depression_broad_GSE125105_significant, Methyl_depression_broad_E_GEOD_72680_significant)
PRF_combined_methyl = PRF_combined_methyl[PRF_combined_methyl$CpG %in% Covered_CpGs_all,]
PRF_Max_LF = max(abs(PRF_combined_methyl$logFC))
PRF_Max_LF = round(PRF_Max_LF, digits = 1)
PRF_Max_LF_seq = seq(from = 0, to = PRF_Max_LF, by = 0.02)

for (i in 1:length(PRF_Max_LF_seq)){
  PRF_current_LF = PRF_Max_LF_seq[i]
  PRF_current_df = PRF_combined_methyl[abs(PRF_combined_methyl$logFC) >= PRF_current_LF,]
  PRF_current_df_PSY_SCR = PRF_current_df[PRF_current_df$Contrast == "PSY-SCR-Dawba:High depr. risk/Low depr. risk",]
  PRF_current_df_GSE125105 = PRF_current_df[PRF_current_df$Contrast == "GSE125105 Depression:Depression/Non-depressed",]
  PRF_current_df_GSE72680 = PRF_current_df[PRF_current_df$Contrast == "E-GEOD-72680 Depression Composite:Depressed/Normal",]
  
  # Counting significant CpGs
  PRF_sign_PSY = nrow(PRF_current_df_PSY_SCR)
  PRF_sign_GSE125105 = nrow(PRF_current_df_GSE125105)
  PRF_sign_GSE72680 = nrow(PRF_current_df_GSE72680)
  PRF_sign_sum = nrow(PRF_current_df)
  PRF_sign_unqiue = length(unique(PRF_current_df$CpG))
  
  #Getting overlapping CpGs
  PRF_sign_PSY_CpGs = PRF_current_df_PSY_SCR$CpG
  PRF_sign_GSE125105_CpGs = PRF_current_df_GSE125105$CpG
  PRF_sign_GSE72680_CpGs = PRF_current_df_GSE72680$CpG
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs))
  
  PRF_intersect_PSY_GSE125105 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs)
  PRF_intersect_PSY_GSE125105 = PRF_intersect_PSY_GSE125105[PRF_intersect_PSY_GSE125105 %!in% PRF_intersect_all]
  
  PRF_intersect_PSY_GSE72680 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_PSY_GSE72680 = PRF_intersect_PSY_GSE72680[PRF_intersect_PSY_GSE72680 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE125105_GSE72680 = intersect(PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_GSE125105_GSE72680 = PRF_intersect_GSE125105_GSE72680[PRF_intersect_GSE125105_GSE72680 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_PSY_GSE125105_count = length(PRF_intersect_PSY_GSE125105)
  PRF_intersect_PSY_GSE72680_count = length(PRF_intersect_PSY_GSE72680)
  PRF_intersect_GSE125105_GSE72680_count = length(PRF_intersect_GSE125105_GSE72680)
  
  # Percent_total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_all = round((PRF_intersect_PSY_GSE125105_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_all = round((PRF_intersect_PSY_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_all = round((PRF_intersect_GSE125105_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  
  # Percent_intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_PSY_GSE125105, PRF_intersect_PSY_GSE72680, PRF_intersect_GSE125105_GSE72680)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_intersect = round((PRF_intersect_PSY_GSE125105_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_intersect = round((PRF_intersect_PSY_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_intersect = round((PRF_intersect_GSE125105_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  
  PRF_output_df = data.frame(
    PRF_sign_PSY = PRF_sign_PSY,
    PRF_sign_GSE125105 = PRF_sign_GSE125105,
    PRF_sign_GSE72680 = PRF_sign_GSE72680,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_PSY_GSE125105_count = PRF_intersect_PSY_GSE125105_count,
    PRF_intersect_PSY_GSE72680_count = PRF_intersect_PSY_GSE72680_count,
    PRF_intersect_GSE125105_GSE72680_count = PRF_intersect_GSE125105_GSE72680_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_PSY_GSE125105_count_percent_all = PRF_intersect_PSY_GSE125105_count_percent_all,
    PRF_intersect_PSY_GSE72680_count_percent_all = PRF_intersect_PSY_GSE72680_count_percent_all,
    PRF_intersect_GSE125105_GSE72680_count_percent_all = PRF_intersect_GSE125105_GSE72680_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_PSY_GSE125105_count_percent_intersect = PRF_intersect_PSY_GSE125105_count_percent_intersect,
    PRF_intersect_PSY_GSE72680_count_percent_intersect = PRF_intersect_PSY_GSE72680_count_percent_intersect,
    PRF_intersect_GSE125105_GSE72680_count_percent_intersect = PRF_intersect_GSE125105_GSE72680_count_percent_intersect
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
                                                 "PSY-SCR & GSE125105",
                                                 "PSY-SCR & GSE72680",
                                                 "GSE125105 & GSE72680",
                                                 "Intersect in any"))
plot = ggplot(data = Dataset_Step_Lf_2, aes(x = LF, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Dark2") +
  scale_y_continuous(n.breaks = 50, limits = c(NA,6), expand = c(0, 0)) +
  scale_x_continuous(n.breaks = 10, limits = c(NA,0.75)) +
  labs(y = "% of overlapping probes", x = "|Log2 fold change|", col = "Overlap type", shape = "Overlap type", linetype = "Overlap type") +
  theme(
    legend.position = c(0.8, 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 13, colour = "black"),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
plot
pdf("Overlap_vs_LF.pdf", width = 10, height = 10)
plot
dev.off()

# Now it is time for steps in p-value
Dataset_Step_P_val = list()
PRF_combined_methyl = rbind(Methyl_depression_broad_PSY_SCR_significant, Methyl_depression_broad_GSE125105_significant, Methyl_depression_broad_E_GEOD_72680_significant)
PRF_combined_methyl = PRF_combined_methyl[PRF_combined_methyl$CpG %in% Covered_CpGs_all,]
PRF_combined_methyl$P.Value.log = -log10(PRF_combined_methyl$P.Value)
PRF_Max_P_val_log = max(PRF_combined_methyl$P.Value.log)
PRF_Max_P_val_log = round(PRF_Max_P_val_log, digits = 1)
PRF_Min_P_val_log = min(PRF_combined_methyl$P.Value.log)
PRF_Min_P_val_log = round(PRF_Min_P_val_log, digits = 1)
PRF_P_val_seq = seq(from = PRF_Min_P_val_log, to = PRF_Max_P_val_log, by = 0.05)

for (i in 1:length(PRF_P_val_seq)){
  PRF_current_P_val_log = PRF_P_val_seq[i]
  PRF_current_df = PRF_combined_methyl[PRF_combined_methyl$P.Value.log >= PRF_current_P_val_log,]
  PRF_current_df_PSY_SCR = PRF_current_df[PRF_current_df$Contrast == "PSY-SCR-Dawba:High depr. risk/Low depr. risk",]
  PRF_current_df_GSE125105 = PRF_current_df[PRF_current_df$Contrast == "GSE125105 Depression:Depression/Non-depressed",]
  PRF_current_df_GSE72680 = PRF_current_df[PRF_current_df$Contrast == "E-GEOD-72680 Depression Composite:Depressed/Normal",]
  
  # Counting significant CpGs
  PRF_sign_PSY = nrow(PRF_current_df_PSY_SCR)
  PRF_sign_GSE125105 = nrow(PRF_current_df_GSE125105)
  PRF_sign_GSE72680 = nrow(PRF_current_df_GSE72680)
  PRF_sign_sum = nrow(PRF_current_df)
  PRF_sign_unqiue = length(unique(PRF_current_df$CpG))
  
  # Getting overlapping CpGs
  PRF_sign_PSY_CpGs = PRF_current_df_PSY_SCR$CpG
  PRF_sign_GSE125105_CpGs = PRF_current_df_GSE125105$CpG
  PRF_sign_GSE72680_CpGs = PRF_current_df_GSE72680$CpG
  PRF_intersect_all = Reduce(intersect, list(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs))
  
  PRF_intersect_PSY_GSE125105 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE125105_CpGs)
  PRF_intersect_PSY_GSE125105 = PRF_intersect_PSY_GSE125105[PRF_intersect_PSY_GSE125105 %!in% PRF_intersect_all]
  
  PRF_intersect_PSY_GSE72680 = intersect(PRF_sign_PSY_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_PSY_GSE72680 = PRF_intersect_PSY_GSE72680[PRF_intersect_PSY_GSE72680 %!in% PRF_intersect_all]
  
  PRF_intersect_GSE125105_GSE72680 = intersect(PRF_sign_GSE125105_CpGs, PRF_sign_GSE72680_CpGs)
  PRF_intersect_GSE125105_GSE72680 = PRF_intersect_GSE125105_GSE72680[PRF_intersect_GSE125105_GSE72680 %!in% PRF_intersect_all]
  
  # Performing counts
  PRF_intersect_all_count = length(PRF_intersect_all)
  PRF_intersect_PSY_GSE125105_count = length(PRF_intersect_PSY_GSE125105)
  PRF_intersect_PSY_GSE72680_count = length(PRF_intersect_PSY_GSE72680)
  PRF_intersect_GSE125105_GSE72680_count = length(PRF_intersect_GSE125105_GSE72680)
  
  # Percent_total
  PRF_intersect_all_count_percent_all = round((PRF_intersect_all_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_all = round((PRF_intersect_PSY_GSE125105_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_all = round((PRF_intersect_PSY_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_all = round((PRF_intersect_GSE125105_GSE72680_count/PRF_sign_unqiue)*100, digits = 2)
  
  # Percent_intersect
  PRF_total_intersect = c(PRF_intersect_all, PRF_intersect_PSY_GSE125105, PRF_intersect_PSY_GSE72680, PRF_intersect_GSE125105_GSE72680)
  PRF_total_intersect = unique(PRF_total_intersect)
  PRF_total_intersect_count = length(PRF_total_intersect)
  PRF_total_intersect_count_percent_all = round((PRF_total_intersect_count/PRF_sign_unqiue)*100, digits = 2)
  
  PRF_intersect_all_count_percent_intersect = round((PRF_intersect_all_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE125105_count_percent_intersect = round((PRF_intersect_PSY_GSE125105_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_PSY_GSE72680_count_percent_intersect = round((PRF_intersect_PSY_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_intersect_GSE125105_GSE72680_count_percent_intersect = round((PRF_intersect_GSE125105_GSE72680_count/PRF_total_intersect_count)*100, digits = 2)
  PRF_output_df = data.frame(
    PRF_sign_PSY = PRF_sign_PSY,
    PRF_sign_GSE125105 = PRF_sign_GSE125105,
    PRF_sign_GSE72680 = PRF_sign_GSE72680,
    PRF_sign_sum = PRF_sign_sum,
    PRF_sign_unqiue = PRF_sign_unqiue,
    PRF_intersect_all_count = PRF_intersect_all_count,
    PRF_intersect_PSY_GSE125105_count = PRF_intersect_PSY_GSE125105_count,
    PRF_intersect_PSY_GSE72680_count = PRF_intersect_PSY_GSE72680_count,
    PRF_intersect_GSE125105_GSE72680_count = PRF_intersect_GSE125105_GSE72680_count,
    PRF_intersect_all_count_percent_all = PRF_intersect_all_count_percent_all,
    PRF_intersect_PSY_GSE125105_count_percent_all = PRF_intersect_PSY_GSE125105_count_percent_all,
    PRF_intersect_PSY_GSE72680_count_percent_all = PRF_intersect_PSY_GSE72680_count_percent_all,
    PRF_intersect_GSE125105_GSE72680_count_percent_all = PRF_intersect_GSE125105_GSE72680_count_percent_all,
    PRF_total_intersect_count_percent_all = PRF_total_intersect_count_percent_all,
    PRF_intersect_all_count_percent_intersect = PRF_intersect_all_count_percent_intersect,
    PRF_intersect_PSY_GSE125105_count_percent_intersect = PRF_intersect_PSY_GSE125105_count_percent_intersect,
    PRF_intersect_PSY_GSE72680_count_percent_intersect = PRF_intersect_PSY_GSE72680_count_percent_intersect,
    PRF_intersect_GSE125105_GSE72680_count_percent_intersect = PRF_intersect_GSE125105_GSE72680_count_percent_intersect
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
                                                    "PSY-SCR & GSE125105",
                                                    "PSY-SCR & GSE72680",
                                                    "GSE125105 & GSE72680",
                                                    "Intersect in any"))
plot = ggplot(data = Dataset_Step_P_val_2, aes(x = P.Val.10.log, y = Score, col = Score_name)) +
  geom_point(aes(shape = Score_name)) +
  geom_line(aes(linetype = Score_name)) + 
  scale_color_brewer(palette="Set2") +
  scale_x_continuous(n.breaks = 10, limits = c(NA,4)) +
  scale_y_continuous(n.breaks = 50, limits = c(NA,6), expand = c(0, 0)) +
  labs(y = "% of overlapping probes", x = "-log10 p-value", col = "Overlap type", shape = "Overlap type", linetype = "Overlap type") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", 
             color = "red", size=0.5) +
  theme(
    legend.position = c(0.8, 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 13, colour = "black"),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
plot
pdf("Overlap_vs_P_val.pdf", width = 10, height = 10)
plot
dev.off()


################### Making Chromosome maps ###################

# Obtaining cytoBand track from UCSC
mySession = browserSession()
genome(mySession) = "hg19"
Full_cytoband_track = list()

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_CytoBand_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_Cytoband_track = getTable(mySession, table = "cytoBand", range =  PRF_CytoBand_grange)
  Full_cytoband_track[[i]] = PRF_Cytoband_track
}

Full_cytoband_track = do.call(rbind, Full_cytoband_track) # last p and first Q are centromers (acen), the end of the last p is the coordinate of a centomere
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Preparing Chromosome table for maps
Chromosome_coord_table_prepared = Chromosome_coord_table
Chromosome_coord_table_prepared$Chromosome = paste0("chr", Chromosome_coord_table_prepared$Chromosome)
Chromosome_coord_table_prepared = Chromosome_coord_table_prepared[,1:2]
Chromosome_coord_table_prepared$start = 1
Centromeres = Full_cytoband_track[Full_cytoband_track$gieStain == "acen",]
Centromere_Starts = Centromeres$chromStart[seq(from = 1, by = 2, to = length(Centromeres$chromStart))]
Centromere_Ends = Centromeres$chromEnd[seq(from = 2, by = 2, to = length(Centromeres$chromEnd))]
Centromeres = Centromeres$chromEnd
Centromeres = Centromeres[seq(from = 1, by = 2, to = length(Centromeres))]
Chromosome_coord_table_prepared = cbind(Chromosome_coord_table_prepared, Centromeres)
colnames(Chromosome_coord_table_prepared) = c("Chromosome", "End", "Start", "Centromere")
Chromosome_coord_table_prepared = Chromosome_coord_table_prepared[c("Chromosome", "Start", "End", "Centromere")]

# Making chromosome df
Chromosome_map_Rideogram = Chromosome_coord_table_prepared
Chromosome_map_Rideogram$Start = 0
Chromosome_map_Rideogram$Centromere = NULL
Chromosome_map_Rideogram$Centromere_St = Centromere_Starts
Chromosome_map_Rideogram$Centromere_End = Centromere_Ends
colnames(Chromosome_map_Rideogram) = c("Chr","Start","End","CE_start","CE_end")

Combined_methyl_df = rbind(Methyl_depression_broad_PSY_SCR,
                           Methyl_depression_broad_GSE125105,
                           Methyl_depression_broad_E_GEOD_72680)
Combined_methyl_df_significant = Combined_methyl_df[Combined_methyl_df$P.Value < 0.05,] # 54905 significant hits
Combined_methyl_df_significant_unique = distinct(Combined_methyl_df_significant, CpG, .keep_all = TRUE) # 52365 unique CpGs

Stats_combined_methyl_df_significant = as.data.frame(table(Combined_methyl_df_significant$CpG))
colnames(Stats_combined_methyl_df_significant) = c("CpG", "Count")
Stats_combined_methyl_df_significant = arrange(Stats_combined_methyl_df_significant, -Count)
Stats_combined_methyl_df_significant_cross_valid = Stats_combined_methyl_df_significant[Stats_combined_methyl_df_significant$Count > 1,]
# 2491 CpGs are valid in more than 1 cohort

# Making heatmap
# Small helper function for the loop
upd_count_gene_helper = function(gene_vector){
  gene_vector = unlist(stri_split_fixed(gene_vector, pattern = ";"))
  gene_vector = unique(gene_vector)
  gene_vector = gene_vector[gene_vector != ""]
  gene_vector = gene_vector[!is.na(gene_vector)]
  if (length(gene_vector) >= 1){
    gene_vector_upd = check_gene_symbol_NIH(PRF_gene_symbols = gene_vector, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded, PRF_replace_NA_with_old = TRUE)
    gene_vector_upd = gene_vector_upd$Suggested.Symbol
    gene_vector_upd_count = length(gene_vector_upd)
    gene_vector_count = length(gene_vector)
    gene_vector = paste0(gene_vector, collapse = ";")
    gene_vector_upd = paste0(gene_vector_upd, collapse = ";")
  } else {
    gene_vector = NA
    gene_vector_upd = NA
    gene_vector_count = 0
    gene_vector_upd_count = 0
  }
  Output = list(
    gene_vector = gene_vector,
    gene_vector_upd = gene_vector_upd,
    gene_vector_count = gene_vector_count,
    gene_vector_upd_count = gene_vector_upd_count
  )
  return(Output)
}

Chrom_heatmap_df_methyl = list()

# MAIN loop to make chromosome maps for DNA methylation
# Outer loop for chromosome
# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
# Note: the loop takes several minuted to run
for (i in 1:nrow(Chromosome_map_Rideogram)){
  PRF_curr_chrom = Chromosome_map_Rideogram$Chr[i]
  writeLines(PRF_curr_chrom)
  PRF_Curr_Methyl_df_annot = annotation[annotation$chr == PRF_curr_chrom,]
  PRF_Curr_Methyl_df_signif = Combined_methyl_df_significant_unique[Combined_methyl_df_significant_unique$Chromosome == PRF_curr_chrom,]
  PRF_curr_chrom_map = Chromosome_map_Rideogram[i,]
  PRF_curr_chrom_vector = seq(from= PRF_curr_chrom_map$Start, by = 1000000, to = PRF_curr_chrom_map$End)
  
  if (PRF_curr_chrom_vector[length(PRF_curr_chrom_vector)] < PRF_curr_chrom_map$End){
    PRF_curr_chrom_vector = c(PRF_curr_chrom_vector, PRF_curr_chrom_map$End)
  }
  
  PRF_Curr_Chrom_Table_Heatmap = list()
  
  # Inner loop for chromosomal intervals
  for (index in 1:(length(PRF_curr_chrom_vector)-1)){
    
    if (index == 1){
      # For the first row
      PRF_Start = PRF_curr_chrom_vector[1]
      PRF_End = PRF_curr_chrom_vector[2]
    } else {
      # For the second and further rows
      PRF_Start = PRF_curr_chrom_vector[index] + 1
      PRF_End = PRF_curr_chrom_vector[index + 1]
    }
    
    RPF_Interval = PRF_Start:PRF_End
    PRF_Curr_Methyl_df_annot_interval = PRF_Curr_Methyl_df_annot[PRF_Curr_Methyl_df_annot$pos %in% RPF_Interval,]
    PRF_Curr_Methyl_df_signif_interval = PRF_Curr_Methyl_df_signif[PRF_Curr_Methyl_df_signif$Position_Hg19 %in% RPF_Interval,]
    
    if (nrow(PRF_Curr_Methyl_df_annot_interval) < 1){
      
      # No CpGs at all in the current interval
      PRF_Illum_Probes_names = NA
      PRF_Signif_CpGs_names = NA
      PRF_Cross_valid_CpGs_names = NA
      PRF_Illum_Probes = 0
      PRF_Signif_CpGs = 0 
      PRF_Cross_valid_CpGs = 0
      
      PRF_Total_CpG_associated_genes = NA
      PRF_CpG_associated_genes_signif = NA
      PRF_CpG_associated_genes_cross_valid = NA
      
      PRF_Total_CpG_associated_genes_upd = NA
      PRF_CpG_associated_genes_signif_upd = NA
      PRF_CpG_associated_genes_cross_valid_upd = NA
      
      PRF_Total_CpG_associated_genes_upd_count = 0
      PRF_CpG_associated_genes_signif_upd_count = 0
      PRF_CpG_associated_genes_cross_valid_upd_count  = 0
      
      PRF_Cross_valid_CpGs_matching_dir_names = NA
      PRF_Cross_valid_CpGs_matching_dir = 0
      PRF_CpG_associated_genes_cross_valid_matching_dir = NA
      PRF_CpG_associated_genes_cross_valid_matching_dir_upd = NA
      PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count = 0
      
    } else {
      
      PRF_Illum_Probes_names = PRF_Curr_Methyl_df_annot_interval$Name
      PRF_Illum_Probes_names = paste0(PRF_Illum_Probes_names, collapse = ";")
      PRF_Illum_Probes = nrow(PRF_Curr_Methyl_df_annot_interval)
      PRF_Total_CpG_associated_genes = PRF_Curr_Methyl_df_annot_interval$UCSC_RefGene_Name
      PRF_Total_CpG_associated_genes_stats = upd_count_gene_helper(PRF_Total_CpG_associated_genes)
      PRF_Total_CpG_associated_genes = PRF_Total_CpG_associated_genes_stats$gene_vector
      PRF_Total_CpG_associated_genes_upd = PRF_Total_CpG_associated_genes_stats$gene_vector_upd
      PRF_Total_CpG_associated_genes_upd_count = PRF_Total_CpG_associated_genes_stats$gene_vector_upd_count
      
      if (nrow(PRF_Curr_Methyl_df_signif_interval) == 0){
        
        # No significant CpGs in the interval
        PRF_Signif_CpGs_names = NA
        PRF_Cross_valid_CpGs_names = NA
        
        PRF_Signif_CpGs = 0 
        PRF_Cross_valid_CpGs = 0
        
        PRF_CpG_associated_genes_signif = NA
        PRF_CpG_associated_genes_cross_valid = NA
        
        PRF_CpG_associated_genes_signif_upd = NA
        PRF_CpG_associated_genes_cross_valid_upd = NA
        
        PRF_CpG_associated_genes_signif_upd_count = 0
        PRF_CpG_associated_genes_cross_valid_upd_count  = 0
        
        PRF_Cross_valid_CpGs_matching_dir_names = NA
        PRF_Cross_valid_CpGs_matching_dir = 0
        PRF_CpG_associated_genes_cross_valid_matching_dir = NA
        PRF_CpG_associated_genes_cross_valid_matching_dir_upd = NA
        PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count = 0
        
      } else {
        
        # Extracting significant CpGs and calculating counts
        PRF_Signif_CpGs_names = PRF_Curr_Methyl_df_signif_interval$CpG
        PRF_Signif_CpGs_names = unique(PRF_Signif_CpGs_names)
        PRF_Cross_valid_CpGs_names = PRF_Signif_CpGs_names[PRF_Signif_CpGs_names %in% Stats_combined_methyl_df_significant_cross_valid$CpG]
        
        PRF_Signif_CpGs = length(PRF_Signif_CpGs_names)
        PRF_Cross_valid_CpGs = length(PRF_Cross_valid_CpGs_names)
        
        PRF_CpG_associated_genes_signif_stats = PRF_Curr_Methyl_df_signif_interval$Gene
        PRF_CpG_associated_genes_signif_stats = upd_count_gene_helper(PRF_CpG_associated_genes_signif_stats)
        PRF_CpG_associated_genes_cross_valid_stats = PRF_Curr_Methyl_df_signif_interval[PRF_Curr_Methyl_df_signif_interval$CpG %in% PRF_Cross_valid_CpGs_names, "Gene"]
        PRF_CpG_associated_genes_cross_valid_stats = upd_count_gene_helper(PRF_CpG_associated_genes_cross_valid_stats)
        
        # Stats for CpGs with a matching direction
        PRF_df_cross_valid_CpGs = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% PRF_Cross_valid_CpGs_names,]
        
        if (nrow(PRF_df_cross_valid_CpGs) < 2){
          PRF_Cross_valid_CpGs_matching_dir_names = NA
          PRF_Cross_valid_CpGs_matching_dir = 0
          PRF_CpG_associated_genes_cross_valid_matching_dir = NA
          PRF_CpG_associated_genes_cross_valid_matching_dir_upd = NA
          PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count = 0
          
        } else {
          
          PRF_matching_dir_index = sapply(unique(PRF_df_cross_valid_CpGs$CpG), function(PRF_x){
            PRF_LFs = PRF_df_cross_valid_CpGs[PRF_df_cross_valid_CpGs$CpG == PRF_x, "logFC"]
            PRF_LFs = PRF_LFs[!is.na(PRF_LFs)]
            PRF_probes_greater_0 = which(PRF_LFs > 0)
            PRF_probes_greater_0 = length(PRF_probes_greater_0)
            PRF_probes_smaller_0 = which(PRF_LFs < 0)
            PRF_probes_smaller_0 = length(PRF_probes_smaller_0)
            PRF_combined = c(PRF_probes_greater_0, PRF_probes_smaller_0)
            if (any(PRF_combined >= 2)){ # At least 2 CpG should match (CpGs signficant in 3 cohorts are counted always)
              return(TRUE)
            }
            return(FALSE)
          })
          
          PRF_Cross_valid_CpGs_matching_dir_names = unique(PRF_df_cross_valid_CpGs$CpG)[PRF_matching_dir_index]
          
          if (length(PRF_Cross_valid_CpGs_matching_dir_names) < 1){
            PRF_Cross_valid_CpGs_matching_dir_names = NA
            PRF_Cross_valid_CpGs_matching_dir = 0
            PRF_CpG_associated_genes_cross_valid_matching_dir = NA
            PRF_CpG_associated_genes_cross_valid_matching_dir_upd = NA
            PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count = 0
            
          } else {
            
            PRF_Cross_valid_CpGs_matching_dir = length(PRF_Cross_valid_CpGs_matching_dir_names)
            PRF_Cross_valid_CpGs_matching_dir_stats = PRF_Curr_Methyl_df_signif_interval[PRF_Curr_Methyl_df_signif_interval$CpG %in% PRF_Cross_valid_CpGs_matching_dir_names, "Gene"]
            PRF_Cross_valid_CpGs_matching_dir_stats = upd_count_gene_helper(PRF_Cross_valid_CpGs_matching_dir_stats)
            PRF_Cross_valid_CpGs_matching_dir_names = paste0(PRF_Cross_valid_CpGs_matching_dir_names, collapse = ";")
            
            if (PRF_Cross_valid_CpGs_matching_dir_names == "NA"){
              PRF_Cross_valid_CpGs_matching_dir_names = NA
            }
            
            PRF_CpG_associated_genes_cross_valid_matching_dir = PRF_Cross_valid_CpGs_matching_dir_stats$gene_vector
            PRF_CpG_associated_genes_cross_valid_matching_dir_upd = PRF_Cross_valid_CpGs_matching_dir_stats$gene_vector_upd
            PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count = PRF_Cross_valid_CpGs_matching_dir_stats$gene_vector_upd_count
          }
        }
        
        PRF_Signif_CpGs_names = paste0(PRF_Signif_CpGs_names, collapse = ";")
        PRF_Cross_valid_CpGs_names = paste0(PRF_Cross_valid_CpGs_names, collapse = ";")
        
        if (PRF_Cross_valid_CpGs_names == "NA"){
          PRF_Cross_valid_CpGs_names = NA
        }
        
        PRF_CpG_associated_genes_signif = PRF_CpG_associated_genes_signif_stats$gene_vector
        PRF_CpG_associated_genes_cross_valid = PRF_CpG_associated_genes_cross_valid_stats$gene_vector
        
        PRF_CpG_associated_genes_signif_upd = PRF_CpG_associated_genes_signif_stats$gene_vector_upd
        PRF_CpG_associated_genes_cross_valid_upd = PRF_CpG_associated_genes_cross_valid_stats$gene_vector_upd
        
        PRF_CpG_associated_genes_signif_upd_count = PRF_CpG_associated_genes_signif_stats$gene_vector_upd_count
        PRF_CpG_associated_genes_cross_valid_upd_count  = PRF_CpG_associated_genes_cross_valid_stats$gene_vector_upd_count
      }
    }
    PRF_Curr_df_stat = data.frame(Chrom = PRF_curr_chrom,Start = PRF_Start, End = PRF_End, 
                                  Illum_Probes_names = PRF_Illum_Probes_names,
                                  Signif_CpGs_names = PRF_Signif_CpGs_names,
                                  Cross_valid_CpGs_names = PRF_Cross_valid_CpGs_names,
                                  Cross_valid_CpGs_matching_dir_names = PRF_Cross_valid_CpGs_matching_dir_names,
                                  Illum_Probes = PRF_Illum_Probes,
                                  Signif_CpGs = PRF_Signif_CpGs,
                                  Cross_valid_CpGs = PRF_Cross_valid_CpGs,
                                  Cross_valid_CpGs_matching_dir = PRF_Cross_valid_CpGs_matching_dir,
                                  Total_CpG_associated_genes = PRF_Total_CpG_associated_genes,
                                  CpG_associated_genes_signif = PRF_CpG_associated_genes_signif,
                                  CpG_associated_genes_cross_valid = PRF_CpG_associated_genes_cross_valid,
                                  CpG_associated_genes_cross_valid_matching_dir = PRF_CpG_associated_genes_cross_valid_matching_dir,
                                  Total_CpG_associated_genes_upd = PRF_Total_CpG_associated_genes_upd,
                                  CpG_associated_genes_signif_upd = PRF_CpG_associated_genes_signif_upd,
                                  CpG_associated_genes_cross_valid_upd = PRF_CpG_associated_genes_cross_valid_upd,
                                  CpG_associated_genes_cross_valid_matching_dir_upd = PRF_CpG_associated_genes_cross_valid_matching_dir_upd,
                                  Total_CpG_associated_genes_upd_count = PRF_Total_CpG_associated_genes_upd_count,
                                  CpG_associated_genes_signif_upd_count = PRF_CpG_associated_genes_signif_upd_count,
                                  CpG_associated_genes_cross_valid_upd_count  = PRF_CpG_associated_genes_cross_valid_upd_count,
                                  CpG_associated_genes_cross_valid_matching_dir_upd_count = PRF_CpG_associated_genes_cross_valid_matching_dir_upd_count)
    PRF_Curr_Chrom_Table_Heatmap[[index]] = PRF_Curr_df_stat
  }
  PRF_Curr_Chrom_Table_Heatmap = list_to_df(PRF_Curr_Chrom_Table_Heatmap)
  Chrom_heatmap_df_methyl[[i]] = PRF_Curr_Chrom_Table_Heatmap
}

rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])
Chrom_heatmap_df_methyl = do.call(rbind, Chrom_heatmap_df_methyl)

# Calculating ratios (Adjusting CpG counts by the number of array probes per sector)
Chrom_heatmap_df_methyl$Signif_CpGs_Ratio = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Chrom_heatmap_df_methyl$Signif_CpGs, Chrom_heatmap_df_methyl$Illum_Probes)
Chrom_heatmap_df_methyl$Signif_CpGs_Ratio_cross_valid = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Chrom_heatmap_df_methyl$Cross_valid_CpGs, Chrom_heatmap_df_methyl$Illum_Probes)
write.csv(Chrom_heatmap_df_methyl, "Chrom_heatmap_df_methyl.csv")

# Test visualization
Chrom_heatmap_plot =  Chrom_heatmap_df_methyl[,c("Chrom",
                                           "Start",
                                           "End",
                                           "Signif_CpGs_Ratio_cross_valid")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#FF0000"), output = "CpG_heatmap_test.svg")
convertSVG("CpG_heatmap_test.svg", device = "png", file = "CpG_heatmap_test.png")


################### Stats for CpGs ###################
# Getting stats for matching directions
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 1,]
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap$Var1
Combined_methyl_df_significant_overlap = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap,]
Combined_methyl_df_significant_overlap = dplyr::arrange(Combined_methyl_df_significant_overlap, CpG)

CpG_table_methyl_depr_Overlap_index = sapply(CpG_table_methyl_depr_Overlap, function(PRF_x){
  PRF_LFs = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG == PRF_x,]
  PRF_LFs = PRF_LFs$logFC
  PRF_LFs_greater_0 = which(PRF_LFs > 0)
  PRF_LFs_greater_0 = length(PRF_LFs_greater_0)
  PRF_LFs_less_0 = which(PRF_LFs < 0)
  PRF_LFs_less_0 = length(PRF_LFs_less_0)
  PRF_LFs_combined_score = c(PRF_LFs_greater_0, PRF_LFs_less_0)
  
  if (any(PRF_LFs_combined_score >= 2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
  
})
table(CpG_table_methyl_depr_Overlap_index)
# 1480 show matching directions, 1011 are not matching

# Checking chromosome map for matching
sum(Chrom_heatmap_df_methyl$Cross_valid_CpGs_matching_dir) # 1480
sum(Chrom_heatmap_df_methyl$Cross_valid_CpGs) # 2491 CpGs are cross-valid
# Everything is matching

# Preparing reporting file
Combined_methyl_df_significant_overlap_matching_dir = CpG_table_methyl_depr_Overlap[CpG_table_methyl_depr_Overlap_index]
Combined_methyl_df_significant_overlap_matching_dir = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% Combined_methyl_df_significant_overlap_matching_dir, ]
Combined_methyl_df_significant_overlap_matching_dir = dplyr::arrange(Combined_methyl_df_significant_overlap_matching_dir, CpG)
tmp_list = list(
  Combined_methyl_df_significant,
  Combined_methyl_df_significant_overlap,
  Combined_methyl_df_significant_overlap_matching_dir
)
tmp_list = lapply(tmp_list, function(x){
  x$Upd_gene_name = sapply(x$Upd_gene_name, function(y){
    if (y == "NA"){
      return(NA)
    }
    y = unlist(stri_split_fixed(y, pattern = ";"))
    y = y[y != ""]
    y = unique(y)
    y = y[!is.na(y)]
    y = y[y != "NA"]
    
    if (length(y) < 1){
      return(NA)
    }
    
    if (length(y) == 1){
      return(y)
    }
    
    y = paste0(y, collapse = ";")
    
    return(y)
  })
  return(x)
})
write.xlsx(x = tmp_list,file = "Methylation_overlaps_results.xlsx", overwrite = TRUE)
rm(tmp_list)

# Getting stats for CpGs significant in all cohorts
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap_3 = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 2,]
CpG_table_methyl_depr_Overlap_3 = CpG_table_methyl_depr_Overlap_3$Var1
CpG_table_methyl_depr_Overlap_3 = as.character(CpG_table_methyl_depr_Overlap_3)
Annot_Overlap_3 = annotation[CpG_table_methyl_depr_Overlap_3, ]
Diff_meth_df_3 = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap_3, ]
Diff_meth_df_3 = dplyr::arrange(Diff_meth_df_3, CpG)
write.xlsx(x = Annot_Overlap_3, file = "Annot_Overlap_3.xlsx", overwrite = TRUE)
write.xlsx(x = Diff_meth_df_3,file = "Diff_meth_df_3.xlsx", overwrite = TRUE)


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

# Making gene universe set for the Illumina 450K array
# Preparing the annotation file
Illumina_450K_anno_filtered = filter_annotation_illumina(annot_df = annotation, mode = "450K")
Illumina_450K_anno_filtered_genes = Illumina_450K_anno_filtered$Updated_gene_names
Illumina_450K_anno_filtered_genes = Illumina_450K_anno_filtered_genes[Illumina_450K_anno_filtered_genes != "NA"]
Illumina_450K_anno_filtered_genes = Illumina_450K_anno_filtered_genes[!is.na(Illumina_450K_anno_filtered_genes)]
Illumina_450K_anno_filtered_genes = unlist(stri_split_fixed(Illumina_450K_anno_filtered_genes, pattern = ";"))
Illumina_450K_anno_filtered_genes = Illumina_450K_anno_filtered_genes[Illumina_450K_anno_filtered_genes != ""]
Illumina_450K_anno_filtered_genes = unique(Illumina_450K_anno_filtered_genes)

# Mapping to Entrez data
ENTREZ_Illumina_450K = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Illumina_450K_anno_filtered_genes){
    return(TRUE)
  }
  return(FALSE)
}),]
ENTREZ_Illumina_450K = ENTREZ_Illumina_450K[ENTREZ_Illumina_450K$chromosome != "chr-",]
write.csv(ENTREZ_Illumina_450K, "ENTREZ_Illumina_450K.csv")
Non_mapped_genes = Illumina_450K_anno_filtered_genes[Illumina_450K_anno_filtered_genes %!in% ENTREZ_Illumina_450K$Symbol]
write(Non_mapped_genes, "Non_mapped_genes_Illumina_450K_anno.txt")

# Mapping to Entrez data by synonyms
ENTREZ_Illumina_450K_idx = vector()

for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(Illumina_450K_anno_filtered$Updated_gene_names, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(Illumina_450K_anno_filtered[Tested_idx, "chr"])){
      ENTREZ_Illumina_450K_idx[i] = TRUE
    } else {
      ENTREZ_Illumina_450K_idx[i] = FALSE
    }
  } else {
    ENTREZ_Illumina_450K_idx[i] = FALSE
  }
}

ENTREZ_Illumina_450K_2 = ENTREZ_genes_Homo_Sapiens[ENTREZ_Illumina_450K_idx,]
write.csv(ENTREZ_Illumina_450K_2, "ENTREZ_Illumina_450K_2.csv")
ENTREZ_Illumina_450K = rbind(ENTREZ_Illumina_450K, ENTREZ_Illumina_450K_2)
write.csv(ENTREZ_Illumina_450K, "ENTREZ_Illumina_450K_full.csv")
ENTREZ_Illumina_450K_ids = unique(ENTREZ_Illumina_450K$GeneID) # Gene IDs for the Illumina 450K array "universe"

# 1) Overlap > 1 similar direction
# Getting CpGs
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 1,]
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap$Var1
CpG_table_methyl_depr_Overlap_index = sapply(CpG_table_methyl_depr_Overlap, function(PRF_x){
  PRF_LFs = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG == PRF_x,]
  PRF_LFs = PRF_LFs$logFC
  PRF_LFs_greater_0 = which(PRF_LFs > 0)
  PRF_LFs_greater_0 = length(PRF_LFs_greater_0)
  PRF_LFs_less_0 = which(PRF_LFs < 0)
  PRF_LFs_less_0 = length(PRF_LFs_less_0)
  PRF_LFs_combined_score = c(PRF_LFs_greater_0, PRF_LFs_less_0)
  if (any(PRF_LFs_combined_score >= 2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
})
table(CpG_table_methyl_depr_Overlap_index)
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap[CpG_table_methyl_depr_Overlap_index]
# 1480 show matching directions, 1011 are not matching

# Getting genes
Genes_Signif_Overlap = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap, "Upd_gene_name"]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unlist(stri_split_fixed(Genes_Signif_Overlap, pattern = ";"))
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unique(Genes_Signif_Overlap)
Genes_Signif_Overlap = Genes_Signif_Overlap[!is.na(Genes_Signif_Overlap)]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != "NA"] # 1100 genes

# Mapping to Entrez
ENTREZ_Genes_Signif_Overlap = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Genes_Signif_Overlap){
    return(TRUE)
  }
  return(FALSE)
}),]
Non_mapped_genes = Genes_Signif_Overlap[Genes_Signif_Overlap %!in% ENTREZ_Genes_Signif_Overlap$Symbol] # 25 genes are not mapped

# Mapping to Entrez through synonyms
ENTREZ_Genes_Signif_overlap_idx = vector()
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(Combined_methyl_df_significant$Upd_gene_name, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(Combined_methyl_df_significant[Tested_idx, "Chromosome"])){
      ENTREZ_Genes_Signif_overlap_idx[i] = TRUE
    } else {
      ENTREZ_Genes_Signif_overlap_idx[i] = FALSE
    }
    
  } else {
    ENTREZ_Genes_Signif_overlap_idx[i] = FALSE
  }
}
table(ENTREZ_Genes_Signif_overlap_idx) # No genes were identified through synonyms
ENTREZ_Genes_Signif_Overlap_ids = unique(ENTREZ_Genes_Signif_Overlap$GeneID) # 1075 IDs
write.csv(ENTREZ_Genes_Signif_Overlap, "ENTREZ_Genes_Signif_Overlap.csv")

# Enrichment
ENRICHMENT_Illumina_Signif_CpGs_Overlap_matching_dit = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Genes_Signif_Overlap_ids,
                                                                                       universe = ENTREZ_Illumina_450K_ids,
                                                                                       categories_to_show = 30,
                                                                                       folder = "Depression_Broad_Enrichment_Signif_Overlap_matching_dir",
                                                                                       plot_name_pref = "Depression_Broad_Enrichment_Signif_Overlap_matching_dir")
writeLines(Genes_Signif_Overlap, "Genes_Signif_Overlap_similar_dir.txt")

# 2) Overlap > 1 (total)
# Getting CpGs and genes
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 1,]
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap$Var1
Genes_Signif_Overlap = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap, "Upd_gene_name"]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unlist(stri_split_fixed(Genes_Signif_Overlap, pattern = ";"))
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unique(Genes_Signif_Overlap)
Genes_Signif_Overlap = Genes_Signif_Overlap[!is.na(Genes_Signif_Overlap)]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != "NA"] # 1838 genes total

# Mapping to Entrez
ENTREZ_Genes_Signif_Overlap = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Genes_Signif_Overlap){
    return(TRUE)
  }
  return(FALSE)
}),]
Non_mapped_genes = Genes_Signif_Overlap[Genes_Signif_Overlap %!in% ENTREZ_Genes_Signif_Overlap$Symbol] # 37 genes are not mapped

# Mapping to Entrez through synonyms
ENTREZ_Genes_Signif_overlap_idx = vector()
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(Combined_methyl_df_significant$Upd_gene_name, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(Combined_methyl_df_significant[Tested_idx, "Chromosome"])){
      ENTREZ_Genes_Signif_overlap_idx[i] = TRUE
    } else {
      ENTREZ_Genes_Signif_overlap_idx[i] = FALSE
    }
    
  } else {
    ENTREZ_Genes_Signif_overlap_idx[i] = FALSE
  }
}
table(ENTREZ_Genes_Signif_overlap_idx) # No genes were identified through synonyms
ENTREZ_Genes_Signif_Overlap_ids = unique(ENTREZ_Genes_Signif_Overlap$GeneID) # 1801 IDs
write.csv(ENTREZ_Genes_Signif_Overlap, "ENTEREZ_Genes_Signif_Overlap_total.csv")

# Enrichment
ENRICHMENT_Illumina_Signif_CpGs_Overlap_matching_dit = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Genes_Signif_Overlap_ids,
                                                                                       universe = ENTREZ_Illumina_450K_ids,
                                                                                       categories_to_show = 30,
                                                                                       folder = "Depression_Broad_Enrichment_Signif_Overlap",
                                                                                       plot_name_pref = "Depression_Broad_Enrichment_Signif_Overlap")

# 3) Overlap > 2 (total)
# Getting CpGs and genes
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 2,]
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap$Var1
Genes_Signif_Overlap = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap, "Upd_gene_name"]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unlist(stri_split_fixed(Genes_Signif_Overlap, pattern = ";"))
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != ""]
Genes_Signif_Overlap = unique(Genes_Signif_Overlap)
Genes_Signif_Overlap = Genes_Signif_Overlap[!is.na(Genes_Signif_Overlap)]
Genes_Signif_Overlap = Genes_Signif_Overlap[Genes_Signif_Overlap != "NA"] # 37 genes total

# Mapping to Entrez
ENTREZ_Genes_Signif_Overlap = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Genes_Signif_Overlap){
    return(TRUE)
  }
  return(FALSE)
}),]
Non_mapped_genes = Genes_Signif_Overlap[Genes_Signif_Overlap %!in% ENTREZ_Genes_Signif_Overlap$Symbol] # 0 genes are not mapped
ENTREZ_Genes_Signif_Overlap_ids = unique(ENTREZ_Genes_Signif_Overlap$GeneID) #37 IDs

# Enrichment
ENRICHMENT_Illumina_Signif_CpGs_Overlap_2 = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Genes_Signif_Overlap_ids,
                                                                            universe = ENTREZ_Illumina_450K_ids,
                                                                            categories_to_show = 30,
                                                                            folder = "Depression_Broad_Enrichment_Signif_Overlap_all",
                                                                            plot_name_pref = "Depression_Broad_Enrichment_Signif_Overlap_all")


################### Characterization of used cohorts ###################

# PSY
PSY_SCR_CHARACT_DF = Phenotypes_SCR
colnames(PSY_SCR_CHARACT_DF) = multiple_stri_replacer(colnames(PSY_SCR_CHARACT_DF), pattern_vector = c("DAWBA_DEPBAND", "Depr_risk", "Is_genotyped"), 
                                                      replacement_vector = c("DAWBA.Depdand", "Depr.risk", "Genotyped"))
colnames(PSY_SCR_CHARACT_DF)

# Exploring variables
PSY_SCR_CHARACT_DF$DAWBA.Depdand #14
PSY_SCR_CHARACT_DF$Depr.risk #15
PSY_SCR_CHARACT_DF$Gender = factor(PSY_SCR_CHARACT_DF$Gender, levels = c("M", "W"), labels = c("Male", "Female")) #8
PSY_SCR_CHARACT_DF$Age #9
PSY_SCR_CHARACT_DF$BMI # 12

PSY_SCR_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = PSY_SCR_CHARACT_DF,
                                                                    study_char = "PSY",
                                                                    contrast_col_number = 15,
                                                                    contrast_vector = c("Low depr. risk", "High depr. risk"),
                                                                    participants_col_number = 6,
                                                                    model_covariates_col_vector = c(8,9,12,5),
                                                                    columns_to_characterise_vector = c(15, 14, 8, 9, 12),
                                                                    Remove_NA_predictors = TRUE,
                                                                    drop_P = TRUE,
                                                                    simplif_P = 3)
openxlsx::write.xlsx(PSY_SCR_CHARACT_DF[["Table"]], file = "Demographics_PSY.xlsx", overwrite = TRUE)

# GSE125105
GSE125105_CHARACT_DF = GSE125105_pheno_full
colnames(GSE125105_CHARACT_DF)
colnames(GSE125105_CHARACT_DF) = multiple_stri_replacer(colnames(GSE125105_CHARACT_DF), 
                                                        pattern_vector = c("Sex",
                                                                           "age",
                                                                           "BMI_fixed",
                                                                           "diagnosis"), 
                                                        replacement_vector = c("Gender",
                                                                               "Age",
                                                                               "BMI",
                                                                               "Diagnosis"))
# Exploring variables
GSE125105_CHARACT_DF$Gender #17
GSE125105_CHARACT_DF$Age #16
GSE125105_CHARACT_DF$Diagnosis #15
GSE125105_CHARACT_DF$BMI #49
# SNP_1_PC 47, SNP_2_PC 48
GSE125105_CHARACT_DF$Gender = factor(GSE125105_CHARACT_DF$Gender, levels = c("F","M"), labels = c("Female", "Male"))

GSE125105_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = GSE125105_CHARACT_DF,
                                                                      study_char = "GSE125105",
                                                                      contrast_col_number = 15,
                                                                      contrast_vector = c("Non-depressed", "Depression"),
                                                                      participants_col_number = 48,
                                                                      model_covariates_col_vector = c(16,17,49,47,48),
                                                                      columns_to_characterise_vector = c(15,17,16,49),
                                                                      Remove_NA_predictors = TRUE,
                                                                      drop_P = TRUE,
                                                                      simplif_P = 3)
openxlsx::write.xlsx(GSE125105_CHARACT_DF[["Table"]], file = "Demographics_GSE125105.xlsx", overwrite = TRUE)

# EGEOD72680
EGEOD72680_CHARACT_DF = E_GEOD_72680_Phenotypes
colnames(EGEOD72680_CHARACT_DF)
EGEOD72680_CHARACT_DF_browsing = EGEOD72680_CHARACT_DF[,c("treatment_for_depression", 
                                                          "beck_depression_inventory_total_score", 
                                                          "Composite_depression_NA_full")]
colnames(EGEOD72680_CHARACT_DF) = multiple_stri_replacer(colnames(EGEOD72680_CHARACT_DF), 
                                                         pattern_vector = c("Sex",
                                                                            "age",
                                                                            "treatment_for_depression", 
                                                                            "treatment_for_bipolar_disorder",
                                                                            "treatment_for_anxiety_disorder",
                                                                            "beck_depression_inventory_total_score",
                                                                            "body_mass_index",
                                                                            "race_ethnicity",
                                                                            "Composite_depression_NA_full"), 
                                                         replacement_vector = c("Gender",
                                                                                "Age",
                                                                                "Depr.treatment",
                                                                                "BP.treatment",
                                                                                "AD.treatment",
                                                                                "BDI.Total",
                                                                                "BMI",
                                                                                "Ethn.background",
                                                                                "Depr.category"))

# Exploring variables
EGEOD72680_CHARACT_DF$Gender #2
EGEOD72680_CHARACT_DF$Age #3
EGEOD72680_CHARACT_DF$Depr.treatment #5
EGEOD72680_CHARACT_DF$BP.treatment #6
EGEOD72680_CHARACT_DF$AD.treatment #8
EGEOD72680_CHARACT_DF$BDI.Total #19
EGEOD72680_CHARACT_DF$BMI #21
EGEOD72680_CHARACT_DF$Ethn.background #23
EGEOD72680_CHARACT_DF$Depr.category #69

EGEOD72680_CHARACT_DF = characterize_dataset_generelized_two_subgroups(dataset = EGEOD72680_CHARACT_DF,
                                                                       study_char = "GSE72680",
                                                                       contrast_col_number = 69,
                                                                       contrast_vector = c("Normal", "Depressed"),
                                                                       participants_col_number = 48,
                                                                       model_covariates_col_vector = c(2,3,23,21,10:15),
                                                                       columns_to_characterise_vector = c(69, 2,3, 23, 21, 5, 19),
                                                                       Remove_NA_predictors = TRUE,
                                                                       drop_P = TRUE,
                                                                       simplif_P = 3)
openxlsx::write.xlsx(EGEOD72680_CHARACT_DF[["Table"]], file = "Demographics_EGEOD72680.xlsx", overwrite = TRUE)






















