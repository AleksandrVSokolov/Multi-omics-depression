################# This script shows analysis for the overlaps, clinical, and drug data ################# 
# Note 1: The equal sign = was used as an assignment operator for typing/productivity reasons
# Note 2: In many cases loops were deliberately used instead of apply functions to enable better control of the variables
# Note 3: Some variables in the loops contain prefixes to enable easy cleanup of the environment once the loop is executed
# Note 4: The script uses data from NIH Drug portal which is now discontinue and available from PubChem only https://www.nlm.nih.gov/pubs/techbull/ja22/ja22_pubchem.html ->
# -> the last downloaded datasets with synonyms and classes are provided here: https://ftp.nlm.nih.gov/projects/chemidlease/chemid-20230222.zip
# Note 5: The script uses data from the DrugBank database (https://go.drugbank.com/). To obtain its dataset, please make a request to https://go.drugbank.com/
# Note 6: The script uses data from the IUPHAR Guide To Pharmacology database (https://www.guidetopharmacology.org/) To obtain its dataset, please use the database download page

Working_directory = "~/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort" # Replace with an appropriate path
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

# A function to check if provided URL exists
valid_url = function(URL, t = 2){
  Connection = url(URL)
  Check = suppressWarnings(try(open.connection(Connection, open = "rt", timeout = t), silent = TRUE)[1])
  suppressWarnings(try(close.connection(Connection), silent = TRUE))
  ifelse(is.null(Check), TRUE, FALSE)
}

# A function to perform searches in ClinicalTrials.gov with two search terms (condition and treatment)
# The function is vectorized and generates all combinations of condtions and treatments
clinical_trial_downloader_two_terms = function(condition_terms = NULL,
                                                treatment_terms = NULL,
                                                folder,
                                                url_suffix = NULL){
  
  if (!(paste0('./', folder) %in% list.dirs())){
    dir.create(folder)
  } else {
    stop("Please specify a new folder")
  }
  
  if (!is.null(condition_terms)){
    condition_terms = stringr::str_trim(condition_terms)
    condition_terms = condition_terms[condition_terms != ""]
    condition_terms = condition_terms[!is.na(condition_terms)]
    condition_terms_initial = condition_terms
    condition_terms = toupper(condition_terms)
    condition_terms = unique(condition_terms)
    condition_terms = stringi::stri_replace_all_fixed(str = condition_terms, pattern = " ", replacement = "+")
    condition_terms = URLencode(condition_terms)
    
    if (length(treatment_terms) >= 1){
      condition_terms_initial = sapply(condition_terms_initial, function(x) paste0(x, rep("", times = length(treatment_terms))))
    }
    
  }
  
  if (!is.null(treatment_terms)){
    treatment_terms = stringr::str_trim(treatment_terms)
    treatment_terms = treatment_terms[treatment_terms != ""]
    treatment_terms = treatment_terms[!is.na(treatment_terms)]
    terms_initial = treatment_terms
    
    if (length(condition_terms) >= 1){
      terms_initial = rep(terms_initial, times = length(condition_terms))
    }
    
    
    treatment_terms = toupper(treatment_terms)
    treatment_terms = unique(treatment_terms)
    treatment_terms = stringi::stri_replace_all_fixed(str = treatment_terms, pattern = " ", replacement = "+")
    treatment_terms = URLencode(treatment_terms)
  }
  
  if (is.null(condition_terms)){
    condition_terms = ""
    condition_terms_initial = ""
  }
  if (is.null(treatment_terms)){
    treatment_terms = ""
    terms_initial = NA
  }
  
  Base_string_1_main = "https://clinicaltrials.gov/ct2/results/download_fields?cond="
  Base_string_2_main = "&intr="
  Base_string_3_main = "&type=Intr&down_count=10000&down_chunk=1&down_fmt=csv&down_flds=all"
  Base_string_1_check = "https://clinicaltrials.gov/ct2/results?cond="
  Base_string_2_check = "&intr="
  Base_string_3_check = "&type=Intr"
  
  if (!is.null(url_suffix)){
    Base_string_3_main = paste0(Base_string_3_main, url_suffix)
    Base_string_3_check = paste0(Base_string_3_check, url_suffix)
  }
  
  
  urls = paste0(Base_string_1_main, condition_terms, Base_string_2_main)
  urls = as.vector(sapply(urls, function(x) paste0(x, treatment_terms, Base_string_3_main)))
  
  urls_check = paste0(Base_string_1_check, condition_terms, Base_string_2_check)
  urls_check = as.vector(sapply(urls_check, function(x) paste0(x, treatment_terms, Base_string_3_check)))
  
  reporting_file = as.vector(sapply(condition_terms, function(x) paste0(x, "&&&", treatment_terms)))
  reporting_file = data.frame(reporting_file)
  colnames(reporting_file) = "Term"
  reporting_file$Index = NA
  reporting_file$URL = NA
  reporting_file$URL_check = NA
  reporting_file$Status = NA
  reporting_file$Trials = NA
  reporting_file$Intervention = NA
  reporting_file$Condition_searched = NA
  
  index = 1
  while(index <= length(urls)){
    
    reporting_file$Index[index] = index
    reporting_file$URL[index] = urls[index]
    reporting_file$URL_check[index] = urls_check[index]
    reporting_file$Intervention[index] = terms_initial[index]
    reporting_file$Condition_searched[index] = condition_terms_initial[index]
    
    if (valid_url(urls_check[index])){
      
      destination = paste0(folder,"/",index, "___",".csv")
      error_flag_ctgov = FALSE
      skip_Flag = FALSE
      error_counter_get_trials = 0
      
      tryCatch({
        download.file(urls[index], destination, quiet = TRUE)
        reporting_file$Status[index] = "Downloaded"
      }, error = function(e){error_flag_ctgov <<- TRUE})
      
      while(error_flag_ctgov){
        Sys.sleep(2**error_counter_get_trials)
        error_counter_get_trials = error_counter_get_trials + 1
        
        warning(paste0('Error detected during getting trials file from   ',
                       urls[index],
                       '  ',
                       Sys.time(),
                       ' Trying to reopen link  #',
                       error_counter_get_trials), immediate. = TRUE)
        
        if (error_counter_get_trials > 10){
          
          warning(paste0("Unable to get trials from ", urls[index]))
          skip_Flag = TRUE
          error_flag_ctgov = FALSE
          
        } else {
          
          tryCatch({
            download.file(urls[index], destination, quiet = TRUE)
            reporting_file$Status[index] = "Downloaded"
            error_flag_ctgov <<- FALSE
          }, error = function(e){error_flag_ctgov <<- TRUE})
          
        }
        
      }
      
      if (skip_Flag){
        suppressWarnings(rm(list = c("error_flag_ctgov", "error_counter_get_trials", "skip_Flag")))
        reporting_file$Status[index] = "Link failed"
        reporting_file$Trials[index] = 0
      }
      
    } else {
      suppressWarnings(rm(list = c("error_flag_ctgov", "error_counter_get_trials", "skip_Flag")))
      reporting_file$Status[index] = "No trials/Invalid"
      reporting_file$Trials[index] = 0
    }
    
    percent = round((index/length(urls))*100, 4)
    writeLines(paste0(index,'   ', percent,' % complete at ', Sys.time()))
    index = index + 1
    Sys.sleep(1)
    
  }
  
  files = list.files(folder)
  if (length(files)<1){
    reporting_file$Trials = 0
    write.csv(reporting_file, paste0(folder, '/', "reporting_file", ".csv"))
    return("0 trials found")
  }
  
  files = files[files != "reporting_file.csv"]
  if (length(files)<1){
    reporting_file$Trials = 0
    write.csv(reporting_file, paste0(folder, '/', "reporting_file", ".csv"))
    return("0 trials found")
  }
  
  paths = paste0(folder,"/", files)
  
  files_list = list()
  for (i in 1:length(paths)){
    df = read.csv(paths[i], header = TRUE, encoding = "UTF-8")
    idx = unlist(stringi::stri_split_fixed(paths[i], pattern = "/"))
    idx = idx[length(idx)]
    idx = unlist(stringi::stri_split_fixed(idx, pattern = "___"))
    idx = idx[1]
    df$Search_term = reporting_file[reporting_file$Index == idx, "Term"]
    df$Requested_inter = reporting_file[reporting_file$Index == idx, "Intervention"]
    df$Requested_cond = reporting_file[reporting_file$Index == idx, "Condition_searched"]
    reporting_file[reporting_file$Index == idx, "Trials"] = nrow(df)
    write.csv(df, paths[i])
    files_list[[i]] = df
  }
  
  write.csv(reporting_file, paste0(folder, "/", "reporting_file", ".csv"))
  output_data = do.call(rbind, files_list)
  return(output_data)
}
# A function to map SNPs from GTEx IDs to dbSNP

map_SNPs_GTEx_dbSNP_b38 = function(PRF_SNP_ID_GTEx, 
                                   PRF_bigBedToBed_location = "/home/aleksandr/UCSC_folder/bigBedToBed",
                                   PRF_db_SNP_location = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb",
                                   PRF_span = 100,
                                   PRF_tmp_file = "tmp_map_file.bed",
                                   PRF_lag_param = 2){
  
  # Parsing SNP info
  PRF_parsed_SNP = unlist(stri_split_fixed(PRF_SNP_ID_GTEx, pattern = "_"))
  
  print(PRF_parsed_SNP)
  
  # Output file
  PRF_out_file = paste0(getwd(), "/", PRF_tmp_file)
  
  # Generate request to UCSC server
  PRF_reqest = paste0(PRF_bigBedToBed_location, 
                      " ",
                      PRF_db_SNP_location,
                      " ",
                      "-chrom=",
                      PRF_parsed_SNP[1],
                      " ",
                      "-start=",
                      as.numeric(PRF_parsed_SNP[2])-PRF_span,
                      " ",
                      "-end=",
                      as.numeric(PRF_parsed_SNP[2])+PRF_span,
                      " ", 
                      PRF_out_file)
  
  # Execute the request to UCSC server
  terminalExecute(PRF_reqest, show = FALSE)
  Sys.sleep(PRF_lag_param)
  
  # Reading obtained file
  PRF_obtained_SNPs = read.table(PRF_out_file, sep = "\t")
  
  # Remove the TMP file
  unlink(PRF_out_file)
  
  # Filtering SNPs
  PRF_filtered_SNPs = PRF_obtained_SNPs[PRF_obtained_SNPs$V3 == as.numeric(PRF_parsed_SNP[2]),]
  if (nrow(PRF_filtered_SNPs) < 1) {
    # No SNPs found
    return(NA)
  }
  
  PRF_filtered_SNPs_ref_alleles = stri_split_fixed(PRF_filtered_SNPs$V5, pattern = ",")
  PRF_filtered_SNPs_alt_alleles = stri_split_fixed(PRF_filtered_SNPs$V7, pattern = ",")
  
  # Checking alleles
  PRF_ref_index = sapply(PRF_filtered_SNPs_ref_alleles, function(x){
    PRF_parsed_SNP[3] %in% x
  })
  PRF_alt_index = sapply(PRF_filtered_SNPs_alt_alleles, function(x){
    PRF_parsed_SNP[4] %in% x
  })
  PRF_ind_matrix = cbind(PRF_ref_index, PRF_alt_index)
  PRF_ind_fin = apply(PRF_ind_matrix, 1, all)
  
  # Required SNP df
  PRF_req_SNP_df = PRF_filtered_SNPs[PRF_ind_fin,]
  
  if (nrow(PRF_req_SNP_df) < 1) {
    # No SNPs found
    return(NA)
  }
  
  if (nrow(PRF_req_SNP_df) > 1) {
    # More than 1 SNPs found
    PRF_SNPs = paste0(PRF_req_SNP_df$V4, collapse = ";")
    return(PRF_SNPs)
  }
  
  return(PRF_req_SNP_df$V4)
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


################### Importing results from the GWAS Catalog analysis ###################

# Getting SNP data
GWAS_Cat_Path = "GWAS_Catalog/"
GWAS_Cat_Filtered_df = smart_fread(paste0(GWAS_Cat_Path, "GWAS_CAT_Depress_specific.csv")) # 2168 records
length(unique(GWAS_Cat_Filtered_df$PUBMEDID)) # 78 pubmed IDs for GWAS 
GWAS_CAT_Depress_specific_SNPs = GWAS_Cat_Filtered_df$SNPS
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[stri_detect_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "rs")]
GWAS_CAT_Depress_specific_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "; ")) # 2254 SNPs
GWAS_CAT_Depress_specific_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = ", ")) # 2281 SNPs
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[!is.na(GWAS_CAT_Depress_specific_SNPs)]
write(GWAS_CAT_Depress_specific_SNPs, "GWAS_CAT_Depress_specific_SNPs.txt")
GWAS_CAT_Depress_specific_SNPs = sapply(GWAS_CAT_Depress_specific_SNPs, function(x){
  if (stri_detect_fixed(x, pattern = "-")){
    x = unlist(stri_split_fixed(x, pattern = "-"))
    x = x[which.max(nchar(x))]
  }
  if (stri_detect_fixed(x, pattern = "_")){
    x = unlist(stri_split_fixed(x, pattern = "_"))
    x = x[which.max(nchar(x))]
  }
  return(x)
})
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[stri_detect_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "rs")] # 2281 SNPs
write(GWAS_CAT_Depress_specific_SNPs)
GWAS_CAT_Depress_specific_SNPs_unique = unique(GWAS_CAT_Depress_specific_SNPs) # 2073 SNPs

# Getting track data
GWAS_CAT_Depress_specific_SNPs_153_track = smart_fread(paste0(GWAS_Cat_Path, "GWAS_CAT_Depress_specific_SNPs_153_track.csv")) # 2052 records
GWAS_CAT_Depress_specific_SNPs_153_mapped = smart_fread(paste0(GWAS_Cat_Path, "GWAS_CAT_Depress_specific_SNPs_153_mapped.csv")) # 2052 records
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = smart_fread(paste0(GWAS_Cat_Path, "Entrez_GWAS_CAT_Depress_specific_SNPs_153_mapped.csv")) # 847 records

# Getting gene data
Frequent_SNP_genes_GWAS_depression = readLines(paste0(GWAS_Cat_Path, "Frequent_SNP_genes_GWAS_depression.txt"))
Frequent_genes_GWAS_depression = readLines(paste0(GWAS_Cat_Path, "Frequent_genes_GWAS_depression.txt"))

# Importing chromosome map
Chrom_map_GWAS_CAT = smart_fread(paste0(GWAS_Cat_Path, "Chrom_map_GWAS_CAT.csv"))
colnames(Chrom_map_GWAS_CAT) = stri_replace_all_fixed(colnames(Chrom_map_GWAS_CAT), pattern = ".", replacement = "_")
colnames(Chrom_map_GWAS_CAT) = paste0("GWAS_", colnames(Chrom_map_GWAS_CAT))


################### Importing results from the DNA methylation analysis ###################

# Importing results for the cohorts
Methyl_data_folder = "DNA_methylation/"
Methyl_depression_broad_PSY_SCR = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_PSY_SCR.csv")) # Replace with another path if needed
Methyl_depression_broad_PSY_SCR_significant = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_PSY_SCR_significant.csv")) # Replace with another path if needed
Methyl_depression_broad_PSY_SCR_significant_LF = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_PSY_SCR_significant_LF.csv")) # Replace with another path if needed

Methyl_depression_broad_E_GEOD_72680 = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_E_GEOD_72680.csv")) # Replace with another path if needed
Methyl_depression_broad_E_GEOD_72680_significant = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_E_GEOD_72680_significant.csv")) # Replace with another path if needed
Methyl_depression_broad_E_GEOD_72680_significant_LF = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_E_GEOD_72680_significant_LF.csv")) # Replace with another path if needed

Methyl_depression_broad_GSE125105 = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_GSE125105.csv")) # Replace with another path if needed
Methyl_depression_broad_GSE125105_significant = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_GSE125105_significant.csv")) # Replace with another path if needed
Methyl_depression_broad_GSE125105_significant_LF = smart_fread(paste0(Methyl_data_folder, "Methyl_depression_broad_GSE125105_significant_LF.csv")) # Replace with another path if needed

Combined_methyl_df = rbind(Methyl_depression_broad_PSY_SCR,
                           Methyl_depression_broad_GSE125105,
                           Methyl_depression_broad_E_GEOD_72680)
Combined_methyl_df_significant = Combined_methyl_df[Combined_methyl_df$P.Value < 0.05,] # 54905 nominally significant findings
Combined_methyl_df_significant_unique = distinct(Combined_methyl_df_significant, CpG, .keep_all = TRUE) # 52365 unique

# Overlap
CpG_table_methyl_depr = as.data.frame(table(Combined_methyl_df_significant$CpG))
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr[CpG_table_methyl_depr$Freq > 1,]
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap$Var1
CpG_table_methyl_depr_Overlap_all = as.character(CpG_table_methyl_depr_Overlap)
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
CpG_table_methyl_depr_Overlap = CpG_table_methyl_depr_Overlap[CpG_table_methyl_depr_Overlap_index]
Df_methyl_Signif_Overlap = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap, ]
Df_methyl_Signif_Overlap_all = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% CpG_table_methyl_depr_Overlap_all, ]

# Chromosome map
Chrom_heatmap_df_methyl = smart_fread(paste0(Methyl_data_folder, "Chrom_heatmap_df_methyl.csv"))
colnames(Chrom_heatmap_df_methyl) = paste0("Meth_",colnames(Chrom_heatmap_df_methyl))


################### Importing results from the transcriptome analysis ###################

# Importing results for the cohorts
Transcr_data_folder = "Transcriptome/"
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif_FDR.csv")) # Replace with another path if needed
Analysis_limma_affymetrix_depr_broad_GSE98793 = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_affymetrix_depr_broad_GSE98793.csv")) # Replace with another path if needed
Analysis_limma_affymetrix_depr_broad_GSE98793_signif = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif.csv")) # Replace with another path if needed
Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_affymetrix_depr_broad_GSE98793_signif_matching_dirs.csv")) # Replace with another path if needed

Analysis_limma_illumina_depr_broad_GSE46743 = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE46743.csv")) # Replace with another path if needed
Analysis_limma_illumina_depr_broad_GSE46743_signif = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE46743_signif.csv")) # Replace with another path if needed
Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE46743_signif_matching_dirs.csv")) # Replace with another path if needed

Analysis_limma_illumina_depr_broad_GSE64930 = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE64930.csv")) # Replace with another path if needed
Analysis_limma_illumina_depr_broad_GSE64930_signif = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE64930_signif.csv")) # Replace with another path if needed
Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs = smart_fread(paste0(Transcr_data_folder, "Analysis_limma_illumina_depr_broad_GSE64930_signif_matching_dirs.csv")) # Replace with another path if needed

# Combined expression dataset
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
Genes_matching_dir_expression = readLines(paste0(Transcr_data_folder, "Genes_matching_dir.txt"))
Combined_expression_df_matching_dir_cohorts = Combined_expression_df[Combined_expression_df$Gene %in% Genes_matching_dir_expression,]

# Chromosome map
Chrom_df_Transcriptome_map = smart_fread(paste0(Transcr_data_folder, "Chrom_df_Transcriptome_map.csv"))
colnames(Chrom_df_Transcriptome_map) = stri_replace_all_fixed(colnames(Chrom_df_Transcriptome_map),
                                                              pattern = "PRF_",
                                                              replacement = "")
colnames(Chrom_df_Transcriptome_map) = paste0("Transcr_",colnames(Chrom_df_Transcriptome_map))


################### Analyzing overlaps between OMICS layers ###################

# GWAS depression
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[!is.na(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)]
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, pattern = ";"))
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unique(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)

# Methylation genes
Genes_matching_CpGs = Df_methyl_Signif_Overlap$Upd_gene_name
Genes_matching_CpGs = Genes_matching_CpGs[!is.na(Genes_matching_CpGs)]
Genes_matching_CpGs = Genes_matching_CpGs[Genes_matching_CpGs != "NA"]
Genes_matching_CpGs = unlist(stri_split_fixed(Genes_matching_CpGs, pattern = ";"))
Genes_matching_CpGs = unique(Genes_matching_CpGs)

# Transcriptome genes
Genes_matching_dir_expression = unique(Genes_matching_dir_expression) # We use only those genes that were analyzable in all 3 cohorts

# Creating Venn Diagrams
Venn_overlaping_genes_3_layers = list(
  "GWAS Frequent Genes" = Frequent_genes_GWAS_depression,
  "DNA methylation (blood)" = Genes_matching_CpGs,
  "Transcriptome (blood)" = Genes_matching_dir_expression
)
Venn_overlaping_genes_3_layers_2 = list(
  "GWAS all genes" = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes,
  "DNA methylation (blood)" = Genes_matching_CpGs,
  "Transcriptome (blood)" = Genes_matching_dir_expression
)
make_Venn_digram_list(named_list = Venn_overlaping_genes_3_layers, palette = 5, plot_full_path = "Venn_overlaping_genes_3_layers.pdf")
make_Venn_digram_list(named_list = Venn_overlaping_genes_3_layers_2, palette = 6, plot_full_path = "Venn_overlaping_genes_3_layers_full_GWAS.pdf")

# Making graphs

# Saving intersect info
Intersect_ALL = Reduce(intersect, list(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes,
                                       Genes_matching_CpGs,
                                       Genes_matching_dir_expression)) # 3 genes "FOXP1" "VPS41" "AKTIP"
Intersect_GWAS_Methyl = intersect(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, Genes_matching_CpGs) 
Intersect_GWAS_Methyl = Intersect_GWAS_Methyl[Intersect_GWAS_Methyl %!in% Intersect_ALL] # 84 genes
Intersect_GWAS_Transcr = intersect(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, Genes_matching_dir_expression) 
Intersect_GWAS_Transcr = Intersect_GWAS_Transcr[Intersect_GWAS_Transcr %!in% Intersect_ALL] # 19 genes
Intersect_Methyl_Transcr = intersect(Genes_matching_CpGs, Genes_matching_dir_expression) 
Intersect_Methyl_Transcr = Intersect_Methyl_Transcr[Intersect_Methyl_Transcr %!in% Intersect_ALL] # 38 genes
write(Intersect_GWAS_Methyl, "Intersect_GWAS_Methyl.txt")
write(Intersect_GWAS_Transcr, "Intersect_GWAS_Transcr.txt")
write(Intersect_Methyl_Transcr, "Intersect_Methyl_Transcr.txt")

GWAS_DS = data_frame(Study = "GWAS all Depr. genes", Gene =  GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)
Meth_DS = data_frame(Study = "DNA methylation (blood)", Gene =  Genes_matching_CpGs)
Trans_DS = data_frame(Study = "Transcriptome (blood)", Gene =  Genes_matching_dir_expression)
Graph_DS = rbind(GWAS_DS, Meth_DS, Trans_DS)
network_data = graph_from_data_frame(d = Graph_DS[,1:2], directed = TRUE)
network_data = toVisNetworkData(network_data)
network_data$nodes$group = sapply(network_data$nodes$id, function(x){
  if (x %in% Graph_DS$Study){
    return("Study")
  }
  if (x %in% Intersect_GWAS_Methyl){
    return("Intersect GWAS-Methyl")
  } else if (x %in% Intersect_GWAS_Transcr){
    return("Intersect GWAS-Transcriptome")
  } else if (x %in% Intersect_Methyl_Transcr){
    return("Intersect Methyl-Transcriptome")
  } else if (x %in% Intersect_ALL){
    return("Triple intersect")
  } else {
    return("Normal gene")
  }
})
net = visNetwork(nodes = network_data$nodes, edges = network_data$edges, height = "5000px", width = "5000px") %>%
  visIgraphLayout(layout = "layout_with_fr", physics = FALSE, randomSeed = 190) %>%
  visGroups(groupname = "Intersect GWAS-Methyl", color = list("background" = "skyblue", border = "black")) %>% 
  visGroups(groupname = "Intersect GWAS-Transcriptome", color = list("background" = "green", border = "black")) %>% 
  visGroups(groupname = "Intersect Methyl-Transcriptome", color = list("background" = "red", border = "black")) %>%
  visGroups(groupname = "Triple intersect", color = list("background" = "yellow", border = "black")) %>% 
  visGroups(groupname = "Normal gene", color = list("background" = "#FBCEB1", border = "black")) %>% 
  visGroups(groupname = "Study", color = list("background" = "orange", border = "black")) %>% 
  visNodes(size = 10, borderWidth = 0.3) %>%
  visEdges(width = 0.5) %>%
  visPhysics(barnesHut = list(gravitationalConstant = -100000, avoidOverlap = 0.8)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)
visSave(net, file = "Findings_3_layers_net.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("Findings_3_layers_net.html", "Findings_3_layers_net.png", vwidth = 5000, vheight = 5000, zoom = 1)

# Small network for nodes in overlaps
int_list = c(Intersect_GWAS_Methyl, Intersect_GWAS_Transcr, Intersect_Methyl_Transcr, Intersect_ALL)
int_list = unique(int_list)
GWAS_DS = data_frame(Study = "GWAS all Depr. genes", Gene =  GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)
Meth_DS = data_frame(Study = "DNA methylation (blood)", Gene =  Genes_matching_CpGs)
Trans_DS = data_frame(Study = "Transcriptome (blood)", Gene =  Genes_matching_dir_expression)
Graph_DS = rbind(GWAS_DS, Meth_DS, Trans_DS)
Graph_DS = Graph_DS[Graph_DS$Gene %in% int_list,]
network_data = graph_from_data_frame(d = Graph_DS[,1:2], directed = TRUE)
network_data = toVisNetworkData(network_data)
network_data$nodes$group = sapply(network_data$nodes$id, function(x){
  if (x %in% Graph_DS$Study){
    return("Study")
  }
  if (x %in% Intersect_GWAS_Methyl){
    return("Intersect GWAS-Methyl")
  } else if (x %in% Intersect_GWAS_Transcr){
    return("Intersect GWAS-Transcriptome")
  } else if (x %in% Intersect_Methyl_Transcr){
    return("Intersect Methyl-Transcriptome")
  } else if (x %in% Intersect_ALL){
    return("Triple intersect")
  } else {
    return("Normal gene")
  }
})

net = visNetwork(nodes = network_data$nodes, edges = network_data$edges, height = "3000px", width = "3000px") %>%
  visIgraphLayout(layout = "layout_with_fr", physics = FALSE, randomSeed = 190) %>%
  visGroups(groupname = "Intersect GWAS-Methyl", color = list("background" = "skyblue", border = "black")) %>% 
  visGroups(groupname = "Intersect GWAS-Transcriptome", color = list("background" = "green", border = "black")) %>% 
  visGroups(groupname = "Intersect Methyl-Transcriptome", color = list("background" = "red", border = "black")) %>% 
  visGroups(groupname = "Triple intersect", color = list("background" = "yellow", border = "black")) %>% 
  visGroups(groupname = "Normal gene", color = list("background" = "#FBCEB1", border = "black")) %>% 
  visGroups(groupname = "Study", color = list("background" = "orange", border = "black")) %>% 
  visNodes(size = 10, borderWidth = 0.3, font = list(size = 25)) %>%
  visEdges(width = 0.5) %>%
  visPhysics(barnesHut = list(gravitationalConstant = -10000, avoidOverlap = 0.8, centralGravity = 0.1)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)
visSave(net, file = "Findings_3_layers_net_small.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("Findings_3_layers_net_small.html", "Findings_3_layers_net_small.png", vwidth = 3000, vheight = 3000)


################### Doing regression of chromosomal regions ###################

# Preparing datasets
Chrom_heatmap_df_methyl_short = Chrom_heatmap_df_methyl[Chrom_heatmap_df_methyl$Meth_Chrom != "chrY",]
Chrom_df_Transcriptome_map_short = Chrom_df_Transcriptome_map[Chrom_df_Transcriptome_map$Transcr_Chrom != "chrY",]
Combined_Chrom_Df = cbind(Chrom_map_GWAS_CAT, Chrom_heatmap_df_methyl_short, Chrom_df_Transcriptome_map_short)

# Removing Centromeres and empty regions
Combined_Chrom_Df_fixed = Combined_Chrom_Df[Combined_Chrom_Df$Transcr_identified_genes != "",]
Combined_Chrom_Df_fixed = Combined_Chrom_Df_fixed[!is.na(Combined_Chrom_Df_fixed$Transcr_identified_genes),]

# Cross-valid CpGs ~ int + Depression SNPs + Illumina probes + Reference genes + error
PRF_fit_lm = lm(Meth_Cross_valid_CpGs ~ GWAS_Ident_SNPs_count + Meth_Illum_Probes + Transcr_identified_genes_count, data = Combined_Chrom_Df_fixed)
PRF_Coefs_Methyl_lm_1 = as.data.frame(summary(PRF_fit_lm)[["coefficients"]])
PRF_Coefs_Methyl_lm_1 = cbind(Dependent.var = "Cross-valid CpGs", Coef = c("Intercept","Ident.SNPs","Illum.Probes", "Ref.Genes"),PRF_Coefs_Methyl_lm_1)
plot(Combined_Chrom_Df_fixed$GWAS_Ident_SNPs_count, Combined_Chrom_Df_fixed$Meth_Cross_valid_CpGs) # No linear relationship

# Cross-valid proteins (matching direction) ~ int + Depression SNPs + Transcriptome array genes + Reference genes + error
PRF_fit_lm = lm(Transcr_Cross_signif_genes_similar_dir_count ~ GWAS_Ident_SNPs_count  + Transcr_identified_transcripts_count + Transcr_identified_genes_count, 
                data = Combined_Chrom_Df_fixed)
PRF_Coefs_Methyl_lm_2 = as.data.frame(summary(PRF_fit_lm)[["coefficients"]])
PRF_Coefs_Methyl_lm_2 = cbind(Dependent.var = "Cross-valid proteins (similar dir.)",Coef = c("Intercept","Ident.SNPs","Array genes", "Ref.Genes"),PRF_Coefs_Methyl_lm_2)
plot(Combined_Chrom_Df_fixed$GWAS_Ident_SNPs_count, Combined_Chrom_Df_fixed$Transcr_Cross_signif_genes_similar_dir_count) # No linear relationship

# Cross-valid proteins (matching direction) ~ int + Cross-valid CpGs + Transcriptome array genes + Reference genes + error
PRF_fit_lm = lm(Transcr_Cross_signif_genes_similar_dir_count ~ Meth_Cross_valid_CpGs  + Transcr_identified_transcripts_count + Transcr_identified_genes_count, 
                data = Combined_Chrom_Df_fixed)
PRF_Coefs_Methyl_lm_3 = as.data.frame(summary(PRF_fit_lm)[["coefficients"]])
PRF_Coefs_Methyl_lm_3 = cbind(Dependent.var = "Cross-valid proteins (similar dir.)", Coef = c("Intercept","Cross-valid CpGs","Array genes", "Ref.Genes"),PRF_Coefs_Methyl_lm_3)
plot(Combined_Chrom_Df_fixed$Meth_Cross_valid_CpGs, Combined_Chrom_Df_fixed$Transcr_Cross_signif_genes_similar_dir_count) # No linear relationship

# Writing the outputs
total_list = list(PRF_Coefs_Methyl_lm_1, PRF_Coefs_Methyl_lm_2, PRF_Coefs_Methyl_lm_3)
openxlsx::write.xlsx(total_list, "genome clustering regresssion_upd.xlsx")

# Spearman correlations
cor.test(Combined_Chrom_Df_fixed$GWAS_Ident_SNPs_count, Combined_Chrom_Df_fixed$Transcr_identified_genes_count, method = "spearman") # p-value = 0.6156, -0.009354482 
cor.test(Combined_Chrom_Df_fixed$Transcr_Cross_signif_genes_similar_dir_count, Combined_Chrom_Df_fixed$Transcr_identified_transcripts_count, method = "spearman") # 0.4430285  p-value < 2.2e-16
cor.test(Combined_Chrom_Df_fixed$Meth_Cross_valid_CpGs, Combined_Chrom_Df_fixed$Meth_Illum_Probes, method = "spearman") # 0.5925004 p-value < 2.2e-16


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

# Enrichment for overlap between methylation and GWAS
Overlapping_genes_gwas_methyl = network_data$nodes
Overlapping_genes_gwas_methyl = Overlapping_genes_gwas_methyl[Overlapping_genes_gwas_methyl$group == "Intersect GWAS-Methyl",]
Overlapping_genes_gwas_methyl = Overlapping_genes_gwas_methyl$id
ENTREZ_Overlapping_genes_gwas_methyl = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Overlapping_genes_gwas_methyl){
    return(TRUE)
  }
  return(FALSE)
}),]
all(Overlapping_genes_gwas_methyl %in% ENTREZ_Overlapping_genes_gwas_methyl$Symbol) # All genes are mapped
ENTREZ_Illumina_450K = smart_fread("DNA_methylation/ENTREZ_Illumina_450K_full.csv") # Replace with another path if needed
ENTREZ_Illumina_450K_ids = unique(ENTREZ_Illumina_450K$GeneID)
debug(run_enrichment_GO_KEGG_gene_set)
ENRICHMENT_Overlapping_genes_gwas_methyl = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Overlapping_genes_gwas_methyl$GeneID,
                                                                           universe = ENTREZ_Illumina_450K_ids,
                                                                           categories_to_show = 30,
                                                                           folder = "Depression_Broad_Enrichment_Overlapping_genes_gwas_methyl",
                                                                           plot_name_pref = "Depression_Broad_Enrichment_Overlapping_genes_gwas_methyl")

# Enrichment for overlapping genes GWAS and expression
Overlapping_genes_gwas_transcript = network_data$nod
Overlapping_genes_gwas_transcript = Overlapping_genes_gwas_transcript[Overlapping_genes_gwas_transcript$group == "Intersect GWAS-Transcriptome",]
Overlapping_genes_gwas_transcript = Overlapping_genes_gwas_transcript$id
ENTREZ_Overlapping_genes_gwas_transcript = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Overlapping_genes_gwas_transcript){
    return(TRUE)
  }
  return(FALSE)
}),]
all(Overlapping_genes_gwas_transcript %in% ENTREZ_Overlapping_genes_gwas_transcript$Symbol) 
ENTREZ_Transcriptome_universe_overlap = smart_fread("Transcriptome/ENTREZ_Transcriptome_universe_overlap.csv") # Replace with another path if needed
all(Overlapping_genes_gwas_transcript %in% ENTREZ_Transcriptome_universe_overlap$Symbol) 
ENRICHMENT_Overlapping_genes_gwas_transcript = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Overlapping_genes_gwas_transcript$GeneID,
                                                                               universe = unique(ENTREZ_Transcriptome_universe_overlap$GeneID),
                                                                               categories_to_show = 30,
                                                                               folder = "Depression_Broad_Enrichment_Overlapping_genes_gwas_transcript",
                                                                               plot_name_pref = "Depression_Broad_Enrichment_Overlapping_genes_gwas_transcript")

# Enrichment for overlap between methylation and expression
Overlapping_genes_methyl_transcript = network_data$nod
Overlapping_genes_methyl_transcript = Overlapping_genes_methyl_transcript[Overlapping_genes_methyl_transcript$group == "Intersect Methyl-Transcriptome",]
Overlapping_genes_methyl_transcript = Overlapping_genes_methyl_transcript$id
ENTREZ_Overlapping_genes_methyl_transcript = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Overlapping_genes_methyl_transcript){
    return(TRUE)
  }
  return(FALSE)
}),]
all(Overlapping_genes_methyl_transcript %in% ENTREZ_Overlapping_genes_methyl_transcript$Symbol) # All genes are mapped
all(Overlapping_genes_methyl_transcript %in% ENTREZ_Illumina_450K$Symbol) # All genes are presented
ENRICHMENT_Overlapping_genes_methyl_transcript = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Overlapping_genes_methyl_transcript$GeneID,
                                                                                 universe = unique(ENTREZ_Transcriptome_universe_overlap$GeneID),
                                                                                 categories_to_show = 30,
                                                                                 folder = "Depression_Broad_Enrichment_Overlapping_genes_methyl_transcript",
                                                                                 plot_name_pref = "Depression_Broad_Enrichment_Overlapping_genes_methyl_transcript")


################### Visualizing chromosomal maps ###################

# Getting the cytoBand track from UCSC
mySession = browserSession("UCSC")
genome(mySession) = "hg19"
Chromosome_coord_table = chromosome_table_getter()
Full_cytoband_track = list()

# Note: the loop is labelled with prefix PRF_ to cleanup the variables
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

# Final chromosome map
# SNP vs Reference
Combined_Chrom_Df$GWAS_SNP_vs_Ref_ratio = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Combined_Chrom_Df$GWAS_Ident_SNPs_count, Combined_Chrom_Df$Transcr_identified_genes_count)
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "GWAS_SNP_vs_Ref_ratio")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#000000"), output = "SNP_heatmap.svg")
convertSVG("SNP_heatmap.svg", device = "png", file = "SNP_heatmap.png")

# Cross-valid CpGs vs ref. array CpGs
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Meth_Signif_CpGs_Ratio_cross_valid")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#FF0000"), output = "CpG_heatmap.svg")
convertSVG("CpG_heatmap.svg", device = "png", file = "CpG_heatmap.png")

# Cross-valid genes vs ref. array genes
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Transcr_Signif_genes_cross_valid_Ratio")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#097969"), output = "Transcr_heatmap.svg")
convertSVG("Transcr_heatmap.svg", device = "png", file = "Transcr_heatmap.png")

# Cross-valid CpGs (matching dir.) vs ref. array CpGs
Combined_Chrom_Df$Meth_Signif_CpGs_Ratio_cross_valid_matching_dir = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Combined_Chrom_Df$Meth_Cross_valid_CpGs_matching_dir, Combined_Chrom_Df$Meth_Illum_Probes)
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Meth_Signif_CpGs_Ratio_cross_valid_matching_dir")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#FF5733"), output = "CpG_match_heatmap.svg") # dark-blue
convertSVG("CpG_match_heatmap.svg", device = "png", file = "CpG_match_heatmap.png")

# Cross-valid genes (matching dir.) vs ref. array genes
Combined_Chrom_Df$Transcr_Signif_genes_cross_valid_matching_dir_Ratio  = mapply(function(x,y){
  if (x == 0){
    return(0)
  }
  result = x/y
  return(result)
}, Combined_Chrom_Df$Transcr_Cross_signif_genes_similar_dir_count, Combined_Chrom_Df$Transcr_identified_transcripts_genes_unique_count)
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Transcr_Signif_genes_cross_valid_matching_dir_Ratio")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#097969"), output = "Transcr_match_heatmap.svg")
convertSVG("Transcr_match_heatmap.svg", device = "png", file = "Transcr_match_heatmap.png")


# Count heatmaps for the main text
# Transcriptome heatmap counts
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Transcr_Cross_signif_genes_similar_dir_count")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#9B26B6"), output = "Transcr_match_count_heatmap.svg")
convertSVG("Transcr_match_count_heatmap.svg", device = "png", file = "Transcr_match_count_heatmap.png")

# SNP heatmap counts
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "GWAS_Ident_SNPs_count")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#000000"), output = "GWAS_heatmap_counts.svg")
convertSVG("GWAS_heatmap_counts.svg", device = "png", file = "GWAS_heatmap_counts.png") # SNPs will be provided as counts

# CpG heatmap counts
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Meth_Cross_valid_CpGs_matching_dir")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#191970"), output = "CpG_match_count_heatmap.svg")
convertSVG("CpG_match_count_heatmap.svg", device = "png", file = "CpG_match_count_heatmap.png")

# Gene reference heatmap
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Transcr_identified_genes_count")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#FF0000"), output = "Gene_ref_heatmap_counts.svg")
convertSVG("Gene_ref_heatmap_counts.svg", device = "png", file = "Gene_ref_heatmap_counts.png")

# Non-zero value at 3 levels
# Columns required: Transcr_Signif_genes_cross_valid_matching_dir_Ratio & GWAS_Ident_SNPs_count & Meth_Cross_valid_CpGs_matching_dir
Combined_Chrom_Df$Non_zero_3 = 0
Combined_Chrom_Df$Non_zero_3  = mapply(function(x,y,z){
  if (all(c(x,y,z) > 0)){
    return(1)
  }
  return(0)
}, Combined_Chrom_Df$Transcr_Cross_signif_genes_similar_dir_count, Combined_Chrom_Df$GWAS_Ident_SNPs_count, Combined_Chrom_Df$Meth_Cross_valid_CpGs_matching_dir)
Chrom_heatmap_plot =  Combined_Chrom_Df[,c("GWAS_Chrom",
                                           "GWAS_Start",
                                           "GWAS_End",
                                           "Non_zero_3")]
colnames(Chrom_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_heatmap_plot, colorset1 = c("#FFFFFF", "#FF0000"), output = "Non_zero_3_heatmap.svg")
convertSVG("Non_zero_3_heatmap.svg", device = "png", file = "Non_zero_3_heatmap.png")

# Saving the map with all stats
write.csv(Combined_Chrom_Df, "Combined_Chrom_Df_final.csv")


################### IUPHAR and Clinical trial analysis ###################

# IUPHAR
IUPHAR_22_Sep_2022 = smart_fread("approved_drug_primary_target_interactions.csv") # Replace with another path if needed
IUPHAR_22_Sep_2022_targets = smart_fread("targets_and_families.csv") # Replace with another path if needed
IUPHAR_22_Sep_2022 = inner_join(IUPHAR_22_Sep_2022, IUPHAR_22_Sep_2022_targets, by = c("Target ID" = "Target id"))
IUPHAR_22_Sep_2022$Ligand = toupper(IUPHAR_22_Sep_2022$Ligand)
IUPHAR_Ligands =  smart_fread("ligand_id_mapping.csv") # Replace with another path if needed
IUPHAR_22_Sep_2022 = inner_join(IUPHAR_22_Sep_2022, IUPHAR_Ligands, by = c("Ligand ID" = "Ligand id"))
IUPHAR_22_Sep_2022$INN = toupper(IUPHAR_22_Sep_2022$INN)
IUPHAR_22_Sep_2022_curated = read.xlsx("IUPHAR_22_Sep_2022_curated.xlsx") # Replace with another path if needed
IUPHAR_22_Sep_2022_curated$Ligand = toupper(IUPHAR_22_Sep_2022_curated$Ligand)

# ChemIDplus Subset Data
# Data comes from parsed XML from https://www.nlm.nih.gov/databases/download/chemidplus.html and https://ftp.nlm.nih.gov/projects/chemidlease
chem_ID_synonyms = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/ChemIDPlusData/synonyms_chemIDplus.txt")
chem_ID_synonyms = as.data.frame(apply(chem_ID_synonyms, 2, toupper))
chem_ID_categories = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/ChemIDPlusData/classes_chemIDplus.txt")
chem_ID_categories = as.data.frame(apply(chem_ID_categories, 2, toupper))
chem_ID_categories_antidepress = chem_ID_categories[stri_detect_fixed(chem_ID_categories$Class, pattern = "ANTIDEPRESS"),]

chem_ID_categories_antidepress$Approvals = sapply(chem_ID_categories_antidepress$Molecule, function(x){
  if (x %in% IUPHAR_22_Sep_2022_curated$Ligand) {
    return(TRUE)
  }
  synonyms = chem_ID_synonyms[chem_ID_synonyms$Molecule == x, "Synonym"]
  
  if (length(synonyms) < 1){
    return(FALSE)
  }
  
  return(any(synonyms %in% IUPHAR_22_Sep_2022_curated$Ligand))
})

# Preparing lists
Overlapping_genes_triple = Intersect_ALL
GWAS_List = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[GWAS_CAT_Depress_specific_SNPs_153_mapped__genes %!in% c(Overlapping_genes_gwas_methyl, 
                                                                                                                      Overlapping_genes_methyl_transcript, 
                                                                                                                      Overlapping_genes_gwas_transcript,
                                                                                                                      Overlapping_genes_triple)]
Methyl_List = Genes_matching_CpGs[Genes_matching_CpGs %!in% c(Overlapping_genes_gwas_methyl, 
                                                              Overlapping_genes_methyl_transcript, 
                                                              Overlapping_genes_gwas_transcript,
                                                              Overlapping_genes_triple)] 
Transcript_List = Genes_matching_dir_expression[Genes_matching_dir_expression %!in% c(Overlapping_genes_gwas_methyl,
                                                                                      Overlapping_genes_methyl_transcript, 
                                                                                      Overlapping_genes_gwas_transcript,
                                                                                      Overlapping_genes_triple)] 
Protein_Lists = list(Overlapping_genes_gwas_methyl,
                     Overlapping_genes_methyl_transcript,
                     Overlapping_genes_gwas_transcript,
                     Overlapping_genes_triple,
                     GWAS_List,
                     Methyl_List,
                     Transcript_List)
names(Protein_Lists) = c(
  "Overlap GWAS-Methylation",
  "Overlap Methylation-Transcriptome",
  "Overlap GWAS-Transcriptome",
  "Triple intersect",
  "GWAS",
  "Methylation",
  "Transcriptome"
)

# Extended drug stats
Proteins_Drugs_Stats = list()

# Note: the loop variables are labelled with prefix PRF_ to enable easy cleanup
for (i in 1:length(Protein_Lists)){
  print(i)
  PRF_Current_List = Protein_Lists[[i]]
  
  # Selecting IUPHAR drugs
  PRF_Drugs_All_df = IUPHAR_22_Sep_2022_curated[IUPHAR_22_Sep_2022_curated$HGNC.symbol %in% PRF_Current_List, ]
  PRF_Drugs_All = PRF_Drugs_All_df$Ligand
  PRF_Drugs_All_Count = length(unique(PRF_Drugs_All))
  PRF_Targets_All = unique(PRF_Drugs_All_df$HGNC.symbol)
  PRF_Targets_All_Count = length(PRF_Targets_All)
  PRF_Targets_All_Percent = round(PRF_Targets_All_Count/length(PRF_Current_List)*100, digits = 1)
  
  # Mapping stats
  PRF_Mapped_Drugs = sapply(unique(PRF_Drugs_All), function(x){
    if (x %in% chem_ID_synonyms$Molecule){
      return(x)
    }
    if (x %in% chem_ID_synonyms$Synonym){
      x = chem_ID_synonyms[chem_ID_synonyms$Synonym == x,]
      x = x$Molecule[1]
      return(x)
    }
    return(NA)
  })
  PRF_Mapped_drugs_names_init = PRF_Mapped_Drugs[!is.na(PRF_Mapped_Drugs)]
  PRF_Mapped_drugs_names = paste0(PRF_Mapped_Drugs[!is.na(PRF_Mapped_Drugs)], collapse = ";")
  PRF_Non_Mapped_drugs_names = PRF_Mapped_Drugs[is.na(PRF_Mapped_Drugs)]
  PRF_Non_Mapped_drugs_names = names(PRF_Non_Mapped_drugs_names)
  PRF_Non_Mapped_drugs_names = paste0(PRF_Non_Mapped_drugs_names, collapse = ";")
  PRF_Mapped_Drugs = PRF_Mapped_Drugs[!is.na(PRF_Mapped_Drugs)]
  PRF_Mapped_Drugs_CID_Categories = chem_ID_categories[chem_ID_categories$Molecule %in% PRF_Mapped_Drugs,]
  PRF_Mapped_Drugs_CID_Categories_ANTIDEPRESS = PRF_Mapped_Drugs_CID_Categories[stri_detect_fixed(PRF_Mapped_Drugs_CID_Categories$Class, pattern = "ANTIDEPRESS"),]
  PRF_Drugs_Antidepress = unique(PRF_Mapped_Drugs_CID_Categories_ANTIDEPRESS$Molecule)
  PRF_Drugs_Antidepress_Count = length(PRF_Drugs_Antidepress)
  PRF_Drugs_Antidepress_df = PRF_Drugs_All_df[PRF_Drugs_All_df$Ligand %in% names(PRF_Mapped_drugs_names_init[PRF_Mapped_drugs_names_init %in% PRF_Drugs_Antidepress]),]
  PRF_Targets_Antidepress = unique(PRF_Drugs_Antidepress_df$HGNC.symbol)
  PRF_Targets_Antidepress_Count = length(PRF_Targets_Antidepress)
  PRF_Targets_Antidepress_Percent = round(PRF_Targets_Antidepress_Count/length(PRF_Current_List)*100, digits = 1)
  
  if (length(PRF_Drugs_Antidepress) < 1){
    PRF_Drugs_Antidepress = NA
  } else {
    PRF_Drugs_Antidepress = paste0(PRF_Drugs_Antidepress, collapse = ";")
  }
  
  if (length(PRF_Targets_Antidepress) < 1){
    PRF_Targets_Antidepress = NA
  } else {
    PRF_Targets_Antidepress = paste0(PRF_Targets_Antidepress, collapse = ";")
  }
  
  PRF_Output_df = data.frame(
    Category = names(Protein_Lists)[i],
    PRF_Drugs_All = paste0(PRF_Drugs_All, collapse = ";"),
    PRF_Drugs_All_Count = PRF_Drugs_All_Count,
    PRF_Targets_All =  paste0(PRF_Targets_All, collapse = ";"),
    PRF_Targets_All_Count = PRF_Targets_All_Count,
    PRF_Targets_All_Percent = PRF_Targets_All_Percent,
    PRF_Mapped_drugs_names = PRF_Mapped_drugs_names,
    PRF_Non_Mapped_drugs_names = PRF_Non_Mapped_drugs_names,
    PRF_Drugs_Antidepress_Count = PRF_Drugs_Antidepress_Count,
    PRF_Targets_Antidepress_Count = PRF_Targets_Antidepress_Count,
    PRF_Targets_Antidepress_Percent = PRF_Targets_Antidepress_Percent,
    PRF_Drugs_Antidepress = PRF_Drugs_Antidepress,
    PRF_Targets_Antidepress = PRF_Targets_Antidepress
  )
  colnames(PRF_Output_df) = stri_replace_all_fixed(colnames(PRF_Output_df), pattern = "PRF_", replacement = "")
  Proteins_Drugs_Stats[[i]] = PRF_Output_df
}
Proteins_Drugs_Stats = list_to_df(Proteins_Drugs_Stats)
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])
colnames(Proteins_Drugs_Stats) = c("Category",
                                     "Drugs_All",
                                     "Drugs_All_Count",
                                     "Targets_All",
                                     "Targets_All_Count",
                                     "Targets_All_Percent",
                                     "ChemID_plus_Mapped_drugs_names",
                                     "ChemID_plus_Non_Mapped_drugs_names",
                                     "Drugs_Antidepress_Count",
                                     "Targets_Antidepress_Count",
                                     "Targets_Antidepress_Percent",
                                     "Drugs_Antidepress",
                                     "Targets_Antidepress")
write.csv(Proteins_Drugs_Stats, "Proteins_Drugs_Stats.csv")

# Adding Clinical Trial searches
Proteins_Drugs_Stats$Drugs_in_depr_trials = NA

# Preparing folders to store sheets
Folders = paste0("Screen_", Proteins_Drugs_Stats$Category)
Folders = multiple_stri_replacer(string = Folders, pattern_vector = c(" ","-"), replacement_vector = c("_","_"))

# Note: the loop variables are labelled with prefix PRF_ to enable easy cleanup
for (i in 1:nrow(Proteins_Drugs_Stats)){
  print(Proteins_Drugs_Stats$Category[i])
  PRF_df = Proteins_Drugs_Stats[i,]
  
  if (PRF_df$Drugs_All_Count == 0){
    Proteins_Drugs_Stats$Drugs_in_depr_trials[i] = 0
    print(Proteins_Drugs_Stats$Drugs_in_depr_trials[i])
    
  } else {
    
    PRF_curr_drugs = PRF_df$Drugs_All
    PRF_curr_drugs = unlist(stri_split_fixed(PRF_curr_drugs, pattern = ";"))
    PRF_curr_drugs = unique(PRF_curr_drugs)
    PRF_curr_screen = clinical_trial_downloader_two_terms(condition_terms = "depression",
                                                          treatment_terms = PRF_curr_drugs,
                                                          folder = Folders[i])
    if (is.character(PRF_curr_screen)){
      Proteins_Drugs_Stats$Drugs_in_depr_trials[i] = 0
      print(Proteins_Drugs_Stats$Drugs_in_depr_trials[i])
    } else {
      Proteins_Drugs_Stats$Drugs_in_depr_trials[i] = length(unique(PRF_curr_screen$Search_term))
      print(Proteins_Drugs_Stats$Drugs_in_depr_trials[i])
    }
  }
}
Proteins_Drugs_Stats_Full = Proteins_Drugs_Stats
Proteins_Drugs_Stats_Full$Drugs_in_depr_trials_names = NA
Proteins_Drugs_Stats_Full$Targets_in_depr_trials_count = 0
Proteins_Drugs_Stats_Full$Targets_in_depr_trials_names = NA

for (i in 1:nrow(Proteins_Drugs_Stats_Full)){
  print(Proteins_Drugs_Stats_Full$Category[i])
  PRF_df = Proteins_Drugs_Stats[i,]
  
  if (PRF_df$Drugs_in_depr_trials == 0){
    Proteins_Drugs_Stats_Full$Drugs_in_depr_trials_names[i] = NA
    Proteins_Drugs_Stats_Full$Targets_in_depr_trials_count[i] = 0
    Proteins_Drugs_Stats_Full$Targets_in_depr_trials_names[i] = NA
  } else {
    PRF_curr_drugs = PRF_df$Drugs_All
    PRF_curr_drugs = unlist(stri_split_fixed(PRF_curr_drugs, pattern = ";"))
    PRF_curr_drugs = unique(PRF_curr_drugs)
    
    PRF_report_file = smart_fread(paste0(Folders[i], "/reporting_file.csv"))
    PRF_Trial_drugs = PRF_curr_drugs[PRF_report_file$Trials != 0]
    Proteins_Drugs_Stats_Full$Drugs_in_depr_trials_names[i] = paste0(PRF_Trial_drugs, collapse = ";")
    PRF_trial_targets = IUPHAR_22_Sep_2022_curated[IUPHAR_22_Sep_2022_curated$Ligand %in% PRF_Trial_drugs, "HGNC.symbol"]
    PRF_trial_targets = PRF_trial_targets[PRF_trial_targets %in% unique(unlist(stri_split_fixed(PRF_df$Targets_All, pattern = ";")))]
    Proteins_Drugs_Stats_Full$Targets_in_depr_trials_count[i] = length(unique(PRF_trial_targets))
    Proteins_Drugs_Stats_Full$Targets_in_depr_trials_names[i] = paste0(unique(PRF_trial_targets), collapse = ";")
  }
}

# Reshaping the dataset for plotting
Proteins_Drugs_Stats_narrow = Proteins_Drugs_Stats_Full[,c("Category",
                                                             "Drugs_All_Count",
                                                             "Targets_All_Count",
                                                             "Drugs_Antidepress_Count",
                                                             "Targets_Antidepress_Count",
                                                             "Drugs_in_depr_trials",
                                                             "Targets_in_depr_trials_count")]
Proteins_Drugs_Stats_narrow = gather(Proteins_Drugs_Stats_narrow, Counts, Numbers, c(2:7))
Proteins_Drugs_Stats_narrow$Counts = factor(Proteins_Drugs_Stats_narrow$Counts, levels = unique(Proteins_Drugs_Stats_narrow$Counts),
                                            labels = c("Drugs (all)", "Targets (all)", "Antidepressant drugs chemID", "Targets (antidepress. chemID)",
                                                       "Drugs depr. trials", "Targets depr. trials"))
Proteins_Drugs_Stats_narrow$Category = factor(Proteins_Drugs_Stats_narrow$Category)

# Generating the plot
Stats_drugs_plot = ggplot(Proteins_Drugs_Stats_narrow, aes(x = Category, y = Numbers, fill = Counts)) +
  geom_col(position = position_dodge(0.8), width = 0.8) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuOr")) +
  geom_text(aes(label = Numbers), position = position_dodge(0.8), vjust = -0.2) +
  labs(x = "Analysis level", y = "Count", fill = "Statistic") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(size = 13, colour = "black"),
    axis.text.x  = element_text(size = 13, colour = "black", angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 13, colour = "black"),
    legend.title = element_text(face = "bold", size = 14),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(size = 1),
    panel.grid.major.y = element_line(size = 0.5, linetype = 2, color =  "black"),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
Stats_drugs_plot

# Saving
write.csv(Proteins_Drugs_Stats_Full, "Proteins_Drugs_Stats_Full_Final.csv")
pdf("Stats_drugs_plot_final.pdf", width = 10, height = 10)
Stats_drugs_plot
dev.off()


################### DrugBank and Clinical trial analysis ###################

# DB path
Drug_Bank_DB_path = "/home/aleksandr/Desktop/WORK/DrugBank_2023/full database.xml" # Replace with another path if needed

# Preparing DB 
Drug_bank_db = dbparser::read_drugbank_xml_db(Drug_Bank_DB_path)

# Running parsers
Drug_bank_db_parsed = dbparser::run_all_parsers()

# Extracting drug datasets
Drug_bank_general = Drug_bank_db_parsed$general_information
Drug_bank_target = Drug_bank_db_parsed$targets
Drug_bank_synonyms = Drug_bank_db_parsed$synonyms
Drug_bank_protein = Drug_bank_db_parsed$targets_polypeptides

# Merging data
Drug_bank_merged = inner_join(Drug_bank_general, Drug_bank_target, by = c("primary_key" = "parent_key"))
Drug_bank_merged = inner_join(Drug_bank_merged, Drug_bank_protein, by = c("id" = "parent_id"))
Drug_bank_merged = Drug_bank_merged[Drug_bank_merged$known_action != "no",]
Drug_bank_merged_curated = Drug_bank_merged[Drug_bank_merged$known_action == "yes",] # Only targets with "known" action

# Converting to the R's native data frame
Drug_bank_merged = as.data.frame(Drug_bank_merged)
Drug_bank_merged_curated = as.data.frame(Drug_bank_merged_curated)
Drug_bank_synonyms = as.data.frame(Drug_bank_synonyms)

# Drug Bank stats
Proteins_Drugs_Stats_Drug_Bank = list()
for (i in 1:length(Protein_Lists)){
  print(i)
  
  # Getting matching drugs
  PRF_matched_drugs = Drug_bank_merged[Drug_bank_merged$gene_name %in% Protein_Lists[[i]], "name.x"]
  PRF_matched_drugs = unique(PRF_matched_drugs)
  PRF_matched_drugs_known_action = Drug_bank_merged_curated[Drug_bank_merged_curated$gene_name %in% Protein_Lists[[i]], "name.x"]
  PRF_matched_drugs_known_action = unique(PRF_matched_drugs_known_action)
  PRF_drugs_known_action_all_ids = Drug_bank_merged[Drug_bank_merged$name.x %in% PRF_matched_drugs_known_action, "primary_key"]
  PRF_drugs_known_action_all_ids = unique(PRF_drugs_known_action_all_ids)
  PRF_matched_drugs_known_action_synonyms = Drug_bank_synonyms[Drug_bank_synonyms$`drugbank-id` %in% PRF_drugs_known_action_all_ids, "synonym"]
  
  # Getting corresponding targets
  PRF_targets_known_action = Drug_bank_merged_curated[Drug_bank_merged_curated$gene_name %in% Protein_Lists[[i]], "gene_name"]
  PRF_targets_known_action = unique(PRF_targets_known_action)
  PRF_Output_df = data.frame(
    category = names(Protein_Lists)[i],
    PRF_matched_drugs = paste0(PRF_matched_drugs, collapse = ";"),
    PRF_matched_drugs_count = length(PRF_matched_drugs),
    PRF_matched_drugs_known_action = paste0(PRF_matched_drugs_known_action, collapse = ";"),
    PRF_matched_drugs_known_action_count = length(PRF_matched_drugs_known_action),
    PRF_matched_drugs_known_action_synonyms = paste0(PRF_matched_drugs_known_action_synonyms, collapse = ";"),
    PRF_matched_targets_known_action = paste0(PRF_targets_known_action, collapse = ";"),
    PRF_matched_targets_known_action_count = length(PRF_targets_known_action)
  )
  colnames(PRF_Output_df) = stri_replace_all_fixed(colnames(PRF_Output_df), pattern = "PRF_", replacement = "")
  Proteins_Drugs_Stats_Drug_Bank[[i]] = PRF_Output_df
}
Proteins_Drugs_Stats_Drug_Bank = list_to_df(Proteins_Drugs_Stats_Drug_Bank)

# Getting clinical trial data
Folders = paste0("DB_Screen_", Proteins_Drugs_Stats_Drug_Bank$category)
Folders = multiple_stri_replacer(string = Folders, pattern_vector = c(" ","-"), replacement_vector = c("_","_"))
Folders = paste0("DB_Screen_folder/", Folders)

# Parsing clinical trial data; the matching is performed for drugs with known action
dir.create("DB_Screen_folder")
# Note: the loop variables are labelled with prefix PRF_ to enable easy cleanup
# Note: the loop takes > 2 hours to run due to capped number of requests per second (change Sys.sleep in clinical_trial_downloader_two_terms to avoid it)
for (i in 1:nrow(Proteins_Drugs_Stats_Drug_Bank)){
  print(Proteins_Drugs_Stats_Drug_Bank$category[i])
  PRF_df = Proteins_Drugs_Stats_Drug_Bank[i,]
  
  if (PRF_df$matched_drugs_known_action_count == 0){
    
    next
    
  } else {
    
    PRF_curr_drugs = PRF_df$matched_drugs_known_action_synonyms
    PRF_curr_drugs = unlist(stri_split_fixed(PRF_curr_drugs, pattern = ";"))
    PRF_curr_drugs = unique(PRF_curr_drugs)
    PRF_curr_screen = clinical_trial_downloader_two_terms(condition_terms = "depression",
                                                          treatment_terms = PRF_curr_drugs,
                                                          folder = Folders[i])
  }
}

# Getting stats for drugs
Proteins_Drugs_Stats_Drug_Bank$matched_drugs_known_action_synonyms_curated = NA
Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_names = NA
Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_count = NA
Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_names_known_action = NA
Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_count = NA

# Note: the loop variables are labelled with prefix PRF_ to enable easy cleanup
for (i in 1:nrow(Proteins_Drugs_Stats_Drug_Bank)){
  print(Proteins_Drugs_Stats_Drug_Bank$category[i])
  PRF_df = Proteins_Drugs_Stats_Drug_Bank[i,]
  
  if (PRF_df$matched_drugs_known_action_count == 0){
    Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_names[i] = NA
    Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_count[i]  = 0
    Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_names_known_action[i]  = NA
    Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_count[i] = 0
    
  } else {
    
    PRF_curr_drugs = PRF_df$matched_drugs_known_action
    PRF_curr_drugs = unlist(stri_split_fixed(PRF_curr_drugs, pattern = ";"))
    PRF_curr_drugs = unique(PRF_curr_drugs)
    PRF_curr_folder = Folders[i]
    PRF_curr_files = list.files(PRF_curr_folder)
    
    if (length(PRF_curr_files) <= 1){
      
      # Only reporting file is there, there are no matching trials
      Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_names[i] = NA
      Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_count[i]  = 0
      Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_names_known_action[i]  = NA
      Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_count[i] = 0
      
    } else {
      
      # Some drugs matched; need to iterate through drugs
      PRF_report_file = paste0(PRF_curr_folder, "/reporting_file.csv")
      PRF_report_file = read.csv(PRF_report_file, header = TRUE, encoding = "UTF-8")
      
      PRF_observed_drugs = vector()
      for (index in 1:length(PRF_curr_drugs)){
        PFR_selected_drug = PRF_curr_drugs[index]
        PFR_DB_ID = Drug_bank_merged[Drug_bank_merged$name.x %in% PFR_selected_drug, "primary_key"]
        PFR_DB_ID = unique(PFR_DB_ID)
        
        # Get all synonyms
        PFR_selected_synonyms = Drug_bank_synonyms[Drug_bank_synonyms$`drugbank-id` %in% PFR_DB_ID, "synonym"]
        
        # Synonyms curation
        PFR_selected_synonyms_nchar = sapply(PFR_selected_synonyms, nchar)
        PFR_selected_synonyms_nchar_capital = str_count(PFR_selected_synonyms, "[A-Z]")
        PFR_selected_synonyms_len_ctrl = PFR_selected_synonyms_nchar > 4
        PFR_selected_synonyms_abbrev_ctrl = PFR_selected_synonyms_nchar_capital != PFR_selected_synonyms_nchar
        PFR_selected_synonyms = intersect(PFR_selected_synonyms[PFR_selected_synonyms_len_ctrl], PFR_selected_synonyms[PFR_selected_synonyms_abbrev_ctrl])
        PRF_curated_synonyms = c(Proteins_Drugs_Stats_Drug_Bank$matched_drugs_known_action_synonyms_curated[i], PFR_selected_synonyms)
        PRF_curated_synonyms = PRF_curated_synonyms[!is.na(PRF_curated_synonyms)]
        Proteins_Drugs_Stats_Drug_Bank$matched_drugs_known_action_synonyms_curated[i] = paste0(PRF_curated_synonyms, collapse = ";")
        
        # Detecting trials
        PRF_selected_trials = PRF_report_file[PRF_report_file$Intervention %in% PFR_selected_synonyms, ]
        
        if (any(PRF_selected_trials$Trials > 0)){
          # at least one synonym worked
          PRF_observed_drugs[index] = PFR_selected_drug
        } else {
          PRF_observed_drugs[index] = NA
        }
      }
      
      # Getting list of drugs
      PRF_observed_drugs = PRF_observed_drugs[!is.na(PRF_observed_drugs)]
      Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_names[i] = paste0(PRF_observed_drugs, collapse = ";")
      Proteins_Drugs_Stats_Drug_Bank$drugs_in_depr_trials_count[i]  = length(PRF_observed_drugs)
      
      # Getting stats for targets
      PRF_targets_known_action = Drug_bank_merged_curated[Drug_bank_merged_curated$name.x %in% PRF_observed_drugs, "gene_name"]
      PRF_targets_list = unlist(stri_split_fixed(PRF_df$matched_targets_known_action, pattern = ";"))
      PRF_targets_known_action = PRF_targets_known_action[PRF_targets_known_action %in% PRF_targets_list]
      PRF_targets_known_action = unique(PRF_targets_known_action)
      Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_names_known_action[i]  = paste0(PRF_targets_known_action, collapse = ";")
      Proteins_Drugs_Stats_Drug_Bank$targets_in_depr_trials_count[i] = length(PRF_targets_known_action)
    }
  }
}

# Reshaping the dataset for plotting
Proteins_Drugs_Stats_Drug_Bank_narrow = Proteins_Drugs_Stats_Drug_Bank[,c("category",
                                                                          "matched_drugs_count",
                                                                          "matched_drugs_known_action_count",
                                                                          "matched_targets_known_action_count",
                                                                          "drugs_in_depr_trials_count",
                                                                          "targets_in_depr_trials_count")]
Proteins_Drugs_Stats_Drug_Bank_narrow = gather(Proteins_Drugs_Stats_Drug_Bank_narrow, Counts, Numbers, c(2:6))
Proteins_Drugs_Stats_Drug_Bank_narrow$Counts = factor(Proteins_Drugs_Stats_Drug_Bank_narrow$Counts, levels = unique(Proteins_Drugs_Stats_Drug_Bank_narrow$Counts),
                                                      labels = c("Drugs (all)", "Drugs (known action target)", "Known action targets",
                                                                 "Drugs depr. trials", "Targets depr. trials"))
Proteins_Drugs_Stats_Drug_Bank_narrow$category = factor(Proteins_Drugs_Stats_Drug_Bank_narrow$category)

# Preparing the plot
Stats_drugs_DB_plot = ggplot(Proteins_Drugs_Stats_Drug_Bank_narrow, aes(x = category, y = Numbers, fill = Counts)) +
  geom_col(position = position_dodge(0.8), width = 0.8) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuOr")) +
  geom_text(aes(label = Numbers), position = position_dodge(0.8), vjust = -0.2) +
  labs(x = "Analysis level", y = "Count", fill = "Statistic") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(size = 13, colour = "black"),
    axis.text.x  = element_text(size = 13, colour = "black", angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 13, colour = "black"),
    legend.title = element_text(face = "bold", size = 14),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(size = 1),
    panel.grid.major.y = element_line(size = 0.5, linetype = 2, color =  "black"),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
Stats_drugs_DB_plot

# Saving the plot
write.csv(Proteins_Drugs_Stats_Drug_Bank, "Proteins_Drugs_Stats_Drug_Bank.csv")
pdf("Stats_drugs_DB_plot.pdf", width = 10, height = 10)
Stats_drugs_DB_plot
dev.off()
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Writing full drug stats
Full_drug_stats = list("IUPHAR" = Proteins_Drugs_Stats_Full, "DrugBank" = Proteins_Drugs_Stats_Drug_Bank)
openxlsx::write.xlsx(Full_drug_stats, "Full_drug_stats.xlsx")


################### DrugBank drug enrichment analysis ###################
DrugBank_drugs_known_target_lists = Proteins_Drugs_Stats_Drug_Bank[, c(
  "category",
  "matched_drugs_known_action",
  "matched_drugs_known_action_synonyms",
  "matched_drugs_known_action_count",
  "drugs_in_depr_trials_names"
)]


DrugBank_drugs_known_target_lists$chemID_names = NA
DrugBank_drugs_known_target_lists$chemID_names_count = NA

# Mapping loop
for (i in 1:nrow(DrugBank_drugs_known_target_lists)){
  PRF_drugs = DrugBank_drugs_known_target_lists$matched_drugs_known_action[i]
  PRF_drugs = unlist(stri_split_fixed(PRF_drugs, pattern = ";"))
  PRF_drugs = PRF_drugs[PRF_drugs!=""]
  PRF_drugs = toupper(PRF_drugs)
  
  # Mapping to chemID names
  PRF_Mapped_Drugs = sapply(PRF_drugs, function(x){
    if (x %in% chem_ID_synonyms$Molecule){
      return(x)
    }
    if (x %in% chem_ID_synonyms$Synonym){
      x = chem_ID_synonyms[chem_ID_synonyms$Synonym == x,]
      x = x$Molecule[1]
      return(x)
    }
    return(NA)
  })
  PRF_Mapped_Drugs = PRF_Mapped_Drugs[!is.na(PRF_Mapped_Drugs)]
  names(PRF_Mapped_Drugs) = NULL
  PRF_Mapped_Drugs = unique(PRF_Mapped_Drugs)
  PRF_Mapped_Drugs_cd = paste0(PRF_Mapped_Drugs, collapse = ";")
  DrugBank_drugs_known_target_lists$chemID_names[i] = PRF_Mapped_Drugs_cd
  DrugBank_drugs_known_target_lists$chemID_names_count[i] = length(PRF_Mapped_Drugs)
}
DrugBank_drugs_known_target_lists$Mapping_fraction = DrugBank_drugs_known_target_lists$chemID_names_count / DrugBank_drugs_known_target_lists$matched_drugs_known_action_count

# Obtaining known action drug bank universe
Drug_Bank_Known_Action_universe = Drug_bank_merged_curated$name.x
Drug_Bank_Known_Action_universe = toupper(Drug_Bank_Known_Action_universe)
Drug_Bank_Known_Action_universe = unique(Drug_Bank_Known_Action_universe)
Drug_Bank_Known_Action_universe_chemID = sapply(Drug_Bank_Known_Action_universe, function(x){
  if (x %in% chem_ID_synonyms$Molecule){
    return(x)
  }
  if (x %in% chem_ID_synonyms$Synonym){
    x = chem_ID_synonyms[chem_ID_synonyms$Synonym == x,]
    x = x$Molecule[1]
    return(x)
  }
  return(NA)
})
Drug_Bank_Known_Action_universe_chemID = Drug_Bank_Known_Action_universe_chemID[!is.na(Drug_Bank_Known_Action_universe_chemID)]
Drug_Bank_Known_Action_universe_classes = chem_ID_categories[chem_ID_categories$Molecule %in% Drug_Bank_Known_Action_universe_chemID,]
Drug_Bank_Known_Action_universe_classes_MESH = Drug_Bank_Known_Action_universe_classes[stri_detect_fixed(Drug_Bank_Known_Action_universe_classes$Class_source, "MESH"),]
Drug_Bank_Known_Action_universe_classes_MESH = Drug_Bank_Known_Action_universe_classes_MESH[, c("Class", "Molecule")]

# Performing enrichments for located drug lists
dir.create("Drug_enrichments")
folders_drug_enrich = paste0("Drug_enrichments/", DrugBank_drugs_known_target_lists$category)

for (i in 1:nrow(DrugBank_drugs_known_target_lists)){
  print(i)
  # Making dir 
  dir.create(folders_drug_enrich[i])
  
  # Extracting drugs
  PRF_drugs = DrugBank_drugs_known_target_lists$chemID_names[i]
  PRF_drugs = unlist(stri_split_fixed(PRF_drugs, pattern = ";"))
  PRF_drugs = PRF_drugs[PRF_drugs!=""]
  
  if (length(PRF_drugs) < 1){
    next
  }
  
  # Performing enrichment
  PRF_drug_enrichment = enricher(gene = PRF_drugs, 
                                 universe = Drug_Bank_Known_Action_universe_classes_MESH$Molecule,
                                 TERM2GENE = Drug_Bank_Known_Action_universe_classes_MESH,
                                 minGSSize = 5, 
                                 maxGSSize = 100)
  PRF_file_name_xlsx = paste0(folders_drug_enrich[i], "/", DrugBank_drugs_known_target_lists$category[i], "_drug_enrich_result.xlsx")
  PRF_erich_df = PRF_drug_enrichment@result
  PRF_erich_df = cbind("Category" = DrugBank_drugs_known_target_lists$category[i], PRF_erich_df)
  colnames(PRF_erich_df) = c("Category",
                             "ID",
                             "Description",
                             "ClassRatio",
                             "BgRatio",
                             "pvalue" ,
                             "p.adjust",
                             "qvalue",
                             "DrugName",
                             "Count")
  openxlsx::write.xlsx(x = PRF_erich_df, file = PRF_file_name_xlsx, overwrite = TRUE)
  
  if (nrow(summary(PRF_drug_enrichment)) < 1){
    writeLines(paste0(DrugBank_drugs_known_target_lists$category[i], " has NO significant results"))
    next
  }
  
  # Making figures
  PRF_file_name = paste0(folders_drug_enrich[i], "/", DrugBank_drugs_known_target_lists$category[i], "_drug_enrich_plot")
  
  # Barplots
  pdf(paste0(PRF_file_name, "_bar.pdf")
      , width = 10, height = 10)
  print(barplot(PRF_drug_enrichment, showCategory = 10, font.size = 10))
  dev.off()
  
  pdf(paste0(PRF_file_name, "_bar_full.pdf")
      , width = 10, height = 10*nrow(summary(PRF_drug_enrichment))/10)
  print(barplot(PRF_drug_enrichment, showCategory = nrow(summary(PRF_drug_enrichment)), font.size = 10))
  dev.off()
  
  # Graph plots
  pdf(paste0(PRF_file_name, "_class_conc.pdf")
      , width = 21, height = 11)
  print(cnetplot(PRF_drug_enrichment, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                 showCategory = 10, cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
  dev.off()
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Making the combined drug enrichment file
Drug_Enrich_reports = list()

for (i in 1:length(folders_drug_enrich)){
  PRF_files = list.files(folders_drug_enrich[i])
  PRF_files = PRF_files[stri_detect_fixed(PRF_files, pattern = ".xlsx")]
  
  if (length(PRF_files) < 1){
    next
  }
  
  PRF_enrich_report = openxlsx::read.xlsx(paste0(folders_drug_enrich[i], "/", PRF_files), sheet = 1)
  Drug_Enrich_reports[[i]] = PRF_enrich_report
}
Drug_Enrich_reports = Drug_Enrich_reports[sapply(Drug_Enrich_reports, is.data.frame)]
openxlsx::write.xlsx(x = Drug_Enrich_reports, 
                     file = "Drug_enrichments/Combined_drug_enrichments.xlsx", 
                     overwrite = TRUE)


################### Additional SNP statistics ###################

HOWARD_DM_SNPS_NAT = GWAS_Cat_Filtered_df[GWAS_Cat_Filtered_df$LINK == "www.ncbi.nlm.nih.gov/pubmed/30718901",]
HOWARD_DM_SNPS_NAT = HOWARD_DM_SNPS_NAT[HOWARD_DM_SNPS_NAT$`P-VALUE` < 5 * 10^-8,] # 92 unique SNPs with GWAS significance
HOWARD_DM_SNPS_NAT = HOWARD_DM_SNPS_NAT$SNPS

# Reproducible SNPs GWAS Cat´
GWAS_CAT_reprod = openxlsx::read.xlsx("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GWAS_Catalog/GWAS_CAT_Depress_specific_genome_sign_SNPs_stats.xlsx")
GWAS_CAT_reprod = GWAS_CAT_reprod[GWAS_CAT_reprod$Count > 1,]
table(HOWARD_DM_SNPS_NAT %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 30 (31%) in reproducible, 66 not in the reproducible


################### Mapping individual OMICs layers and overlaps to GTEx database ###################
# 1sr post-mortem data from GTEx portal (We primarily use the whole blood since it was used as the main sample source)
# GTEx_Analysis_v8_trans_eGenes_fdr05.txt trans-eQTLs mapped individually in each tissue, at FDR < 0.05

Trans_eQTLs_GTEx = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GTEx_portal/GTEx_Analysis_v8_trans_eGenes_fdr05.txt")
Cis_eQTLs_GTEx = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GTEx_portal/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt")

# Checking IDs 
Trans_eQTLs_GTEx$variant_id %in% Cis_eQTLs_GTEx$variant_id # Only 5 SNPs overlap in both lists...
table(Trans_eQTLs_GTEx$variant_id %in% Cis_eQTLs_GTEx$variant_id)

# Mapping Trans_eQTLs_GTEx to dbSNP ids, 9 variants were not mapped to dbSNP153
Trans_eQTLs_GTEx$dbSNP153id = sapply(Trans_eQTLs_GTEx$variant_id, function(x){
  
  x = map_SNPs_GTEx_dbSNP_b38(PRF_SNP_ID_GTEx = x, PRF_lag_param = 5)
  return(x)
  
})

# Mappings SNPs to depression GWAS results
Trans_eQTLs_GTEx_depr = Trans_eQTLs_GTEx[Trans_eQTLs_GTEx$dbSNP153id %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name,] # None were mapped to depression SNPs from GWAS Cat
Cis_eQTLs_GTEx_depr = Cis_eQTLs_GTEx[Cis_eQTLs_GTEx$rs_id_dbSNP151_GRCh38p7 %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name, ] # 6 SNPs were mapped
Cis_eQTLs_GTEx_depr$Gene_SNP = sapply(Cis_eQTLs_GTEx_depr$rs_id_dbSNP151_GRCh38p7, function(x){
  x = GWAS_CAT_Depress_specific_SNPs_153_mapped[GWAS_CAT_Depress_specific_SNPs_153_mapped$name == x, "Mapped_gene_fixed"]
  return(x)
})
check_gene_symbol_NIH(Cis_eQTLs_GTEx_depr$gene_name, PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded)
c("RP5-1115A15","PRKAR2A","MST1R","BTN2A2","AKR1C3","DNASE1L2")

Cis_eQTLs_GTEx_depr$gene_name %in% Genes_matching_dir_expression # No overlap with genes matching direction expression

# Cis eQTLs in the brain
files = list.files("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GTEx_portal/GTEx_Analysis_v8_eQTL")
files = files[stri_detect_fixed(files, pattern = "Brain")]
files = files[stri_detect_fixed(files, pattern = "egenes")]
brain_names = sapply(files, function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[1]
  return(x)
})
files_paths = paste0("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GTEx_portal/GTEx_Analysis_v8_eQTL/", files)

brain_eQTLS_GTEx = list()
for ( i in 1:length(brain_names)){
  PRF_df = smart_fread(files_paths[i])
  PRF_df = cbind(Tissue = brain_names[i], PRF_df)
  brain_eQTLS_GTEx[[i]] = PRF_df
}
brain_eQTLS_GTEx = do.call(rbind, brain_eQTLS_GTEx)

brain_eQTLS_GTEx_depr = brain_eQTLS_GTEx[brain_eQTLS_GTEx$rs_id_dbSNP151_GRCh38p7 %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name, ]
table(brain_eQTLS_GTEx_depr$gene_name %in% Genes_matching_dir_expression) # 2 genes were mapped with matching dir in blood
repeated_brain_eQTLs_GTEx = as.data.frame(table(brain_eQTLS_GTEx_depr$rs_id_dbSNP151_GRCh38p7))
repeated_brain_eQTLs_GTEx = repeated_brain_eQTLs_GTEx[repeated_brain_eQTLs_GTEx$Freq > 1,]

# Exploring mapped eQTL and gene associations
brain_eQTLS_GTEx_depr_mapped = brain_eQTLS_GTEx_depr[brain_eQTLS_GTEx_depr$gene_name %in% Genes_matching_dir_expression, ]
Combined_expression_brain_eQTLS_GTEx_depr_mapped = Combined_expression_df[Combined_expression_df$Gene %in% brain_eQTLS_GTEx_depr_mapped$gene_name,]

# Stats blood
table(Cis_eQTLs_GTEx_depr$rs_id_dbSNP151_GRCh38p7 %in% HOWARD_DM_SNPS_NAT) # 4 not in Howard DM, 2 in Howard DM, 6 total "rs301799"   "rs13084037"
table(Cis_eQTLs_GTEx_depr$rs_id_dbSNP151_GRCh38p7 %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 6 not in reprod

# Stats brain
table(repeated_brain_eQTLs_GTEx$Var1 %in% HOWARD_DM_SNPS_NAT) # 18 not in Howard DM, 3 in Howard DM, 21 total rs12624433 rs13084037 rs1933802 
table(repeated_brain_eQTLs_GTEx$Var1 %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 17 not in reprod, 4 in reprod, 21 total rs12624433 rs1933802  rs2721811  rs9517313 
table(unique(brain_eQTLS_GTEx_depr$rs_id_dbSNP151_GRCh38p7) %in% HOWARD_DM_SNPS_NAT) # 7 in Howard DM, 74 not
table(unique(brain_eQTLS_GTEx_depr$rs_id_dbSNP151_GRCh38p7) %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 8 in reprod, 73 not


################### Mapping individual OMICs layers and overlaps to BIOS QTL browser ###################

# Importing BIOS QTL browser data
Cis_meQTL_BIOS = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/BIOS_qtl/2015_09_02_cis_meQTLsFDR0.05-CpGLevel.txt")
Trans_meQTL_BIOS = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/BIOS_qtl/2015_09_02_trans_meQTLsFDR0.05-CpGLevel.txt")
Cis_eQTM_BIOS =  smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/BIOS_qtl/2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt")
gene_eQTLs_BIOS =  smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/BIOS_qtl/gene_level_eQTLs.txt", sep = "\t")

# Start with eQTLs
gene_eQTLs_BIOS_depr = gene_eQTLs_BIOS[gene_eQTLs_BIOS$SNPName %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name, ]
gene_eQTLs_BIOS_depr_expanded = multiple_expander(df = gene_eQTLs_BIOS_depr, cols_to_expand = c(5,17), pattern = ",")
gene_eQTLs_BIOS_depr_snps = gene_eQTLs_BIOS_depr$SNPName
gene_eQTLs_BIOS_depr_snps = unique(gene_eQTLs_BIOS_depr_snps) 
length(gene_eQTLs_BIOS_depr_snps) # 599 unique SNPs are eQTLs

associated_genes_eQTLs_BIOS_depr = gene_eQTLs_BIOS_depr$HGNCName
associated_genes_eQTLs_BIOS_depr = unlist(stri_split_fixed(associated_genes_eQTLs_BIOS_depr, pattern = ","))
associated_genes_eQTLs_BIOS_depr = unique(associated_genes_eQTLs_BIOS_depr)
length(associated_genes_eQTLs_BIOS_depr) # 1142 unique genes

# Matching to transcriptome
associated_genes_eQTLs_BIOS_depr[associated_genes_eQTLs_BIOS_depr %in% Genes_matching_dir_expression]
table(associated_genes_eQTLs_BIOS_depr %in% Genes_matching_dir_expression) # 35 genes are among the ones with matching transcriptome direction
eQTLs_BIOS_depr_matched = gene_eQTLs_BIOS_depr[gene_eQTLs_BIOS_depr$HGNCName %in% Genes_matching_dir_expression,]
Combined_expression_brain_eQTLs_BIOS_depr_matched = Combined_expression_df[Combined_expression_df$Gene %in% eQTLs_BIOS_depr_matched$HGNCName,]

associated_genes_eQTLs_BIOS_depr[associated_genes_eQTLs_BIOS_depr %in% Overlapping_genes_gwas_transcript] # 7 genes are in the overlap "ACTR5"  "METTL9" "RNF24"  "HIP1"   "PRRC2A" "RTN4"   "MRPL37"

# Matching to reprod SNPs and  Howard DM
table(gene_eQTLs_BIOS_depr_snps %in% HOWARD_DM_SNPS_NAT) # 561 not in Howard DM, 38 in Howard DM
table(gene_eQTLs_BIOS_depr_snps %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 564 not in reproducible, 35 in reproducible

# Reviewing DNA methylation associations

# Trans meQTL
Trans_meQTL_BIOS_depr = Trans_meQTL_BIOS[Trans_meQTL_BIOS$SNPName %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name, ] # 99 in total
Trans_meQTL_BIOS_depr_SNPs = Trans_meQTL_BIOS_depr$SNPName
length(unique(Trans_meQTL_BIOS_depr_SNPs)) # 13 unique SNPs affecting DNA methylation at FDR
Trans_meQTL_BIOS_depr_CpGs = Trans_meQTL_BIOS_depr$ProbeName
length(unique(Trans_meQTL_BIOS_depr_CpGs))  # 96 unique CpGs

# Matching to reproducible CpGs
table(Trans_meQTL_BIOS_depr_CpGs %in% CpG_table_methyl_depr_Overlap_all) # No CpGs were detected in the overlap

# Matching to reprod SNPs and  Howard DM
table(unique(Trans_meQTL_BIOS_depr_SNPs) %in% HOWARD_DM_SNPS_NAT) # 12 not in Howard DM, 1 in Howard DM
table(unique(Trans_meQTL_BIOS_depr_SNPs) %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 10 not in reproducible, 3 in reproducible

# Cis meQTL
Cis_meQTL_BIOS_depr = Cis_meQTL_BIOS[Cis_meQTL_BIOS$SNPName %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name, ] # 211 in total
Cis_meQTL_BIOS_depr_SNPs = Cis_meQTL_BIOS_depr$SNPName
length(unique(Cis_meQTL_BIOS_depr_SNPs)) # 121 unique SNPs affecting DNA methylation at FDR
Cis_meQTL_BIOS_depr_CpGs = Cis_meQTL_BIOS_depr$ProbeName
length(unique(Cis_meQTL_BIOS_depr_CpGs))  # 210 unique CpGs

# Matching to reproducible CpGs
table(Cis_meQTL_BIOS_depr_CpGs %in% CpG_table_methyl_depr_Overlap_all) # 209 not in reproducible, 2 in reproducible "cg26688911" "cg07325168"
Cis_meQTL_BIOS_depr_CpGs[Cis_meQTL_BIOS_depr_CpGs %in% CpG_table_methyl_depr_Overlap_all]
Cis_meQTL_BIOS_depr_matched = Cis_meQTL_BIOS_depr[Cis_meQTL_BIOS_depr$ProbeName %in% CpG_table_methyl_depr_Overlap_all,]
Cis_meQTL_BIOS_depr_matched$SNP_gene = sapply(Cis_meQTL_BIOS_depr_matched$SNPName, function(x){
  x = GWAS_CAT_Depress_specific_SNPs_153_mapped[GWAS_CAT_Depress_specific_SNPs_153_mapped$name == x, "Mapped_gene_fixed"]
  return(x)
})
Cis_meQTL_BIOS_depr_matched$CpG_gene = sapply(Cis_meQTL_BIOS_depr_matched$ProbeName, function(x){
  x = Combined_methyl_df_significant_unique[Combined_methyl_df_significant_unique$CpG == x, "Upd_gene_name"]
  return(x)
})
Cis_meQTL_BIOS_depr_matched$SNP_gene == Cis_meQTL_BIOS_depr_matched$CpG_gene # Both match "LPIN3" "VPS41"
Cis_meQTL_BIOS_depr_matched$SNP_gene %in% Intersect_GWAS_Methyl # No overlaps
"LPIN3" %in% Df_methyl_Signif_Overlap_all$Upd_gene_name
"VPS41" %in% Intersect_ALL # This gene is in triple intersect

# Matching to reproducible SNPs and Howard DM
table(unique(Cis_meQTL_BIOS_depr_SNPs) %in% HOWARD_DM_SNPS_NAT) # 113 not in Howard DM, 8 in Howard DM
table(unique(Cis_meQTL_BIOS_depr_SNPs) %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 111 not in reproducible, 10 in reprod
table(unique(Cis_meQTL_BIOS_depr_matched$SNPName) %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`) # 1 not in reproducible, 1 in reprod

# Cis eQTM mapping
Cis_eQTM_BIOS_depr = Cis_eQTM_BIOS[Cis_eQTM_BIOS$SNPName %in% CpG_table_methyl_depr_Overlap_all,] # it is not an SNP, 112 pairs, 85 unique SNPs
Cis_eQTM_BIOS_depr_genes = Cis_eQTM_BIOS_depr$HGNCName
Cis_eQTM_BIOS_depr_genes = unique(Cis_eQTM_BIOS_depr_genes)

# meQTL relation
table(Cis_meQTL_BIOS_depr_CpGs %in% Cis_eQTM_BIOS_depr$SNPName) # No overlaps between meQTL and Cis eQTM
c("cg26688911", "cg07325168") %in% Cis_eQTM_BIOS_depr$SNPName

# Matching to transcriptome
table(Cis_eQTM_BIOS_depr_genes %in% Genes_matching_dir_expression) # One gene overlaps with transcriptome
Cis_eQTM_BIOS_depr_genes[Cis_eQTM_BIOS_depr_genes %in% Genes_matching_dir_expression] # "PCID2"
Cis_eQTM_BIOS_depr_genes[Cis_eQTM_BIOS_depr_genes %in% Overlapping_genes_methyl_transcript] # None of them overlap with methyl-transcriptome
Cis_eQTM_BIOS_depr_genes_mapped = Cis_eQTM_BIOS_depr[Cis_eQTM_BIOS_depr$HGNCName %in% Genes_matching_dir_expression, ]
Cis_eQTM_BIOS_depr_genes_mapped_methylation = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% Cis_eQTM_BIOS_depr_genes_mapped$SNPName,]

################### cis-eQTM in GSE56046 ###################
all_eQTMs = smart_fread("eQTM_GSE56046/all_eQTMs.csv")

all_eQTMs$FDR = p.adjust(all_eQTMs$`Pr(>|t|)`, method = "fdr")
all_eQTMs_FDR = all_eQTMs[all_eQTMs$FDR < 0.05,]
all_eQTMs_FDR_depress =  all_eQTMs_FDR[all_eQTMs_FDR$CpG %in% CpG_table_methyl_depr_Overlap_all,] 
Cis_eQTM_GSE56046_depr_genes = unique(all_eQTMs_FDR_depress$PRF_Gene_symbol)

# Matching to transcriptome
table(Cis_eQTM_GSE56046_depr_genes %in% Genes_matching_dir_expression) # 12 genes overlap with transcriptome
Cis_eQTM_GSE56046_depr_genes[Cis_eQTM_GSE56046_depr_genes %in% Genes_matching_dir_expression]
# "PPFIA1"   "CERT1"    "DENND5A"  "TNFRSF1B" "LIN7A"    "ARHGAP27" "NUP43"    "RAB37"    "SLC16A5"  "ELANE"    "LBR"      "SLC11A1" 
Cis_eQTM_GSE56046_depr_genes[Cis_eQTM_GSE56046_depr_genes %in% Overlapping_genes_methyl_transcript] # 2 genes overlap with methyl-transcriptome
# "ARHGAP27" "RAB37"
Cis_eQTM_GSE56046_depr_genes_mapped = all_eQTMs_FDR_depress[all_eQTMs_FDR_depress$PRF_Gene_symbol %in% Genes_matching_dir_expression, ]
Cis_eQTM_GSE56046_depr_genes_mapped_methylation = Combined_methyl_df_significant[Combined_methyl_df_significant$CpG %in% all_eQTMs_FDR_depress$CpG,]

################### Network visualization ###################

# Preparing datasets

# eQTL GTEx
brain_eQTLS_GTEx_depr_mapped

# eQTL BIOS
gene_eQTLs_BIOS_depr_matched = gene_eQTLs_BIOS_depr_expanded[gene_eQTLs_BIOS_depr_expanded$HGNCName %in% Genes_matching_dir_expression,]
all(gene_eQTLs_BIOS_depr_matched$SNPName %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name) # All are among depression SNPs

# cis meQTL BIOS
Cis_meQTL_BIOS_depr_matched

# cis eQTM BIOS
Cis_eQTM_BIOS_depr_genes_mapped

# cis eQTM GSE56046
Cis_eQTM_GSE56046_depr_genes_mapped

# Making graphs
eQTL_GTEx_graph = brain_eQTLS_GTEx_depr_mapped[,c("rs_id_dbSNP151_GRCh38p7", "gene_name")]
colnames(eQTL_GTEx_graph) = c("Start", "End")

eQTL_BIOS_graph = gene_eQTLs_BIOS_depr_matched[,c("SNPName", "HGNCName")]
colnames(eQTL_BIOS_graph) = c("Start", "End")

meQTL_graph = Cis_meQTL_BIOS_depr_matched[,c("SNPName", "ProbeName")]
colnames(meQTL_graph) = c("Start", "End")

eQTM_graph = Cis_eQTM_BIOS_depr_genes_mapped[,c("SNPName", "HGNCName")]
colnames(eQTM_graph) = c("Start", "End")

eQTM_graph_GSE56046 = Cis_eQTM_GSE56046_depr_genes_mapped[,c("CpG", "PRF_Gene_symbol")]
colnames(eQTM_graph_GSE56046) = c("Start", "End")

Graph_DS = rbind(eQTL_GTEx_graph, 
                 eQTL_BIOS_graph,
                 meQTL_graph,
                 eQTM_graph,
                 eQTM_graph_GSE56046)
network_data = graph_from_data_frame(d = Graph_DS[,1:2], directed = TRUE)
network_data = toVisNetworkData(network_data)

# Setting nodes
network_data$nodes$group = sapply(network_data$nodes$id, function(x){
  if (x %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name & x %!in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`){
    return("SNP")
    
  } else if (x %in% GWAS_CAT_Depress_specific_SNPs_153_mapped$name & x %in% GWAS_CAT_reprod$`Genome-wide.significant.SNP`){
    return("Reprod.SNP")
    
  } else if (x %in% CpG_table_methyl_depr_Overlap_all & x %!in% CpG_table_methyl_depr_Overlap){
    return("CpG")
    
  } else if (x %in% CpG_table_methyl_depr_Overlap_all & x %in% CpG_table_methyl_depr_Overlap){
    return("CpG.Matching.dir")
    
  } else if (x %in% Intersect_ALL){
      return("Gene intersect")
    
  } else {
    return("Transcript")
  }
  
})

# Setting edges
network_data$edges$color = NA
network_data$edges$label = NA
network_data$edges$arrows = "to"
network_data$edges$font.color = NA
network_data$edges$font.size = 15

for (i in 1:nrow(network_data$edges)){
  PRF_Start = network_data$edges$from[i]
  PRF_End = network_data$edges$to[i]
  
  PRF_gr_Start = network_data$nodes[network_data$nodes$id == network_data$edges$from[i], "group" ]
  PRF_gr_End = network_data$nodes[network_data$nodes$id == network_data$edges$to[i], "group" ]
  
  if (PRF_gr_Start %in% c("SNP", "Reprod.SNP") & PRF_gr_End %in% c("CpG", "CpG.Matching.dir")){
    network_data$edges$color[i] = "blue"
    network_data$edges$font.color[i] = "blue"
    network_data$edges$label[i] = "meQTL"
  } else if (PRF_gr_Start %in% c("SNP", "Reprod.SNP") & PRF_gr_End %in% c("Gene intersect", "Transcript")){
    network_data$edges$color[i] = "red"
    network_data$edges$font.color[i] = "red"
    network_data$edges$label[i] = "eQTL"
  } else if (PRF_gr_Start %in% c("CpG", "CpG.Matching.dir") & PRF_gr_End %in% c("Gene intersect", "Transcript")){
    network_data$edges$color[i] = "brown"
    network_data$edges$font.color[i] = "brown"
    network_data$edges$label[i] = "eQTM"
    
    # hilhight GSE56046
    if (PRF_Start %in% eQTM_graph_GSE56046$Start & PRF_End %in% eQTM_graph_GSE56046$End){
      network_data$edges$color[i] = "green"
      network_data$edges$font.color[i] = "green"
      network_data$edges$label[i] = "eQTM GSE56046"
    }
    
  }
}

net = visNetwork(nodes = network_data$nodes, edges = network_data$edges, height = "2000px", width = "2000px") %>%
  visIgraphLayout(layout = "layout_with_fr", physics = FALSE, randomSeed = 190) %>%
  visGroups(groupname = "SNP", color = list("background" = "skyblue", border = "black")) %>% 
  visGroups(groupname = "Reprod.SNP", color = list("background" = "green", border = "black")) %>% 
  visGroups(groupname = "CpG", color = list("background" = "#FBCEB1", border = "black")) %>% 
  visGroups(groupname = "CpG.Matching.dir", color = list("background" = "violet", border = "black")) %>% 
  visGroups(groupname = "Gene intersect", color = list("background" = "red", border = "black")) %>% 
  visGroups(groupname = "Transcript", color = list("background" = "orange", border = "black")) %>% 
  visNodes(size = 15, borderWidth = 0.3, font = list(size = 29)) %>%
  visEdges(width = 1.5, smooth = list(type = "curvedCW", roundness = 0.1)) %>%
  visPhysics(barnesHut = list(gravitationalConstant = -100000, avoidOverlap = 0.9)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)
visSave(net, file = "mapped_hits.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("mapped_hits.html", "mapped_hits.png", vwidth = 2000, vheight = 2000, zoom = 1)

# Writing outputs

QTL_analysis_list = list(
  "blood_eQTLS_GTEx" = Cis_eQTLs_GTEx_depr,
  "brain_eQTLS_GTEx" = brain_eQTLS_GTEx_depr,
  "brain_eQTLS_GTEx_matched" = brain_eQTLS_GTEx_depr_mapped,
  "eQTLs_BIOS" = gene_eQTLs_BIOS_depr,
  "eQTLs_BIOS_matched" = gene_eQTLs_BIOS_depr_matched,
  "cis_meQTL_BIOS_depr" = Cis_meQTL_BIOS_depr,
  "cis_meQTL_BIOS_depr_matched" = Cis_meQTL_BIOS_depr_matched,
  "blood_eQTMs_BIOS" = Cis_eQTM_BIOS_depr,
  "blood_eQTMs_BIOS_matched" = Cis_eQTM_BIOS_depr_genes_mapped,
  "CD14_eQTMs_GSE56046" = all_eQTMs_FDR_depress,
  "CD14_eQTMs_GSE56046_matched" = Cis_eQTM_GSE56046_depr_genes_mapped
)
write.xlsx(QTL_analysis_list,"QTL_analysis_list.xlsx", overwrite = TRUE)


################### Generating aggreagted data for matched genes ###################
Intersect_ALL

# GWAS
selection = sapply(GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed, function(x){
  x = str_trim(unlist(stri_split_fixed(x, pattern = ";")))
  index = x %in% Intersect_ALL
  index = any(index)
  return(index)
})
GWAS_intersect_all = GWAS_CAT_Depress_specific_SNPs_153_mapped[selection,]
selection_studies = sapply(GWAS_Cat_Filtered_df$SNPS, function(x){
  x = unlist(stri_split_fixed(x, pattern = "; "))
  x = unlist(stri_split_fixed(x, pattern = ", "))
  x = x[!is.na(x)]
  
  if (length(x) < 1){
    return(FALSE)
  }
  
  x = sapply(x, function(z){
    
    if (stri_detect_fixed(z, pattern = "-")){
      z = unlist(stri_split_fixed(z, pattern = "-"))
      z = z[which.max(nchar(z))]
    }
    
    if (stri_detect_fixed(z, pattern = "_")){
      z = unlist(stri_split_fixed(z, pattern = "_"))
      z = z[which.max(nchar(z))]
    }
    return(z)
  })
  index = x %in% GWAS_intersect_all$name
  index = any(index)
  return(index)
})
GWAS_intersect_studies = GWAS_Cat_Filtered_df[selection_studies,]
GWAS_intersect_studies_filtered = GWAS_intersect_studies[,c("PUBMEDID",
                                                            "FIRST AUTHOR",
                                                            "DISEASE/TRAIT",
                                                            "INITIAL SAMPLE SIZE",
                                                            "REPLICATION SAMPLE SIZE",
                                                            "SNPS",
                                                            "MAPPED_GENE",
                                                            "P-VALUE",
                                                            "OR or BETA")]
colnames(GWAS_intersect_studies_filtered) = c(
  "PMID",
  "First author",
  "DISEASE/TRAIT",
  "Initial sample size",
  "Replication sample size",
  "SNP",
  "Mapped Gene",
  "P-value",
  "OR or BETA"
)

# DNA methylation
DNAm_index = sapply(Df_methyl_Signif_Overlap$Upd_gene_name, function(x){
  x = x[!is.na(x)]
  x = x[x != "NA"]
  x = unlist(stri_split_fixed(x, pattern = ";"))
  x = unique(x)
  x = x[!is.na(x)]
  if (length(x) < 1){
    return(FALSE)
  }
  index = x %in% Intersect_ALL
  index = any(index)
  return(index)
})

Triple_intersect_DNA_methyl = Df_methyl_Signif_Overlap[DNAm_index,]
Triple_intersect_DNA_methyl = arrange(Triple_intersect_DNA_methyl, CpG)
Triple_intersect_DNA_methyl_selected = Triple_intersect_DNA_methyl[,c(
  "Contrast",
  "Upd_gene_name",
  "CpG",
  "Relation_to_Gene",
  "logFC",
  "P.Value"
)]
colnames(Triple_intersect_DNA_methyl_selected) = c(
  "Contrast",
  "Gene",
  "CpG",
  "Relation_to_Gene",
  "logFC",
  "P.Value.nominal"
)

# Transcriptome
Triple_transcriptome = Combined_expression_df[Combined_expression_df$Gene %in% Intersect_ALL, ]
Triple_transcriptome = arrange(Triple_transcriptome, ID)
Triple_transcriptome_selected = Triple_transcriptome[, c(
  "Contrast",
  "Gene",
  "ID",
  "logFC",
  "P.Value"
)]
colnames(Triple_transcriptome_selected) = c(
  "Contrast",
  "Gene",
  "Probe ID",
  "logFC",
  "P.Value.nominal"
)

# QTL data
Triple_QTL_data =  gene_eQTLs_BIOS_depr_matched[gene_eQTLs_BIOS_depr_matched$HGNCName %in% Intersect_ALL,]
Triple_QTL_data = Triple_QTL_data[,c(
  "SNPName",
  "ProbeName",
  "HGNCName",
  "CisTrans",
  "PValue",
  "Beta (SE)",
  "FDR"
)]


# Writing files
openxlsx::write.xlsx(GWAS_intersect_studies_filtered, "GWAS_intersect_studies_filtered.xlsx", overwrite = TRUE)
openxlsx::write.xlsx(Triple_intersect_DNA_methyl_selected, "Triple_intersect_DNA_methyl_selected.xlsx", overwrite = TRUE)
openxlsx::write.xlsx(Triple_transcriptome_selected, "Triple_transcriptome_selected.xlsx", overwrite = TRUE)
openxlsx::write.xlsx(Triple_QTL_data, "Triple_QTL_data.xlsx", overwrite = TRUE)

