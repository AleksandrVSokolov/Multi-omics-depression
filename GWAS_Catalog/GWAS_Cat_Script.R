
################# This script shows analysis for the GWAS Catalog database ################# 
# Note 1: The equal sign = was used as an assignment operator as authors don't buy the idea of using <- for typing/productivity reasons
# Note 2: In many cases loops were deliberately used instead of apply functions to enable better control of the variables (even though loops in R are slow and computationally inefficient)
# Note 3: Some variables in the loops contain prefixes to enable easy cleanup of the environment once the loop is executed
# Note 4: To execute the script, several raw files and data from UCSC and NIH should be downloaded 

Working_directory = "~/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GWAS_Catalog" # Replace with an appropriate path
setwd(Working_directory)

getOption("scipen") # default number notation is 0
options(scipen=999)


################### Package import ###################
# importing packages that may be nessesary
library("fun")
library("stringr")
library("dplyr")
library("ggplot2")
library("openxlsx")
library("grid")
library("gdata")
library("RColorBrewer")
library("networkD3")
library("webshot")
library("htmlwidgets")
library("magrittr")
library("igraph")
library("visNetwork")
library("data.table")
library("XML")
library("rvest")
library("RCurl")
library("HGNChelper")
library("stringi")
library("httr")
library("lubridate")
library("rjson")
library("rtracklayer")
library("rstudioapi")
library("tidyr")
library("Gviz")
library("limma")
library("FactoMineR")
library("ggthemes")
library("igraph")
library("RSelenium")
library("lumi")
library("outliers")
library("svglite")
library("scatterplot3d")
library("sva")
library("jsonlite")
library("ggrepel")
library("parallel")
library("bacon")
library("gridExtra")
library("grid")
library("ggplotify")
library("clusterProfiler")
library("DOSE")
library("enrichplot")
library("chromoMap")
library("UniprotR")
library("RIdeogram")


################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 14, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

# function to detect at least one pattern in a string
multiple_stri_detector = function(string, pattern_vector){
  output_list = list()
  for (i in 1:length(pattern_vector)){
    output_list[[i]] = stri_detect_fixed(str = string, pattern = pattern_vector[i])
  }
  output_list = do.call(rbind, output_list)
  apply(output_list, 2, any)
}


# function to expand a data frame where several columns contain condensed cells with a specified separator
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

# function to map SNPs to dbSNP153. It uses bigBedNamedItems from UCSC Genome Browser. The function was not tested on Windows-based OS
# bigBedNamedItems has to be downloaded and activated on a computer
# paths specified in the function should be replaced with appropriate paths for the current computer
lin_req_SNP_db_153 = function(SNPs, File_name, lag_param){
  writeLines(text = SNPs, con = "/home/aleksandr/Desktop/WORK/SNP_CHECK/myIds.txt")
  Info_string = "If the first use, please type the following:\nchmod +x /home/aleksandr/UCSC_folder/bigBedNamedItems\n\n"
  writeLines(Info_string)
  Base_string = "/home/aleksandr/UCSC_folder/bigBedNamedItems -nameFile /home/aleksandr/Desktop/WORK/SNP_CHECK/dbSnp153.bb /home/aleksandr/Desktop/WORK/SNP_CHECK/myIds.txt /home/aleksandr/Desktop/WORK/SNP_CHECK/"
  Request_linux = paste0(Base_string, File_name)
  terminalExecute(Request_linux)
  Sys.sleep(lag_param)
  File_name_2 = paste0("/home/aleksandr/Desktop/WORK/SNP_CHECK/",File_name)
  Data = fread(File_name_2)
  colnames(Data) = c("chrom", "chromStart", "chromEnd", "name",
                     "ref", "altCount", "alts", "shiftBases", "freqSourceCount",
                     "minorAlleleFreq", "majorAllele", "minorAllele", "maxFuncImpact",
                     "class", "ucscNotes", "_dataOffset", "_dataLen")
  return(Data)
}

# function to get coordinate table for human chromosomes (hg19 genome build)
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

# this function add mapped genes to UCSC SNP153 track
# mapping is performed directly by gene coordinates obtained from UCSC tracks 
gene_mapper = function(SNP_153_Track, gene_track){
  
  SNP_153_Track$Mapped_transcript = NA
  SNP_153_Track$Mapped_gene = NA
  
  for (i in 1:nrow(SNP_153_Track)){
    
    Chr = SNP_153_Track$chrom[i]
    Curr_track_genes = gene_track[gene_track$chrom == Chr,]
    Curr_mapped_transcript = vector()
    Curr_mapped_genes = vector()
    
    for (j in 1:nrow(Curr_track_genes)){
      Curr_gene_coord = Curr_track_genes$txStart[j]:Curr_track_genes$txEnd[j]
      if (SNP_153_Track$chromEnd[i] %in% Curr_gene_coord){
        Curr_mapped_genes = c(Curr_mapped_genes, Curr_track_genes$name2[j])
        Curr_mapped_transcript = c(Curr_mapped_transcript, Curr_track_genes$name[j])
      }
    }
    
    if (length(Curr_mapped_genes) >=1){
      Curr_mapped_genes = paste0(Curr_mapped_genes, collapse = ";")
      Curr_mapped_transcript = paste0(Curr_mapped_transcript, collapse = ";")
      SNP_153_Track$Mapped_transcript[i] = Curr_mapped_transcript
      SNP_153_Track$Mapped_gene[i] = Curr_mapped_genes
    }
    
    Percent = round(i/nrow(SNP_153_Track), digits = 4)*100
    writeLines(paste0(i, "   ", Percent, "_%_complete_at_", Sys.time()))
  }
  return(SNP_153_Track)
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


################### Curating of the GWAS Catalog database ###################

# Importing files
GWAS_CAT_file = "gwas_catalog_v1.0-associations_e105_r2022-04-07.tsv" # Replace with another path if needed
GWAS_CAT_INIT = smart_fread(GWAS_CAT_file)
GWAS_CAT_INIT$`DISEASE/TRAIT` = toupper(GWAS_CAT_INIT$`DISEASE/TRAIT`)
GWAS_CAT_Depress = GWAS_CAT_INIT[stri_detect_fixed(GWAS_CAT_INIT$`DISEASE/TRAIT`, pattern = "DEPRESS"),]
GWAS_CAT_Depress = GWAS_CAT_Depress[!stri_detect_fixed(GWAS_CAT_Depress$`DISEASE/TRAIT`, pattern = "BIPOLAR"),]

# Depression-specific dataset
GWAS_CAT_Depress_specific = GWAS_CAT_Depress[!sapply(GWAS_CAT_Depress$`DISEASE/TRAIT`, function(x){
  x = multiple_stri_detector(string = x, pattern_vector = c("DEPRESSIVE SYMPTOMS X DEPENDENT STRESSFUL LIFE EVENTS INTERACTION",
                                                            "DEPRESSIVE SYMPTOMS X INDEPENDENT STRESSFUL LIFE EVENTS INTERACTION",
                                                            "SEVERE DEPRESSIVE DISORDERS IN CORONARY ARTERY DISEASE",
                                                            "BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                                            "DEPRESSIVE SYMPTOMS X HERPES SIMPLEX 1 INFECTION INTERACTION",
                                                            "COGNITIVE IMPAIRMENT IN DEPRESSION",
                                                            "COGNITIVE DECLINE IN DEPRESSION",
                                                            "TRAUMA EXPOSURE IN MAJOR DEPRESSIVE DISORDER-NEGATIVE INDIVIDUALS",
                                                            "ANOREXIA NERVOSA, ATTENTION-DEFICIT/HYPERACTIVITY DISORDER, AUTISM SPECTRUM DISORDER, B",
                                                            "BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER AND DEPRESSION (COMBINED)",
                                                            "ANXIETY IN MAJOR DEPRESSIVE DISORDER",
                                                            "ANXIETY AND MAJOR DEPRESSIVE DISORDER",
                                                            "AUTISM SPECTRUM DISORDER, ATTENTION DEFICIT-HYPERACTIVITY DISORDER, BIPOLAR DISORDER, M",
                                                            "STRESS SENSITIVITY (NEUROTICISM SCORE X MAJOR DEPRESSIVE DISORDER STATUS INTERACTION)",
                                                            "RESPONSE TO KETAMINE IN BIPOLAR DISORDER OR MAJOR DEPRESSION",
                                                            "SUICIDE IDEATION SCORE IN MAJOR DEPRESSIVE DISORDER",
                                                            "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER",
                                                            "MAJOR DEPRESSIVE DISORDER (STRESSFUL LIFE EVENTS INTERACTION)",
                                                            "DEPRESSION IN RESPONSE TO INTERFERON-BASED THERAPY IN CHRONIC HEPATITIS C",
                                                            "BROAD DEPRESSION OR BIPOLAR DISORDER",
                                                            "BROAD DEPRESSION OR SCHIZOPHRENIA",
                                                            "SUICIDAL IDEATION IN DEPRESSION OR BIPOLAR DISORDER",
                                                            "SUICIDE ATTEMPTS IN DEPRESSION OR BIPOLAR DISORDER",
                                                            "FUNCTIONAL IMPAIRMENT IN MAJOR DEPRESSIVE DISORDER, BIPOLAR DISORDER AND SCHIZOPHRENIA",
                                                            "PHARMACOKINETICS OF ANTIDEPRESSANT DRUGS IN SEVERE MENTAL DISORDER (CONCENTRATION DO",
                                                            "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER",
                                                            "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER OR BIPOLAR DISORDER OR SCHIZOPHRENIA",
                                                            "DEPRESSIVE SYMPTOMS (STRESSFUL LIFE EVENTS INTERACTION)",
                                                            "QT INTERVAL (TRICYCLIC/TETRACYCLIC ANTIDEPRESSANT USE INTERACTION)",
                                                            "RR INTERVAL (TRICYCLIC/TETRACYCLIC ANTIDEPRESSANT USE INTERACTION)",
                                                            "DEPRESSIVE EPISODES IN BIPOLAR DISORDER",
                                                            "DEPRESSIVE AND MANIC EPISODES IN BIPOLAR DISORDER",
                                                            "SLEEP DISTURBANCE (INCREASED SLEEP) IN DEPRESSIVE DISORDERS",
                                                            "SLEEP DISTURBANCE (DECREASED SLEEP) IN DEPRESSIVE DISORDERS",
                                                            "DEPRESSION X VITAMIN D PRS INTERACTION (COJO ADJUSTED)",
                                                            "DEPRESSION SCORE X VITAMIN D INTERACTION",
                                                            "RESPONSE TO COGNITIVE-BEHAVIOURAL THERAPY IN ANXIETY AND MAJOR DEPRESSIVE DISORDERS",
                                                            "RESPONSE TO ANTIDEPRESSANTS (PAIN RELIEF EFFICACY)",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER X SEX INTERACTION",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER X SEX INTERACTION (3DF)",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER OR RECURRENT MAJOR DEPRESSIVE DISORDER X SEX INTERACTI",
                                                            "ENDOMETRIOSIS OR DEPRESSION (PLEIOTROPY)",
                                                            "BIPOLAR DISORDER VS MAJOR DEPRESSIVE DISORDER (ORDINARY LEAST SQUARES (OLS))",
                                                            "MAJOR DEPRESSIVE DISORDER VS ADHD (ORDINARY LEAST SQUARES (OLS))",
                                                            "MAJOR DEPRESSIVE DISORDER VS ANOREXIA NERVOSA (ORDINARY LEAST SQUARES (OLS)",
                                                            "MAJOR DEPRESSIVE DISORDER VS AUTISM SPECTRUM DISORDER (ORDINARY LEAST SQUARES",
                                                            "MAJOR DEPRESSIVE DISORDER OR STRESS-RELATED DISORDER",
                                                            "DEPRESSION X MATERNAL SMOKING DURING PREGNANCY INTERACTION",
                                                            "SCHIZOPHRENIA, BIPOLAR DISORDER OR RECURRENT MAJOR DEPRESSIVE DISORDER",
                                                            "MAJOR DEPRESSIVE DISORDER OR STROKE (PLEIOTROPY)",
                                                            "COGNITIVE FUNCTION (PROCESSING SPEED) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                                            "COGNITIVE FUNCTION (EXECUTIVE FUNCTION) X MAJOR DEPRESSIVE DISORDER INTERACTION (1DF)",
                                                            "COGNITIVE FUNCTION (EXECUTIVE FUNCTION) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                                            "COGNITIVE FUNCTION (DELAYED MEMORY) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                                            "COGNITIVE FUNCTION (PROCESSING SPEED) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                                            "COGNITIVE FUNCTION (GLOBAL COGNITION) X MAJOR DEPRESSIVE DISORDER INTERACTION (1DF)",
                                                            "DEPRESSION X HERPES SIMPLEX 1 INFECTION ",
                                                            "DEPRESSION IN MULTIPLE SCLEROSIS",
                                                            "MENSTRUATION QUALITY OF LIFE IMPACT (DEPRESS",
                                                            "HYPERSOMNIA DURING A MAJOR DEPRESSIVE EPIS",
                                                            "(MTAG)",
                                                            "VASCULAR ENDOTHELIAL GROWTH FACTOR LEVELS"
  ))
  return(x)
}),]
Exclusion_terms_depression =  c("DEPRESSIVE SYMPTOMS X DEPENDENT STRESSFUL LIFE EVENTS INTERACTION",
                                "DEPRESSIVE SYMPTOMS X INDEPENDENT STRESSFUL LIFE EVENTS INTERACTION",
                                "SEVERE DEPRESSIVE DISORDERS IN CORONARY ARTERY DISEASE",
                                "BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                "DEPRESSIVE SYMPTOMS X HERPES SIMPLEX 1 INFECTION INTERACTION",
                                "COGNITIVE IMPAIRMENT IN DEPRESSION",
                                "COGNITIVE DECLINE IN DEPRESSION",
                                "TRAUMA EXPOSURE IN MAJOR DEPRESSIVE DISORDER-NEGATIVE INDIVIDUALS",
                                "ANOREXIA NERVOSA, ATTENTION-DEFICIT/HYPERACTIVITY DISORDER, AUTISM SPECTRUM DISORDER, B",
                                "BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER AND DEPRESSION (COMBINED)",
                                "ANXIETY IN MAJOR DEPRESSIVE DISORDER",
                                "ANXIETY AND MAJOR DEPRESSIVE DISORDER",
                                "AUTISM SPECTRUM DISORDER, ATTENTION DEFICIT-HYPERACTIVITY DISORDER, BIPOLAR DISORDER, M",
                                "STRESS SENSITIVITY (NEUROTICISM SCORE X MAJOR DEPRESSIVE DISORDER STATUS INTERACTION)",
                                "RESPONSE TO KETAMINE IN BIPOLAR DISORDER OR MAJOR DEPRESSION",
                                "SUICIDE IDEATION SCORE IN MAJOR DEPRESSIVE DISORDER",
                                "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER",
                                "MAJOR DEPRESSIVE DISORDER (STRESSFUL LIFE EVENTS INTERACTION)",
                                "DEPRESSION IN RESPONSE TO INTERFERON-BASED THERAPY IN CHRONIC HEPATITIS C",
                                "BROAD DEPRESSION OR BIPOLAR DISORDER",
                                "BROAD DEPRESSION OR SCHIZOPHRENIA",
                                "SUICIDAL IDEATION IN DEPRESSION OR BIPOLAR DISORDER",
                                "SUICIDE ATTEMPTS IN DEPRESSION OR BIPOLAR DISORDER",
                                "FUNCTIONAL IMPAIRMENT IN MAJOR DEPRESSIVE DISORDER, BIPOLAR DISORDER AND SCHIZOPHRENIA",
                                "PHARMACOKINETICS OF ANTIDEPRESSANT DRUGS IN SEVERE MENTAL DISORDER (CONCENTRATION DO",
                                "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER",
                                "SUICIDE ATTEMPTS IN MAJOR DEPRESSIVE DISORDER OR BIPOLAR DISORDER OR SCHIZOPHRENIA",
                                "DEPRESSIVE SYMPTOMS (STRESSFUL LIFE EVENTS INTERACTION)",
                                "QT INTERVAL (TRICYCLIC/TETRACYCLIC ANTIDEPRESSANT USE INTERACTION)",
                                "RR INTERVAL (TRICYCLIC/TETRACYCLIC ANTIDEPRESSANT USE INTERACTION)",
                                "DEPRESSIVE EPISODES IN BIPOLAR DISORDER",
                                "DEPRESSIVE AND MANIC EPISODES IN BIPOLAR DISORDER",
                                "SLEEP DISTURBANCE (INCREASED SLEEP) IN DEPRESSIVE DISORDERS",
                                "SLEEP DISTURBANCE (DECREASED SLEEP) IN DEPRESSIVE DISORDERS",
                                "DEPRESSION X VITAMIN D PRS INTERACTION (COJO ADJUSTED)",
                                "DEPRESSION SCORE X VITAMIN D INTERACTION",
                                "RESPONSE TO COGNITIVE-BEHAVIOURAL THERAPY IN ANXIETY AND MAJOR DEPRESSIVE DISORDERS",
                                "RESPONSE TO ANTIDEPRESSANTS (PAIN RELIEF EFFICACY)",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER X SEX INTERACTION",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER OR MAJOR DEPRESSIVE DISORDER X SEX INTERACTION (3DF)",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER OR RECURRENT MAJOR DEPRESSIVE DISORDER X SEX INTERACTI",
                                "ENDOMETRIOSIS OR DEPRESSION (PLEIOTROPY)",
                                "BIPOLAR DISORDER VS MAJOR DEPRESSIVE DISORDER (ORDINARY LEAST SQUARES (OLS))",
                                "MAJOR DEPRESSIVE DISORDER VS ADHD (ORDINARY LEAST SQUARES (OLS))",
                                "MAJOR DEPRESSIVE DISORDER VS ANOREXIA NERVOSA (ORDINARY LEAST SQUARES (OLS)",
                                "MAJOR DEPRESSIVE DISORDER VS AUTISM SPECTRUM DISORDER (ORDINARY LEAST SQUARES",
                                "MAJOR DEPRESSIVE DISORDER OR STRESS-RELATED DISORDER",
                                "DEPRESSION X MATERNAL SMOKING DURING PREGNANCY INTERACTION",
                                "SCHIZOPHRENIA, BIPOLAR DISORDER OR RECURRENT MAJOR DEPRESSIVE DISORDER",
                                "MAJOR DEPRESSIVE DISORDER OR STROKE (PLEIOTROPY)",
                                "COGNITIVE FUNCTION (PROCESSING SPEED) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                "COGNITIVE FUNCTION (EXECUTIVE FUNCTION) X MAJOR DEPRESSIVE DISORDER INTERACTION (1DF)",
                                "COGNITIVE FUNCTION (EXECUTIVE FUNCTION) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                "COGNITIVE FUNCTION (DELAYED MEMORY) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                "COGNITIVE FUNCTION (PROCESSING SPEED) X MAJOR DEPRESSIVE DISORDER INTERACTION (2DF)",
                                "COGNITIVE FUNCTION (GLOBAL COGNITION) X MAJOR DEPRESSIVE DISORDER INTERACTION (1DF)",
                                "DEPRESSION X HERPES SIMPLEX 1 INFECTION ",
                                "DEPRESSION IN MULTIPLE SCLEROSIS",
                                "MENSTRUATION QUALITY OF LIFE IMPACT (DEPRESS",
                                "HYPERSOMNIA DURING A MAJOR DEPRESSIVE EPIS",
                                "(MTAG)",
                                "VASCULAR ENDOTHELIAL GROWTH FACTOR LEVELS"
)
Exclusion_terms_depression = unique(Exclusion_terms_depression)
# Saving exclusion terms
write(Exclusion_terms_depression, "Exclusion_terms_depression.txt")

# Extracting SNPs, curating names, performing filtering
write.csv(GWAS_CAT_Depress_specific, "GWAS_CAT_Depress_specific.csv")
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific$SNPS
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[stri_detect_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "rs")]
GWAS_CAT_Depress_specific_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "; ")) #2254 SNPs
GWAS_CAT_Depress_specific_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = ", ")) #2281 SNPs
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[!is.na(GWAS_CAT_Depress_specific_SNPs)]
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
GWAS_CAT_Depress_specific_SNPs = GWAS_CAT_Depress_specific_SNPs[stri_detect_fixed(GWAS_CAT_Depress_specific_SNPs, pattern = "rs")] #2281 SNPs
GWAS_CAT_Depress_specific_SNPs_count = as.data.frame(table(GWAS_CAT_Depress_specific_SNPs))
GWAS_CAT_Depress_specific_SNPs_count = arrange(GWAS_CAT_Depress_specific_SNPs_count, -Freq)
colnames(GWAS_CAT_Depress_specific_SNPs_count) = c("SNP", "Count")
write.xlsx(GWAS_CAT_Depress_specific_SNPs_count, "GWAS_CAT_Depress_specific_SNPs_count.xlsx")
GWAS_CAT_Depress_specific_SNPs_unique = unique(GWAS_CAT_Depress_specific_SNPs) #2073 SNPs

# mapping to dbSNP153
GWAS_CAT_Depress_specific_SNPs_153_track = lin_req_SNP_db_153(SNPs = GWAS_CAT_Depress_specific_SNPs_unique, File_name = "GWAS_CAT_Depress_specific_SNPs_unique.txt", lag_param = 10)
GWAS_CAT_Depress_specific_SNPs_153_track = GWAS_CAT_Depress_specific_SNPs_153_track[GWAS_CAT_Depress_specific_SNPs_153_track$name %in% GWAS_CAT_Depress_specific_SNPs_unique,]
GWAS_CAT_Depress_specific_SNPs_153_track = GWAS_CAT_Depress_specific_SNPs_153_track[!stri_detect_fixed(GWAS_CAT_Depress_specific_SNPs_153_track$chrom, pattern = "_")]
all(GWAS_CAT_Depress_specific_SNPs_153_track$name %in% GWAS_CAT_Depress_specific_SNPs_unique)
write.csv(GWAS_CAT_Depress_specific_SNPs_153_track, "GWAS_CAT_Depress_specific_SNPs_153_track.csv")
# 2052 SNPs identified


################### Mapping to genes ###################

Chromosome_coord_table = chromosome_table_getter()

# Importing USCS Genome browser tracks
# UCSC IDs to gene symbols
UCSC_Mapping = smart_fread("UCSC_ID_MAP/kgXref_hg19.txt") # https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1380809257_xAl2AOtRghxobFCyaRjWX2SAZXr0 # replace with an appropriate path
# UCSC known gene
UCSC_known_gene_track = list()
mySession = browserSession()
genome(mySession) = "hg19"

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_knownGene_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_knownGene_track = getTable(mySession, table = "knownGene", range =  PRF_knownGene_grange)
  
  PRF_knownGene_track$name2 = sapply(PRF_knownGene_track$name, function(x){
    x = UCSC_Mapping[UCSC_Mapping$`#kgID` == x, "geneSymbol"]
    return(x)
  })
  
  UCSC_known_gene_track[[i]] = PRF_knownGene_track
}

UCSC_known_gene_track = do.call(rbind, UCSC_known_gene_track)
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")]) # removes all variables with prefix


# complete wgEncodeGencodeCompV38lift37
wgEncodeGencodeCompV38lift37_gene_track = list()
wgEncodeGencodePseudoGeneV38lift37_gene_track = list()
wgEncodeGencodeV38lift37_combined_gene_track = list()
mySession = browserSession()
genome(mySession) = "hg19"

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_track = getTable(mySession, table = "wgEncodeGencodeCompV38lift37", range =  PRF_grange)
  wgEncodeGencodeCompV38lift37_gene_track[[i]] = PRF_track
  
  # Pseudogenes
  PRF_track_2 = getTable(mySession, table = "wgEncodeGencodePseudoGeneV38lift37", range =  PRF_grange)
  wgEncodeGencodePseudoGeneV38lift37_gene_track[[i]] = PRF_track_2
  
  # Combined
  PRF_track_combined = rbind(PRF_track, PRF_track_2)
  wgEncodeGencodeV38lift37_combined_gene_track[[i]] = PRF_track_combined
  
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")]) # removes all variables with prefix
wgEncodeGencodeCompV38lift37_gene_track = do.call(rbind, wgEncodeGencodeCompV38lift37_gene_track)
wgEncodeGencodePseudoGeneV38lift37_gene_track = do.call(rbind, wgEncodeGencodePseudoGeneV38lift37_gene_track)
wgEncodeGencodeV38lift37_combined_gene_track = do.call(rbind, wgEncodeGencodeV38lift37_combined_gene_track)



# ncbiRefSeq 
ncbiRefSeq_tracks_full = list()
mySession = browserSession()
genome(mySession) = "hg19"

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_track = getTable(mySession, table = "ncbiRefSeq", range =  PRF_grange)
  ncbiRefSeq_tracks_full[[i]] = PRF_track
}

rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")]) # removes all variables with prefix
ncbiRefSeq_tracks_full = do.call(rbind, ncbiRefSeq_tracks_full)


# mapping may be slow since loops in R are slow to execute
GWAS_CAT_Depress_specific_SNPs_153_mapped = gene_mapper(SNP_153_Track = GWAS_CAT_Depress_specific_SNPs_153_track, 
                                                        gene_track = wgEncodeGencodeV38lift37_combined_gene_track)
GWAS_CAT_Depress_specific_SNPs_153_mapped_UCSC = gene_mapper(SNP_153_Track = GWAS_CAT_Depress_specific_SNPs_153_track, 
                                                             gene_track = UCSC_known_gene_track)
GWAS_CAT_Depress_specific_SNPs_153_mapped_refSeq = gene_mapper(SNP_153_Track = GWAS_CAT_Depress_specific_SNPs_153_track, 
                                                               gene_track = ncbiRefSeq_tracks_full)
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_UCSC =  GWAS_CAT_Depress_specific_SNPs_153_mapped_UCSC$Mapped_gene
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_refSeq = GWAS_CAT_Depress_specific_SNPs_153_mapped_refSeq$Mapped_gene

# aggregate information from 3 mappings
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed = mapply(function(x,y,z){
  Combined_vector = c(x, y, z)
  Combined_vector = Combined_vector[!is.na(Combined_vector)]
  if (length(Combined_vector) > 0){
    Combined_vector = unlist(stri_split_fixed(Combined_vector, pattern = ";"))
    Combined_vector = paste0(Combined_vector, collapse = ";")
    return(Combined_vector)
  } else {
    return(NA)
  }
}, GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene, 
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_UCSC, 
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_refSeq)

# remove duplicates
GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed = sapply(GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed, function(x){
  genes = unlist(stri_split_fixed(x, pattern = ";"))
  genes = unique(genes)
  if (length(genes) > 1){
    genes = paste0(genes, collapse = ";")
  }
  return(genes)
})

# classify SNPs as Intergenic and Intragenic
GWAS_CAT_Depress_specific_SNPs_153_mapped$Intragenic = sapply(GWAS_CAT_Depress_specific_SNPs_153_mapped$Mapped_gene_fixed, function(x){
  if(is.na(x)){
    return("Intergenic")
  }
  return("Intragenic")
})
write.csv(GWAS_CAT_Depress_specific_SNPs_153_mapped, "GWAS_CAT_Depress_specific_SNPs_153_mapped.csv")

################### Generating chromosome and gene stats ###################

Chromosome_order = c(1:22, "X")
Chromosome_order = paste0("chr", Chromosome_order)
GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered = GWAS_CAT_Depress_specific_SNPs_153_mapped %>% arrange(factor(chrom, levels = Chromosome_order))
Stats_chromosome_full = GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered %>%
  group_by(., chrom, Intragenic) %>%
  summarise(Unique_SNPs = length(unique(name)))  %>%
  as.data.frame(.)
Stats_chromosome_full = ungroup(Stats_chromosome_full)
Stats_chromosome_full = arrange(Stats_chromosome_full, factor(chrom, levels = rev(Chromosome_order)))
Stats_chromosome_full$chrom = factor(Stats_chromosome_full$chrom, levels = rev(Chromosome_order))

Stats_chromosome_plot = ggplot(Stats_chromosome_full, aes(x = chrom, y = Unique_SNPs, fill = Intragenic)) +
  geom_col(position = position_dodge(1)) +
  scale_fill_manual(values = c("#412C84", "#FF8D40")) +
  geom_text(aes(label = Unique_SNPs), position = position_dodge(1), hjust = -0.2) +
  labs(x = "Chromosome", y = "Number of unique depression-related SNPs", fill = "Genetic context") +
  coord_flip() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 13, colour = "black"),
    legend.title = element_text(face = "bold", size = 14),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(size = 1),
    panel.grid.major.y = element_line(size = 0.5, linetype = 3, color =  "black"),
    strip.text = element_text(face = "bold", size = 14, colour = "blue")
  )
Stats_chromosome_plot

pdf("Stats_chromosome_plot.pdf", width = 10, height = 10)
Stats_chromosome_plot
dev.off()

# Top genes per chromosome
Top_genes_per_chromosome = list()
Top_genes_per_chromosome_top5 = list()
GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded = GWAS_CAT_Depress_specific_SNPs_153_mapped[GWAS_CAT_Depress_specific_SNPs_153_mapped$Intragenic == "Intragenic", ]
GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded = as.data.frame(GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded)
GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded = multiple_expander(df = GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded, cols_to_expand = 22, pattern = ";")

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:length(Chromosome_order)){
  PRF_Curr_chr = Chromosome_order[i]
  PRF_Curr_df = GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded[GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded$chrom == PRF_Curr_chr,]
  PRF_Curr_Stat = as.data.frame(table(PRF_Curr_df$Mapped_gene_fixed))
  PRF_Curr_Stat = cbind(PRF_Curr_chr, PRF_Curr_Stat)
  colnames(PRF_Curr_Stat) = c("Chrom","Gene", "Count")
  PRF_Curr_Stat = arrange(PRF_Curr_Stat, -Count)
  Top_genes_per_chromosome[[i]] = PRF_Curr_Stat
  
  # selecting top 5
  PRF_Curr_Stat_5 = PRF_Curr_Stat[1:5,]
  PRF_Curr_Stat_5 = PRF_Curr_Stat_5[!is.na(PRF_Curr_Stat_5$Chrom),]
  Top_genes_per_chromosome_top5[[i]] = PRF_Curr_Stat_5
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

Top_genes_per_chromosome = do.call(rbind, Top_genes_per_chromosome)
Top_genes_per_chromosome_top5 = do.call(rbind, Top_genes_per_chromosome_top5)
openxlsx::write.xlsx(Top_genes_per_chromosome, "Top_genes_per_chromosome.xlsx", overwrite = TRUE)
openxlsx::write.xlsx(Top_genes_per_chromosome_top5, "Top_genes_per_chromosome_top5.xlsx", overwrite = TRUE)

# Additional stats for genes
STATS_genes_all = as.data.frame(table(GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded$Mapped_gene_fixed))
STATS_genes_all = arrange(STATS_genes_all, -Freq)
colnames(STATS_genes_all) = c("Gene", "Count")
write.xlsx(STATS_genes_all, "STATS_genes_all.xlsx")

# Stats for "high" confidence SNPs
GWAS_CAT_Depress_specific_genome_sign = GWAS_CAT_Depress_specific
GWAS_CAT_Depress_specific_genome_sign$`P-VALUE` = as.numeric(GWAS_CAT_Depress_specific_genome_sign$`P-VALUE`)
GWAS_CAT_Depress_specific_genome_sign = GWAS_CAT_Depress_specific_genome_sign[GWAS_CAT_Depress_specific_genome_sign$`P-VALUE` < (5 * 10^-8),]
GWAS_CAT_Depress_specific_genome_sign_SNPs = GWAS_CAT_Depress_specific_genome_sign$SNPS
GWAS_CAT_Depress_specific_genome_sign_SNPs = GWAS_CAT_Depress_specific_genome_sign_SNPs[stri_detect_fixed(GWAS_CAT_Depress_specific_genome_sign_SNPs, pattern = "rs")]
GWAS_CAT_Depress_specific_genome_sign_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_genome_sign_SNPs, pattern = "; ")) # 1029 SNPs
GWAS_CAT_Depress_specific_genome_sign_SNPs = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_genome_sign_SNPs, pattern = ", ")) # 1029 SNPs
GWAS_CAT_Depress_specific_genome_sign_SNPs = GWAS_CAT_Depress_specific_genome_sign_SNPs[!is.na(GWAS_CAT_Depress_specific_genome_sign_SNPs)] # 1029 SNPs, 879 unique
GWAS_CAT_Depress_specific_genome_sign_SNPs_stats = as.data.frame(table(GWAS_CAT_Depress_specific_genome_sign_SNPs))
GWAS_CAT_Depress_specific_genome_sign_SNPs_stats = arrange(GWAS_CAT_Depress_specific_genome_sign_SNPs_stats, -Freq) # 114 SNPs reported more than once, 0.1296928 reproducibility

colnames(GWAS_CAT_Depress_specific_genome_sign_SNPs_stats) = c("Genome-wide significant SNP", "Count")
write.xlsx(GWAS_CAT_Depress_specific_genome_sign_SNPs_stats, "GWAS_CAT_Depress_specific_genome_sign_SNPs_stats.xlsx", overwrite = TRUE)
GWAS_CAT_Depress_specific_genome_sign_SNPs_mapped = GWAS_CAT_Depress_specific_SNPs_153_mapped[GWAS_CAT_Depress_specific_SNPs_153_mapped$name %in% GWAS_CAT_Depress_specific_genome_sign_SNPs,] # 865 unique SNPs mapped to dbSNP153 track


################### Enrichment analysis ###################

# Preparing reference Entrez dataset and gene universe
ENTREZ_genes_Homo_Sapiens = smart_fread("/home/aleksandr/Desktop/WORK/Broad_Depression_Paper_folder/Depression_omics_multi_cohort/GWAS_Catalog/Entrez Homo sapiens/data") # Replace with an appropriate path
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
Gene_universe = ENTREZ_genes_Homo_Sapiens$GeneID
Gene_universe = as.character(Gene_universe)
Gene_universe = unique(Gene_universe)


##### 1) Enrichment for all intragenic SNPs ##### 
# Preparing reference universe and mapped gene list
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$Mapped_gene_fixed
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[!is.na(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)]
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, pattern = ";"))
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unique(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)

# Mapping
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% GWAS_CAT_Depress_specific_SNPs_153_mapped__genes){
    return(TRUE)
  }
  return(FALSE)
}),] #843 genes are mapped directly
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped[ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$chromosome != "chr-",]
# 843 genes are mapped still
Non_mapped_genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[GWAS_CAT_Depress_specific_SNPs_153_mapped__genes %!in% ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$Symbol]
# 249 genes are not mapped
write(Non_mapped_genes, "Non_mapped_genes_GWAS.txt")

# mapping through synonyms
ENTREZ_Genes_GWAS_idx = vector()
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$Mapped_gene_fixed, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered[Tested_idx, "chrom"])){
      ENTREZ_Genes_GWAS_idx[i] = TRUE
    } else {
      ENTREZ_Genes_GWAS_idx[i] = FALSE
    }
    
  } else {
    ENTREZ_Genes_GWAS_idx[i] = FALSE
  }
}
table(ENTREZ_Genes_GWAS_idx) # 14 genes are mapped by synonyms
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped_2 = ENTREZ_genes_Homo_Sapiens[ENTREZ_Genes_GWAS_idx,]
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = rbind(
  ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped,
  ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped_2
)
any(duplicated(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$GeneID)) #duplicates
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = distinct(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped, GeneID, .keep_all = TRUE) #847 genes are mapped to Entrez
write.csv(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped, "Entrez_GWAS_CAT_Depress_specific_SNPs_153_mapped.csv")

ENRICHMENT_GWAS_CAT_Depress_intragenic = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$GeneID,
                                                                         universe = Gene_universe,
                                                                         categories_to_show = 30,
                                                                         folder = "Enrichment_Depression_GWAS_CAT",
                                                                         plot_name_pref = "Depr_GWAS_intragen")

##### 2) Enrichment for overlapping SNPs ##### 
Depr_SNP_freq = as.data.frame(table(GWAS_CAT_Depress_specific_SNPs))
Depr_SNP_freq_2 = Depr_SNP_freq[Depr_SNP_freq$Freq > 1,]
GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered_more_1 = GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered[GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$name %in% Depr_SNP_freq_2$GWAS_CAT_Depress_specific_SNPs,]
# 166 SNPs appeared more than once
Depression_frequent_genes = GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered_more_1$Mapped_gene_fixed
Depression_frequent_genes = Depression_frequent_genes[!is.na(Depression_frequent_genes)]
Depression_frequent_genes = unlist(stri_split_fixed(Depression_frequent_genes, pattern = ";"))
Depression_frequent_SNP_genes_stats = as.data.frame(table(Depression_frequent_genes))
Depression_frequent_SNP_genes_stats = arrange(Depression_frequent_SNP_genes_stats, -Freq)
colnames(Depression_frequent_SNP_genes_stats) = c("Gene", "Count")
write.xlsx(Depression_frequent_SNP_genes_stats, "Depression_frequent_SNP_genes_stats.xlsx")
Depression_frequent_genes = unique(Depression_frequent_genes) # 121 gene with frequent SNPs
write(Depression_frequent_genes, "Frequent_SNP_genes_GWAS_depression.txt")

ENTREZ_Depression_frequent_genes = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% Depression_frequent_genes){
    return(TRUE)
  }
  return(FALSE)
}),] # 98 genes are mapped directly

Non_mapped_genes = Depression_frequent_genes[Depression_frequent_genes %!in% ENTREZ_Depression_frequent_genes$Symbol] #23 genes are not mapped
ENTREZ_Genes_GWAS_idx = vector()

# Mapping through synonyms
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$Mapped_gene_fixed, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered[Tested_idx, "chrom"])){
      ENTREZ_Genes_GWAS_idx[i] = TRUE
    } else {
      ENTREZ_Genes_GWAS_idx[i] = FALSE
    }
    
  } else {
    ENTREZ_Genes_GWAS_idx[i] = FALSE
  }
  
}
table(ENTREZ_Genes_GWAS_idx) # 1 gene is mapped by a synonym
ENTREZ_Depression_frequent_genes_2 = ENTREZ_genes_Homo_Sapiens[ENTREZ_Genes_GWAS_idx,]
ENTREZ_Depression_frequent_genes = rbind(
  ENTREZ_Depression_frequent_genes,
  ENTREZ_Depression_frequent_genes_2
)
any(duplicated(ENTREZ_Depression_frequent_genes$GeneID)) # NO duplicates
ENRICHMENT_GWAS_CAT_Depress_intragenic_freq_SNPs = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_Depression_frequent_genes$GeneID,
                                                                                   universe = Gene_universe,
                                                                                   categories_to_show = 30,
                                                                                   folder = "Enrichment_Depression_freq_SNPs",
                                                                                   plot_name_pref = "Enrichment_Depression_freq_SNPs")

##### 3) Enrichment for frequent genes ##### 
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$Mapped_gene_fixed
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[!is.na(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)]
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unlist(stri_split_fixed(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, pattern = ";"))
GWAS_CAT_gene_freq = as.data.frame(table(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)) # 1093 unqiue genes
GWAS_CAT_gene_freq_2 = GWAS_CAT_gene_freq[GWAS_CAT_gene_freq$Freq >1,] # 238 frequent genes
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[GWAS_CAT_Depress_specific_SNPs_153_mapped__genes %in% GWAS_CAT_gene_freq_2$GWAS_CAT_Depress_specific_SNPs_153_mapped__genes]
# 809 genes with duplicates
GWAS_CAT_Depress_specific_SNPs_153_mapped__genes = unique(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes)
# 238 genes without duplicates
write(GWAS_CAT_Depress_specific_SNPs_153_mapped__genes, "Frequent_genes_GWAS_depression.txt")

ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% GWAS_CAT_Depress_specific_SNPs_153_mapped__genes){
    return(TRUE)
  }
  return(FALSE)
}),] # 192 genes are mapped directly

ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped[ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$chromosome != "chr-",]
# 192 genes are mapped still
Non_mapped_genes = GWAS_CAT_Depress_specific_SNPs_153_mapped__genes[GWAS_CAT_Depress_specific_SNPs_153_mapped__genes %!in% ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$Symbol]
# 46 genes are not mapped

# Mapping through synonyms
ENTREZ_Genes_GWAS_idx = vector()
for (i in 1:nrow(ENTREZ_genes_Homo_Sapiens)){
  print(i)
  y = ENTREZ_genes_Homo_Sapiens$Synonyms[i]
  Splitted_synonyms = unlist(stri_split_fixed(str = y, pattern = "|"))
  Splitted_synonyms = c(Splitted_synonyms, ENTREZ_genes_Homo_Sapiens$Symbol_from_nomenclature_authority[i])
  Splitted_synonyms = str_trim(Splitted_synonyms)
  
  if (any(Splitted_synonyms %in% Non_mapped_genes)){
    Tested_idx = lapply(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered$Mapped_gene_fixed, function(x) unlist(stri_split_fixed(x, pattern = ";")))
    Tested_idx = sapply(Tested_idx, function(x){
      x = x[x != "NA"]
      x = x[!is.na(x)]
      if (any(Splitted_synonyms %in% x)){
        return(TRUE)
      }
      return(FALSE)
    })
    
    if (ENTREZ_genes_Homo_Sapiens$chromosome[i] %in% unique(GWAS_CAT_Depress_specific_SNPs_153_mapped_ordered[Tested_idx, "chrom"])){
      ENTREZ_Genes_GWAS_idx[i] = TRUE
    } else {
      ENTREZ_Genes_GWAS_idx[i] = FALSE
    }
    
  } else {
    ENTREZ_Genes_GWAS_idx[i] = FALSE
  }
}
table(ENTREZ_Genes_GWAS_idx) # 3 genes are mapped by synonyms

ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped_2 = ENTREZ_genes_Homo_Sapiens[ENTREZ_Genes_GWAS_idx,]
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = rbind(
  ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped,
  ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped_2
)
any(duplicated(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$GeneID)) # duplicates
ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped = distinct(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped, GeneID, .keep_all = TRUE) # 193 genes are mapped to Entrez
write.csv(ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped, "Entrez_frequent_genes_GWAS_depression.csv")

ENRICHMENT_GWAS_CAT_Depress_intragenic_freq_genes = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_GWAS_CAT_Depress_specific_SNPs_153_mapped$GeneID,
                                                                                    universe = Gene_universe,
                                                                                    categories_to_show = 30,
                                                                                    folder = "Enrichment_Depression_freq_genes",
                                                                                    plot_name_pref = "Enrichment_Depression_freq_genes")

################### Chromosome maps ###################
# Getting cytoBand track from UCSC
Chromosome_coord_table = chromosome_table_getter()
Full_cytoband_track = list()

# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
for (i in 1:nrow(Chromosome_coord_table)){
  PRF_chromosome = paste0("chr", Chromosome_coord_table$Chromosome[i])
  print(PRF_chromosome)
  PRF_CytoBand_grange = GRanges(PRF_chromosome, IRanges(1, Chromosome_coord_table$`Total length (bp)`[i]), strand = "*")
  PRF_Cytoband_track = getTable(mySession, table = "cytoBand", range =  PRF_CytoBand_grange)
  Full_cytoband_track[[i]] = PRF_Cytoband_track
}

Full_cytoband_track = do.call(rbind, Full_cytoband_track) #last p and first Q are centromers (acen), the end of the last p is the coordinate of a centomere
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

# Making chromosome df
Chromosome_map_Rideogram = Chromosome_coord_table_prepared
Chromosome_map_Rideogram$Start = 0
Chromosome_map_Rideogram$Centromere = NULL
Chromosome_map_Rideogram$Centromere_St = Centromere_Starts
Chromosome_map_Rideogram$Centromere_End = Centromere_Ends
colnames(Chromosome_map_Rideogram) = c("Chr","Start","End","CE_start","CE_end")

# Preparing Complete SNP dataset
GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded_all = rbind(GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded,
                                                               GWAS_CAT_Depress_specific_SNPs_153_mapped[GWAS_CAT_Depress_specific_SNPs_153_mapped$Intragenic == "Intergenic", ])
write.csv(GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded_all, "GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded_all.csv")

# Obtaining Chromosome map in a loop
# Note: the loop variables contain prefix PRF_ to enable easy cleanup after the loop is executed
Chrom_depr_heatmap = list()

# Outer loop walking through chromosomes
for (i in 1:nrow(Chromosome_map_Rideogram)){
  PRF_curr_chrom = Chromosome_map_Rideogram$Chr[i]
  writeLines(PRF_curr_chrom)
  PRF_Curr_SNP_df = GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded_all[GWAS_CAT_Depress_specific_SNPs_153_mapped_expanded_all$chrom == PRF_curr_chrom,]
  PRF_curr_chrom_map = Chromosome_map_Rideogram[i,]
  
  if (nrow(PRF_Curr_SNP_df) < 1){
    break
  }
  
  PRF_curr_chrom_vector = seq(from= PRF_curr_chrom_map$Start, by = 1000000, to = PRF_curr_chrom_map$End)
  
  if (PRF_curr_chrom_vector[length(PRF_curr_chrom_vector)] < PRF_curr_chrom_map$End){
    PRF_curr_chrom_vector = c(PRF_curr_chrom_vector, PRF_curr_chrom_map$End)
  }
  
  PRF_Curr_Chrom_Table_Heatmap = list()
  
  # Inner loop for a chromosome
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
    PRF_Curr_SNP_df_interval = PRF_Curr_SNP_df[PRF_Curr_SNP_df$chromStart %in% RPF_Interval,]
    
    if (nrow(PRF_Curr_SNP_df_interval) < 1){
      PRF_identified_SNP = NA
      PRF_identified_Genes = NA
      PRF_SNP_Count_unique = 0
      PRF_Gene_Count_unique = 0
    } else {
      PRF_identified_SNP = PRF_Curr_SNP_df_interval$name
      PRF_identified_SNP = unique(PRF_identified_SNP)
      PRF_SNP_Count_unique = length(PRF_identified_SNP)
      PRF_identified_SNP = paste0(PRF_identified_SNP, collapse = ";")
      
      # identifying genes
      PRF_identified_Genes = PRF_Curr_SNP_df_interval$Mapped_gene_fixed
      PRF_identified_Genes = PRF_identified_Genes[!is.na(PRF_identified_Genes)]
      PRF_identified_Genes = unlist(stri_split_fixed(PRF_identified_Genes, pattern = ";"))
      PRF_identified_Genes = unique(PRF_identified_Genes)
      PRF_identified_Genes = PRF_identified_Genes[!is.na(PRF_identified_Genes)]
      PRF_Gene_Count_unique = length(PRF_identified_Genes)
      
      if (length(PRF_identified_Genes) > 0){
        PRF_identified_Genes = paste0(PRF_identified_Genes, collapse = ";")
      } else {
        PRF_identified_Genes = NA
      }
      
    }
    
    PRF_Curr_df_stat = data.frame(Chrom = PRF_curr_chrom,Start = PRF_Start, End = PRF_End, Ident.SNPs.count = PRF_SNP_Count_unique, 
                                  Ident.gene.count = PRF_Gene_Count_unique, Ident.SNPs = PRF_identified_SNP, Ident.Genes = PRF_identified_Genes)
    PRF_Curr_Chrom_Table_Heatmap[[index]] = PRF_Curr_df_stat
  }
  
  if (length(PRF_Curr_Chrom_Table_Heatmap) > 1){
    PRF_Curr_Chrom_Table_Heatmap = do.call(rbind, PRF_Curr_Chrom_Table_Heatmap)
  } else {
    PRF_Curr_Chrom_Table_Heatmap = PRF_Curr_Chrom_Table_Heatmap[[1]]
  }
  
  Chrom_depr_heatmap[[i]] = PRF_Curr_Chrom_Table_Heatmap
}
rm(list = ls()[stri_detect_fixed(ls(), pattern = "PRF_")])

# Checking chromosomal map
sum(Chrom_depr_heatmap$Ident.SNPs.count) # 2052 SNPs in total that matches to the number of SNPs mapped to dbSNP153

Chrom_depr_heatmap = do.call(rbind, Chrom_depr_heatmap)
write.csv(Chrom_depr_heatmap, "Chrom_map_GWAS_CAT.csv")

# Test visualization
Chrom_depr_heatmap_plot = Chrom_depr_heatmap[,1:4]
colnames(Chrom_depr_heatmap_plot) = c("Chr", "Start", "End", "Value")
ideogram(karyotype = Chromosome_map_Rideogram, overlaid = Chrom_depr_heatmap_plot, colorset1 = c("#FFFFFF", "#000000"))
convertSVG("chromosome.svg", device = "png")




