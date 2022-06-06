#'
#' @name ceRNAFunction
#' @title Functional enrichment analysis
#' @description A function to conduct Functional enrichment analysis for
#' biological interpretation. The databases supported in this function include
#' Gene Ontology (GO; http://geneontology.org/) and Kyoto Encyclopedia of
#' Genes and Genomes (KEGG; https://www.genome.jp/kegg/).
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param window_size the number of samples for each window and usually about
#' one third of total samples
#' @param cor_method selection of correlation methods, including pearson and
#' spearman (default: pearson)
#' @param cor_threshold_peak peak threshold of correlation value between 0 and 1
#' (default: 0.85)
#'
#' @importfrom SPONGE sponge_gene_miRNA_interaction_filter
#' @examples
#' ceRNAFunction(
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' pairs_cutoff = 1
#' )
#'
#' @export


ceRNAFunction <- function(path_prefix = NULL,
                          project_name,
                          disease_name,
                          pairs_cutoff){

  if (is.null(path_prefix)){
    path_prefix <- getwd()
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }else{
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }

  time1 <- Sys.time()
  #setwd(paste0(project_name,'-',disease_name))

  if(!dir.exists(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/'))){
    dir.create(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/'))
  }

  if(!dir.exists(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/'))){
    dir.create(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/'))
  }

  message('\u25CF Step5: Dowstream Analyses - Functional analysis')
  Res <- readRDS(paste0(project_name,'-',disease_name,'/03_identifiedPairs/', project_name, '-', disease_name,'_finalpairs.rds'))
  Res_dataframe <- Reduce(rbind, Res)
  mir_unique <- unique(Res_dataframe[,1])
  mir_df_final <- c()
  skip_count <- 0
  remain_count <- 0
  for (i in mir_unique){
    mir_df <- Res_dataframe[Res_dataframe[,1]== i,]
    if (dim(mir_df)[1] <= pairs_cutoff || is.null(dim(mir_df)[1])){
      #print(paste0(i, ' is skipped.'))
      skip_count <- skip_count + 1
      next
    }else{
      # select ceRNAs for pathway analysis using specific miRNA
      #print(paste0(i, ' will be involved.'))
      mir_df_final <- rbind(mir_df_final, mir_df)
      remain_count <- remain_count + 1
      }
  }

  message(paste0('\u2605 The ceRNAs related to ', skip_count, ' miRNAs are skipped.'))
  message(paste0('\u2605 The ceRNAs related to ', remain_count, ' miRNAs are involved.'))

  genes <- as.data.frame(Reduce(rbind,mir_df_final[,2]))
  each_gene <- Reduce(rbind,stringr::str_split(genes[,1], " "))
  each_gene <- unique(c(each_gene[,1], each_gene[,2]))
  message('\u2605 Data preprocessing ...')

  # prematch to OrgDb
  db <- GOSemSim::load_OrgDb(OrgDb = "org.Hs.eg.db")
  gene.df <- suppressWarnings(AnnotationDbi::select(db,
                                                keys = each_gene,
                                                keytype = "SYMBOL",
                                                columns=c("SYMBOL", "ENSEMBL", "ENTREZID")))

  # kegg ora
  message('\u2605 Running KEGG ORA analysis ...')
  kk <- clusterProfiler::enrichKEGG(gene = gene.df[,3],
                                    organism = 'hsa',
                                    pvalueCutoff = 0.01)

  kk_df <- kk@result
  write.csv(kk_df,paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/',project_name,'-',disease_name,'_kegg_ora.csv'), row.names = F)
  kk_babble <- enrichplot::dotplot(kk, showCategory=10,orderBy = "x")+ ggplot2::ggtitle('Dotplot for ORA based on KEGG')
  kk_bar <- barplot(kk,showCategory = 10)+ ggplot2::ggtitle('Barplot for ORA based on KEGG')

  # go ora
  message('\u2605 Running GO ORA analysis ...')
  go_level <- c('CC','MF','BP')

  runGOora <- function(go_level){
    #go_level <- go_level[1]
    ego <- clusterProfiler::enrichGO(gene          = gene.df[,3],
                                     OrgDb         = org.Hs.eg.db,
                                     ont           = go_level,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01,
                                     qvalueCutoff  = 0.05,
                                     readable      = TRUE)
    clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  }
  plotGOora_babble <- function(go_level, enrichResult){
    enrichplot::dotplot(enrichResult, showCategory=10,orderBy = "x")+ ggplot2::ggtitle(paste0('Dotplot for ORA based on GO:', go_level))
  }

  go <- purrr::map(go_level, runGOora)
  names(go) <- go_level
  go_babble <- purrr::map2(go_level, go, plotGOora_babble)

  go_cc <- go[[1]]@result
  go_cc$GOLevel <- 'CC'
  go_mf <- go[[2]]@result
  go_mf$GOLevel <- 'MF'
  go_bp <- go[[3]]@result
  go_bp$GOLevel <- 'BP'
  go_all <- rbind(go_cc, go_mf, go_bp)
  write.csv(go_all,paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/',project_name, '-',disease_name, '_go_ora.csv'), row.names = F)
  go_tobarplot <- rbind(head(go_cc, n=10), head(go_mf, n=10), head(go_bp, n=10))
  go_tobarplot$p.adjust.convert <- -log10(go_tobarplot$p.adjust)
  go_bar <- ggpubr::ggbarplot(go_tobarplot, x = "Description", y = "p.adjust.convert",
                              fill = "GOLevel",
                              color = "white",
                              palette = "jco",
                              sort.val = "asc",
                              sort.by.groups = TRUE,
                              x.text.angle = 0,
                              y.text.angle = 0) + ggpubr::rotate()
  # merge kegg and GO results
  gg_babble <- cowplot::plot_grid(kk_babble, go_babble[[1]], go_babble[[2]], go_babble[[3]], labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol=2)
  ggplot2::ggsave(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/', project_name,'-',disease_name,'_function_babble.png'), height = 10, width = 20,dpi = 300)
  gg_bar <- cowplot::plot_grid(kk_bar, go_bar, labels = c('A', 'B'), label_size = 12, ncol = 2)
  ggplot2::ggsave(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/', project_name,'-',disease_name,'_function_bar.png'), height = 8, width = 20,dpi = 300)

  closeAllConnections()

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))

  message('\u2605\u2605\u2605 All analyses has completed! \u2605\u2605\u2605')
}

