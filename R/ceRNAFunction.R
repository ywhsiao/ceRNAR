#'
#' @name ceRNAFunction
#' @title Functional enrichment analysis
#' @description A function to conduct Functional enrichment analysis for
#' biological interpretation. The databases supported in this function include
#' Gene Ontology (GO; http://geneontology.org/) and Kyoto Encyclopedia of
#' Genes and Genomes (KEGG; https://www.genome.jp/kegg/).
#'
#' @import org.Hs.eg.db
#' @import grDevices
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default: demo)
#' @param disease_name the abbreviation of disease that users are interested in
#' (default: DLBC)
#' @param pairs_cutoff at least the number of ceRNA pairs that a mirna must have
#'  (default: 1)
#'
#' @returns a list object
#' @export
#'
#' @examples
#' ceRNAFunction(
#' path_prefix = NULL,
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' pairs_cutoff = 1
#' )
#'
#'


ceRNAFunction <- function(path_prefix = NULL,
                          project_name = 'demo',
                          disease_name = 'DLBC',
                          pairs_cutoff = 1){

  if (is.null(path_prefix)){
    path_prefix <- fs::path_home()
  }else{
    path_prefix <- path_prefix
  }

  if (!stringr::str_detect(path_prefix, '/$')){
    path_prefix <- paste0(path_prefix, '/')
  }

  time1 <- Sys.time()
  #setwd(paste0(project_name,'-',disease_name))

  if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/'))){
    dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/'))
  }

  if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/'))){
    dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/'))
  }

  message('\u25CF Step5: Dowstream Analyses - Functional analysis')
  Res <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/03_identifiedPairs/', project_name, '-', disease_name,'_finalpairs.rds'))
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

  kk_df <- as.data.frame(kk@result)
  utils::write.csv(kk_df,paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/',project_name,'-',disease_name,'_kegg_ora.csv'), row.names = FALSE)
  kk_babble <- enrichplot::dotplot(kk, showCategory=10,orderBy = "x")+ ggplot2::ggtitle('Dotplot for ORA based on KEGG')
  kk_bar <- graphics::barplot(kk,showCategory = 10)+ ggplot2::ggtitle('Barplot for ORA based on KEGG')

  # go ora
  message('\u2605 Running GO ORA analysis ...')

  go_level <- c('CC','MF','BP')

  runGOora <- function(go_level){
    #go_level <- go_level[1]
    ego <- clusterProfiler::enrichGO(gene          = gene.df[,3],
                                     OrgDb         = db,
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

  if (dim(go[[1]]@result)[1] != 0){
    GO_cc <- as.data.frame(go[[1]]@result)
    GO_cc$GOLevel <- 'CC'
  }
  if (dim(go[[2]]@result)[1] != 0){
    GO_mf <- as.data.frame(go[[2]]@result)
    GO_mf$GOLevel <- 'MF'
  }
  if (dim(go[[3]]@result)[1] != 0){
    GO_bp <- as.data.frame(go[[3]]@result)
    GO_bp$GOLevel <- 'BP'
  }
  GloEnv <- ls(pattern = '^GO_')
  go_lst <- list()
  go_lst_tobar <- list()
  for (i in 1:length(GloEnv)){
      go_lst[[i]] <- get(GloEnv[i])
      go_lst_tobar[[i]] <- utils::head(get(GloEnv[i]), n=10)
  }
  go_all <- as.data.frame(Reduce(rbind, go_lst))
  utils::write.csv(go_all,paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/',project_name, '-',disease_name, '_go_ora.csv'), row.names = FALSE)
  go_tobarplot <- as.data.frame(Reduce(rbind, go_lst_tobar))
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
  ggplot2::ggsave(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/', project_name,'-',disease_name,'_function_babble.png'), height = 10, width = 20,dpi = 300)
  gg_bar <- cowplot::plot_grid(kk_bar, go_bar, labels = c('A', 'B'), label_size = 12, ncol = 2)
  ggplot2::ggsave(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/functionResults/', project_name,'-',disease_name,'_function_bar.png'), height = 8, width = 20,dpi = 300)

  # return as a list object
  function_plots <- list(gg_babble, gg_bar)

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
  message('\u2605\u2605\u2605 All analyses has completed! \u2605\u2605\u2605')
  function_plots
}


