#'
#' @name ceRNASurvival
#' @title Survival analysis and visualization
#' @description A function to analyze the survival outcome when people carry the
#' identified ceRNAs and visualize the results using Kaplan-Meier plot
#'
#' @import utils
#' @import grDevices
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param mirnas a list of mirna name
#'
#' @export
#'
#' @examples
#' ceRNASurvival(
#' path_prefix = '~/',
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' mirnas = 'hsa-miR-101-3p'
#' )
#'

ceRNASurvival <- function(path_prefix,
                          project_name = 'demo',
                          disease_name = 'DLBC',
                          mirnas = 'hsa-miR-101-3p'){

  if (!stringr::str_detect(path_prefix, '/')){
    path_prefix <- paste0(path_prefix, '/')
  }

  time1 <- Sys.time()
  #setwd(paste0(project_name,'-',disease_name))
  if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/'))){
    dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/'))
  }

  if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/survivalResults/'))){
    dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/survivalResults/'))
  }

  message('\u25CF Step5: Dowstream Analyses - Survival analysis')

  # expression data
  geneExp <- as.data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mrna.csv'),header = T, stringsAsFactors = F))
  row.names(geneExp) <- geneExp[,1]
  geneExp <- geneExp[,-1]
  names(geneExp) <- substring(names(geneExp),1,12)

  # survival data
  survivalData <- as.data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_survival.csv'),header = T, stringsAsFactors = F))
  row.names(survivalData) <- survivalData[,1]
  survivalData <- survivalData[,-1]
  survivalData <- survivalData[,colnames(survivalData)%in%c('OS', 'OS.time')]

  # mirna-gene pairs results
  Res <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/03_identifiedPairs/', project_name, '-', disease_name,'_finalpairs.rds'))
  Res_dataframe <- Reduce(rbind, Res)

  mirna <- as.data.frame(Reduce(rbind, Res_dataframe[,1]))
  gene <- as.data.frame(Reduce(rbind, Res_dataframe[,2]))
  mirna_gene_pairs <- cbind(mirna, gene)
  names(mirna_gene_pairs) <- c('mirna', 'gene')
  mirna_uni <- mirnas
  message('\u2605 Total number of identified miRNA in ',project_name, '-', disease_name,' cohort: ', length(mirna_uni), '.')
  message('\u2605 Total number of identified miRNA-ceRNAs in ',project_name, '-', disease_name,' cohort: ', dim(mirna_gene_pairs)[1], '.')

  runforeachmirna <- function(each_mirna){
    if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/survivalResults/', each_mirna))){
      dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/survivalResults/', each_mirna))
    }
    df <- mirna_gene_pairs[mirna_gene_pairs$mirna == each_mirna,]
    if (dim(df)[1] !=1) {
      each_gene <- as.data.frame(Reduce(rbind,stringr::str_split(df[,2], " ")))
    }else{
      each_gene <- t(as.data.frame(stringr::str_split(df[,2], " ")[[1]]))
    }

    get_ExpSurData <- function(which_gene){
      #which_gene <- as.character(each_gene[1,1])
      which_geneExp <- as.data.frame(t(geneExp[row.names(geneExp) == which_gene,]))
      which_geneExp_sur <- merge(which_geneExp, survivalData, by = 'row.names')
      row.names(which_geneExp_sur) <- which_geneExp_sur[,1]
      which_geneExp_sur$level <- ifelse(which_geneExp_sur[,2] > stats::median(which_geneExp_sur[,2]), 1, 0)
      which_geneExp_sur <- which_geneExp_sur[,-1:-2]
      names(which_geneExp_sur)[1:2] <- c('status','time')
      return(which_geneExp_sur)
    }
    draw_surCurve <- function(which_gene){
      which_geneExp_sur <- get_ExpSurData(which_gene)
      fit<- survival::survfit(survival::Surv(time, status)~level, data = which_geneExp_sur)
      survplot <- survminer::ggsurvplot(fit, data = which_geneExp_sur, pval = TRUE, pval.coord = c(100, 0.03),
                            conf.int = TRUE, legend = c(0.9,0.9),
                            palette = c("#E7B800", "#2E9FDF"),
                            title = which_gene)
    }
    plot_lst <- list()
    for (j in 1:dim(each_gene)[1]){
      #j =1
      which_gene <- as.character(each_gene[j,])
      surCurve_list <- purrr::map(which_gene, draw_surCurve)
      grDevices::png(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/survivalResults/', each_mirna, '/', each_mirna,'_', which_gene[1],'-', which_gene[2], '_triplets.png'),height = 1500, width = 2800, res = 300)
      survminer::arrange_ggsurvplots(surCurve_list, print = TRUE, ncol=2, nrow=1,
                                     title = paste0(each_mirna,'_', which_gene[1],'-', which_gene[2], ' triplets'),
                                     surv.plot.height = 0.6)
      grDevices::dev.off()
      plot_lst[[j]] <- survminer::arrange_ggsurvplots(surCurve_list, print = TRUE, ncol=2, nrow=1,
                                                  title = paste0(each_mirna,'_', which_gene[1],'-', which_gene[2], ' triplets'),
                                                  surv.plot.height = 0.6)
    }
    plot_lst
  }
  if (length(mirna_uni)==1){
    all_plots <- runforeachmirna(mirna_uni)
  }else{
    tmp <- purrr::map(mirna_uni, runforeachmirna)
  }

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
  message('\u2605\u2605\u2605 Survival analysis has completed! \u2605\u2605\u2605')

  all_plots

}


