#'
#' @name ceRNAIntergate
#' @title Integration of the possible ceRNA pairs among published tools
#' @description A function to integrate the possible ceRNA pairs that are found
#' by ceRNAR algorithm with those from other tools, such as SPONGE (List et al.,
#' 2019) and RJAMI (Hornakova et al.,)2018.
#'
#' @importfrom SPONGE sponge_gene_miRNA_interaction_filter
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#'
#' @examples
#' ceRNAIntegrate(
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' )
#'
#' @export
#'
ceRNAIntegrate <- function(path_prefix = NULL,
                           project_name = 'TCGA',
                           disease_name){

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

  if(!dir.exists(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/integration/'))){
    dir.create(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/integration/'))
  }

  message('\u25CF Step5: Dowstream Analyses - Integration')

  dict <- readRDS(paste0(project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_MirnaTarget_dictionary.rds'))
  mirna <- data.frame(data.table::fread(paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mirna.csv')),row.names = 1)
  mrna <- data.frame(data.table::fread(paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mrna.csv')),row.names = 1)
  d <- as.data.frame(matrix(0,nrow = dim(mrna)[1], ncol = dim(mirna)[1]))
  names(d) <- row.names(mirna)
  row.names(d) <- row.names(mrna)
  for (i in 1:dim(dict)[1]){
    #i=1
    gene_pair <- dict[i,][[2]]
    d[gene_pair,i] <- 1
  }

  # SPONGE
  doParallel::registerDoParallel(cores=4)
  mir_expr <- t(mirna)
  gene_expr <- t(mrna)
  genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
    gene_expr = gene_expr,
    mir_expr = mir_expr,
    mir_predicted_targets = as.matrix(d))

  ceRNA_interactions <- sponge(gene_expr = gene_expr,
                               mir_expr = mir_expr,
                               mir_interactions = genes_miRNA_candidates)
  mscor_null_model <- sponge_build_null_model(number_of_datasets = 100,
                                              number_of_samples = dim(gene_expr)[1])
  sponge_result <- sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                                   null_model = mscor_null_model)
  sponge_result_sig <- sponge_result[sponge_result$p.adj<=0.05,]
  sponge_result_sig$genepairs_1 <- paste0(sponge_result_sig$geneA,'|',sponge_result_sig$geneB)
  sponge_result_sig$genepairs_2 <- paste0(sponge_result_sig$geneB,'|',sponge_result_sig$geneA)
  write.csv(sponge_result_sig, paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/integration/',project_name,'-',disease_name,'_sponge.csv'), row.names = F)

  #JAMI
  mir_exp <- mirna
  gene_exp <- mrna
  df_lst <- list()
  for (i in 1:dim(dict)[1]){
    #i=1
    gene_pair <- dict[i,][[2]]
    tmp <- as.data.frame(t(combn(gene_pair,2)))
    tmp$V3<- dict[i,][[1]]
    df_lst[[i]] <- tmp
  }
  gene_mir_interactions_triplets <- Reduce(rbind, df_lst)
  names(gene_mir_interactions_triplets) <- c('geneA','geneB','mirnas')
  RJAMI::test_jvm()
  RJAMI::jami_settings(pvalueCutOff = 0.05)
  RJAMI::jami_settings(tripleFormat = FALSE)
  result <- RJAMI::jami(gene_miRNA_interactions = gene_mir_interactions_triplets,
                        gene_expr = gene_exp,
                        mir_expr = mir_exp)
  rjami_result <- result$result[,1:5]
  rjami_result_sig <- rjami_result[rjami_result$p.value <=0.05,]
  rjami_result_sig$triplets <- paste0(rjami_result_sig$miRNA,'|',rjami_result_sig$Source, '|', rjami_result_sig$Target)
  write.csv(sponge_result_sig, paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/integration/',project_name,'-',disease_name,'_jami.csv'), row.names = F)

  # our results
  our_result <- as.data.frame(read.csv(paste0(project_name,'-',disease_name,'/',project_name,'-',disease_name,'_finalpairs.csv')))
  cand_pair <- Reduce(rbind,stringr::str_split(our_result$cand.ceRNA,' '))
  our_result <- cbind(our_result[,1:2],cand_pair)
  our_result <- our_result[,-2]
  names(our_result)[2:3] <- c("geneA","geneB")
  our_result$triplets <- paste0(our_result$miRNA,'|', our_result$geneA, '|', our_result$geneB)
  our_result$genepairs <- paste0(our_result$geneA, '|', our_result$geneB)

  # integrate
  sponge_integrate <- c(intersect(our_result$genepairs,sponge_result_sig$genepairs_1),intersect(our_result$genepairs,sponge_result_sig$genepairs_2))
  rjami_integrate <- intersect(our_result$triplets,rjami_result_sig$triplets)
  our_result$sponge <- '-'
  our_result$rjami <- '-'
  our_result$sponge[our_result$genepairs%in%sponge_integrate] <- 'yes'
  our_result$rjami[our_result$triplets%in%rjami_integrate] <- 'yes'
  write.csv(our_result, paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/integration/',project_name,'-',disease_name,'_integrate.csv'), row.names = F)

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))

  message('\u2605\u2605\u2605 All analyses has completed! \u2605\u2605\u2605')
}
