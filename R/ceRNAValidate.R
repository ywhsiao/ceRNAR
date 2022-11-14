#'
#' @name ceRNAValidate
#' @title Externally experimental validation for the potential ceRNA pairs
#' @description A function to validate the potential ceRNA pairs based on the miRSponge database
#' (http://www.bio-bigdata.net/miRSponge)
#'
#' @import utils
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#'
#' @returns a dataframe output
#' @export
#'
#' @examples
#' ceRNAValidate(
#' path_prefix = NULL,
#' project_name = 'demo',
#' disease_name = 'DLBC'
#' )
#'
#'

ceRNAValidate <- function(path_prefix = NULL,
                          project_name = 'demo',
                          disease_name = 'DLBC'){

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

  if(!dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/external_validation/'))){
    dir.create(paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/external_validation/'))
  }

  message('\u25CF Step5: Dowstream Analyses - External validation')

  datapreparing <- function(file){
    candidate_ceRNA <- as.data.frame(utils::read.csv(file))
    candidate_ceRNA_list <- candidate_ceRNA[,c(1,2)]
    candidate_ceRNA_list$miRNA <- unlist(candidate_ceRNA_list$miRNA)
    candidate_ceRNA_list <- tidyr::separate(candidate_ceRNA_list,cand.ceRNA, c("ceRNA1", "ceRNA2"), " ")
    candidate_ceRNA_list$miRNA <- gsub('hsa-','',candidate_ceRNA_list$miRNA)
    candidate_ceRNA_list$miRNA_short <- gsub('-5p','',candidate_ceRNA_list$miRNA)
    candidate_ceRNA_list$miRNA_short <- gsub('-3p','',candidate_ceRNA_list$miRNA_short)
    candidate_ceRNA_list$miRNA_short1 <- gsub('[[:alpha:]]{1}$', '', candidate_ceRNA_list$miRNA_short)
    candidate_ceRNA_list
  }
  pair <-datapreparing(paste0(path_prefix, project_name,'-',disease_name,'/',project_name,'-',disease_name,'_finalpairs.csv'))
  total_findings <- dim(pair)[1]
  pair_lst <-  dplyr::group_split(pair,miRNA_short)
  # miRSponge clean data
  miRSponge_exp <- get0("ext_val", envir = asNamespace("ceRNAR"))
  miRSponge_exp$all <- paste0(miRSponge_exp$miRNA.Name, '|', miRSponge_exp$Target.Name, '|', miRSponge_exp$Experimental.method)
  miRSponge_exp$miRNA_short <- gsub('-3p|-5p','',miRSponge_exp$miRNA)
  miRSponge_exp <- miRSponge_exp[,c(1,3,6,7)]

  # search based on miRNA
  cancer_df = pair_lst
  search_based_on_miRNA <- function(cancer_df){
    miRNA_lst <- list()
    for (i in 1:length(cancer_df)){
      #i=1
      if (dim(miRSponge_exp[miRSponge_exp$miRNA_short==unique(unlist(cancer_df[[i]]['miRNA_short'])),])[1]!=0 | dim(miRSponge_exp[miRSponge_exp$miRNA_short==unique(unlist(cancer_df[[i]]['miRNA_short1'])),])[1]!=0){
        miRNA_lst[[i]]<- c(i,unique(unlist(cancer_df[[i]]['miRNA_short'])), unique(unlist(cancer_df[[i]]['miRNA_short1'])))
      }
    }
    miRNA_df <- as.data.frame(Reduce(rbind, miRNA_lst))
    names(miRNA_df) <- c('index', 'miRNA_short', 'miRNA_short1')
    row.names(miRNA_df) <- NULL
    miRNA_df
  }
  miRNA <- search_based_on_miRNA(pair_lst)
  matched1 <- merge(miRNA, miRSponge_exp, by='miRNA_short', all.x = TRUE)
  matched2 <- merge(miRNA, miRSponge_exp, by.x='miRNA_short1', by.y = 'miRNA_short', all.x = TRUE)
  final <- unique(rbind(matched1, matched2))
  final_lst <- dplyr::group_split(final, index)
  # search for each gene
  idx <- sort(unique(final$index))
  proof_lst <- list()
  for (j in 1:length(idx)){
    #j=1
    df <- as.data.frame(cancer_df[[as.numeric(idx[j])]])
    exp <- final[final$index==as.numeric(idx[j]),]
    for (i in 1:dim(df)[1]){
      if (df[i,'ceRNA2'] %in% exp$Target.Name & df[i,'ceRNA1'] %in% exp$Target.Name){
        df$evidence[i] <- 'exp_proof_inboth'
        df$evidence1[i] <- 'ceRNA1|ceRNA2'
      }else if (df[i,'ceRNA1'] %in% exp$Target.Name){
        df$evidence[i] <- 'exp_proof_inceRNA1'
        df$evidence1[i] <- 'ceRNA1'
      }else if (df[i,'ceRNA2'] %in% exp$Target.Name){
        df$evidence[i] <- 'exp_proof_inceRNA2'
        df$evidence1[i] <- 'ceRNA2'
      }else{
        df$evidence[i] <- '-'
        df$evidence1[i] <- '-'
      }

    }
    proof_lst[[j]] <- df
  }
  proof_df <- Reduce(rbind, proof_lst)
  with_evidence <- proof_df[proof_df$evidence!='-',]
  with_evidence$exp_both <- '-'
  # link to evidence
  for (i in 1:dim(with_evidence)[1]){
    #print(i)
    #i =8
    `%!in%` <- Negate(`%in%`)
    if (with_evidence$evidence1[i] == 'ceRNA1|ceRNA2'){
      ceRNA1_tmp1 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA1[i],]
      ceRNA1_tmp2 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short1[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA1[i],]
      ceRNA1_tmp <- rbind(ceRNA1_tmp1,ceRNA1_tmp2)
      ceRNA2_tmp1 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA2[i],]
      ceRNA2_tmp2 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short1[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA2[i],]
      ceRNA2_tmp <- rbind(ceRNA2_tmp1,ceRNA2_tmp2)
      with_evidence$exp_both[i] <- paste0(ceRNA1_tmp$all[1], '|||', ceRNA2_tmp$all[1])
    } else if(with_evidence$evidence1[i] == 'ceRNA1'){
      ceRNA1_tmp1 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA1[i],]
      ceRNA1_tmp2 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short1[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA1[i],]
      ceRNA1_tmp <- rbind(ceRNA1_tmp1,ceRNA1_tmp2)
      for (j in 1:dim(ceRNA1_tmp)[1]){
        if (paste0("exp_ceRNA1_", j) %!in% colnames(with_evidence)){
          with_evidence[i, ncol(with_evidence) + 1] <- ceRNA1_tmp[j,'all']           # Append new column
          colnames(with_evidence)[ncol(with_evidence)] <- paste0("exp_ceRNA1_", j)
        }else{
          with_evidence[i, paste0("exp_ceRNA1_", j)] <- ceRNA1_tmp[j,'all']
        }
      }
    }else if(with_evidence$evidence1[i] == 'ceRNA2'){
      ceRNA2_tmp1 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA2[i],]
      ceRNA2_tmp2 <- miRSponge_exp[miRSponge_exp$miRNA_short == with_evidence$miRNA_short1[i] & miRSponge_exp$Target.Name == with_evidence$ceRNA2[i],]
      ceRNA2_tmp <- rbind(ceRNA2_tmp1,ceRNA2_tmp2)
      for (j in 1:dim(ceRNA2_tmp)[1]){
        if (paste0("exp_ceRNA2_", j) %!in% colnames(with_evidence)){
          with_evidence[i, ncol(with_evidence) + 1] <- ceRNA2_tmp[j,'all']           # Append new column
          colnames(with_evidence)[ncol(with_evidence)] <- paste0("exp_ceRNA2_", j)
        }else{
          with_evidence[i, paste0("exp_ceRNA2_", j)] <- ceRNA2_tmp[j,'all']
        }
      }
    }


  }
  with_evidence[is.na(with_evidence)] <- '-'
  common_pairs <- dim(with_evidence)[1]
  utils::write.csv(with_evidence, paste0(path_prefix, project_name,'-',disease_name,'/04_downstreamAnalyses/external_validation/',project_name,'-',disease_name,'_with_target_exp_evidence.csv'), row.names = FALSE)
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
  message('\u2605\u2605\u2605 All analyses has completed! \u2605\u2605\u2605')

  return(with_evidence)
}
