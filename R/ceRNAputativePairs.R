#'
#' @name ceRNAputativePairs
#' @title Extraction of putative mRNA-miRNA pairs
#' @description A function to obtain putative mRNA-miRNA pairs from several databases
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param filtering three different filtering criteria, including strict,
#' moderate and less. (Default: less)
#'
#' @export
#'
#' @examples
#' ceRNAputativePairs(
#' path_prefix = '~/',
#' project_name ='demo',
#' disease_name = 'DLBC',
#' filtering = 'less'
#' )
#'

ceRNAputativePairs <- function(path_prefix,
                               project_name = 'demo',
                               disease_name = 'DLBC',
                               filtering = 'less'){

  if (!stringr::str_detect(path_prefix, '/')){
    path_prefix <- paste0(path_prefix, '/')
  }

  time1 <- Sys.time()

  ## match to interaction database
  message('\u25CF Step 2: Obtaining putative mRNA-miRNA pairs')
  # load processed mRNA and miRNA data
  annot_cdRNA_unique <- as.data.frame(data.table::fread(paste0(path_prefix, '/', project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mrna.csv')))
  row.names(annot_cdRNA_unique) <- annot_cdRNA_unique[,1]
  annot_cdRNA_unique <- annot_cdRNA_unique[,-1]
  miRNA_with_precurer <- as.data.frame(data.table::fread(paste0(path_prefix, '/', project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mirna.csv')))
  row.names(miRNA_with_precurer) <- miRNA_with_precurer[,1]
  miRNA_with_precurer <- miRNA_with_precurer[,-1]

  # miRNA-mRNA validation
  target <- get0("mirna_mrna_pairsdb", envir = asNamespace("ceRNAR"))
  if (filtering == 'strict'){
    target.t.val <- target[target$evidence_levels == "Strong" & target$total_counts == 7,]
    message('\u2605 Filtering: strict')
  }else if (filtering == 'moderate') {
    target.t.val <- target[target$evidence_levels == "Strong" | target$total_counts == 7,]
    message('\u2605 Filtering: moderate')
  }else if (filtering == 'less'){
    target.t.val <- target[target$evidence_levels == "Strong" | target$total_counts >= 6,]
    message('\u2605 Filtering: less')
  }


  miRNA_f <- intersect(unique(target.t.val[,c('miRNA_names')]),row.names(miRNA_with_precurer))
  target.t.val <- target.t.val[target.t.val$miRNA_names %in% miRNA_f,]  # overlap mirna between target database and gse miRNA profiles
  target.t.val <- target.t.val[target.t.val$gene_names %in% rownames(annot_cdRNA_unique),]   # overlap gene between target database and gse mRNA profiles

  message('\u2605 Total putative mRNA-miRNA pairs are ', sum(rownames(annot_cdRNA_unique) %in% target.t.val$gene_names), '!')

  if (dir.exists(paste0(path_prefix, '/', project_name,'-',disease_name,'/02_potentialPairs')) == FALSE){
    dir.create(paste0(path_prefix, '/', project_name,'-',disease_name,'/02_potentialPairs'))
  }

  # input for ceRNA identification step
  dictionary <- c()
  for (mirna in miRNA_f) {
    #mirna = 'hsa-let-7a-3p'
    tmp_list <- list(target.t.val[target.t.val$miRNA_names == mirna,2])
    if(length(unlist(tmp_list[[1]])) >1){
      dictionary <- rbind(dictionary,c(mirna,tmp_list))
    }
  }
  saveRDS(dictionary, paste0(path_prefix, '/', project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_MirnaTarget_dictionary.rds'))
  message('(\u2714) Putative results have been created and stored!')
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')

}
