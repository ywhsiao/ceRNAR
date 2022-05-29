#'
#' @name ceRNAIntergate
#' @title Integration of the possible ceRNA pairs among published tools
#' @description A function to integrate the possible ceRNA pairs that are found
#' by ceRNAR algorithm with those from other tools, such as SPONGE (List et al.,
#' 2019) and RJAMI (Hornakova et al.,)2018.
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
  genes_miRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
    gene_expr = gene_expr,
    mir_expr = mir_expr,
    mir_predicted_targets = as.matrix(d))

  ## minor revision of sponge function from SPONGE: rm() issue on original code ##
  sponge <- function(gene_expr,
                     mir_expr,
                     mir_interactions = NULL,
                     log.level = "ERROR",
                     log.every.n = 1e5,
                     log.file = NULL,
                     selected.genes = NULL,
                     gene.combinations = NULL,
                     each.miRNA = FALSE,
                     min.cor = 0.1,
                     parallel.chunks = 1e3,
                     random_seed = NULL,
                     result_as_dt = FALSE){
    check_and_convert_expression_data <- function(expr_data){

      if(is(expr_data, "big.matrix.descriptor")){
        expr_data <- attach.big.matrix(expr_data)
        if(length(mwhich(expr_data, 1:ncol(expr_data), NA, 'eq', 'OR')) > 0){
          stop("NA values found in expression data. Can not proceed")
        }
        else return(expr_data)
      }

      if(!is(expr_data, "matrix") && !is(expr_data, "ExpressionSet")){
        stop("expression matrix has to be of class matrix, ExpressionSet or big.matrix.descriptor")
      }

      if(is(expr_data, "ExpressionSet")){
        expr_data <- t(exprs(expr_data))
      }

      #check for NA values that make elasticnet crash
      if(anyNA(expr_data)) stop("NA values found in expression data. Can not proceed")

      return(expr_data)
    }

    genes_pairwise_combinations <- function(number.of.genes){
      #t(combnPrim(number.of.genes, 2))
      t(combn_prim(number.of.genes, 2))
    }

    split_rows <- function(x, ...){
      it <- idiv(nrow(x), ...)
      i <- 1L
      nextEl <- function() {
        n <- as.integer(nextElem(it))
        j <- i
        i <<- i + n
        x[seq(j, length = n), , drop = FALSE]
      }
      object <- list(nextElem = nextEl)
      class(object) <- c("abstractiter", "iter")
      object
    }

    processChunk <- function(gene_combis, attached_gene_expr, attached_mir_expr, mir_interactions,
                             all_mirs, each.miRNA, min.cor){
      if(is.null(mir_interactions))
        mir_intersect <- all_mirs

      foreach(geneA_idx = gene_combis$geneA_idx,
              geneB_idx = gene_combis$geneB_idx,
              geneA = gene_combis$geneA,
              geneB = gene_combis$geneB,
              .export = c("compute_pcor", "fn_get_shared_miRNAs"),
              .combine=function(...) rbindlist(list(...))) %do% {

                source_expr <- attached_gene_expr[,geneA_idx]
                target_expr <- attached_gene_expr[,geneB_idx]

                #check correlation
                dcor <- cor(source_expr, target_expr)

                if(is.na(dcor)) return(NULL)

                if(!is.null(min.cor)){
                  if(dcor < min.cor)
                    return(NULL)
                }

                #check if miRNA interaction information is provided, otherwise we
                #consider ALL miRNAs in each comparison
                if(!is.null(mir_interactions)){
                  mir_intersect <- fn_get_shared_miRNAs(geneA, geneB,
                                                        mir_interactions)

                  #check if shared miRNAs are in expression matrix
                  if(length(setdiff(mir_intersect, all_mirs)) > 0){
                    mir_intersect <- intersect(mir_intersect, all_mirs)
                  }

                  #check if there are actually any shared mirnas
                  if(length(mir_intersect) == 0){
                    return(NULL)
                  }
                }

                if(each.miRNA){
                  result <- foreach(mirna = mir_intersect,
                                    .export = c("compute_pcor"),
                                    .combine = function(...) rbindlist(list(...)),
                                    .inorder = TRUE) %do%{
                                      m_expr <- attached_mir_expr[,which(all_mirs == mirna)]
                                      compute_pcor(source_expr, target_expr, m_expr,
                                                   geneA, geneB, dcor)
                                    }
                  result$miRNA <- mir_intersect
                  return(result)
                }
                else{
                  m_expr <- attached_mir_expr[,which(all_mirs %in% mir_intersect)]
                  compute_pcor(source_expr, target_expr, m_expr,
                               geneA, geneB, dcor)
                }
              }
    }

    compute_pcor <- function(source_expr, target_expr, m_expr,
                             geneA, geneB, dcor){

      pcor <- tryCatch({
        pcor.test(source_expr, target_expr, m_expr)
      }, warning = function(w) {
        logdebug(w)
        suppressWarnings(pcor.test(source_expr, target_expr, m_expr))
      }, error = function(e) {
        logerror(e)
        return(NULL)
      })

      if(is.null(pcor)) return(NULL)

      list(geneA = geneA,
           geneB = geneB,
           df = pcor$gp,
           cor =  dcor,
           pcor = pcor$estimate,
           mscor = dcor - pcor$estimate
      )
    }

    fn_get_shared_miRNAs <- function(geneA, geneB, mir_interactions){

      source_sign_mirs <- mir_interactions[[geneA]]

      if(!is.null(source_sign_mirs)){
        source_sign_mirs <- as.character(source_sign_mirs$mirna)
      }
      target_sign_mirs <- mir_interactions[[geneB]]

      if(!is.null(target_sign_mirs)){
        target_sign_mirs <- as.character(target_sign_mirs$mirna)
      }
      mir_intersect <- intersect(source_sign_mirs, target_sign_mirs)

      return(mir_intersect)
    }

    if(!is.null(log.file))
      addHandler(writeToFile, file=log.file, level=log.level)
    else{
      basicConfig(level = log.level)
    }

    #handle bigmemory objects and check expression matrices
    foreach_packages <- c("logging", "ppcor", "foreach",
                          "iterators", "data.table")

    if(is(gene_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
      loginfo("Detected gene expression big matrix descriptor")
      gene_expr_big_memory <- TRUE
      gene_expr_description <- gene_expr
      gene_expr <- check_and_convert_expression_data(gene_expr)
      foreach_packages <- union(foreach_packages, "bigmemory")
    }
    else{
      gene_expr_big_memory <- FALSE
      gene_expr <- check_and_convert_expression_data(gene_expr)
    }

    if(is(mir_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
      loginfo("Detected miRNA expression big matrix descriptor")
      mir_expr_big_memory <- TRUE
      mir_expr_description <- mir_expr
      mir_expr <- check_and_convert_expression_data(mir_expr)
      foreach_packages <- union(foreach_packages, "bigmemory")
    }
    else{
      mir_expr_big_memory <- FALSE
      mir_expr <- check_and_convert_expression_data(mir_expr)
    }

    if(is.null(mir_interactions)){
      logwarn("No information on miRNA gene interactions was provided,
                all miRNAs will be considered and runtime will likely explode.")
      genes <- colnames(gene_expr)
    }
    else{
      #filter out genes without miR interactions
      mir_interactions <- Filter(Negate(is.null), mir_interactions)

      #names of genes for which we have expr and miRNA data
      genes <- intersect(colnames(gene_expr), names(mir_interactions))
    }

    #now only compute for selected genes
    if(is.null(selected.genes)){
      sel.genes <- genes
    }
    else{
      available.selected.genes <- intersect(selected.genes, genes)
      if(length(available.selected.genes) == 0){
        stop("None of the selected genes is found in the data")
      }
      else if(length(available.selected.genes) < length(selected.genes)){
        warning(paste("Some genes are not found in the data:",
                      paste(setdiff(selected.genes, genes), collapse=","), sep=""))
      }
      sel.genes <- available.selected.genes
    }

    #use indices for faster access
    genes.as.indices <- FALSE

    #make sure sel.genes order is used
    gene_expr <- gene_expr[,sel.genes]

    #for getting indices of things
    all_mirs <- colnames(mir_expr)
    all_genes <- colnames(gene_expr)

    #all pairwise combinations of selected genes
    if(is.null(gene.combinations)){
      loginfo("Computing all pairwise combinations of genes")
      #consider only genes that have miRNA interactions

      gene.combinations <-
        genes_pairwise_combinations(length(sel.genes))

      gene.combinations <- data.frame(gene.combinations,
                                      all_genes[gene.combinations[,1]],
                                      all_genes[gene.combinations[,2]],
                                      stringsAsFactors = FALSE)

    } else{
      if(ncol(gene.combinations) > 2)
        stop("gene.combinations is expected to have two columns with gene identifiers")

      colnames(gene.combinations) <- c("geneA", "geneB")

      gene.combinations$geneA <- as.character(gene.combinations$geneA)
      gene.combinations$geneB <- as.character(gene.combinations$geneB)

      valid_genes_A <- which(gene.combinations$geneA %in% all_genes)
      valid_genes_B <- which(gene.combinations$geneB %in% all_genes)
      valid_genes <- intersect(valid_genes_A, valid_genes_B)

      if(length(valid_genes) > 0){
        gene.combinations <- gene.combinations[valid_genes,]

        gene.combinations <- data.frame(
          which(all_genes %in% gene.combinations$geneA),
          which(all_genes %in% gene.combinations$geneB),
          gene.combinations,
          stringsAsFactors = FALSE)
      } else stop("No valid gene combinations selected")
    }
    colnames(gene.combinations) <- c("geneA_idx", "geneB_idx", "geneA", "geneB")

    loginfo("Beginning SPONGE run...")

    if(!gene_expr_big_memory) gene_expr_description <- gene_expr
    if(!mir_expr_big_memory) mir_expr_description <- mir_expr

    # rm(gene_expr)
    # rm(mir_expr)

    num_of_samples <- nrow(gene_expr)
    num_of_tasks <- min(max(1, ceiling(nrow(gene.combinations) / 1000)),
                        parallel.chunks)

    SPONGE_result <- foreach(
      gene_combis =
        split_rows(gene.combinations,
                   chunks = num_of_tasks),
      i = iterators::icount(),
      .combine = function(...)
        rbindlist(list(...)),
      .multicombine = TRUE,
      .packages = foreach_packages,
      .export = c("fn_get_shared_miRNAs", "processChunk", "compute_pcor"),
      .options.RNG = random_seed
    ) %dorng% {
      if (!is.null(log.file))
        addHandler(writeToFile, file = log.file, level = log.level)
      else{
        basicConfig(level = log.level)
      }

      loginfo(paste("SPONGE: worker is processing chunk: ", i, sep = ""))

      #attach bigmemory objects if necessary. avoid using names of actual
      #matrix objects because they would then be exported to the workers
      if (gene_expr_big_memory)
        attached_gene_expr <- attach.big.matrix(gene_expr_description)
      else
        attached_gene_expr <- gene_expr_description

      if (mir_expr_big_memory)
        attached_mir_expr <- attach.big.matrix(mir_expr_description)
      else
        attached_mir_expr <- mir_expr_description

      #if(require(pryr)) logdebug(paste("current memory used by worker:", pryr::mem_used()))

      result <-
        processChunk(
          gene_combis,
          attached_gene_expr,
          attached_mir_expr,
          mir_interactions,
          all_mirs,
          each.miRNA,
          min.cor
        )

      loginfo(paste("SPONGE finished chunk:", i, "of", num_of_tasks))
      if(is.null(result)) return(list())
      else return(result)
    }
    loginfo("SPONGE completed successfully. Returning results.")

    if(result_as_dt) return(SPONGE_result)
    else return(as.data.frame(SPONGE_result))
  }
  ceRNA_interactions <- sponge(gene_expr = gene_expr,
                               mir_expr = mir_expr,
                               mir_interactions = genes_miRNA_candidates)
  precomputed_cov_matrices <- SPONGE::precomputed_cov_matrices
  mscor_null_model <- SPONGE::sponge_build_null_model(number_of_datasets = 100,
                                                      number_of_samples = dim(gene_expr)[1])
  sponge_result <- SPONGE::sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                                   null_model = mscor_null_model)
  sponge_result_sig <- sponge_result[sponge_result$p.adj<=0.05,]
  sponge_result_sig$genepairs_1 <- paste0(sponge_result_sig$geneA,'|',sponge_result_sig$geneB)
  sponge_result_sig$genepairs_2 <- paste0(sponge_result_sig$geneB,'|',sponge_result_sig$geneA)

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
