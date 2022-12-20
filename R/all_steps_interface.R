#'
#' @name All_steps_interface
#' @title A function for conducting the algorithm through the toy data
#' @description A function to allow users to conducting the
#'
#' @import foreach
#' @import utils
#' @import tidyverse
#' @importFrom stats pnorm
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default: demo)
#' @param disease_name the abbreviation of disease that users are interested in
#' (default: DLBC)
#' @param gene_exp location of gene expression data (default: gene_exp)
#' @param mirna_exp location of miRNA expression data (default: mirna_exp)
#' @param surv_data location of survival data (default: surv_data)
#' @param filtering three different filtering criteria, including strict,
#' moderate and less. (Default: less)
#' @param window_size the number of samples for each window (default:10)
#' @param cor_method selection of correlation methods, including pearson and
#' spearman (default: pearson)
#' @param cor_threshold_peak peak threshold of correlation value between 0 and 1
#' (default: 0.85)
#'
#' @returns a dataframe object
#' @export
#' @export
#'
#' @examples
#' data(gene_exp)
#' data(mirna_exp)
#' data(surv_data)
#' All_steps_interface(
#' path_prefix = NULL,
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' gene_exp = gene_exp,
#' mirna_exp = mirna_exp,
#' surv_data = surv_data,
#' filtering = 'less',
#' window_size = 10,
#' cor_method = 'pearson',
#' cor_threshold_peak = 0.85)

All_steps_interface <- function(path_prefix = NULL,
                                project_name = 'demo',
                                disease_name = 'DLBC',
                                gene_exp = gene_exp,
                                mirna_exp = mirna_exp,
                                surv_data = surv_data,
                                filtering = 'less',
                                window_size = 10,
                                cor_method = 'pearson',
                                cor_threshold_peak = 0.85){
  if (is.null(path_prefix)){
    path_prefix <- fs::path_home()
  }else{
    path_prefix <- path_prefix
  }

  if (!stringr::str_detect(path_prefix, '/$')){
    path_prefix <- paste0(path_prefix, '/')
  }
  # import data
  ceRNACustomize <- function(path_prefix = NULL,
                             project_name = 'demo',
                             disease_name = 'DLBC',
                             gene_exp = gene_exp,
                             mirna_exp = mirna_exp,
                             surv_data = surv_data){
    if (is.null(path_prefix)){
      path_prefix <- fs::path_home()
    }else{
      path_prefix <- path_prefix
    }

    if (!stringr::str_detect(path_prefix, '/$')){
      path_prefix <- paste0(path_prefix, '/')
    }

    if (dir.exists(paste0(path_prefix, project_name,'-',disease_name)) == FALSE){
      dir.create(paste0(path_prefix, project_name,'-',disease_name))
    }

    if (dir.exists(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata')) == FALSE){
      dir.create(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata'))
    }

    time1 <- Sys.time()
    message('\u25CF Step 1 for Customized data: Checking the data...')
    # check input
    if (is.null(project_name)|is.null(disease_name)|is.null(gene_exp)|is.null(mirna_exp)){
      stop("project_name/disease_name/gene_exp/mirna_exp has not provided!")
    }else{
      if (!is(gene_exp, "character")){
        exp <- gene_exp
      }else{
        exp <- as.data.frame(data.table::fread(gene_exp,header = TRUE))
        row.names(exp) <- exp[,1]
        exp <- exp[,-1]
      }
      if (!is(mirna_exp, "character")){
        mirna <- mirna_exp
      }else{
        mirna <- as.data.frame(data.table::fread(mirna_exp,header = TRUE))
        row.names(mirna) <- mirna[,1]
        mirna <- mirna[,-1]
      }

    }
    if (!is.null(surv_data)){
      if (!is(surv_data, "character")){
        surv <- surv_data
      }else{
        surv <- as.data.frame(data.table::fread(surv_data,header = TRUE))
      }

      message('(\u2714) Input data involve mRNA, miRNA and survival data.')
    }else{
      message('(\u2714) Input data ONLY involve mRNA and miRNA.')
    }

    # check geneid
    #gtf_df <- ceRNAR:::gencode_v22_annot
    ensem2symbol <- get0("ensem2symbol", envir = asNamespace("ceRNAR"))
    genecode_v22 <- unique(ensem2symbol$gene_name)
    exp_gene_name <- row.names(exp)
    if (sum(exp_gene_name %in% genecode_v22) == 0){
      stop("Gene names in this dataset ate not matched!")
    }else{
      message('(\u2714) All the gene names are matched!')
    }

    # check mirnaid
    # ID_converter <- ceRNAR:::hsa_pre_mature_matching
    ID_converter <- get0("hsa_pre_mature_matching", envir = asNamespace("ceRNAR"))
    mature_mirna <- unique(ID_converter$Mature_ID)
    mirna_name <- row.names(mirna)
    if (sum(mirna_name %in% mature_mirna) == 0){
      stop("Mature miRNA names in this dataset ate not matched!")
    }else{
      message('(\u2714) All the mature miRNA names are matched!')
    }


    # check survival data
    if (sum(colnames(surv) %in% c('OS','OS.time')) == 0){
      stop("Survival data did not contain 'OS' and 'OS.time' columns!")
    }else{
      message("(\u2714) Survival data contains 'OS' and 'OS.time' columns!")
    }
    # check whether the sample id is matched
    if (length(row.names(surv)) != length(colnames(exp)) && length(row.names(surv)) != length(colnames(mirna))){
      if (sum(row.names(surv) %in% colnames(exp)) == 0 && sum(row.names(surv) %in% colnames(mirna)) ==0 ){
        stop("(\u2714) Sample ids are not matched!")
      }else{
        stop("The number of samples are not euqal!")
      }
    }else{
      message("(\u2714) Sample ids are matched!")
    }

    data.table::fwrite(as.data.frame(exp),paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mrna.csv'), row.names = TRUE)
    data.table::fwrite(as.data.frame(mirna),paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mirna.csv'), row.names = TRUE)
    data.table::fwrite(surv, paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_survival.csv'), row.names = TRUE)

    time2 <- Sys.time()
    diftime <- difftime(time2, time1, units = 'min')
    message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
    message('\u2605\u2605\u2605 All the files are checked for further next step! \u2605\u2605\u2605')
  }

  ceRNACustomize(path_prefix = path_prefix,
                 project_name = 'demo',
                 disease_name = 'DLBC',
                 gene_exp = gene_exp,
                 mirna_exp = mirna_exp,
                 surv_data = surv_data)

  # putative pairs
  ceRNAputativePairs <- function(path_prefix = NULL,
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

  ceRNAputativePairs(path_prefix = path_prefix,
                     project_name = project_name,
                     disease_name = disease_name,
                     filtering = filtering)

  # method: pair filtering +segment clustering
  ceRNAMethod <- function(path_prefix = NULL,
                          project_name = 'demo',
                          disease_name = 'DLBC',
                          window_size = 10,
                          cor_method = 'pearson',
                          cor_threshold_peak = 0.85){

    if (is.null(path_prefix)){
      path_prefix <- fs::path_home()
    }else{
      path_prefix <- path_prefix
    }

    if (!stringr::str_detect(path_prefix, '/$')){
      path_prefix <- paste0(path_prefix, '/')
    }

    # ceRNApairfiltering
    ceRNApairFilering <- function(path_prefix = NULL,
                                  project_name = 'demo',
                                  disease_name = 'DLBC',
                                  window_size = 10,
                                  cor_method = 'pearson'){

      if (is.null(path_prefix)){
        path_prefix <- fs::path_home()
      }else{
        path_prefix <- path_prefix
      }

      if (!stringr::str_detect(path_prefix, '/$')){
        path_prefix <- paste0(path_prefix, '/')
      }

      time1 <- Sys.time()
      message('\u25CF Step 3: Filtering putative mRNA-miRNA pairs using sliding window approach')

      # setwd(paste0(project_name,'-',disease_name))
      # import example data & putative pairs
      dict <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_MirnaTarget_dictionary.rds'))
      mirna <- data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mirna.csv')),row.names = 1)
      mrna <- data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mrna.csv')),row.names = 1)

      mirna_total <- unlist(dict[,1])
      message(paste0('\u2605 total miRNA: ', length(mirna_total)))


      slidingWindow <- function(window_size, mirna_total, cor_method){
        num_core <- parallel::detectCores()
        if (num_core <=2){
          cl <- parallel::makeCluster(num_core)
          message(paste0('Number of CPU used: ', num_core))
        }else{
          cl <- parallel::makeCluster(num_core-1)
          message(paste0('Number of CPU used: ', num_core-1))
        }
        doParallel::registerDoParallel(cl)

        # reate a cluster
        #message('\u2605 Number of computational cores: ', num_workers, '.')
        parallel_d <- foreach(mir=1:length(mirna_total), .export = c('dict','mirna', 'mrna'))  %dopar%  {
          #mir = 50
          mir = mirna_total[mir]
          gene <- as.character(data.frame(dict[dict[,1]==mir,][[2]])[,1])
          gene <- intersect(gene,rownames(mrna))
          gene_mir <- data.frame("miRNA"=t(mirna[rownames(mirna)==mir,]), t(mrna[gene,]))

          names(gene_mir) <- c("miRNA",names(gene_mir)[-1])
          names(gene_mir) <- gsub("\\.","\\-",names(gene_mir))

          w <- window_size
          N <- dim(gene_mir)[1] # total samples
          gene_pair <- utils::combn(gene,2)
          gene_pair_index <- rbind(match(gene_pair[1,], names(gene_mir)),match(gene_pair[2,], names(gene_mir)))
          total_pairs <- choose(length(gene),2)
          #message(paste0("\u2605 ",mir,"'s total potential ceRNA-miRNA pairs: ",total_pairs))

          getcorr <- function(r,s){
            #r=2
            #s=3
            y <- gene_mir[,c(1,r,s)]
            y <- y[order(y$miRNA),]
            corr <- zoo::rollapply(y, width=w, function(x) stats::cor(x[,2],x[,3],method = cor_method), by.column=FALSE)
            miRNA <- zoo::rollapply(y, width=w, function(x) mean(x[,1]), by.column=FALSE)
            data <- data.frame(cbind(miRNA,corr))
            data <- data[order(data$miRNA),]
            data$corr
          }
          cor_all <- purrr::map2_dfc(gene_pair_index[1,1:total_pairs],gene_pair_index[2,1:total_pairs],getcorr) %>%
            suppressMessages()
          getordermiRNA <- function(){
            y <- gene_mir[,c(1,1,1)]
            y <- y[order(y$miRNA),]
            #corr <- zoo::rollapply(y, width=w, function(x) cor(x[,2],x[,3],method = cor_method), by.column=FALSE)
            miRNA <- zoo::rollapply(y, width=w, function(x) mean(x[,1]), by.column=FALSE)
            #data <- data.frame(cbind(miRNA,corr))
            #data <- data[order(data$miRNA),]
            #data$miRNA
          }
          win_miRNA <- getordermiRNA()
          triplet = data.frame("miRNA"=win_miRNA,"corr"=cor_all)
          if (dim(triplet)[2]==2){
            names(triplet)[2] <- gsub('\\.\\.\\.','\\.init',names(triplet)[2])
          }else{
            names(triplet)[2] <- gsub('\\.\\.\\.\\.','\\.init',names(triplet)[2])
            names(triplet) <- gsub('\\.\\.\\.\\.','\\.V',names(triplet))
          }
          triplet

        }
        parallel::stopCluster(cl)
        parallel_d
      }
      Realdata <- slidingWindow(window_size,mirna_total, 'pearson')
      saveRDS(Realdata,paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_pairfiltering.rds'))

      time2 <- Sys.time()
      diftime <- difftime(time2, time1, units = 'min')
      message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))

      message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')

    }

    ceRNApairFilering(path_prefix = path_prefix,
                      project_name = project_name,
                      disease_name = disease_name,
                      window_size = window_size,
                      cor_method = cor_method)

    # SegmentClustering + PeakMerging
    SegmentClusteringPlusPeakMerging <- function(path_prefix = NULL,
                                                 project_name = 'demo',
                                                 disease_name = 'DLBC',
                                                 cor_threshold_peak = 0.85,
                                                 window_size = 10){

      if (is.null(path_prefix)){
        path_prefix <- fs::path_home()
      }else{
        path_prefix <- path_prefix
      }

      if (!stringr::str_detect(path_prefix, '/$')){
        path_prefix <- paste0(path_prefix, '/')
      }

      time1 <- Sys.time()

      message('\u25CF Step 4: Clustering segments using CBS algorithm plus Mearging peaks')

      dict <- readRDS(paste0(path_prefix, project_name, '-', disease_name, '/02_potentialPairs/', project_name,'-', disease_name, '_MirnaTarget_dictionary.rds'))
      mirna <- data.frame(data.table::fread(paste0(path_prefix, project_name, '-', disease_name, '/01_rawdata/', project_name, '-', disease_name, '_mirna.csv')), row.names = 1)
      mrna <- data.frame(data.table::fread(paste0(path_prefix, project_name, '-', disease_name, '/01_rawdata/', project_name, '-', disease_name, '_mrna.csv')), row.names = 1)
      mirna_total <- unlist(dict[,1])
      d <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/', project_name,'-',disease_name,'_pairfiltering.rds'))

      sigCernaPeak <- function(index, d, cor_threshold_peak, window_size) {
        w <- window_size
        mir = mirna_total[index]
        gene <- as.character(data.frame(dict[dict[, 1] == mir,][[2]])[, 1])
        gene <- intersect(gene, rownames(mrna))
        gene_pair <- combn(gene, 2)
        total_pairs <- choose(length(gene), 2)

        num_core <- parallel::detectCores()
        if (num_core <=2){
          cl <- parallel::makeCluster(num_core)
          message(paste0('Number of CPU used: ', num_core))
        }else{
          cl <- parallel::makeCluster(num_core-1)
        }
        doParallel::registerDoParallel(cl)

        tmp <- foreach(p = 1:total_pairs, .combine = "rbind") %dopar%{
          #print(paste0("no_of_index:", index, "|", "no_of_pairs:",p))
          cand.ceRNA = c()
          location = list()
          r = gene_pair[1, p]
          s = gene_pair[2, p]
          triplet <- d[[index]][, c(1, p + 1)]
          names(triplet) <- c("miRNA", "corr")
          if (!any(duplicated(triplet$miRNA))){
            if (sum(is.na(triplet$corr)) == 0) {
              CNA.object <- DNAcopy::CNA(triplet$corr, rep(1, dim(triplet)[1]), triplet$miRNA)
              names(CNA.object) <- c("chrom", "maploc", paste("gene", r, "and", s))
              sink("/dev/null")
              result <- DNAcopy::segment(CNA.object)
              sink()
              if (sum(result$output$num.mark <= 3) >= 1) {
                tooshort <- which(result$output$num.mark <= 3)
                num.mark <- c(0, cumsum(result$output$num.mark),
                              data.table::last(cumsum(result$output$num.mark)))
                if (1 %in% diff(tooshort)) {
                  cc = 1
                  lag = 0
                  for (q in 1:(length(tooshort) - 1)) {
                    if (tooshort[q + 1] - tooshort[q] == 1) {
                      result$output[tooshort[q], "loc.end"] <- result$output[tooshort[q + 1], "loc.end"]
                      result$output[tooshort[q], "seg.mean"] <- t(matrix(result$output[tooshort[q]:tooshort[q + 1], "num.mark"])) %*% matrix(result$output[tooshort[q]:tooshort[q + 1], "seg.mean"])/sum(result$output[tooshort[q]:tooshort[q + 1], "num.mark"])
                      result$output[tooshort[q], "num.mark"] <- sum(result$output[tooshort[q]:tooshort[q + 1], "num.mark"])
                      result$output[tooshort[q + 1], ] <- result$output[tooshort[q],
                      ]
                      lag[cc] <- tooshort[q]
                      cc <- cc + 1
                    }
                  }
                  result$output <- result$output[-lag, ]
                  row.names(result$output) <- 1:dim(result$output)[1]
                  tooshort <- which(result$output$num.mark <= 3)
                }
                if (length(tooshort) >= 1) {
                  cc = 1
                  lag = c()
                  for (t in seq_along(tooshort)) {
                    long_seg <- which(result$output$num.mark > 3)
                    diff = abs(tooshort[t] - long_seg)
                    closest_seg <- long_seg[which(diff == min(diff))]
                    if (length(closest_seg) >= 2) {
                      b <- abs(result$output$seg.mean[closest_seg] - result$output$seg.mean[tooshort[t]])
                      closest_seg <- closest_seg[b == min(b)][1]
                    }
                    result$output[tooshort[t], "loc.start"] <- min(result$output[tooshort[t],  "loc.start"], result$output[closest_seg, "loc.start"])
                    result$output[tooshort[t], "loc.end"] <- max(result$output[tooshort[t], "loc.end"], result$output[closest_seg, "loc.end"])
                    result$output[tooshort[t], "seg.mean"] <- t(matrix(result$output[tooshort[t]:closest_seg, "num.mark"])) %*% matrix(result$output[tooshort[t]:closest_seg, "seg.mean"])/sum(result$output[tooshort[t]:closest_seg, "num.mark"])
                    result$output[tooshort[t], "num.mark"] <- sum(result$output[tooshort[t]:closest_seg, "num.mark"])
                    result$output[closest_seg, ] <- result$output[tooshort[t],]
                    lag[cc] <- tooshort[t]
                    cc <- cc + 1
                  }
                  result$output <- result$output[-lag, ]
                  row.names(result$output) <- 1:dim(result$output)[1]
                }
              }
              cand.corr <- c(-1, result$output$seg.mean, -1)
              peak.loc <- quantmod::findPeaks(cand.corr) - 2
              no_merg_loc <- c()
              no_merg_count <- 1
              if (sum(cand.corr[peak.loc + 1] > cor_threshold_peak) >= 2) {
                for (i in 1:(length(peak.loc) - 1)) {
                  if (is.na(sum(result$output[(peak.loc[i] + 1):(peak.loc[i + 1] - 1), "num.mark"]))) {
                    break
                  }
                  if (sum(result$output[(peak.loc[i] + 1):(peak.loc[i + 1] - 1), "num.mark"]) > w) {
                    no_merg_loc[no_merg_count] <- peak.loc[i]
                  }
                }
                peak.loc <- peak.loc[-which(peak.loc == no_merg_loc)]
                while (sum(cand.corr[peak.loc + 1] > cor_threshold_peak) >= 2) {
                  num.mark <- c(0, cumsum(result$output$num.mark), data.table::last(cumsum(result$output$num.mark)))
                  TestPeak.pval <- c()
                  if (length(peak.loc) > 2) {
                    for (i in 1:(length(peak.loc) - 1)) {
                      z1 <- psych::fisherz(mean(triplet$corr[(num.mark[peak.loc[i]] + 1):num.mark[peak.loc[i] + 1]], na.rm = TRUE))
                      z2 <- psych::fisherz(mean(triplet$corr[(num.mark[peak.loc[i]] + 1):num.mark[peak.loc[i + 1] + 1]], na.rm = TRUE))
                      N1 <- length(triplet$corr[(num.mark[peak.loc[i]] + 1):num.mark[peak.loc[i] + 1]])
                      N2 <- length(triplet$corr[(num.mark[peak.loc[i]] + 1):num.mark[peak.loc[i + 1] + 1]])
                      TestPeak.pval[i] <- 2 * pnorm(abs(z1 - z2)/sqrt(1/(N1 - 3) + 1/(N2 - 3)), lower.tail = FALSE)
                    }
                  }
                  if (sum(TestPeak.pval > 0.05) != 0) {
                    TestPeak.p <- TestPeak.pval[TestPeak.pval > 0.05]
                    mergp.loc <- which(TestPeak.pval %in% TestPeak.p)
                    distance <- c()
                    if (length(peak.loc) > 2) {
                      for (i in 1:(length(peak.loc) - 1)) {
                        distance[i] <- sum(result$output[(peak.loc[i] + 1):(peak.loc[i + 1] - 1), "num.mark"])
                      }
                    }
                    peak_min <- mergp.loc[distance[mergp.loc] == min(distance[mergp.loc])]
                    p_merg <- intersect(mergp.loc, peak_min)
                    if (length(peak_min) >= 2) {
                      peak_min <- mergp.loc[TestPeak.pval[p_merg] == min(TestPeak.pval[p_merg])]
                    }
                    peak_min <- peak_min[1]
                    result$output[peak.loc[peak_min], "loc.end"] <- result$output[peak.loc[peak_min + 1], "loc.end"]
                    result$output[peak.loc[peak_min], "seg.mean"] <- t(matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min + 1], "num.mark"])) %*% matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min +
                                                                                                                                                                                                            1], "seg.mean"])/sum(result$output[peak.loc[peak_min]:peak.loc[peak_min +
                                                                                                                                                                                                                                                                             1], "num.mark"])
                    result$output[peak.loc[peak_min], "num.mark"] <- sum(result$output[peak.loc[peak_min]:peak.loc[peak_min + 1], "num.mark"])
                    result$output <- result$output[-c((peak.loc[peak_min] + 1):peak.loc[peak_min + 1]), ]
                    row.names(result$output) <- 1:dim(result$output)[1]
                    cand.corr.new <- c(-1, result$output$seg.mean, -1)
                    peak.loc.new <- quantmod::findPeaks(cand.corr.new) - 2
                    no_merg_loc <- c()
                    no_merg_count <- 1
                    if (length(peak.loc.new) >= 2) {
                      for (i in 1:(length(peak.loc.new) - 1)) {
                        if (is.na(sum(result$output[(peak.loc.new[i] + 1):(peak.loc.new[i + 1] - 1), "num.mark"]))) {
                          break
                        }
                        if (sum(result$output[(peak.loc.new[i] + 1):(peak.loc.new[i + 1] - 1), "num.mark"]) > w) {
                          no_merg_loc[no_merg_count] <- peak.loc.new[i]
                        }
                      }
                    }
                    peak.loc.new <- peak.loc.new[-as.numeric(no_merg_loc)]
                    if (length(peak.loc.new) == length(peak.loc))
                      break
                    peak.loc <- peak.loc.new
                  }
                  else break
                }
              }
              num.mark <- c(0, cumsum(result$output$num.mark), data.table::last(cumsum(result$output$num.mark)))
              max_seg <- which(result$output$seg.mean == max(result$output$seg.mean))
              min_seg <- which(result$output$seg.mean == min(result$output$seg.mean))
              if (length(max_seg) != 1) {
                max_seg <- max_seg[2]
              }
              #print(paste0("min_len:", length(max_seg)))
              if (length(min_seg) != 1) {
                min_seg <- min_seg[1]
              }
              z1 <- psych::fisherz(result$output$seg.mean[max_seg])
              z2 <- psych::fisherz(result$output$seg.mean[min_seg])
              N1 <- result$output[max_seg, "num.mark"]
              N2 <- result$output[min_seg, "num.mark"]
              Test <- 2 * pnorm(abs(z1 - z2)/sqrt(1/(N1 - 3) + 1/(N2 - 3)), lower.tail = FALSE)
              if (!is.na(Test) && Test < 0.05) {
                if (sum(cand.corr[peak.loc + 1] > cor_threshold_peak) > 0 && sum(cand.corr[peak.loc + 1] > cor_threshold_peak) <= 2) {
                  cand.ceRNA = paste(r, s)
                  peak.loc = sort(c(peak.loc, no_merg_loc))
                  True_peak <- peak.loc[cand.corr[peak.loc + 1] > cor_threshold_peak]
                  location = result$output[True_peak, c("loc.start", "loc.end")]
                  if (!is.null(cand.ceRNA)) {
                    tmp <- list(miRNA = mir, cand.ceRNA = cand.ceRNA, location = location, numOfseg = result$output$num.mark[True_peak])
                  }
                }
              }
            }
          }

        }

        parallel::stopCluster(cl)
        tmp
      }

      testfunction <- purrr::map(1:length(mirna_total), sigCernaPeak,readRDS(paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_pairfiltering.rds')),cor_threshold_peak,window_size)
      FinalResult <- purrr::compact(testfunction)

      if (dir.exists(paste0(path_prefix, project_name, '-', disease_name,'/03_identifiedPairs')) == FALSE){
        dir.create(paste0(path_prefix, project_name, '-', disease_name,'/03_identifiedPairs'))
      }
      saveRDS(FinalResult,paste0(path_prefix, project_name,'-',disease_name,'/03_identifiedPairs/',project_name,'-',disease_name,'_finalpairs.rds'))

      final_df <- as.data.frame(Reduce(rbind, FinalResult))
      flat_df <-  final_df %>%
        tidyr::unnest(location) %>%
        tidyr::unnest(numOfseg)
      flat_df <- as.data.frame(flat_df)
      data.table::fwrite(flat_df, paste0(path_prefix, project_name,'-', disease_name,'/',project_name,'-', disease_name, '_finalpairs.csv'), row.names = FALSE)

      time2 <- Sys.time()
      diftime <- difftime(time2, time1, units = 'min')

      message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
      message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')
      flat_df
    }

    final_results <- SegmentClusteringPlusPeakMerging(path_prefix = path_prefix,
                                                      project_name = project_name,
                                                      disease_name = disease_name,
                                                      cor_threshold_peak = cor_threshold_peak,
                                                      window_size = window_size)


    as.data.frame(final_results)

  }

  final_results <- ceRNAMethod(path_prefix = path_prefix,
                               project_name = project_name,
                               disease_name = disease_name,
                               window_size = window_size,
                               cor_method = cor_method,
                               cor_threshold_peak = 0.85)

  as.data.frame(final_results)
}
