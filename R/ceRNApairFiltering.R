#'
#' @name ceRNApairFiltering
#' @title one of three steps in main ceRNAR algorithm
#' @description A function to conduct one of three steps in algorithm,
#' that is, pairs filtering
#'
#' @import foreach
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param window_size the number of samples for each window
#' @param cor_method selection of correlation methods, including pearson and
#' spearman (default: pearson)
#'
#' @examples
#' ceRNApairFilering(
#' project_name = 'demo',
#' disease_name= 'DLBC',
#' window_size = 45/5
#' )
#'
#' @export


ceRNApairFilering <- function(path_prefix = NULL,
                              project_name,
                              disease_name,
                              window_size,
                              cor_method='pearson'){
  if (is.null(path_prefix)){
    path_prefix <- getwd()
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }else{
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }

  time1 <- Sys.time()
  message('\u25CF Step3: Filtering putative mRNA-miRNA pairs using sliding window approach')

  # create a cluster
  n <- parallel::detectCores()
  message('\u2605 Number of cores: ', n-2, '/', n, '.')
  cl <- parallel::makeCluster(n-3, outfile="")
  doSNOW::registerDoSNOW(cl)

  # setwd(paste0(project_name,'-',disease_name))
  # import example data & putative pairs
  dict <- readRDS(paste0(project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_MirnaTarget_dictionary.rds'))
  mirna <- data.frame(data.table::fread(paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mirna.csv')),row.names = 1)
  mrna <- data.frame(data.table::fread(paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mrna.csv')),row.names = 1)
  mirna_total <- unlist(dict[,1])
  message(paste0('\u2605 total miRNA: ', length(mirna_total)))
  slidingWindow <- function(window_size, mirna_total, cor_method){
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
      gene_pair <- combn(gene,2)
      gene_pair_index <- rbind(match(gene_pair[1,], names(gene_mir)),match(gene_pair[2,], names(gene_mir)))
      total_pairs <- choose(length(gene),2)
      message(paste0("\u2605 ",mir,"'s total potential ceRNA-miRNA pairs: ",total_pairs))

      getcorr <- function(r,s){
        #r=2
        #s=3
        y <- gene_mir[,c(1,r,s)]
        y <- y[order(y$miRNA),]
        corr <- zoo::rollapply(y, width=w, function(x) cor(x[,2],x[,3],method = cor_method), by.column=FALSE)
        miRNA <- zoo::rollapply(y, width=w, function(x) mean(x[,1]), by.column=FALSE)
        data <- data.frame(cbind(miRNA,corr))
        data <- data[order(data$miRNA),]
        data$corr
      }
      cor_all <- purrr::map2_dfc(gene_pair_index[1,1:total_pairs],gene_pair_index[2,1:total_pairs],getcorr)
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
    parallel_d
  }

  Realdata <- slidingWindow(window_size,mirna_total, 'pearson')
  saveRDS(Realdata,paste0(project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_pairfiltering.rds'))

  # close a cluster
  #closeAllConnections()
  CatchupPause <- function(Secs){
    Sys.sleep(Secs) #pause to let connection work
    future:::ClusterRegistry("stop")
  }
  CatchupPause(3)
  #parallel::stopCluster(cl)

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))

  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')

}

