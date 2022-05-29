#'
#' @name ceRNAModule
#' @title Network analysis and visualization
#' @description A function to analyze and visualize potential network of
#' identified ceRNAs
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param window_size the number of samples for each window and usually about
#' one third of total samples
#' @param pairs_cutoff at least the number of ceRNA pairs that a mirna must have
#' @param column_sum the number of ceRNAs
#'
#' @examples
#' ceRNAModule(
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' pairs_cutoff = 5,
#' column_sum = 1
#' )
#'
#' @export


ceRNAModule <- function(path_prefix=NULL,
                        project_name,
                        disease_name,
                        pairs_cutoff,
                        column_sum){
  if (is.null(path_prefix)){
    path_prefix <- getwd()
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }else{
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }
  time1 <- Sys.time()
  # setwd(paste0(project_name,'-',disease_name))
  if(!dir.exists(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/'))){
    dir.create(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/'))
  }

  if(!dir.exists(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults'))){
    dir.create(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults'))
  }
  message('\u25CF Step4: Dowstream Analyses - Network analysis')
  Dict <- readRDS(paste0(project_name,'-',disease_name,'/02_potentialPairs/', project_name, '-', disease_name, '_MirnaTarget_dictionary.rds'))
  Res <- readRDS(paste0(project_name,'-',disease_name,'/03_identifiedPairs/', project_name, '-', disease_name, '_finalpairs.rds'))
  # total number of putative genes (targeted gene not data gene)
  all_putative_p <- unlist(lapply(Dict[,2], function(x) unlist(c(x))))
  message('\u2605 Number of targeted genes: ', length(unique(all_putative_p)), '.')
  message('\u2605 Number of putative triplets: ', sum(unlist(lapply(Dict[,2], function(x) choose(length(unlist(c(x))),2) ))), '.')

  # get numbers of ceRNA pair in each miRNA
  ceRNA_pair_count <- unlist(lapply(Res,function(x) dim(x)[1]))

  # get those miRNA names
  rank_count <- sort(ceRNA_pair_count,decreasing = TRUE) # get all
  rank_count <- unique(rank_count)
  summary_miR=c()
  for(i in 1:length(rank_count)){
    # i=1
    rank_miRNA <- which(lapply(Res,function(x) dim(x)[1]==rank_count[i])==TRUE)   # max ceRNA pairs
    miR_name <- unlist(lapply(Res[rank_miRNA], function(x) unique(x[,1])))
    summary_miR <- rbind(summary_miR,c(list(miR_name),rank_count[i]))
  }

  # save miR summary table
  write.csv(summary_miR, paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults/miR_summary_', project_name, '-', disease_name, '.csv'), row.names = F)

  # Average ceRNA pairs of overall miRNAs
  message('\u2605 Average ceRNA pairs of overall miRNAs: ', sum(unlist(summary_miR[,2]))/dim(summary_miR)[1], '.')

  # get hubgene in each miRNA
  rank_count <- unique(rank_count)
  hubgene_name=list()
  for(i in 1:length(rank_count)){
    rank_miRNA <- which(lapply(Res,function(x) dim(x)[1]==rank_count[i])==TRUE)  # max ceRNA pairs
    tmp=c()
    for(j in rank_miRNA){
      a <- names(sort(table(unlist(strsplit(unlist(Res[[j]][,2]), " "))),decreasing=TRUE)[1])
      count <- sort(table(unlist(strsplit(unlist(Res[[j]][,2]), " "))),decreasing=TRUE)[1]

      tmp=rbind(tmp,c(a,count))
    }
    hubgene_name[[i]] <- tmp
  }

  write.csv(matrix(unlist(hubgene_name),ncol=2,byrow=TRUE),paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults/hubgene_miR_', project_name, '-', disease_name, '.csv'), row.names = F)
  hubg_onlyname <- unlist(sapply(hubgene_name, function(x) x[,1]))

  # get all ceRNAs (Res list => Res dataframe)
  Res <- as.data.frame(Reduce(rbind,purrr::compact(Res)))

  hubgene_all <- sort(table(unlist(strsplit(as.character(Res[,2])," "))), decreasing = TRUE)
  write.csv(hubgene_all, paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults/hubgene_all_', project_name, '-', disease_name, '.csv'), row.names = F)

  # average number of coexpressed gene
  message('\u2605 Average number of coexpressed gene: ', sum(hubgene_all)/length(hubgene_all), '.')
  message('\u2605 Number of genes being found: ', length(hubgene_all), '.')

  # general network
  mir_unique <- unique(Res[,1])
  network_plot <- function(a, column_sum){
    node <- unique(unlist(strsplit(unlist(a[,1]), " ")))
    m <- matrix(0, nrow = length(node), ncol = length(node))
    colnames(m) <- node
    rownames(m) <- node
    pair_list <- strsplit(unlist(a[,1]), " ")
    for(p in 1:length(pair_list)){
      m[pair_list[[p]][1],pair_list[[p]][2]] <- 1
      m[pair_list[[p]][2],pair_list[[p]][1]] <- 1

    }
    net <- network::network(m, directed = TRUE)

    g <- 1:length(node)

    ggnet::ggnet(m, label = ifelse(colSums(m)>=column_sum, colnames(m), NA), alpha = 1, color="black",
          node.group = g, node.color =rep('orange',length(node)),  segment.color = "grey50",
          legend.position="none",  weight = colSums(m) ,mode="circle") # columnsum to tune
  }

  for (i in mir_unique){
    mir_df <- Res[Res[,1]== i,]
    if (dim(mir_df)[1] <= pairs_cutoff || is.null(dim(mir_df)[1])){
      #print(paste0(i, ' is skipped.'))
      next
    }else{
      # general network
      #print(paste0(i, ' is processed.'))
      network_plot(mir_df[,-1], column_sum)
      ggplot2::ggsave(paste0(project_name,'-',disease_name,'/04_downstreamAnalyses/moduleResults/',i,'_ceRNAs_network.png'), height = 10, width = 10, dpi = 300)
    }
  }
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
  message('\u2605\u2605\u2605 Network analysis has completed! \u2605\u2605\u2605')
}


