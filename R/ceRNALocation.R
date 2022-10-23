#'
#' @name ceRNALocation
#' @title Visualization for peak location
#' @description A function to visualize the peak location at certain miRNA level
#'
#' @import utils
#' @import grDevices
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param mirna a specific mirna name (such as hsa-miR-101-3p) assigned by user
#' @param window_size the number of samples for each window and usually about
#' one third of total samples
#'
#' @export
#'
#' @examples
#' ceRNALocation(
#' path_prefix = '~/',
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' mirna='hsa-miR-101-3p',
#' window_size = 45/5
#' )
#'

ceRNALocation <- function(path_prefix,
                          project_name,
                          disease_name,
                          mirna,
                          window_size){

  if (!stringr::str_detect(path_prefix, '/')){
    path_prefix <- paste0(path_prefix, '/')
  }

  if(!dir.exists(paste0(path_prefix, '/', project_name,'-',disease_name,'/04_downstreamAnalyses/'))){
    dir.create(paste0(path_prefix, '/', project_name,'-',disease_name,'/04_downstreamAnalyses/'))
  }

  if(!dir.exists(paste0(path_prefix, '/', project_name,'-',disease_name,'/04_downstreamAnalyses/peakLocationResults/'))){
    dir.create(paste0(path_prefix, '/', project_name,'-',disease_name,'/04_downstreamAnalyses/peakLocationResults/'))
  }

  message('\u25CF Step4: Dowstream Analyses - Peak Location analysis for ', mirna, '.')
  Res <- readRDS(paste0(path_prefix, '/', project_name,'-',disease_name,'/03_identifiedPairs/', project_name, '-', disease_name,'_finalpairs.rds'))

  # get all ceRNAs
  Res_dataframe <- as.data.frame(Reduce(rbind,purrr::compact(Res)))
  a=Res_dataframe[Res_dataframe$miRNA==mirna ,2:4]

  # get mirna expression data
  mirExp <- as.data.frame(data.table::fread(paste0(path_prefix,'/',project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mirna.csv'),header = T,stringsAsFactors = F))
  row.names(mirExp) <- mirExp[,1]
  mirExp <- mirExp[,-1]
  names(mirExp) <- substring(names(mirExp),1,12)

  # get location
  a_order <- a[order(unlist(lapply(a[,2], function(x) dplyr::last(x[,1]))), decreasing = TRUE),]
  start_ls <- lapply(a_order[,2], function(x) stats::na.omit(x[,1]))
  end_ls <- lapply(a_order[,2], function(x) stats::na.omit(x[,2]))

  start <- unlist(start_ls)
  end <- unlist(end_ls)
  location <- end_ls

  ystart=rep(0:(dim(a_order)[1]-1),unlist(lapply(location, length)))
  yend=ystart+1

  d=data.frame(x1=start, x2=end, y1=ystart, y2=yend,
               ylab=unlist(rep(a_order[,1],unlist(lapply(location, length)))),
               yloc=rep(seq(0.5,dim(a_order)[1],by=1),unlist(lapply(location, length))),
               Count=as.matrix(stats::na.omit(unlist(a_order[,3]))))
  N=dim(mirExp)[2]
  w=window_size
  tmp <- as.data.frame(t(mirExp[row.names(mirExp)==mirna,]))
  tmp$sample <- row.names(tmp)
  mirExp_subset <- tmp[order(tmp[,1], decreasing = T),]

  #x_range <- c(mean(mirExp_subset[1:w,1]),mean(mirExp_subset[N:(N-w+1),1]))
  # to customize based on the total pairs
  p <- ggplot2::ggplot(data=d,ggplot2::aes(fill=Count/N*100))+
    ggplot2::geom_rect(mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
              colour="steelblue", alpha=0.9)+
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),axis.ticks.y = ggplot2::element_blank())+ # Remove axis ticks and tick mark labels
    ggplot2::scale_fill_gradient2(low="black", high="steelblue", name="Sample (%)")
  p1 <- p + ggplot2::scale_y_continuous(breaks = d$yloc, labels = d$ylab)+
    ggplot2::theme_bw()+
    ggplot2::xlab("miRNA expression (in terms of RPM)")+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.text.y = ggplot2::element_text(size=6),
                   axis.title=ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())+
    ggplot2::ggtitle(mirna)

  p2 <- ggplot2::qplot(mirExp_subset[,1],geom = 'histogram', fill=I("steelblue"),col=I("white"),binwidth=0.1,xlab='miRNA expression (in terms of RPM)')+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))+
    ggplot2::ylab('Count')

  grDevices::png(paste0(path_prefix, '/', project_name,'-',disease_name,'/04_downstreamAnalyses/peakLocationResults/',project_name,'-',disease_name, '_', mirna, '_peakLocation.png'), width = 2000, height = 3500, res = 300)
  merged_p <- egg::ggarrange(p1,p2,
            ncol=2, nrow=2, widths=c(4, 1), heights=c(4, 0.5))
  merged_p
  grDevices::dev.off()

  return(merged_p)

  message('\u2605\u2605\u2605 Peak location analysis has completed! \u2605\u2605\u2605')

}


