#'
#' @name SegmentClusteringPlusPeakMerging
#' @title two of three steps in main ceRNAR algorithm
#' @description A function to conduct two of three steps in algorithm, that is, segment clustering and peak merging
#'
#' @import foreach
#' @import parallel
#' @import utils
#' @rawNamespace import(magrittr, except = set_names)
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default: demo)
#' @param disease_name the abbreviation of disease that users are interested in
#' (default: DLBC)
#' @param cor_threshold_peak peak threshold of correlation value between 0 and 1
#' (default: 0.85)
#' @param window_size the number of samples for each window (default:10)
#'
#' @returns a tabular output
#' @export
#'
#' @examples
#' SegmentClusteringPlusPeakMerging(
#' path_prefix =  '~/',
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' cor_threshold_peak = 0.85,
#' window_size = 10
#' )
#'

SegmentClusteringPlusPeakMerging <- function(path_prefix,
                                             project_name = 'demo',
                                             disease_name = 'DLBC',
                                             cor_threshold_peak = 0.85,
                                             window_size = 10){

  time1 <- Sys.time()

  message('\u25CF Step4: Clustering segments using CBS algorithm plus Mearging peaks')

  # setwd(paste0(project_name,'-',disease_name))
  dict <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_MirnaTarget_dictionary.rds'))
  mirna <- data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mirna.csv')),row.names = 1)
  mrna <- data.frame(data.table::fread(paste0(path_prefix, project_name,'-',disease_name,'/01_rawdata/',project_name,'-',disease_name,'_mrna.csv')),row.names = 1)
  mirna_total <- unlist(dict[,1])
  d <- readRDS(paste0(path_prefix, project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_pairfiltering.rds'))

  ## create a cluster
  message('\u2605 Number of computational cores: ',parallel::detectCores()-3,'/',parallel::detectCores(), '.')
  sigCernaPeak <- function(index,d, cor_threshold_peak, window_size){
    #index=1
    #print(paste0('microRNA:', index))
    w <- window_size
    mir = mirna_total[index]
    gene <- as.character(data.frame(dict[dict[,1]==mir,][[2]])[,1])
    gene <- intersect(gene,rownames(mrna))

    gene_pair <- combn(gene,2)
    total_pairs <- choose(length(gene),2)
    #tmp <- NULL
    doParallel::registerDoParallel(parallel::detectCores()-3)
    #tmp <- tryCatch({
    tmp <- foreach(p=1:total_pairs, .combine = "rbind")  %dopar%  {
      #lst <- list()
      #for (p in 1:total_pairs){ # test foreach
      #p=1
      print(paste0("no_of_index:", index,"|", "no_of_pairs:",p))
      cand.ceRNA=c()
      location=list()
      r=gene_pair[1,p]
      s=gene_pair[2,p]
      triplet <- d[[index]][,c(1,p+1)]
      names(triplet) <- c("miRNA","corr")

      if(sum(is.na(triplet$corr)) ==0){
        # 01. SegmentClustering method: using CBS ("DNAcopy")
        SegmentClustering <- function(triplet){
          CNA.object <- DNAcopy::CNA(triplet$corr,rep(1,dim(triplet)[1]),triplet$miRNA)
          names(CNA.object) <- c("chrom","maploc",paste("gene",r,"and",s))  ### gene name
          result <- DNAcopy::segment(CNA.object)
          result
        }
        result <- SegmentClustering(triplet)

        # 02. peakMerging method
        # 02-1. merge too short segment
        if(sum(result$output$num.mark<=3)>=1){ ### 2
          tooshort <- which(result$output$num.mark<=3)
          num.mark <- c(0,cumsum(result$output$num.mark),data.table::last(cumsum(result$output$num.mark)))
          ### merge short neighbor segment first
          if(1 %in% diff(tooshort)){
            cc=1
            lag=0
            for(q in 1:(length(tooshort)-1)){
              if(tooshort[q+1]-tooshort[q]==1){
                result$output[tooshort[q],"loc.end"] <- result$output[tooshort[q+1],"loc.end"]
                result$output[tooshort[q],"seg.mean"] <- t(matrix(result$output[tooshort[q]:tooshort[q+1],"num.mark"]))%*%matrix(result$output[tooshort[q]:tooshort[q+1],"seg.mean"])/sum(result$output[tooshort[q]:tooshort[q+1],"num.mark"])
                result$output[tooshort[q],"num.mark"] <- sum(result$output[tooshort[q]:tooshort[q+1],"num.mark"])
                result$output[tooshort[q+1],] <- result$output[tooshort[q],]
                lag[cc] <- tooshort[q]
                cc <- cc+1
              }
            }
            result$output <- result$output[-lag,]
            row.names(result$output) <- 1:dim(result$output)[1]
            tooshort <- which(result$output$num.mark<=3)
          }
          ### merge short segment to long neighbor segment
          if(length(tooshort)>=1){
            cc=1
            lag=c()
            for(t in seq_along(tooshort)){
              long_seg <- which(result$output$num.mark>3)
              diff=abs(tooshort[t]-long_seg)
              closest_seg <- long_seg[which(diff==min(diff))]
              if(length(closest_seg)>=2){
                b <- abs(result$output$seg.mean[closest_seg]-result$output$seg.mean[tooshort[t]])
                closest_seg <- closest_seg[b==min(b)][1]
              }
              result$output[tooshort[t],"loc.start"] <- min(result$output[tooshort[t],"loc.start"],result$output[closest_seg,"loc.start"])
              result$output[tooshort[t],"loc.end"] <- max(result$output[tooshort[t],"loc.end"],result$output[closest_seg,"loc.end"])
              result$output[tooshort[t],"seg.mean"] <- t(matrix(result$output[tooshort[t]:closest_seg,"num.mark"]))%*%matrix(result$output[tooshort[t]:closest_seg,"seg.mean"])/sum(result$output[tooshort[t]:closest_seg,"num.mark"])
              result$output[tooshort[t],"num.mark"] <- sum(result$output[tooshort[t]:closest_seg,"num.mark"])
              result$output[closest_seg,] <- result$output[tooshort[t],]
              lag[cc] <- tooshort[t]
              cc <- cc+1
            }
            result$output <- result$output[-lag,]
            row.names(result$output) <- 1:dim(result$output)[1]
          }

        }

        cand.corr <- c(-1,result$output$seg.mean,-1)  #add 0 to detect peaks happened at head and tail
        peak.loc <- quantmod::findPeaks(cand.corr)-2

        # 02-2. While loop for merging peak
        no_merg_loc <- c()
        no_merg_count <- 1
        if(sum(cand.corr[peak.loc+1] > cor_threshold_peak) >=2){ ### para 0.5
          for(i in 1:(length(peak.loc)-1)){
            if (is.na(sum(result$output[(peak.loc[i]+1):(peak.loc[i+1]-1),"num.mark"]))) {break}
            if(sum(result$output[(peak.loc[i]+1):(peak.loc[i+1]-1),"num.mark"]) > w){
              no_merg_loc[no_merg_count] <- peak.loc[i]
            }
          }
          #tryCatch({
          peak.loc <- peak.loc[-which(peak.loc==no_merg_loc)]
          #},error=function(e){})
          while(sum(cand.corr[peak.loc+1] > cor_threshold_peak) >=2){  ### para 0.5
            num.mark <- c(0,cumsum(result$output$num.mark),data.table::last(cumsum(result$output$num.mark)))
            TestPeak.pval <- c()
            if (length(peak.loc)>2){
              for(i in 1:(length(peak.loc)-1)){
                #TestPeak.pval[i] <- t.test(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]],triplet$corr[(num.mark[peak.loc[i+1]]+1):num.mark[peak.loc[i+1]+1]])$p.value
                z1 <- psych::fisherz(mean(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]],na.rm=T))
                z2 <- psych::fisherz(mean(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i+1]+1]],na.rm=T))
                N1 <- length(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]])
                N2 <- length(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i+1]+1]])
                TestPeak.pval[i] <- 2*pnorm(abs(z1-z2)/sqrt(1/(N1-3)+1/(N2-3)),lower.tail = FALSE)
                #TestPeak.pval[i] <- t.test(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]],triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i+1]+1]])$p.value

              }
            }
            if(sum(TestPeak.pval>0.05)!=0){  ### para:alpha 0.05
              TestPeak.p <- TestPeak.pval[TestPeak.pval>0.05]
              mergp.loc <- which(TestPeak.pval%in%TestPeak.p)
              #peak_min <- which(TestPeak.p==max(TestPeak.p))
              distance <- c()
              if (length(peak.loc)>2){
                for(i in 1:(length(peak.loc)-1)){
                  distance[i] <- sum(result$output[(peak.loc[i]+1):(peak.loc[i+1]-1),"num.mark"])
                }
              }
              peak_min <- mergp.loc[distance[mergp.loc]==min(distance[mergp.loc])]
              p_merg <- intersect(mergp.loc,peak_min)
              if(length(peak_min)>=2){
                peak_min <- mergp.loc[TestPeak.pval[p_merg]==min(TestPeak.pval[p_merg])]
              }
              peak_min <- peak_min[1]

              result$output[peak.loc[peak_min],"loc.end"] <- result$output[peak.loc[peak_min+1],"loc.end"]
              result$output[peak.loc[peak_min],"seg.mean"] <- t(matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"]))%*%matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"seg.mean"])/sum(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"])
              result$output[peak.loc[peak_min],"num.mark"] <- sum(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"])
              result$output <- result$output[-c((peak.loc[peak_min]+1):peak.loc[peak_min+1]),]
              row.names(result$output) <- 1:dim(result$output)[1]

              cand.corr.new <- c(-1,result$output$seg.mean,-1)  #add 0 to detect peaks happened at head and tail
              peak.loc.new <- quantmod::findPeaks(cand.corr.new)-2

              #tryCatch({
              no_merg_loc <- c()
              no_merg_count <- 1
              if(length(peak.loc.new)>=2){
                for(i in 1:(length(peak.loc.new)-1)){
                  if (is.na(sum(result$output[(peak.loc.new[i]+1):(peak.loc.new[i+1]-1),"num.mark"]))) {break}
                  if(sum(result$output[(peak.loc.new[i]+1):(peak.loc.new[i+1]-1),"num.mark"])> w){
                    no_merg_loc[no_merg_count] <- peak.loc.new[i]
                  }
                }
              }
              peak.loc.new <- peak.loc.new[-as.numeric(no_merg_loc)]
              #},error=function(e){})

              if(length(peak.loc.new)==length(peak.loc)) break
              peak.loc <- peak.loc.new


            }else break
          }

        }

        # 03. test significance of the highest peak vs the lowest
        num.mark <- c(0,cumsum(result$output$num.mark),data.table::last(cumsum(result$output$num.mark)))
        max_seg <- which(result$output$seg.mean==max(result$output$seg.mean))
        min_seg <- which(result$output$seg.mean==min(result$output$seg.mean))

        if (length(max_seg)!=1){
          max_seg <- max_seg[2]
        }
        print(paste0('min_len:',length(max_seg)))
        if (length(min_seg)!=1){
          min_seg <- min_seg[1]
        }

        z1 <- psych::fisherz(result$output$seg.mean[max_seg])
        z2 <- psych::fisherz(result$output$seg.mean[min_seg])
        N1 <- result$output[max_seg,"num.mark"]
        N2 <- result$output[min_seg,"num.mark"]
        Test <- 2*pnorm(abs(z1-z2)/sqrt(1/(N1-3)+1/(N2-3)),lower.tail = FALSE)
        # generate final output
        if(!is.na(Test) && Test < 0.05){
          if(sum(cand.corr[peak.loc+1] > cor_threshold_peak) >0 && sum(cand.corr[peak.loc+1] > cor_threshold_peak) <=2){  ### para 0.5
            cand.ceRNA=paste(r,s)

            #tryCatch({
            peak.loc=sort(c(peak.loc,no_merg_loc)) #put back the peak that can't be merged
            #},error=function(e){})

            True_peak <- peak.loc[cand.corr[peak.loc+1] > cor_threshold_peak]
            location=result$output[True_peak,c("loc.start","loc.end")]

            if(!is.null(cand.ceRNA)){
              lst <- list(miRNA=mir,cand.ceRNA=cand.ceRNA,location=location,numOfseg=result$output$num.mark[True_peak])
              lst
            }

          }
        }
      }
    }
    #}
    #},error=function(e){e})
    #tmp <- do.call(rbind,lst)
    tmp
  }
  # seed=NULL
  # txt=list()
  # for (i in 1:length(mirna_total)){
  #   print(i)
  #   txt[[i]] <- sigCernaPeak(i,d, cor_threshold_peak, window_size)
  # }
  # txt_final <- do.call(rbind, txt)

  testfunction <- purrr::map(1:length(mirna_total), sigCernaPeak,readRDS(paste0(path_prefix,project_name,'-',disease_name,'/02_potentialPairs/',project_name,'-',disease_name,'_pairfiltering.rds')),cor_threshold_peak,window_size)
  FinalResult <- purrr::compact(testfunction)

  if (dir.exists(paste0(path_prefix, project_name, '-', disease_name,'/03_identifiedPairs')) == FALSE){
    dir.create(paste0(path_prefix, project_name, '-', disease_name,'/03_identifiedPairs'))
  }
  saveRDS(FinalResult,paste0(path_prefix, project_name,'-',disease_name,'/03_identifiedPairs/',project_name,'-',disease_name,'_finalpairs.rds'))

  final_df <- as.data.frame(Reduce(rbind, FinalResult))
  flat_df <-  final_df %>%
    tidyr::unnest(location) %>%
    tidyr::unnest(numOfseg)
  data.table::fwrite(flat_df, paste0(path_prefix, project_name,'-', disease_name,'/',project_name,'-', disease_name, '_finalpairs.csv'), row.names = FALSE)

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')

  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' min.'))
  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')
  flat_df
}

