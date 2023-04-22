#'
#' @name ceRNATCGA
#' @title Retrieval of public TCGA data from GDC Xena Hub
#' @description A function to retrieve TCGA data
#'
#' @import utils
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default = 'TCGA')
#' @param disease_name the abbreviation of disease that users are interested in
#' (default = 'DLBC')
#'
#' @returns file
#' @export
#'
#' @examples
#' ceRNATCGA(
#' path_prefix = NULL,
#' project_name = 'TCGA',
#' disease_name = 'DLBC'
#' )
#'

ceRNATCGA <- function(path_prefix = NULL,
                      project_name = 'TCGA',
                      disease_name = 'DLBC'){

  if (is.null(path_prefix)){
    path_prefix <- tempdir()
  }else{
    path_prefix <- path_prefix
  }

  if (!stringr::str_detect(path_prefix, '/$')){
    path_prefix <- paste0(path_prefix, '/')
  }

  if (dir.exists(paste0(path_prefix, project_name,'-', disease_name)) == FALSE){
    dir.create(paste0(path_prefix, project_name,'-', disease_name))
  }

  if (dir.exists(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata')) == FALSE){
    dir.create(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'))
  }

  time1 <- Sys.time()

  (TCGA <- curatedTCGAData::curatedTCGAData(
    disease_name, c("RNASeq*", "miRNA*"), version = "2.0.1", dry.run = FALSE
  ))
  exportClass(TCGA, dir = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'), fmt = "csv", ext = ".csv")
  file.rename(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'_', disease_name,'_miRNASeqGene-20160128.csv'), paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mirna.csv'))
  file.rename(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'_', disease_name,'_RNASeq2Gene-20160128.csv'), paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mrna.csv'))

  file.rename(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'_colData.csv'), paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_phenotype.csv'))

  junk <- dir(path=paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'),  pattern="TCGA_") # ?dir
  file.remove(paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/',junk))

  temp <-  list.files(path = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'), pattern="*.csv")
  temp_name <- gsub(paste0(project_name,'-', disease_name,"_"),"",gsub(".csv", "", temp))

  for (i in 1:length(temp)) assign(temp_name[i], data.frame(data.table::fread(paste0(path_prefix,project_name,'-', disease_name, '/01_rawdata/', temp[i]),header = TRUE, fill = TRUE), row.names = 1))

  # extract survival information from
  survival <- phenotype[,c("vital_status", "days_to_death", "days_to_last_followup")]
  for (i in 1:dim(survival)[1]){
    #i =1
    if (is.na(survival$days_to_death[i])){
      survival$days_to_death[i] <- survival$days_to_last_followup[i]
    }
  }
  survival <- survival[,-3]
  names(survival) <- c('OS', 'OS.time')

  # id in names()
  htseq_fpkm <- mrna
  names(htseq_fpkm) <- gsub("\\.","-", names(htseq_fpkm))
  htseq_fpkm <- htseq_fpkm[,substring(names(htseq_fpkm),16,17)=='A-']
  htseq_fpkm <- htseq_fpkm[,!(substring(names(htseq_fpkm),14,15) %in% seq(10,29,1))]
  htseq_fpkm <- htseq_fpkm[,!duplicated(substring(names(htseq_fpkm),1,12))]
  names(htseq_fpkm) <-substring(names(htseq_fpkm),1,12)
  names(mirna) <- gsub("\\.","-", names(mirna))
  mirna <- mirna[,substring(names(mirna),16,17)=='A-']
  mirna <- mirna[,!(substring(names(mirna),14,15) %in% seq(10,29,1))]
  mirna <- mirna[,!duplicated(substring(names(mirna),1,12))]
  names(mirna) <-substring(names(mirna),1,12)
  # get union sampleID
  g1 <- list(names(mirna), names(htseq_fpkm), row.names(survival))

  union_sampleID <- sort(Reduce(intersect, g1))
  message('\u2605 TCGA-', disease_name, ' cohort contains ',length(union_sampleID), ' tumor samples!')

  # id in row.names()
  clinicDataSort <- function(df){
    df <- df[row.names(df) %in% union_sampleID,]
    df$rowname <- row.names(df)
    df <-  df[sort(df$rowname),]
    df[ , -which(names(df) %in% c("rowname"))]
  }

  survival <- clinicDataSort(survival)
  # id in names()
  htseq_fpkm <- htseq_fpkm[,names(htseq_fpkm) %in% union_sampleID]
  htseq_fpkm <- htseq_fpkm[ , order(names(htseq_fpkm))]
  htseq_fpkm <- as.data.frame(scale(log(htseq_fpkm+1,2), center = TRUE, scale = TRUE))
  mirna <- mirna[,names(mirna) %in% union_sampleID]
  mirna <- mirna[ , order(names(mirna))]
  mirna <- as.data.frame(scale(log(mirna+1,2), center = TRUE, scale = TRUE))

  # mRNA:log2(fpkm+1)，miRNA:log2(RPM+1)
  # mirna <- (2^mirna-1)*1000  #RPKM: X1000
  # htseq_fpkm <- 2^htseq_fpkm -1

  # focus on protein coding RNA because downloaded mRNA expression matrix includes both codingRNA and lncRNA
  ensem2symbol <- get0("ensem2symbol", envir = asNamespace("ceRNAR"))
  rownames(ensem2symbol) <- ensem2symbol$gene_id
  cdRNA <- htseq_fpkm[rownames(htseq_fpkm) %in% ensem2symbol$gene_name[ensem2symbol$gene_type=="protein_coding"],]

  #miRNA id conversion
  #ID_converter <- ceRNAR:::hsa_pre_mature_matching
  ID_converter <- get0("hsa_pre_mature_matching", envir = asNamespace("ceRNAR"))
  mirna$ID <- row.names(mirna)
  miRNA_with_precurer <- merge(ID_converter, mirna, by='ID')[,-1]
  miRNA_with_precurer <- stats::aggregate(. ~ Mature_ID, data = miRNA_with_precurer, mean)
  row.names(miRNA_with_precurer) <- miRNA_with_precurer[,1]
  miRNA_with_precurer <- miRNA_with_precurer[,-1]

  # store processed data
  data.table::fwrite(as.data.frame(cdRNA),paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mrna.csv'), row.names = TRUE)
  data.table::fwrite(as.data.frame(miRNA_with_precurer),paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mirna.csv'), row.names = TRUE)
  data.table::fwrite(as.data.frame(survival), paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_survival.csv'), row.names = TRUE)
  message('(\u2714) All files have been preprocessed!')
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')
}
