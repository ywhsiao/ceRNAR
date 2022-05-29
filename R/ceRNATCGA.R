#'
#' @name ceRNATCGA
#' @title Retrieval of public TCGA data from GDC Xena Hub
#' @description A function to retrieve TCGA data from GDC Xena Hub
#' (https://xenabrowser.net/datapages/)
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param timout the allowance time for downloading TCGA data
#'
#' @examples
#' ceRNATCGA(
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' )
#'
#' @export

ceRNATCGA <- function(path_prefix = NULL,
                      project_name = 'TCGA',
                      disease_name,
                      timeout = 5000000){

  if (is.null(path_prefix)){
    path_prefix <- getwd()
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }else{
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }

  if (dir.exists(paste0(project_name,'-',disease_name)) == FALSE){
    dir.create(paste0(project_name,'-',disease_name))
  }
  setwd(paste0(project_name,'-',disease_name))

  if (dir.exists('01_rawdata') == FALSE){
    dir.create('01_rawdata')
  }
  setwd('01_rawdata')

  time1 <- Sys.time()

  # download cancer files (phenotype, survival, miRNA and mRNA) from gdc resource
  downloadFromGDC <- function(project,cancer,timeout=5000000){
    message('\u25CF Step1 for TCGA data: Downloading the data ...')
    options(timeout=timeout) # to test the largest data
    HelpersMG::wget(paste0('https://gdc.xenahubs.net/download/',project_name,'-', disease_name,'.GDC_phenotype.tsv.gz'))
    HelpersMG::wget(paste0('https://gdc.xenahubs.net/download/',project_name,'-', disease_name,'.survival.tsv'))
    HelpersMG::wget(paste0('https://gdc.xenahubs.net/download/',project_name,'-', disease_name,'.mirna.tsv.gz'))
    HelpersMG::wget(paste0('https://gdc.xenahubs.net/download/',project_name,'-', disease_name,'.htseq_fpkm.tsv.gz'))
  }
  downloadFromGDC(project,cancer)

  # unzip
  temp <-  list.files(pattern="*.gz")
  head(temp)
  for (i in 1:length(temp)) R.utils::gunzip(temp[i], remove=TRUE)
  message('(\u2714) All files have been and downloaded and uncompressed!')

  # match sample id
  temp <-  list.files(pattern="*.tsv")
  temp_name <- gsub(".*\\.","",gsub(".tsv", "", temp))
  for (i in 1:length(temp)) assign(temp_name[i], data.frame(data.table::fread(temp[i], sep = '\t',header = TRUE, dec = ".", fill = TRUE), row.names = 1))
  # id in a column
  GDC_phenotype <- as.data.frame(dplyr::filter(GDC_phenotype, grepl("01A", row.names(GDC_phenotype))))
  row.names(GDC_phenotype) <- gsub('-01A','',row.names(GDC_phenotype))
  survival <- as.data.frame(dplyr::filter(survival, grepl("01A", row.names(survival))))
  row.names(survival) <- gsub('-01A','',row.names(survival))
  survival <- survival[,colnames(survival)%in%c('OS', 'OS.time')]
  # id in names()
  htseq_fpkm <- dplyr::select(as.data.frame(htseq_fpkm), ends_with('01A'))
  names(htseq_fpkm) <- gsub("\\.","-", names(htseq_fpkm))
  names(htseq_fpkm) <- gsub('-01A','', names(htseq_fpkm))
  mirna <- dplyr::select(as.data.frame(mirna), ends_with('01A'))
  names(mirna) <- gsub("\\.","-", names(mirna))
  names(mirna) <- gsub('-01A','', names(mirna))
  # get union sampleID
  g1 <- list(names(mirna), names(htseq_fpkm), row.names(GDC_phenotype), row.names(survival))
  union_sampleID <- sort(Reduce(intersect, g1))
  message('\u2605 TCGA-', disease_name, ' cohort contains ',length(union_sampleID), ' tumor samples!')

  # id in row.names()
  clinicDataSort <- function(df){
    df <- df[row.names(df) %in% union_sampleID,]
    df$rowname <- row.names(df)
    df <-  df[sort(df$rowname),]
    df[ , -which(names(df) %in% c("rowname"))]
  }

  GDC_phenotype <- clinicDataSort(GDC_phenotype)
  survival <- clinicDataSort(survival)
  # id in names()
  htseq_fpkm <- htseq_fpkm[,names(htseq_fpkm) %in% union_sampleID]
  htseq_fpkm <- htseq_fpkm[ , order(names(htseq_fpkm))]
  htseq_fpkm <- scale(log(htseq_fpkm+1,2), center = T, scale = T)
  mirna <- mirna[,names(mirna) %in% union_sampleID]
  mirna <- mirna[ , order(names(mirna))]
  mirna <- scale(log(mirna+1,2), center = T, scale = T)

  # mRNA:log2(fpkm+1)，miRNA:log2(RPM+1)
  #mirna <- (2^mirna-1)*1000  #RPKM: X1000
  #htseq_fpkm <- 2^htseq_fpkm -1

  # focus on protein coding RNA because downloaded mRNA expression matrix includes both codingRNA and lncRNA
  gtf_df <- ceRNAR:::gencode_v22_annot
  ensem2symbol <- gtf_df[gtf_df$type == 'gene',c('gene_id', 'gene_type', 'gene_name')]
  rownames(ensem2symbol) <- ensem2symbol$gene_id
  cdRNA <- htseq_fpkm[rownames(htseq_fpkm) %in% ensem2symbol$gene_id[ensem2symbol$gene_type=="protein_coding"],]
  annot_cdRNA <- merge(ensem2symbol, cdRNA, by = 'row.names')
  rownames(annot_cdRNA) <- annot_cdRNA[,1]
  annot_cdRNA <- annot_cdRNA[,-1:-3]
  annot_cdRNA_unique <- aggregate(. ~ gene_name, data = annot_cdRNA, mean) # longer time
  row.names(annot_cdRNA_unique) <- annot_cdRNA_unique$gene_name
  annot_cdRNA_unique <- annot_cdRNA_unique[,-1]

  #miRNA id conversion
  ID_converter <- ceRNAR:::hsa_pre_mature_matching
  mirna$ID <- row.names(mirna)
  miRNA_with_precurer <- merge(ID_converter, mirna, by='ID')[,-1]
  miRNA_with_precurer <- aggregate(. ~ Mature_ID, data = miRNA_with_precurer, mean)
  row.names(miRNA_with_precurer) <- miRNA_with_precurer[,1]
  miRNA_with_precurer <- miRNA_with_precurer[,-1]

  # store processed data
  data.table::fwrite(annot_cdRNA_unique,paste0(project_name,'-', disease_name,'_mrna.csv'), row.names = T)
  data.table::fwrite(miRNA_with_precurer,paste0(project_name,'-', disease_name,'_mirna.csv'), row.names = T)
  data.table::fwrite(GDC_phenotype,paste0(project_name,'-', disease_name,'_phenotype.csv'), row.names = T)
  data.table::fwrite(survival, paste0(project_name,'-', disease_name,'_survival.csv'), row.names = T)
  message('(\u2714) All files have been preprocessed!')
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')
}



