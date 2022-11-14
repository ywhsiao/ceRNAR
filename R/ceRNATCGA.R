#'
#' @name ceRNATCGA
#' @title Retrieval of public TCGA data from GDC Xena Hub
#' @description A function to retrieve TCGA data from GDC Xena Hub
#' (https://xenabrowser.net/datapages/)
#'
#' @import utils
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default = 'TCGA')
#' @param disease_name the abbreviation of disease that users are interested in
#' (default = 'DLBC')
#' @param timeout the allowance time for downloading TCGA data
#' (default = 1000)
#'
#' @returns file
#' @export
#'
#' @examples
#' ceRNATCGA(
#' path_prefix = NULL,
#' project_name = 'TCGA',
#' disease_name = 'DLBC',
#' timeout = 500000
#' )
#'



ceRNATCGA <- function(path_prefix = NULL,
                      project_name = 'TCGA',
                      disease_name = 'DLBC',
                      timeout = 500000){

  if (is.null(path_prefix)){
    path_prefix <- fs::path_home()
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

  # download cancer files (phenotype, survival, miRNA and mRNA) from gdc resource
  downloadFromGDC <- function(project, cancer, timeout=5000000){
    message('\u25CF TCGA data: Downloading the data ...')
    options(timeout=timeout) # to test the largest data
    HelpersMG::wget(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project_name,'-', disease_name,'.GDC_phenotype.tsv.gz'), destfile = paste0(path_prefix,project_name,'-', disease_name, '/01_rawdata/',project_name,'-', disease_name,'.GDC_phenotype.tsv.gz'))
    HelpersMG::wget(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project_name,'-', disease_name,'.GDC_phenotype.tsv.gz'), destfile = paste0(path_prefix,project_name,'-', disease_name, '/01_rawdata/',project_name,'-', disease_name,'.GDC_phenotype.tsv.gz'))
    HelpersMG::wget(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project_name,'-', disease_name,'.survival.tsv'), destfile = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/',project_name,'-', disease_name,'.survival.tsv'))
    HelpersMG::wget(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project_name,'-', disease_name,'.mirna.tsv.gz'), destfile = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/',project_name,'-', disease_name,'.mirna.tsv.gz'))
    HelpersMG::wget(paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project_name,'-', disease_name,'.htseq_fpkm.tsv.gz'), destfile = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/',project_name,'-', disease_name,'.htseq_fpkm.tsv.gz'))
  }

  downloadFromGDC(project,cancer)

  # unzip
  temp <-  list.files(path = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'), pattern=".gz")
  for (i in 1:length(temp)) R.utils::gunzip(paste0(path_prefix,project_name,'-', disease_name, '/01_rawdata/',temp[i]), remove=TRUE, overwrite = TRUE)
  message('(\u2714) All files have been and downloaded and uncompressed!')

  # match sample id
  # Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.
  # ignore GBM because its mirna data only contain 5 samples
  temp <-  list.files(path = paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata'), pattern="*.tsv")
  temp_name <- gsub(".*\\.","",gsub(".tsv", "", temp))
  for (i in 1:length(temp)) assign(temp_name[i], data.frame(data.table::fread(paste0(path_prefix,project_name,'-', disease_name, '/01_rawdata/', temp[i]), sep = '\t',header = TRUE, dec = ".", fill = TRUE), row.names = 1))
  # id in a column
  GDC_phenotype <- GDC_phenotype[substring(row.names(GDC_phenotype),16)=='A',]
  GDC_phenotype <- GDC_phenotype[!(substring(row.names(GDC_phenotype),14,15) %in% seq(10,29,1)),]
  GDC_phenotype <- GDC_phenotype[!duplicated(GDC_phenotype$submitter_id),]
  row.names(GDC_phenotype) <-substring(row.names(GDC_phenotype),1,12)
  survival <- survival[substring(row.names(survival),16)=='A',]
  survival <- survival[!(substring(row.names(GDC_phenotype),14,15) %in% seq(10,29,1)),]
  survival <- survival[,colnames(survival)%in%c('OS', 'OS.time')]
  survival <- survival[!duplicated(substring(row.names(survival),1,12)),]
  row.names(survival) <-substring(row.names(survival),1,12)
  # id in names()
  names(htseq_fpkm) <- gsub("\\.","-", names(htseq_fpkm))
  htseq_fpkm <- htseq_fpkm[,substring(names(htseq_fpkm),16)=='A']
  htseq_fpkm <- htseq_fpkm[,!(substring(names(htseq_fpkm),14,15) %in% seq(10,29,1))]
  htseq_fpkm <- htseq_fpkm[,!duplicated(substring(names(htseq_fpkm),1,12))]
  names(htseq_fpkm) <-substring(names(htseq_fpkm),1,12)
  names(mirna) <- gsub("\\.","-", names(mirna))
  mirna <- mirna[,substring(names(mirna),16)=='A']
  mirna <- mirna[,!(substring(names(mirna),14,15) %in% seq(10,29,1))]
  mirna <- mirna[,!duplicated(substring(names(mirna),1,12))]
  names(mirna) <-substring(names(mirna),1,12)
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
  cdRNA <- htseq_fpkm[rownames(htseq_fpkm) %in% ensem2symbol$gene_id[ensem2symbol$gene_type=="protein_coding"],]
  annot_cdRNA <- merge(ensem2symbol, cdRNA, by = 'row.names')
  rownames(annot_cdRNA) <- annot_cdRNA[,1]
  annot_cdRNA <- annot_cdRNA[,-1:-3]
  annot_cdRNA_unique <- stats::aggregate(. ~ gene_name, data = annot_cdRNA, mean) # longer time
  row.names(annot_cdRNA_unique) <- annot_cdRNA_unique$gene_name
  annot_cdRNA_unique <- annot_cdRNA_unique[,-1]

  #miRNA id conversion
  #ID_converter <- ceRNAR:::hsa_pre_mature_matching
  ID_converter <- get0("hsa_pre_mature_matching", envir = asNamespace("ceRNAR"))
  mirna$ID <- row.names(mirna)
  miRNA_with_precurer <- merge(ID_converter, mirna, by='ID')[,-1]
  miRNA_with_precurer <- stats::aggregate(. ~ Mature_ID, data = miRNA_with_precurer, mean)
  row.names(miRNA_with_precurer) <- miRNA_with_precurer[,1]
  miRNA_with_precurer <- miRNA_with_precurer[,-1]

  # store processed data
  data.table::fwrite(as.data.frame(annot_cdRNA_unique),paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mrna.csv'), row.names = TRUE)
  data.table::fwrite(as.data.frame(miRNA_with_precurer),paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_mirna.csv'), row.names = TRUE)
  data.table::fwrite(as.data.frame(GDC_phenotype),paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_phenotype.csv'), row.names = TRUE)
  data.table::fwrite(as.data.frame(survival), paste0(path_prefix, project_name,'-', disease_name, '/01_rawdata/', project_name,'-', disease_name,'_survival.csv'), row.names = TRUE)
  message('(\u2714) All files have been preprocessed!')
  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
  message('\u2605\u2605\u2605 Ready to next step! \u2605\u2605\u2605')
}
