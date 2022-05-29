#'
#' @name ceRNACustomize
#' @title A function for uploading customized data
#' @description A function to allow users to upload their own data
#'
#' @import data.table
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign
#' @param disease_name the abbreviation of disease that users are interested in
#' @param gene_exp location of gene expression data
#' @param mirna_exp location of miRNA expression data
#' @param surv_data location of survival data
#'
#' @examples
#' load('data/ExampleData.Rdata')
#' ceRNACustomize(
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' gene_exp = gene_exp,
#' mirna_exp = mirna_exp,
#' surv_data = surv_data
#' )
#'
#' @export


ceRNACustomize <- function(path_prefix = NULL,
                           project_name = NULL,
                           disease_name = NULL,
                           gene_exp = NULL,
                           mirna_exp = NULL,
                           surv_data){
  if (!is.null(path_prefix)){
    setwd(path_prefix)
    message('Your current directory: ', getwd())
  }else{
    path_prefix <- getwd()
    message('Your current directory: ', getwd())
  }

  time1 <- Sys.time()
  message('\u25CF Step 1 for Customized data: Checking the data...')
  # check input
  if (is.null(project_name)|is.null(disease_name)|is.null(gene_exp)|is.null(mirna_exp)){
    stop("project_name/disease_name/gene_exp/mirna_exp has not provided!")
  }else{
    if ('gene_exp' %in% ls()){
      exp <- gene_exp
    }else{
      exp <- as.data.frame(fread(gene_exp,header = T))
    }
    if ('mirna_exp' %in% ls()){
      mirna <- mirna_exp
    }else{
      mirna <- as.data.frame(fread(mirna_exp,header = T))
    }

  }
  if (!is.null(surv_data)){
    if ('surv_data' %in% ls()){
      surv <- surv_data
    }else{
      surv <- as.data.frame(fread(surv_data,header = T))
    }

    message('(\u2714) Input data involve mRNA, miRNA and survival data.')
  }else{
    message('(\u2714) Input data ONLY involve mRNA and miRNA.')
  }

  # check geneid
  gtf_df <- ceRNAR:::gencode_v22_annot
  genecode_v22 <- unique(gtf_df$gene_name)
  exp_gene_name <- row.names(exp)
  if (sum(exp_gene_name %in% genecode_v22) == 0){
    stop("Gene names in this dataset ate not matched!")
  }else{
    message('(\u2714) All the gene names are matched!')
  }

  # check mirnaid
  ID_converter <- ceRNAR:::hsa_pre_mature_matching
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
  dirname <- paste0(project_name,'-',disease_name)
  if (dir.exists(dirname) == FALSE){
    dir.create(dirname)
  }
  if (dir.exists(paste0(project_name,'-',disease_name,'/01_rawdata')) == FALSE){
    dir.create(paste0(project_name,'-',disease_name,'/01_rawdata'))
  }
  fwrite(exp,paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mrna.csv'), row.names = T)
  fwrite(mirna,paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_mirna.csv'), row.names = T)
  fwrite(surv, paste0(project_name,'-',disease_name,'/01_rawdata/',project_name,'-', disease_name,'_survival.csv'), row.names = T)

  CatchupPause <- function(Secs){
    Sys.sleep(Secs) #pause to let connection work
    closeAllConnections()
  }
  CatchupPause(3)

  time2 <- Sys.time()
  diftime <- difftime(time2, time1, units = 'min')
  message(paste0('\u2605 Consuming time: ',round(as.numeric(diftime)), ' minutes.'))
  message('\u2605\u2605\u2605 All the files are checked for further next step! \u2605\u2605\u2605')
}

