#'
#' @name all_steps_interface
#' @title A function for conducting the algorithm through the toy data
#' @description A function to allow users to conducting the
#'
#' @import utils
#'
#' @param path_prefix user's working directory
#' @param project_name the project name that users can assign (default: demo)
#' @param disease_name the abbreviation of disease that users are interested in
#' (default: DLBC)
#' @param gene_exp location of gene expression data (default: gene_exp)
#' @param mirna_exp location of miRNA expression data (default: mirna_exp)
#' @param surv_data location of survival data (default: surv_data)
#' @param window_size the number of samples for each window (default:10)
#' @param cor_method selection of correlation methods, including pearson and
#' spearman (default: pearson)
#' @param cor_threshold_peak peak threshold of correlation value between 0 and 1
#' (default: 0.85)
#'
#' @returns a dataframe object
#' @export
#'
#' @examples
#' data(gene_exp)
#' data(mirna_exp)
#' data(surv_data)
#' all_steps_interface(
#' path_prefix = path_prefix,
#' project_name = 'demo',
#' disease_name = 'DLBC',
#' gene_exp = gene_exp,
#' mirna_exp = mirna_exp,
#' surv_data = surv_data,
#' filtering = 'less',
#' window_size = 10,
#' cor_method = 'pearson',
#' cor_threshold_peak = 0.85)

all_steps_interface <- function(path_prefix,
                                project_name = 'demo',
                                disease_name = 'DLBC',
                                gene_exp = gene_exp,
                                mirna_exp = mirna_exp,
                                surv_data = surv_data,
                                filtering = 'less',
                                window_size = 10,
                                cor_method = 'pearson',
                                cor_threshold_peak = 0.85){
  # import data
  ceRNAR::ceRNACustomize(path_prefix = path_prefix,
                         project_name = 'demo',
                         disease_name = 'DLBC',
                         gene_exp = gene_exp,
                         mirna_exp = mirna_exp,
                         surv_data = surv_data)
  # putative pairs
  ceRNAR::ceRNAputativePairs(path_prefix = path_prefix,
                             project_name = 'demo',
                             disease_name = 'DLBC',
                             filtering = 'less')
  # method: pair filtering +segment clustering
  ceRNAR::ceRNAMethod(path_prefix = path_prefix,
                      project_name = 'demo',
                      disease_name = 'DLBC',
                      window_size = 10,
                      cor_method = 'pearson',
                      cor_threshold_peak = 0.85)


}
