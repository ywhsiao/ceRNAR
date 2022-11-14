## ----style, echo=FALSE, results="asis", message=FALSE-------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ---- warning=FALSE, message=FALSE, echo=FALSE--------------------------------
library(ceRNAR)

## -----------------------------------------------------------------------------
data(gene_exp)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(gene_exp[1:5,1:5])

## -----------------------------------------------------------------------------
data(mirna_exp)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(mirna_exp[1:5,1:5])

## -----------------------------------------------------------------------------
data(surv_data)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(surv_data[1:5,1:2])

## -----------------------------------------------------------------------------
ceRNACustomize(project_name = 'demo', disease_name = 'DLBC', gene_exp = gene_exp,  mirna_exp = mirna_exp, surv_data = surv_data)

## -----------------------------------------------------------------------------
ceRNATCGA(project_name = 'TCGA', disease_name = 'DLBC', timeout = 5000000)

## -----------------------------------------------------------------------------
ceRNAputativePairs(project_name = 'demo', disease_name = 'DLBC', filtering = 'less')

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
ceRNApairFilering(project_name = 'demo', disease_name = "DLBC", window_size = 10)
cernar_result <- SegmentClusteringPlusPeakMerging(project_name = 'demo', disease_name = "DLBC", window_size = 10)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
cernar_result <- ceRNAMethod(project_name = 'demo', disease_name = 'DLBC', window_size = 10)

## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::kable(cernar_result[1:5,])

## ----message=FALSE, results='hide'--------------------------------------------
function_results <- ceRNAFunction(project_name = 'demo', disease_name = 'DLBC', pairs_cutoff = 1)

## ----message=FALSE, fig.dim = c(20, 10), out.width="100%"---------------------
function_results[[1]]

## ----message=FALSE, fig.dim = c(20, 10), out.width="100%"---------------------
function_results[[2]]

## ----message=FALSE, results='hide', fig.dim = c(6, 4), out.width='100%'-------
ceRNALocation(project_name = 'demo', disease_name = 'DLBC', mirna = 'hsa-miR-101-3p', window_size = 10)

## ----message=FALSE, results='hide'--------------------------------------------
survplot_results <- ceRNASurvival(project_name = 'demo', disease_name = 'DLBC', mirnas = 'hsa-miR-101-3p')

## ----message=FALSE, fig.dim = c(8, 4), out.width='100%'-----------------------
survplot_results[[1]]

## ----message=FALSE, results='hide'--------------------------------------------
network_results <- ceRNAModule(project_name = 'demo', disease_name = 'DLBC', pairs_cutoff = 5, column_sum = 1)

## ----message=FALSE, fig.dim = c(8, 6), out.width= '100%'----------------------
network_results[[1]]

## ----message=FALSE, results='hide'--------------------------------------------
external_val_result <- ceRNAValidate(project_name = 'demo', disease_name = 'DLBC')

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(external_val_result[1:5,1:8], row.names = FALSE)

## ----message=FALSE, results='hide'--------------------------------------------
library(SPONGE)
integrated_result <- ceRNAIntegrate(project_name = 'demo', disease_name = 'DLBC')

## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::kable(integrated_result[1:5,], row.names = FALSE)

## -----------------------------------------------------------------------------
sessionInfo()

