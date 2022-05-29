
## ceRNAR package

We present a novel method for identification of ceRNA-miRNA triplets at
specific miRNA expression level. Building on correlation-based approach,
we present a computationally efficient approach to be applied following
identification of ceRNA-miRNA triplets that is related to the etiology
of diseases.

Our approach leverages the genom-wide expression profiles and the
correlation-based approach to identify ceRNA-miRNA triplets.It is
particularly well suited for the studies having both mRNA and miRNA gene
expression data. When applied to such studies, ceRNAR provides a scheme
for unveiling the novel ceRNAs that regulate the biological systems and
analyzing these novel findings for further biological interpretation.

This package also provides several downstream analyses, including
functional, network, survival and peak location analyses, by which to
further investigate the identified ceRNA-miRNA triplets and visualize
the analytical results to aid in the understanding of the role of such
triplets in biological mechanism.

### Installation

``` r
install.packages("devtools")
library(devtools)
install_github("ywhsiao/ceRNAR")
library(ceRNAR)
```

### Example code

#### I. using an example data to conduct all functions in this package

1.  import example data and check the data

``` r
data(gene_exp)
data(mirna_exp)
data(surv_data)

ceRNACustomize(project_name = 'demo', disease_name = 'DLBC', gene_exp = gene_exp, 
               mirna_exp = mirna_exp, surv_data = surv_data)
```

2.  obtain putative mRNA-ceRNA pairs

``` r
ceRNAputativePairs(project_name = 'demo', disease_name = 'DLBC', filtering = 'less')
```

3.  conduct main algorithm through one of following ways

-   through `ceRNAMethod()`

``` r
ceRNAMethod(project_name = 'demo', disease_name = "DLBC", window_size = 45/5)
```

-   through `ceRNApairFilering()` and
    `SegmentClusteringPlusPeakMerging()`

``` r
library(foreach)
ceRNApairFilering(project_name = 'demo', disease_name = "DLBC", window_size = 45/5)
SegmentClusteringPlusPeakMerging(project_name = 'demo', disease_name = "DLBC", 
                                 window_size = 45/5)
```

4.  conduct downstream analysis

``` r
ceRNAFunction(project_name = 'demo', disease_name = 'DLBC', pairs_cutoff = 1)
ceRNALocation(project_name = 'demo', disease_name = 'DLBC', mirna = 'hsa-miR-101-3p', 
              window_size = 45/5)
ceRNAModule(project_name = 'demo', disease_name = 'DLBC', pairs_cutoff = 5, 
            column_sum = 1)
ceRNASurvival(project_name = 'demo', disease_name = 'DLBC', mirnas = 'hsa-miR-101-3p')
ceRNAValidate(project_name = 'demo', disease_name = 'DLBC')
library(SPONGE)
ceRNAIntegrate(project_name = 'demo', disease_name = 'DLBC')

```

#### II. (Alternative) retrieving TCGA data from GDC Xena Hub

1.  retrieve TCGA data, and check the data

``` r
ceRNATCGA(project_name = 'TCGA',disease_name = 'DLBC', timeout = 5000000)
```

2.  change project name and repeat above-mentioned steps from step 2 to
    step 4
