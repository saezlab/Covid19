SARS-CoV-2 dataset: Differential expression analysis
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
14/04/2020

### License Info

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Please check <http://www.gnu.org/licenses/>.

## Introduction

The present script takes the RNAseq data from the study *"SARS-CoV-2
launches* *a unique transcriptional signature from in vitro, ex vivo,
and in vivo systems"*

<https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1>

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>

It performs a differential expression analysis to compare:

  - Human lung epithelial cells (NHBE): mock treated vs infected with
    SARS-CoV-2.

  - A549 alveolar cancer cell line: mock treated vs infected with
    SARS-CoV-2.

  - A549 cell line does not express ACE2, the receptor used by
    SARS-CoV-2 to penetrate into human cells. Therefore A549 were also
    transduced with ACE2 and then mock treated or infected with
    SARS-CoV-2

  - Calu-3 human lung epithelial cancer cell line: mock treated vs
    infected with SARS-CoV-2.

We used the *DESeq2* R package.

## Getting Started

We first load the required libraries.

``` r
library(dplyr)
library(DESeq2)
```

We also read the raw counts from the original experiment (GSE147507)

``` r
## Raw counts table
GSE147507_raw_counts <- 
    read.csv("RawData/GSE147507_RawReadCounts_Human.tsv", sep = "\t")
```

## NHBE mock treated vs infected with SARS-CoV-2

We first select the Series 1, which corresponds to independent
biological triplicates of primary human lung epithelium (NHBE) that were
either mock treated or infected with SARS-CoV-2.

``` r
## We select series 1 as described above.
count_NHBEvsCOV2_df <- GSE147507_raw_counts[,c(2:7)]
row.names(count_NHBEvsCOV2_df) <- GSE147507_raw_counts$X
```

We create a dataframe describing the experimental design.

``` r
targets_NHBEvsCOV2 <- 
    as.data.frame(matrix(NA,length(names(count_NHBEvsCOV2_df)),1))
names(targets_NHBEvsCOV2) <- c("condition")
row.names(targets_NHBEvsCOV2) <- names(count_NHBEvsCOV2_df)
targets_NHBEvsCOV2$condition <- 
    gsub("Series1_", "", row.names(targets_NHBEvsCOV2))
targets_NHBEvsCOV2$condition <- 
    factor(gsub("_[1-3]$", "", targets_NHBEvsCOV2$condition))
targets_NHBEvsCOV2
```

    ##                                 condition
    ## Series1_NHBE_Mock_1             NHBE_Mock
    ## Series1_NHBE_Mock_2             NHBE_Mock
    ## Series1_NHBE_Mock_3             NHBE_Mock
    ## Series1_NHBE_SARS.CoV.2_1 NHBE_SARS.CoV.2
    ## Series1_NHBE_SARS.CoV.2_2 NHBE_SARS.CoV.2
    ## Series1_NHBE_SARS.CoV.2_3 NHBE_SARS.CoV.2

Then, we perform the differential expression analysis with *DESeq2*

``` r
## Create deseq2 object
dds_NHBEvsCOV2 <- 
    DESeqDataSetFromMatrix(countData = as.matrix(count_NHBEvsCOV2_df), 
    colData = targets_NHBEvsCOV2, design = ~ condition)

## Set control
dds_NHBEvsCOV2$condition <- relevel(dds_NHBEvsCOV2$condition, 
    ref = levels(targets_NHBEvsCOV2$condition)[1])

## Carry out diff exp
dds_NHBEvsCOV2 <- DESeq(dds_NHBEvsCOV2)
```

We finally save the table with the results of the analysis and the
normalised counts

``` r
## See the comparisons carried out
comparison_NHBEvsCOV2 <- resultsNames(dds_NHBEvsCOV2)

## Get results table
results_NHBEvsCOV2 <- 
    results(dds_NHBEvsCOV2, name=comparison_NHBEvsCOV2[2])

## Save the object
saveRDS(results_NHBEvsCOV2, file="IntermediateFiles/dds_results_NHBEvsCOV2.rds")

## Extract normalised counts data
dds_NHBEvsCOV2 <- estimateSizeFactors(dds_NHBEvsCOV2)
norm_counts_NHBEvsCOV2 <- counts(dds_NHBEvsCOV2, normalized = TRUE)
saveRDS(norm_counts_NHBEvsCOV2, 
    file="IntermediateFiles/counts_norm_NHBEvsCOV2.rds")
```

## A549 mock treated vs infected with SARS-CoV-2

Then, we select the Series 2, which corresponds to independent
biological triplicates of the A549 alveolar cancer cell line that were
either mock treated or infected with SARS-CoV-2.

``` r
## We select series 2 as described above.
count_A549vsCOV2_df <- GSE147507_raw_counts[,c(8:13)]
row.names(count_A549vsCOV2_df) <- GSE147507_raw_counts$X
```

We create a dataframe describing the experimental design.

``` r
targets_A549vsCOV2 <- 
    as.data.frame(matrix(NA,length(names(count_A549vsCOV2_df)),1))
names(targets_A549vsCOV2) <- c("condition")
row.names(targets_A549vsCOV2) <- names(count_A549vsCOV2_df)
targets_A549vsCOV2$condition <- 
    gsub("Series2_", "", row.names(targets_A549vsCOV2))
targets_A549vsCOV2$condition <- 
    factor(gsub("_[1-3]$", "", targets_A549vsCOV2$condition))
targets_A549vsCOV2
```

    ##                                 condition
    ## Series2_A549_Mock_1             A549_Mock
    ## Series2_A549_Mock_2             A549_Mock
    ## Series2_A549_Mock_3             A549_Mock
    ## Series2_A549_SARS.CoV.2_1 A549_SARS.CoV.2
    ## Series2_A549_SARS.CoV.2_2 A549_SARS.CoV.2
    ## Series2_A549_SARS.CoV.2_3 A549_SARS.CoV.2

Then, we perform the differential expression analysis with *DESeq2*

``` r
## Create deseq2 object
dds_A549vsCOV2 <- 
    DESeqDataSetFromMatrix(countData = as.matrix(count_A549vsCOV2_df), 
    colData = targets_A549vsCOV2, design = ~ condition)

## Set control
dds_A549vsCOV2$condition <- relevel(dds_A549vsCOV2$condition, 
    ref = levels(targets_A549vsCOV2$condition)[1])

## Carry out diff exp
dds_A549vsCOV2 <- DESeq(dds_A549vsCOV2)
```

We finally save the table with the results of the analysis and the
normalised counts:

``` r
## See the comparisons carried out
comparison_A549vsCOV2 <- resultsNames(dds_A549vsCOV2)

## Get results table
results_A549vsCOV2 <- 
    results(dds_A549vsCOV2, name=comparison_A549vsCOV2[2])

## Save the object
saveRDS(results_A549vsCOV2, file="IntermediateFiles/dds_results_A549vsCOV2.rds")

## Extract normalised counts data
dds_A549vsCOV2 <- estimateSizeFactors(dds_A549vsCOV2)
norm_counts_A549vsCOV2 <- counts(dds_A549vsCOV2, normalized = TRUE)
saveRDS(norm_counts_A549vsCOV2, 
    file="IntermediateFiles/counts_norm_A549vsCOV2.rds")
```

## A549 transduced with ACE2 mock treated vs infected with SARS-CoV-2

A549 does not express the ACE2 receptor used by SARS-CoV-2 to infect
cells. Therefore, A549 were also transduced with a vector expressing
human ACE2, were also mock treated or infected with SARS-CoV-2.

``` r
## We select series 6 as described above.
count_A549ACE2vsCOV2_df <- GSE147507_raw_counts[,c(28:33)]
row.names(count_A549ACE2vsCOV2_df) <- GSE147507_raw_counts$X
```

We create a dataframe describing the experimental design.

``` r
targets_A549ACE2vsCOV2 <- 
    as.data.frame(matrix(NA,length(names(count_A549ACE2vsCOV2_df)),1))
names(targets_A549ACE2vsCOV2) <- c("condition")
row.names(targets_A549ACE2vsCOV2) <- names(count_A549ACE2vsCOV2_df)
targets_A549ACE2vsCOV2$condition <- 
    gsub("Series2_", "", row.names(targets_A549ACE2vsCOV2))
targets_A549ACE2vsCOV2$condition <- 
    factor(gsub("_[1-3]$", "", targets_A549ACE2vsCOV2$condition))
targets_A549ACE2vsCOV2
```

    ##                                                   condition
    ## Series6_A549.ACE2_Mock_1             Series6_A549.ACE2_Mock
    ## Series6_A549.ACE2_Mock_2             Series6_A549.ACE2_Mock
    ## Series6_A549.ACE2_Mock_3             Series6_A549.ACE2_Mock
    ## Series6_A549.ACE2_SARS.CoV.2_1 Series6_A549.ACE2_SARS.CoV.2
    ## Series6_A549.ACE2_SARS.CoV.2_2 Series6_A549.ACE2_SARS.CoV.2
    ## Series6_A549.ACE2_SARS.CoV.2_3 Series6_A549.ACE2_SARS.CoV.2

Then, we perform the differential expression analysis with *DESeq2*

``` r
## Create deseq2 object
dds_A549ACE2vsCOV2 <- 
    DESeqDataSetFromMatrix(countData = as.matrix(count_A549ACE2vsCOV2_df), 
    colData = targets_A549ACE2vsCOV2, design = ~ condition)

## Set control
dds_A549ACE2vsCOV2$condition <- relevel(dds_A549ACE2vsCOV2$condition, 
    ref = levels(targets_A549ACE2vsCOV2$condition)[1])

## Carry out diff exp
dds_A549ACE2vsCOV2 <- DESeq(dds_A549ACE2vsCOV2)
```

We finally save the table with the results of the analysis and the
normalised counts:

``` r
## See the comparisons carried out
comparison_A549ACE2vsCOV2 <- resultsNames(dds_A549ACE2vsCOV2)

## Get results table
results_A549ACE2vsCOV2 <- 
    results(dds_A549ACE2vsCOV2, name=comparison_A549ACE2vsCOV2[2])

## Save the object
saveRDS(results_A549ACE2vsCOV2, file="IntermediateFiles/dds_results_A549ACE2vsCOV2.rds")

## Extract normalised counts data
dds_A549ACE2vsCOV2 <- estimateSizeFactors(dds_A549ACE2vsCOV2)
norm_counts_A549ACE2vsCOV2 <- counts(dds_A549ACE2vsCOV2, normalized = TRUE)
saveRDS(norm_counts_A549ACE2vsCOV2, 
    file="IntermediateFiles/counts_norm_A549ACE2vsCOV2.rds")
```

## CALU-3 mock treated vs infected with SARS-CoV-2

We then select the Series 7, which corresponds to independent biological
triplicates of a human lung epithelium cancer cell line named Calu-3,
that were either mock treated or infected with SARS-CoV-2.

``` r
## We select series 7 as described above.
count_CALU3vsCOV2_df <- GSE147507_raw_counts[,c(34:39)]
row.names(count_CALU3vsCOV2_df) <- GSE147507_raw_counts$X
```

We create a dataframe describing the experimental design.

``` r
targets_CALU3vsCOV2 <- 
    as.data.frame(matrix(NA,length(names(count_CALU3vsCOV2_df)),1))
names(targets_CALU3vsCOV2) <- c("condition")
row.names(targets_CALU3vsCOV2) <- names(count_CALU3vsCOV2_df)
targets_CALU3vsCOV2$condition <- 
    gsub("Series1_", "", row.names(targets_CALU3vsCOV2))
targets_CALU3vsCOV2$condition <- 
    factor(gsub("_[1-3]$", "", targets_CALU3vsCOV2$condition))
targets_CALU3vsCOV2
```

    ##                                           condition
    ## Series7_Calu3_Mock_1             Series7_Calu3_Mock
    ## Series7_Calu3_Mock_2             Series7_Calu3_Mock
    ## Series7_Calu3_Mock_3             Series7_Calu3_Mock
    ## Series7_Calu3_SARS.CoV.2_1 Series7_Calu3_SARS.CoV.2
    ## Series7_Calu3_SARS.CoV.2_2 Series7_Calu3_SARS.CoV.2
    ## Series7_Calu3_SARS.CoV.2_3 Series7_Calu3_SARS.CoV.2

Then, we perform the differential expression analysis with *DESeq2*

``` r
## Create deseq2 object
dds_CALU3vsCOV2 <- 
    DESeqDataSetFromMatrix(countData = as.matrix(count_CALU3vsCOV2_df), 
    colData = targets_CALU3vsCOV2, design = ~ condition)

## Set control
dds_CALU3vsCOV2$condition <- relevel(dds_CALU3vsCOV2$condition, 
    ref = levels(targets_CALU3vsCOV2$condition)[1])

## Carry out diff exp
dds_CALU3vsCOV2 <- DESeq(dds_CALU3vsCOV2)
```

We finally save the table with the results of the analysis and the
normalised counts

``` r
## See the comparisons carried out
comparison_CALU3vsCOV2 <- resultsNames(dds_CALU3vsCOV2)

## Get results table
results_CALU3vsCOV2 <- 
    results(dds_CALU3vsCOV2, name=comparison_CALU3vsCOV2[2])

## Save the object
saveRDS(results_CALU3vsCOV2, file="IntermediateFiles/dds_results_CALU3vsCOV2.rds")

## Extract normalised counts data
dds_CALU3vsCOV2 <- estimateSizeFactors(dds_CALU3vsCOV2)
norm_counts_CALU3vsCOV2 <- counts(dds_CALU3vsCOV2, normalized = TRUE)
saveRDS(norm_counts_CALU3vsCOV2, 
    file="IntermediateFiles/counts_norm_CALU3vsCOV2.rds")
```

## Session Info Details

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/local/lib/R/lib/libRblas.so
    ## LAPACK: /usr/local/lib/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DESeq2_1.26.0               SummarizedExperiment_1.16.0
    ##  [3] DelayedArray_0.12.0         BiocParallel_1.20.0        
    ##  [5] matrixStats_0.56.0          Biobase_2.46.0             
    ##  [7] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
    ##  [9] IRanges_2.20.1              S4Vectors_0.24.1           
    ## [11] BiocGenerics_0.32.0         dplyr_0.8.5                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7            splines_3.6.3          Formula_1.2-3         
    ##  [4] assertthat_0.2.1       latticeExtra_0.6-29    blob_1.2.1            
    ##  [7] GenomeInfoDbData_1.2.2 yaml_2.2.1             RSQLite_2.2.0         
    ## [10] pillar_1.4.3           backports_1.1.5        lattice_0.20-41       
    ## [13] glue_1.4.0             digest_0.6.25          RColorBrewer_1.1-2    
    ## [16] XVector_0.26.0         checkmate_2.0.0        colorspace_1.4-1      
    ## [19] htmltools_0.4.0        Matrix_1.2-18          XML_3.99-0.3          
    ## [22] pkgconfig_2.0.3        genefilter_1.68.0      zlibbioc_1.32.0       
    ## [25] xtable_1.8-4           purrr_0.3.3            scales_1.1.0          
    ## [28] jpeg_0.1-8.1           tibble_3.0.0           htmlTable_1.13.3      
    ## [31] annotate_1.64.0        ggplot2_3.3.0          ellipsis_0.3.0        
    ## [34] nnet_7.3-13            cli_2.0.2              survival_3.1-11       
    ## [37] magrittr_1.5           crayon_1.3.4           memoise_1.1.0         
    ## [40] evaluate_0.14          fansi_0.4.1            foreign_0.8-76        
    ## [43] tools_3.6.3            data.table_1.12.8      lifecycle_0.2.0       
    ## [46] stringr_1.4.0          locfit_1.5-9.4         munsell_0.5.0         
    ## [49] cluster_2.1.0          AnnotationDbi_1.48.0   compiler_3.6.3        
    ## [52] rlang_0.4.5            grid_3.6.3             RCurl_1.98-1.1        
    ## [55] rstudioapi_0.11        htmlwidgets_1.5.1      bitops_1.0-6          
    ## [58] base64enc_0.1-3        rmarkdown_2.1          gtable_0.3.0          
    ## [61] DBI_1.1.0              R6_2.4.1               gridExtra_2.3         
    ## [64] knitr_1.28             bit_1.1-15.2           Hmisc_4.4-0           
    ## [67] stringi_1.4.6          Rcpp_1.0.4             geneplotter_1.64.0    
    ## [70] vctrs_0.2.4            rpart_4.1-15           acepack_1.4.1         
    ## [73] png_0.1-7              tidyselect_1.0.0       xfun_0.12
