SARS-CoV-2 dataset: Identifying Cell cycle Phase
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
30/04/2020

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

The present script deals with the RNAseq data from the study
*"SARS-CoV-2 launches* *a unique transcriptional signature from in
vitro, ex vivo, and in vivo systems"*

<https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1>

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>

It uses a **Progeny** based approach to compute the activity of the
different cell cycles phases in the samples. In order to do so, we used
the [cyclebase](https://cyclebase.org/CyclebaseSearch) database, which
contains summarized information of genome-wide cell-cycle-related
experiments. Our approach is applied to the four conditions under study.

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

## Getting Started

We first load the required libraries.

``` r
library(dplyr)
library(tidyr)
library(tibble)
library(progeny)
library(ggplot2)
```

## Prior Knowledge information about the cell cycle phase

From the [cyclebase](https://cyclebase.org/CyclebaseSearch) database, we
downloaded the list of genes that show a periodic peak in a given
cellular cycle stage. The information looks as follows:

<br><br> ![](RawData/CyclePhasePhoto.png) <br><br>

``` r
periodicGenes <- 
    read.csv(file = "RawData/human_periodic.tsv", header = TRUE, 
        sep= "\t", stringsAsFactors = FALSE) %>%
    dplyr::filter(organism == "9606")  
    
GenesPhase <- read.csv(file = "RawData/human_periodic_cellphase.csv", 
    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
```

We are now going to use this information to build a signature matrix. To
do so, we perform a -log10 transformation of the pvalues. In addition,
we select the top 35 genes (according to their pvalue) for each cell
phase. We chose 35, because one of the phases is only related to 36
significant genes. Proceeding this way, we start the analysis in the
same conditions for all the phases.

``` r
## Parameter to decide the number of significant genes per cell phase
top <- 35 # Because one phase has only 36 related genes. 

periodicGenes_HGNC <- periodicGenes %>% 
    dplyr::left_join(GenesPhase,  by = c("gene" = "Identifier")) %>%
    dplyr::select(Matched.name, periodicity_pvalue, Peaktime) %>% 
    dplyr::mutate(LogPvalue = -log10(periodicity_pvalue)) %>%
    dplyr::distinct(Matched.name, Peaktime,  .keep_all = TRUE) %>% 
    dplyr::group_by(Peaktime) %>% 
    dplyr::top_n(top, wt = LogPvalue) %>%
    dplyr::ungroup(Peaktime) %>%
    dplyr::select(-periodicity_pvalue) %>% 
    tidyr::spread(Peaktime, LogPvalue, fill=0) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

## This step is not needed if we run Progeny with permutaitons. It is a 
## normalization by column to assure that all phases are in the same conditions.
periodicGenes_HGNC_Norm <- apply(periodicGenes_HGNC, 2, function(x) x/sum(x))
```

## Running Progeny with permutations to determine Cell Phase activity

We are now going to take the results from the differential expression
analysis for our four cell lines under study. We use the statistic to
run Progeny with permutations with the above-defined signature matrix.

### NHBE cell line

First for the NHBE cell line:

``` r
dds_NHBEvsCOV2 <- readRDS("IntermediateFiles/dds_results_NHBEvsCOV2.rds")

## We prepare the dataframe to Run progeny
dds_NHBEvsCOV2_df <- as.data.frame(dds_NHBEvsCOV2) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") 

expr <- data.frame(names = row.names(dds_NHBEvsCOV2_df), row.names = NULL, 
    dds_NHBEvsCOV2_df)
model <- data.frame(names = row.names(periodicGenes_HGNC_Norm), 
    row.names = NULL, periodicGenes_HGNC_Norm)

results_NHBE <- progenyPerm(expr, model, k = 10000, z_scores = TRUE) %>% t()
colnames(results_NHBE) <- "NES"

## We prepare the data for the plot. 
Cycle_NHBEvsCOV2_zscore_df <- as.data.frame(results_NHBE) %>% 
    rownames_to_column(var = "Stage") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Stage = factor(Stage))
```

### A549 cell line

Now for the A549 cell line:

``` r
dds_A549vsCOV2 <- readRDS("IntermediateFiles/dds_results_A549vsCOV2.rds")

## We prepare the dataframe to Run progeny
dds_A549vsCOV2_df <- as.data.frame(dds_A549vsCOV2) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") 

expr <- data.frame(names = row.names(dds_A549vsCOV2_df), row.names = NULL, 
    dds_A549vsCOV2_df)
model <- data.frame(names = row.names(periodicGenes_HGNC_Norm), 
    row.names = NULL, periodicGenes_HGNC_Norm)

results_A549 <- progenyPerm(expr, model, k = 10000, z_scores = TRUE) %>% t()
colnames(results_A549) <- "NES"

## We prepare the data for the plot. 
Cycle_A549vsCOV2_zscore_df <- as.data.frame(results_A549) %>% 
    rownames_to_column(var = "Stage") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Stage = factor(Stage))
```

### A549 cell line transduced with ACE2 expression

We repeat the process for the A549 cell line transduced with ACE2
expression

``` r
dds_A549ACE2vsCOV2 <- 
    readRDS("IntermediateFiles/dds_results_A549ACE2vsCOV2.rds")

## We prepare the dataframe to Run progeny
dds_A549ACE2vsCOV2_df <- as.data.frame(dds_A549ACE2vsCOV2) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") 

expr <- data.frame(names = row.names(dds_A549ACE2vsCOV2_df), row.names = NULL, 
    dds_A549ACE2vsCOV2_df)
model <- 
    data.frame(names = row.names(periodicGenes_HGNC_Norm), row.names = NULL, 
    periodicGenes_HGNC_Norm)

results_A549ACE2 <- progenyPerm(expr, model, k = 10000, z_scores = TRUE) %>% t()
colnames(results_A549ACE2) <- "NES"

## We prepare the data for the plot. 
Cycle_A549ACE2vsCOV2_zscore_df <- as.data.frame(results_A549ACE2) %>% 
    rownames_to_column(var = "Stage") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Stage = factor(Stage))
```

### CALU-3 cell line

We finally repeat the process for the CALU-3 cell line

``` r
dds_CALU3vsCOV2 <- readRDS("IntermediateFiles/dds_results_CALU3vsCOV2.rds")

## We prepare the dataframe to Run progeny
dds_CALU3vsCOV2_df <- as.data.frame(dds_CALU3vsCOV2) %>% 
    rownames_to_column(var = "GeneID") %>% 
    dplyr::select(GeneID, stat) %>% 
    dplyr::filter(!is.na(stat)) %>% 
    column_to_rownames(var = "GeneID") 

expr <- data.frame(names = row.names(dds_CALU3vsCOV2_df), row.names = NULL, 
    dds_CALU3vsCOV2_df)
model <- data.frame(names = row.names(periodicGenes_HGNC_Norm), 
    row.names = NULL, periodicGenes_HGNC_Norm)

results_CALU3 <- progenyPerm(expr, model, k = 10000, z_scores = TRUE) %>% t()
colnames(results_CALU3) <- "NES"

## We prepare the data for the plot. 
Cycle_CALU3vsCOV2_zscore_df <- as.data.frame(results_CALU3) %>% 
    rownames_to_column(var = "Stage") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Stage = factor(Stage))
```

## Plotting the results

We display the results in a single plot

``` r
Cycle_NHBEvsCOV2_zscore_df<- Cycle_NHBEvsCOV2_zscore_df %>% 
    add_column(cellLine = "NHBE")

Cycle_A549vsCOV2_zscore_df<- Cycle_A549vsCOV2_zscore_df %>% 
    add_column(cellLine = "A549")

Cycle_A549ACE2vsCOV2_zscore_df<- Cycle_A549ACE2vsCOV2_zscore_df %>% 
    add_column(cellLine = "A549_ACE2")

Cycle_CALU3vsCOV2_zscore_df<- Cycle_CALU3vsCOV2_zscore_df %>% 
    add_column(cellLine = "CALU-3")

All_Lines <- bind_rows(Cycle_NHBEvsCOV2_zscore_df, Cycle_A549vsCOV2_zscore_df,
    Cycle_A549ACE2vsCOV2_zscore_df,Cycle_CALU3vsCOV2_zscore_df)

p <- ggplot(All_Lines,aes(x = Stage, y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    facet_wrap(~cellLine) +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Cell Cycle Stage") 
```

![](CellCycleStage_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

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
    ## [11] BiocGenerics_0.32.0         ggplot2_3.3.0              
    ## [13] progeny_1.9.7               tibble_3.0.0               
    ## [15] tidyr_1.0.2                 dplyr_0.8.5                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7            splines_3.6.3          Formula_1.2-3         
    ##  [4] assertthat_0.2.1       latticeExtra_0.6-29    blob_1.2.1            
    ##  [7] GenomeInfoDbData_1.2.2 yaml_2.2.1             ggrepel_0.8.2         
    ## [10] RSQLite_2.2.0          pillar_1.4.3           backports_1.1.5       
    ## [13] lattice_0.20-41        glue_1.4.0             digest_0.6.25         
    ## [16] RColorBrewer_1.1-2     XVector_0.26.0         checkmate_2.0.0       
    ## [19] colorspace_1.4-1       htmltools_0.4.0        Matrix_1.2-18         
    ## [22] XML_3.99-0.3           pkgconfig_2.0.3        genefilter_1.68.0     
    ## [25] zlibbioc_1.32.0        xtable_1.8-4           purrr_0.3.3           
    ## [28] scales_1.1.0           jpeg_0.1-8.1           annotate_1.64.0       
    ## [31] htmlTable_1.13.3       farver_2.0.3           ellipsis_0.3.0        
    ## [34] withr_2.1.2            nnet_7.3-13            cli_2.0.2             
    ## [37] survival_3.1-11        magrittr_1.5           crayon_1.3.4          
    ## [40] memoise_1.1.0          evaluate_0.14          fansi_0.4.1           
    ## [43] foreign_0.8-76         tools_3.6.3            data.table_1.12.8     
    ## [46] lifecycle_0.2.0        stringr_1.4.0          locfit_1.5-9.4        
    ## [49] munsell_0.5.0          cluster_2.1.0          AnnotationDbi_1.48.0  
    ## [52] compiler_3.6.3         rlang_0.4.5            grid_3.6.3            
    ## [55] RCurl_1.98-1.1         rstudioapi_0.11        htmlwidgets_1.5.1     
    ## [58] labeling_0.3           bitops_1.0-6           base64enc_0.1-3       
    ## [61] rmarkdown_2.1          gtable_0.3.0           DBI_1.1.0             
    ## [64] R6_2.4.1               gridExtra_2.3          knitr_1.28            
    ## [67] bit_1.1-15.2           Hmisc_4.4-0            stringi_1.4.6         
    ## [70] Rcpp_1.0.4             geneplotter_1.64.0     vctrs_0.2.4           
    ## [73] rpart_4.1-15           acepack_1.4.1          png_0.1-7             
    ## [76] tidyselect_1.0.0       xfun_0.12
