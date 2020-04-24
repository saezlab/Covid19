SARS-CoV-2 dataset: Enrichment of CARNIVAL results
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
15/04/2020

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

The present script deal with the RNAseq data from the study *"SARS-CoV-2
launches* *a unique transcriptional signature from in vitro, ex vivo,
and in vivo systems"*

<https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1>

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>

It uses the nodes from the **CARNIVAL** output to run an enrichment
analysis. The enrichment is performed for the different situations in
which **CARNIVAL** was executed for the two condiitons under study:

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

## Reading input data for Enrichment Analysis

To perform the enrichment analysis, we need to read the following input
files:

  - Output from CARNIVAL: to obtain the significant genes and the
    background genes

  - Datasets from MSigDB: describing the pathways in which our
    significant genes are known to be involved in.

  - Differential expression analysis: To evaluate the activity of the
    genes involved in the different enriched pathaways.

We first load the required packages and we define some functions.

``` r
library(readr)
library(piano)
library(dplyr)
library(ggplot2)
library(omicToolsTest)

## Function to extract the nodes that appear in CARNIVAL network and the 
## background genes (all genes present in the prior knowledge network).
## It returns a list with two objects: the success and the background genes.
extractCARNIVALnodes <- function(CarnivalResults){

    CarnivalNetwork <- 
        as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE)
    
    colnames(CarnivalNetwork) <- c("source", "sign", "target", "Weight")

    ## We define the set of nodes interesting for our condition
    sucesses <- unique(c(gsub("_.*","",CarnivalNetwork$source), 
        gsub("_.*","",CarnivalNetwork$target)))

    CarnivalAttributes <- as.data.frame(CarnivalResults$nodesAttributes, 
        stringsAsFactors = FALSE)

    ## We define the background as all the genes in our prior knowledge network.
    bg <- unique(gsub("_.*","",CarnivalAttributes$Node))     
    
    return(list(sucesses = sucesses, bg= bg))
}

### Function to print a barplot with the enriched pathways.
BarplotEnrichment <- function(PathwaysSelect, Interesting_pathways){ 
    
    p <- ggplot(PathwaysSelect, aes(x = reorder(pathway, pvalue), 
            y = -log10(pvalue))) + 
        geom_bar(aes(fill = sign), stat = "identity") +
        scale_fill_gradient2(low = "darkblue", high = "indianred", 
            mid = "whitesmoke", midpoint = 0) + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
            colour = ifelse(levels(reorder(PathwaysSelect$pathway, 
                PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                "red", "grey40"),
            face = ifelse(levels(reorder(PathwaysSelect$pathway, 
                PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                "bold", "plain")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
        xlab("")
    return(p)
}    
```

### Reading and formatting CARNIVAL output

We read the CARNIVAL results generated in the previous script. We define
two different gene sets in order tor conduct the enrichment. The first
set contains the nodes that appear in the CARNIVAL output and are
therefore relevant in the context of our input transcriptomic data. The
second set contains all the genes in our prior knowledge network which
are used as the backgroud.

``` r
## NHBE 
CarnivalResultsNHBE_noinput <- readRDS("ResultsCARNIVAL/NHBE_noinput.rds")
NodesNHBE_noinput <- extractCARNIVALnodes(CarnivalResultsNHBE_noinput)    
CarnivalResultsNHBE_ACE2input <- readRDS("ResultsCARNIVAL/NHBE_ACE2input.rds")
NodesNHBE_ACE2input <- extractCARNIVALnodes(CarnivalResultsNHBE_ACE2input)  
CarnivalResultsNHBE_HostVirusinput <- 
   readRDS("ResultsCARNIVAL/NHBE_HostVirusinput.rds")
NodesNHBE_HostVirusinput <- 
    extractCARNIVALnodes(CarnivalResultsNHBE_HostVirusinput)  

## A549
CarnivalResultsA549_noinput <- readRDS("ResultsCARNIVAL/A549_noinput.rds")
NodesA549_noinput <- extractCARNIVALnodes(CarnivalResultsA549_noinput)    
CarnivalResultsA549_ACE2input <- readRDS("ResultsCARNIVAL/A549_ACE2input.rds")
NodesA549_ACE2input <- extractCARNIVALnodes(CarnivalResultsA549_ACE2input)  
CarnivalResultsA549_HostVirusinput <- 
    readRDS("ResultsCARNIVAL/A549_HostVirusinput.rds")
NodesA549_HostVirusinput <- 
    extractCARNIVALnodes(CarnivalResultsA549_HostVirusinput)  

## A549 transfected with ACE2
CarnivalResultsA549ACE2_noinput <- 
    readRDS("ResultsCARNIVAL/A549ACE2_noinput.rds")
NodesA549ACE2_noinput <- extractCARNIVALnodes(CarnivalResultsA549ACE2_noinput)    
CarnivalResultsA549ACE2_ACE2input <- 
    readRDS("ResultsCARNIVAL/A549ACE2_ACE2input.rds")
NodesA549ACE2_ACE2input <- 
    extractCARNIVALnodes(CarnivalResultsA549ACE2_ACE2input)  
CarnivalResultsA549ACE2_HostVirusinput <- 
    readRDS("ResultsCARNIVAL/A549ACE2_HostVirusinput.rds")
NodesA549ACE2_HostVirusinput <- 
    extractCARNIVALnodes(CarnivalResultsA549ACE2_HostVirusinput) 

## CALU-3
CarnivalResultsCALU3_noinput <- readRDS("ResultsCARNIVAL/CALU3_noinput.rds")
NodesCALU3_noinput <- extractCARNIVALnodes(CarnivalResultsCALU3_noinput)    
CarnivalResultsCALU3_ACE2input <- readRDS("ResultsCARNIVAL/CALU3_ACE2input.rds")
NodesCALU3_ACE2input <- extractCARNIVALnodes(CarnivalResultsCALU3_ACE2input)  
CarnivalResultsCALU3_HostVirusinput <- 
   readRDS("ResultsCARNIVAL/CALU3_HostVirusinput.rds")
NodesCALU3_HostVirusinput <- 
    extractCARNIVALnodes(CarnivalResultsCALU3_HostVirusinput)  
```

### Reading Pathway data sets from MSigDB

We downloaded from MSigDB <https://www.gsea-msigdb.org/> the following
dataset: c2.cp.v7.0.symbols.gmt. It contains several pathways from
different resources and the genes that are known to be involved in those
pathways.

``` r
pathways <- gmt_to_csv("RawData/c2.cp.v7.0.symbols.gmt")
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

### Reading and formatting statistic from DEG

We read the results from the differential expression analysis. The
statistic of the genes will be mapped later on in the different
significant pathways.

``` r
## Differential expression table
dds_NHBEvsCOV2 <- readRDS("IntermediateFiles/dds_results_NHBEvsCOV2.rds") %>%
    as.data.frame() %>% 
    select(stat)
dds_A549vsCOV2 <- readRDS("IntermediateFiles/dds_results_A549vsCOV2.rds") %>%
    as.data.frame() %>% 
    select(stat)
dds_A549ACE2vsCOV2 <- 
    readRDS("IntermediateFiles/dds_results_A549ACE2vsCOV2.rds") %>%
    as.data.frame() %>% 
    select(stat)
dds_CALU3vsCOV2 <- readRDS("IntermediateFiles/dds_results_CALU3vsCOV2.rds") %>%
    as.data.frame() %>% 
    select(stat)
```

## Performing Enrichment Analysis and plotting the Results

Using the **Piano** R package, we run a gene set analysis (GSA) based on
a list of significant genes (CARNIVAL nodes) and a gene set collection
(background). It uses Fisher’s exact test.

### CARNIVAL output with no perturbation

#### NHBE cell line

First for the NHBE cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_NHBE_noinput <- runGSAhyper(NodesNHBE_noinput$sucesses, 
    universe = NodesNHBE_noinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_NHBE_noinput <- as.data.frame(sig_pathways_NHBE_noinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_NHBE_noinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_NHBE_noinput), 
    function(x, kinases, pathways){
        return(mean(dds_NHBEvsCOV2[row.names(dds_NHBEvsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_NHBE_noinput <- 
    sig_pathways_df_NHBE_noinput[!is.nan(sig_pathways_df_NHBE_noinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.0001.

``` r
PathwaysSelect_NHBE_noinput <- sig_pathways_df_NHBE_noinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.0001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_NHBE_noinput <- c("PID_IFNG_PATHWAY",
                          "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                          "PID_IL6_7_PATHWAY")

p_NHBE_noinput <- BarplotEnrichment(PathwaysSelect_NHBE_noinput, 
    Interesting_pathways_NHBE_noinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

#### A549 cell line

Then, for the A549 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_A549_noinput <- runGSAhyper(NodesA549_noinput$sucesses, 
    universe = NodesA549_noinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549_noinput <- as.data.frame(sig_pathways_A549_noinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549_noinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549_noinput), 
    function(x, kinases, pathways){
        return(mean(dds_A549vsCOV2[row.names(dds_A549vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549_noinput <- 
    sig_pathways_df_A549_noinput[!is.nan(sig_pathways_df_A549_noinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549_noinput <- sig_pathways_df_A549_noinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549_noinput <- c("PID_IFNG_PATHWAY", 
    "BIOCARTA_IFNA_PATHWAY", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")

p_A549_noinput <- BarplotEnrichment(PathwaysSelect_A549_noinput, 
    Interesting_pathways_A549_noinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

#### A549 cell line transfected with ACE2

Then, for the A549 cell line transfected with ACE2

``` r
## We run GSA hyper Geometric test
sig_pathways_A549ACE2_noinput <- runGSAhyper(NodesA549ACE2_noinput$sucesses, 
    universe = NodesA549ACE2_noinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549ACE2_noinput <- as.data.frame(sig_pathways_A549ACE2_noinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549ACE2_noinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549ACE2_noinput), 
    function(x, kinases, pathways){
        return(mean(dds_A549ACE2vsCOV2[row.names(dds_A549ACE2vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549ACE2_noinput <- 
    sig_pathways_df_A549ACE2_noinput[!is.nan(sig_pathways_df_A549ACE2_noinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549ACE2_noinput <- sig_pathways_df_A549ACE2_noinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.05) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549ACE2_noinput <- c("")

p_A549ACE2_noinput <- BarplotEnrichment(PathwaysSelect_A549ACE2_noinput, 
    Interesting_pathways_A549ACE2_noinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

#### CALU-3 cell line

Now for the CALU-3 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_CALU3_noinput <- runGSAhyper(NodesCALU3_noinput$sucesses, 
    universe = NodesCALU3_noinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_CALU3_noinput <- as.data.frame(sig_pathways_CALU3_noinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_CALU3_noinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_CALU3_noinput), 
    function(x, kinases, pathways){
        return(mean(dds_CALU3vsCOV2[row.names(dds_CALU3vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_CALU3_noinput <- 
    sig_pathways_df_CALU3_noinput[!is.nan(sig_pathways_df_CALU3_noinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.0001.

``` r
PathwaysSelect_CALU3_noinput <- sig_pathways_df_CALU3_noinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.0001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_CALU3_noinput <- c("PID_IL6_7_PATHWAY", 
    "PID_IFNG_PATHWAY", "PID_IL4_2PATHWAY", 
    "REACTOME_INTERLEUKIN_17_SIGNALING")

p_CALU3_noinput <- BarplotEnrichment(PathwaysSelect_CALU3_noinput, 
    Interesting_pathways_CALU3_noinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### CARNIVAL output with ACE2 as perturbation

#### NHBE cell line

First for the NHBE cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_NHBE_ACE2input <- runGSAhyper(NodesNHBE_ACE2input$sucesses, 
    universe = NodesNHBE_ACE2input$bg, gsc = loadGSC(pathways))
sig_pathways_df_NHBE_ACE2input <- 
    as.data.frame(sig_pathways_NHBE_ACE2input$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_NHBE_ACE2input$sign <- 
    unlist(lapply(row.names(sig_pathways_df_NHBE_ACE2input), 
    function(x, kinases, pathways){
        return(mean(dds_NHBEvsCOV2[row.names(dds_NHBEvsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_NHBE_ACE2input <- 
    sig_pathways_df_NHBE_ACE2input[!is.nan(sig_pathways_df_NHBE_ACE2input$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.0001.

``` r
PathwaysSelect_NHBE_ACE2input <- sig_pathways_df_NHBE_ACE2input %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.0001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_NHBE_ACE2input <- c("PID_IL6_7_PATHWAY", 
    "PID_IFNG_PATHWAY", "PID_IL12_2PATHWAY", 
    "REACTOME_INTERLEUKIN_1_PROCESSING")

p_NHBE_ACE2input <- BarplotEnrichment(PathwaysSelect_NHBE_ACE2input, 
    Interesting_pathways_NHBE_ACE2input)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

#### A549 cell line

Then, for the A549 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_A549_ACE2input <- runGSAhyper(NodesA549_ACE2input$sucesses, 
    universe = NodesA549_ACE2input$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549_ACE2input <- as.data.frame(sig_pathways_A549_ACE2input$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549_ACE2input$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549_ACE2input), 
    function(x, kinases, pathways){
        return(mean(dds_A549vsCOV2[row.names(dds_A549vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549_ACE2input <- 
    sig_pathways_df_A549_ACE2input[!is.nan(sig_pathways_df_A549_ACE2input$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549_ACE2input <- sig_pathways_df_A549_ACE2input %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549_ACE2input <- c("PID_IFNG_PATHWAY", 
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
    "REACTOME_INTERFERON_GAMMA_SIGNALING")

p_A549_ACE2input <- BarplotEnrichment(PathwaysSelect_A549_ACE2input, 
    Interesting_pathways_A549_ACE2input)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

#### A549 cell line transfected with ACE2

Then, for the A549 cell line transfected with ACE2

``` r
## We run GSA hyper Geometric test
sig_pathways_A549ACE2_ACE2input <- runGSAhyper(NodesA549ACE2_ACE2input$sucesses, 
    universe = NodesA549ACE2_ACE2input$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549ACE2_ACE2input <- 
    as.data.frame(sig_pathways_A549ACE2_ACE2input$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549ACE2_ACE2input$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549ACE2_ACE2input), 
    function(x, kinases, pathways){
        return(mean(dds_A549ACE2vsCOV2[row.names(dds_A549ACE2vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549ACE2_ACE2input <- 
    sig_pathways_df_A549ACE2_ACE2input[!is.nan(sig_pathways_df_A549ACE2_ACE2input$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549ACE2_ACE2input <- sig_pathways_df_A549ACE2_ACE2input %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.0001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549ACE2_ACE2input <- c("BIOCARTA_CXCR4_PATHWAY")

p_A549ACE2_ACE2input <- BarplotEnrichment(PathwaysSelect_A549ACE2_ACE2input, 
    Interesting_pathways_A549ACE2_ACE2input)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

#### CALU-3 cell line

First for the CALU-3 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_CALU3_ACE2input <- runGSAhyper(NodesCALU3_ACE2input$sucesses, 
    universe = NodesCALU3_ACE2input$bg, gsc = loadGSC(pathways))
sig_pathways_df_CALU3_ACE2input <- 
    as.data.frame(sig_pathways_CALU3_ACE2input$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_CALU3_ACE2input$sign <- 
    unlist(lapply(row.names(sig_pathways_df_CALU3_ACE2input), 
    function(x, kinases, pathways){
        return(mean(dds_CALU3vsCOV2[row.names(dds_CALU3vsCOV2) %in% pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_CALU3_ACE2input <- 
    sig_pathways_df_CALU3_ACE2input[!is.nan(sig_pathways_df_CALU3_ACE2input$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.0001.

``` r
PathwaysSelect_CALU3_ACE2input <- sig_pathways_df_CALU3_ACE2input %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_CALU3_ACE2input <- c("PID_IFNG_PATHWAY")

p_CALU3_ACE2input <- BarplotEnrichment(PathwaysSelect_CALU3_ACE2input, 
    Interesting_pathways_CALU3_ACE2input)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### CARNIVAL output with host-virus interaction as Perturbation

#### NHBE cell line

First for the NHBE cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_NHBE_HostVirusinput <- 
    runGSAhyper(NodesNHBE_HostVirusinput$sucesses, 
    universe = NodesNHBE_HostVirusinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_NHBE_HostVirusinput <- 
    as.data.frame(sig_pathways_NHBE_HostVirusinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_NHBE_HostVirusinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_NHBE_HostVirusinput), 
    function(x, kinases, pathways){
        return(mean(dds_NHBEvsCOV2[row.names(dds_NHBEvsCOV2) %in% 
            pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_NHBE_HostVirusinput <- 
    sig_pathways_df_NHBE_HostVirusinput[!is.nan(sig_pathways_df_NHBE_HostVirusinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.00001.

``` r
PathwaysSelect_NHBE_HostVirusinput <- sig_pathways_df_NHBE_HostVirusinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.00001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_NHBE_HostVirusinput <- c("PID_IL2_1PATHWAY", 
    "PID_IL12_2PATHWAY","PID_IL6_7_PATHWAY","PID_IFNG_PATHWAY", 
    "BIOCARTA_IL6_PATHWAY")
    
p_NHBE_HostVirusinput <- BarplotEnrichment(PathwaysSelect_NHBE_HostVirusinput, 
    Interesting_pathways_NHBE_HostVirusinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

#### A549 cell line

Then, for the A549 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_A549_HostVirusinput <- 
    runGSAhyper(NodesA549_HostVirusinput$sucesses, 
    universe = NodesA549_HostVirusinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549_HostVirusinput <- as.data.frame(sig_pathways_A549_HostVirusinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549_HostVirusinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549_HostVirusinput), 
    function(x, kinases, pathways){
        return(mean(dds_A549vsCOV2[row.names(dds_A549vsCOV2) %in% 
            pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549_HostVirusinput <- 
    sig_pathways_df_A549_HostVirusinput[!is.nan(sig_pathways_df_A549_HostVirusinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549_HostVirusinput <- sig_pathways_df_A549_HostVirusinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549_HostVirusinput <- c("PID_IFNG_PATHWAY")
    
p_A549_HostVirusinput <- BarplotEnrichment(PathwaysSelect_A549_HostVirusinput, 
    Interesting_pathways_A549_HostVirusinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

#### A549 cell line transfected with ACE2

Then, for the A549 cell line transfected with ACE2

``` r
## We run GSA hyper Geometric test
sig_pathways_A549ACE2_HostVirusinput <- 
    runGSAhyper(NodesA549ACE2_HostVirusinput$sucesses, 
    universe = NodesA549ACE2_HostVirusinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_A549ACE2_HostVirusinput <- 
    as.data.frame(sig_pathways_A549ACE2_HostVirusinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_A549ACE2_HostVirusinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_A549ACE2_HostVirusinput), 
    function(x, kinases, pathways){
        return(mean(dds_A549ACE2vsCOV2[row.names(dds_A549ACE2vsCOV2) %in% 
            pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_A549ACE2_HostVirusinput <- 
    sig_pathways_df_A549ACE2_HostVirusinput[!is.nan(sig_pathways_df_A549ACE2_HostVirusinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.001.

``` r
PathwaysSelect_A549ACE2_HostVirusinput <- sig_pathways_df_A549ACE2_HostVirusinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_A549ACE2_HostVirusinput <- c("BIOCARTA_CXCR4_PATHWAY",
    "PID_IL2_1PATHWAY", "PID_IL23_PATHWAY")
    
p_A549ACE2_HostVirusinput <- BarplotEnrichment(PathwaysSelect_A549ACE2_HostVirusinput, 
    Interesting_pathways_A549ACE2_HostVirusinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

#### CALU-3 cell line

Finally for the CALU-3 cell line

``` r
## We run GSA hyper Geometric test
sig_pathways_CALU3_HostVirusinput <- 
    runGSAhyper(NodesCALU3_HostVirusinput$sucesses, 
    universe = NodesCALU3_HostVirusinput$bg, gsc = loadGSC(pathways))
sig_pathways_df_CALU3_HostVirusinput <- 
    as.data.frame(sig_pathways_CALU3_HostVirusinput$resTab)

## We map the t-stastic into the resulted enriched pathways.
sig_pathways_df_CALU3_HostVirusinput$sign <- 
    unlist(lapply(row.names(sig_pathways_df_CALU3_HostVirusinput), 
    function(x, kinases, pathways){
        return(mean(dds_CALU3vsCOV2[row.names(dds_CALU3vsCOV2) %in% 
            pathways[pathways$term == x,1],1], na.rm = TRUE))
    },kinases = kinases, pathways = pathways))

sig_pathways_df_CALU3_HostVirusinput <- 
    sig_pathways_df_CALU3_HostVirusinput[!is.nan(sig_pathways_df_CALU3_HostVirusinput$sign),]
```

We format the results and we prepare them to be plotted. For
visualization purposes, we just select pathways with adjusted p-values
lower than 0.00001.

``` r
PathwaysSelect_CALU3_HostVirusinput <- sig_pathways_df_CALU3_HostVirusinput %>%
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, sign) %>%
    dplyr::filter(`Adjusted p-value` <= 0.05) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

We finally plot the results highlighting the most relevant pathways.

``` r
Interesting_pathways_CALU3_HostVirusinput <- 
    c("SIF_IL4RECEPTOR_IN_B_LYPHOCYTES", "KEGG_TTGF_BETA_SIGNALING_PATHWAY")
    
p_CALU3_HostVirusinput <- BarplotEnrichment(PathwaysSelect_CALU3_HostVirusinput, 
    Interesting_pathways_CALU3_HostVirusinput)
```

![](CarnivalEnrichment_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

## R session Info

``` r
sessionInfo()
```

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
    ## [11] BiocGenerics_0.32.0         omicToolsTest_0.1.0        
    ## [13] ggplot2_3.3.0               dplyr_0.8.5                
    ## [15] piano_2.2.0                 readr_1.3.1                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fgsea_1.12.0           colorspace_1.4-1       ellipsis_0.3.0        
    ##   [4] htmlTable_1.13.3       XVector_0.26.0         base64enc_0.1-3       
    ##   [7] rstudioapi_0.11        farver_2.0.3           ggrepel_0.8.2         
    ##  [10] DT_0.13                bit64_0.9-7            AnnotationDbi_1.48.0  
    ##  [13] fansi_0.4.1            splines_3.6.3          geneplotter_1.64.0    
    ##  [16] knitr_1.28             Formula_1.2-3          jsonlite_1.6.1        
    ##  [19] annotate_1.64.0        cluster_2.1.0          dbplyr_1.4.2          
    ##  [22] png_0.1-7              pheatmap_1.0.12        shinydashboard_0.7.1  
    ##  [25] snowfall_1.84-6.1      graph_1.64.0           shiny_1.4.0.2         
    ##  [28] compiler_3.6.3         httr_1.4.1             backports_1.1.5       
    ##  [31] assertthat_0.2.1       Matrix_1.2-18          fastmap_1.0.1         
    ##  [34] limma_3.42.0           cli_2.0.2              later_1.0.0           
    ##  [37] acepack_1.4.1          visNetwork_2.0.9       htmltools_0.4.0       
    ##  [40] tools_3.6.3            igraph_1.2.5           gtable_0.3.0          
    ##  [43] glue_1.4.0             GenomeInfoDbData_1.2.2 rappdirs_0.3.1        
    ##  [46] fastmatch_1.1-0        Rcpp_1.0.4             slam_0.1-47           
    ##  [49] vctrs_0.2.4            gdata_2.18.0           xfun_0.12             
    ##  [52] stringr_1.4.0          mime_0.9               lifecycle_0.2.0       
    ##  [55] gtools_3.8.2           XML_3.99-0.3           zlibbioc_1.32.0       
    ##  [58] scales_1.1.0           hms_0.5.3              promises_1.1.0        
    ##  [61] relations_0.6-9        RColorBrewer_1.1-2     sets_1.0-18           
    ##  [64] yaml_2.2.1             curl_4.3               memoise_1.1.0         
    ##  [67] gridExtra_2.3          rpart_4.1-15           latticeExtra_0.6-29   
    ##  [70] reshape_0.8.8          stringi_1.4.6          RSQLite_2.2.0         
    ##  [73] genefilter_1.68.0      checkmate_2.0.0        caTools_1.18.0        
    ##  [76] rlang_0.4.5            pkgconfig_2.0.3        bitops_1.0-6          
    ##  [79] evaluate_0.14          lattice_0.20-41        purrr_0.3.3           
    ##  [82] labeling_0.3           htmlwidgets_1.5.1      cowplot_1.0.0         
    ##  [85] bit_1.1-15.2           tidyselect_1.0.0       GSEABase_1.48.0       
    ##  [88] plyr_1.8.6             magrittr_1.5           R6_2.4.1              
    ##  [91] gplots_3.0.3           UniProt.ws_2.26.0      Hmisc_4.4-0           
    ##  [94] DBI_1.1.0              foreign_0.8-76         pillar_1.4.3          
    ##  [97] withr_2.1.2            nnet_7.3-13            survival_3.1-11       
    ## [100] RCurl_1.98-1.1         tibble_3.0.0           crayon_1.3.4          
    ## [103] KernSmooth_2.23-16     BiocFileCache_1.10.2   rmarkdown_2.1         
    ## [106] jpeg_0.1-8.1           locfit_1.5-9.4         grid_3.6.3            
    ## [109] data.table_1.12.8      marray_1.64.0          blob_1.2.1            
    ## [112] digest_0.6.25          xtable_1.8-4           httpuv_1.5.2          
    ## [115] munsell_0.5.0          shinyjs_1.1
