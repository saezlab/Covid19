## Applying Saezlab Tools to Covid-19 related Datasets

This repository contains the scripts used to apply some of our tools to Covid-19
related datasets. In particular, we take the RNAseq data from the study: 

*"SARS-CoV-2 launches a unique transcriptional signature from in vitro, ex vivo, and in vivo systems"* 

<https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1>

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>

We focused in the experiments containing these conditions:

+ Human lung epithelial cells (**NHBE**): mock treated vs infected with 
SARS-CoV-2. 

+ **A549** alveolar cancer cell line: mock treated vs infected with SARS-CoV-2.   

We detail below the different scripts and analysis performed:

+ We first applied a traditional differential expression analysis using the 
**DESeq2** R package. This is aligned with the script of Agatha Treveil from 
Tamas' group. 

+ Then, we used the normalised counts and the stastics generated in the previous 
script to run **Dorothea** and **Progeny**. Doing so, we estimated Transcription
factors activity and pathway activity in the SARS-CoV-2 infected lines in 
comparison with the mock treated. 

+ We finally run **CARNIVAL** using different perturbations conditions related 
to the virus action. In order to do so, we used TF and Pathway activity scores 
generated in the previous script. CARNIVAL also requires a prior knowledge 
network that was extracted from **Omnipath** using **OmnipathR**. 

### License Info

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Please check http://www.gnu.org/licenses/.