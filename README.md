# CNApp ##
CNApp is a user-friendly web tool that offers sample- and cohort-level association analyses, allowing a comprehensive and integrative exploration of CNAs with clinical and molecular variables. CNApp generates genome-wide profiles with tumor-purity correctiion, calculates CNA levels by computing broad, focal and global CNA scores (BCS, FCS and GCS, respectively), identifies CNAs associated with clinical features and molecular annotated variables and uses machine learning-based predictions to classify samples by using segmented data from either microarrays or next-generation sequencing.
CNApp provides a unique scenario to comprehensively analyze CNAs and integrate them with molecular and clinical features.

Functions of CNApp comprise three main sections: 1- Re-Seg & Score: re-segmentation, CNA scores computation, variable association and survival analysis, 2- Region profile: genome-wide CNA profiling, descriptive regions assessment and 3- Classifier model: machine learning classification model predictions.

Please give us credit and cite CNApp when you use it for your integrative CNA analysis:

CNApp: a web-based tool for integrative analysis of genomic copy number alterations in cancer
Sebastia Franch-Exposito, Laia Bassaganyas, Maria Vila-Casadesus, Eva Hernandez-Illan, Roger Esteban-Fabro, Marcos Diaz-Gay, Juan Jose Lozano, Antoni Castells, Josep M. Llovet, Sergi Castellvi-Bel, Jordi Camps
bioRxiv 479667; doi: https://doi.org/10.1101/479667 .

## Running the app ##

First __check dependencies needed to run CNApp__ by copying the following code in R (or RStudio):


```R

if(!require(shiny)) install.packages("shiny")
if(!require(shinyBS)) install.packages("shinyBS")
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(shinyWidgets)) install.packages("shinyWidgets")
if(!require(shinydashboard)) install.packages("shinydashboard")

if(!require(V8)) install.packages("V8")
#if any issue here like 'ERROR: configuration failed for package ‘curl’' or 'ERROR: configuration failed for package ‘V8’' follow printed instructions on your command-line window

if(!require(httr)) install.packages("httr")
#if any issue here like 'ERROR: configuration failed for package ‘openssl’' follow printed instructions on your command-line window

if(!require(plotly)) install.packages("plotly")
if(!require(randomcoloR)) install.packages("randomcoloR")
if(!require(heatmaply)) install.packages("heatmaply")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggsignif)) install.packages("ggsignif")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(randomForest)) install.packages("randomForest")
if(!require(doBy)) install.packages("doBy")
if(!require(parallel)) install.packages("parallel")
if(!require(caret)) install.packages("caret")
if(!require(XML)) install.packages("XML")
#if any issue here like 'ERROR: configuration failed for package ‘XML’' ty installing 'libxml2' into your OS

if (!require(devtools)) install.packages("devtools")

if(!require(limma)) {source("http://www.bioconductor.org/biocLite.R");biocLite("limma")}
if(!require(GenomicFeatures)) {source("http://www.bioconductor.org/biocLite.R");biocLite("GenomicFeatures")}
#if any issue here like 'ERROR: configuration failed for package ‘RMySQL’' ty installing 'libmysqlclient' into your OS
if(!require(GenomicAlignments)) {source("http://www.bioconductor.org/biocLite.R");biocLite("GenomicAlignments")}
if(!require(GenVisR)) {source("http://www.bioconductor.org/biocLite.R");biocLite("GenVisR")}

    ############## to use BiocManager ####################
    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #    install.packages("BiocManager")
    #BiocManager::install("limma")
    #BiocManager::install("GenomicFeatures")
    #BiocManager::install("GenomicAlignments")
    #BiocManager::install("GenVisR")
    ######################################################

if(!require(shinysky)) {library(devtools); devtools::install_github("AnalytixWare/ShinySky")}
#if(!require(GenVisR)) {library(devtools); devtools::install_github("griffithlab/GenVisR")}

if (!require(webshot)) install.packages("webshot")

if(!require(phantomjs)) {webshot::install_phantomjs()}

```
**If warning outputs are printed, try again!**
**(or try line by line...)**


There are many ways to download and run CNApp:

Run from GitHub repository:

```R
# Easiest way is to use runGitHub to run CNApp from GitHub:
library(shiny)
runGitHub("CNApp", "ait5")
```

By loading CNApp from compressed app url:

```R
# Run a tar or zip file directly
runUrl("https://github.com/ait5/CNApp/archive/master.tar.gz")
runUrl("https://github.com/ait5/CNApp/archive/master.zip")
```
By cloning or downloading repository:
```
# Using runApp(),  first clone the repository with git. If you have cloned it into
# ~/CNApp, first go to that directory, then use runApp().
setwd("~/CNApp")
runApp()
```

