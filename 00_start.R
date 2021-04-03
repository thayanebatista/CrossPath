# --- execute o arquivo por inteiro

# --- Instalacao das bibliotecas necessárias

# --- instalação do biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager",  update = F)
}
# --- instalação do devtools para instalação deo ggcorrplot
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools",  update = F)
  
}
# --- instalação biblioteca KEGGREST
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  BiocManager::install("KEGGREST", update = F)
}
# --- instalação tidyverse
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse", update = F)
}
# --- intalação TCGA
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks", update = F)
}
# --- instalação DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = F)
}
# --- instalação dplyr
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", update = F)
}
# --- instalação qdapregex
if (!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman", update = F)
} 
if (!requireNamespace("qdapRegex", quietly = TRUE)){
  pacman::p_load_gh("trinker/qdapRegex")
  install.packages("qdapRegex", update = F)
}
# --- instalação ggcorrplot
if(!requireNamespace("ggcorrplot", quietly = TRUE)){
  devtools::install_github("kassambara/ggcorrplot")
  install.packages("ggcorrplot", update = F)
}
# --- instalação do EnhancedVolcano
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)){
  BiocManager::install("EnhancedVolcano", update = F)
}
# --- instalçao ggplot2   
if (!requireNamespace("ggplot2", quietly = TRUE)){
  install.packages("ggplot2", update = F)
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)){
  install.packages("RColorBrewer", update = F)
}

# --- inicialização das libs

library(KEGGREST)
library(tidyverse)
library(DESeq2)
library(TCGAbiolinks)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)
library(qdapRegex)
library(readr)

setwd("~/CrossPath/processing")
rstudioapi::navigateToFile("data-processing-01.R")
rstudioapi::navigateToFile("data-processing-02.R")
