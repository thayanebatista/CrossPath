# --- execute o arquivo por inteiro

# --- Instalacao das bibliotecas necess�rias

# --- instala��o do biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager",  update = F)
}
# --- instala��o do devtools para instala��o deo ggcorrplot
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools",  update = F)
  
}
# --- instala��o biblioteca KEGGREST
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  BiocManager::install("KEGGREST", update = F)
}
# --- instala��o tidyverse
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse", update = F)
}
# --- intala��o TCGA
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks", update = F)
}
# --- instala��o DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = F)
}
# --- instala��o dplyr
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", update = F)
}
# --- instala��o qdapregex
if (!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman", update = F)
} 
if (!requireNamespace("qdapRegex", quietly = TRUE)){
  pacman::p_load_gh("trinker/qdapRegex")
  install.packages("qdapRegex", update = F)
}
# --- instala��o ggcorrplot
if(!requireNamespace("ggcorrplot", quietly = TRUE)){
  devtools::install_github("kassambara/ggcorrplot")
  install.packages("ggcorrplot", update = F)
}
# --- instala��o do EnhancedVolcano
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)){
  BiocManager::install("EnhancedVolcano", update = F)
}
# --- instal�ao ggplot2   
if (!requireNamespace("ggplot2", quietly = TRUE)){
  install.packages("ggplot2", update = F)
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)){
  install.packages("RColorBrewer", update = F)
}

# --- inicializa��o das libs

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
