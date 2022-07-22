
install.packages("gtExtras")
BiocManager::install('tidyverse')
# or if wanting the dev version
# if needed install.packages("remotes")
# remotes::install_github("jthomasmock/gtExtras")

# Tau and Amyloid Associated genes
library(readxl)
library(dplyr)
library(tidyverse)
library(gt)
library(gtExtras)
X43587_2021_138_MOESM3_ESM <- read_excel("G:/Projects/Common/43587_2021_138_MOESM3_ESM.xlsx", 
                                         sheet = "Table S13")
AD.SeedGenes <- X43587_2021_138_MOESM3_ESM$`AD seed genes`[2:length(X43587_2021_138_MOESM3_ESM$`AD seed genes`)]

TAU.SeedGenes <- X43587_2021_138_MOESM3_ESM$`Tauopathy seed genes`[2:length(X43587_2021_138_MOESM3_ESM$`Tauopathy seed genes`)]
AMYLOID.SeedGenes <- X43587_2021_138_MOESM3_ESM$`Amyloid seed genes`[2:length(X43587_2021_138_MOESM3_ESM$`Amyloid seed genes`)]
AMYLOID.SeedGenes <- AMYLOID.SeedGenes[AMYLOID.SeedGenes!='NA']
AMYLOID.SeedGenes <- AMYLOID.SeedGenes[grepl('NA',as.character(AMYLOID.SeedGenes))==F]
library(stringr)
AMYLOID.SeedGenes <- stringr::str_replace_na(AMYLOID.SeedGenes)
AMYLOID.SeedGenes <- str_remove(string = AMYLOID.SeedGenes, pattern = '^NA')
AMYLOID.SeedGenes <- stringr::str_remove_all(AMYLOID.SeedGenes,'')

library(stringi)
AMYLOID.SeedGenes <- stringi::stri_remove_empty(AMYLOID.SeedGenes)
AMYLOID.SeedGenes <- stringi::stri_remove_empty(AMYLOID.SeedGenes)
# basic use
df <- X43587_2021_138_MOESM3_ESM[-1,-c(3:5,12:19)][,c(3:5)]
df <- df[complete.cases(df),]
colnames(df) <- c('Symbol','Entrez','Name')
df %>%
  # head() %>%
  gt() 
