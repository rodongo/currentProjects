# ATAC-seq Data Analysis Tutorial
# install devtools if not previously installed
setwd('G:/Thesis Data/ATAC-Tutorial/')
if (!is.element('devtools', installed.packages()[,"Package"])) install.packages('devtools')

# install dependencies
devtools::install_github("demuellae/muLogR")
devtools::install_github("demuellae/muRtools")
devtools::install_github("demuellae/muReportR")

# install ChrAccR
devtools::install_github("GreenleafLab/ChrAccR", dependencies=TRUE)


# hg38 annotation package
install.packages("https://muellerf.s3.amazonaws.com/data/ChrAccR/data/annotation/ChrAccRAnnotationHg38_0.0.1.tar.gz")
library(ChrAccR)

# Example Data
devtools::install_github("GreenleafLab/ChrAccRex")

library(ggplot2)
# use a grid-less theme
theme_set(muRtools::theme_nogrid())

# download and unzip the dataset
datasetUrl <- "https://s3.amazonaws.com/muellerf/data/ChrAccR/data/tutorial/tcells.zip"
downFn <- "tcells.zip"
download.file(datasetUrl, downFn)
unzip(downFn, exdir=".")

# prepare the sample annotation table
sampleAnnotFn <- file.path("tcells", "samples.tsv")
bamDir <- file.path("tcells", "bam")
sampleAnnot <- read.table(sampleAnnotFn, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# add a column that ChrAccR can use to find the correct bam file for each sample
sampleAnnot[,"bamFilenameFull"] <- file.path(bamDir, sampleAnnot[,"bamFilename"])

setConfigElement("annotationColumns", c("cellType", "stimulus", "donor"))
setConfigElement("colorSchemes", c(
  getConfigElement("colorSchemes"), 
  list(
    "stimulus"=c("U"="royalblue", "S"="tomato"),
    "cellType"=c("TCD8naive"="purple3", "TCD8EM"="purple4", "TeffNaive"="seagreen3", "TeffMem"="seagreen4")
  )
))
setConfigElement("filteringCovgCount", 1L)
setConfigElement("filteringCovgReqSamples", 0.25)
setConfigElement("filteringSexChroms", TRUE)
setConfigElement("normalizationMethod", "quantile")
diffCompNames <- c(
  "TeffNaive_U vs TeffMem_U [sampleGroup]",
  "TCD8naive_U vs TCD8EM_U [sampleGroup]"
)
setConfigElement("differentialCompNames", diffCompNames)
# ChrAccR uses sets of predefined genomic regions to summarize count data. By default, this includes 
# 200bp-tiling windows. The following code generates a list of custom region sets that we can 
# optionally use in the analysis. Here, we use 500-bp tiling regions a consensus set of ATAC-seq peaks 
# that were identified in a large cancer dataset (Corces et al. 2018).
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("HDF5Array")
BiocManager::install("DelayedArray")
library(BSgenome.Hsapiens.UCSC.hg38)
library(HDF5Array)
library(DelayedArray)
# Download an ATAC peak set from a Pan-cancer dataset [doi:10.1126/science.aav1898]
cancerPeaks <- read.table("https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec", sep="\t", header=TRUE, comment.char="", quote="", stringsAsFactors=FALSE)
cancerPeaks.1 <- GenomicRanges::GRanges(cancerPeaks[,"seqnames"], IRanges::IRanges(start=cancerPeaks[,"start"], end=cancerPeaks[,"end"], names=cancerPeaks[,"name"]))
cancerPeaks.1 <- sort(cancerPeaks.1) # sort
# a list of custom region sets
regionSetList <- list(
  tiling500bp = muRtools::getTilingRegions("hg38", width=500L, onlyMainChrs=TRUE),
  cancerPeaks = cancerPeaks.1
)

run_atac("ChrAccR_analysis", "bamFilenameFull", sampleAnnot, genome="hg38", sampleIdCol="sampleId", regionSets=regionSetList)


dsa <- ChrAccRex::loadExample('dsAtac_ia_example')
