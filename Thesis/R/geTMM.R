geTMM <- function(
  countsMat,#Raw counts data matrix
  sampleGrps,#A factor object with sample-group information
  geneLength,#Length of each gene in the raw count matrix
  preNorm
){
  require(edgeR)
  countsMat <- as.matrix(countsMat)
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  geneLength <- geneLength[keep]
  
  # if not rpk, convert and normalize with gene length
  if(preNorm=='FPKM'){
    rpk <- countsMat
  }
  else{
    rpk <- countsMat/(geneLength/10^3)
  }
  # Create a DGEList object using edgeR
  DGE <- DGEList(counts=rpk, group=sampleGrps)
  # normalize for library size by cacluating scaling factor using TMM (default method)
  DGE <- calcNormFactors(DGE, method = 'TMM')
  # count per million read (normalized count)
  normCounts <- cpm(DGE)
  return(as.data.frame(normCounts))
}