# This function extracts important features (column variables) from principal components
# following PCA.
plotLoadings3 <- function(
  pcaobj, # Provide a PCA object input
  components, # The components you wish to include
  rangeRetain = 0.05 # The percentage  to be considered
)
{
  x <- pcaobj$loadings[, components, drop = FALSE]
  retain <- c()
  for (i in seq_along(components)) {
    offset <- (max(x[, i]) - min(x[, i])) * rangeRetain
    uppercutoff <- max(x[, i]) - offset
    lowercutoff <- min(x[, i]) + offset
    retain <- unique(c(retain, which(x[, i] >= uppercutoff), 
                       which(x[, i] <= lowercutoff)))
  }
  x <- x[retain, , drop = FALSE]
  PC <- Loading <- NULL
  x <- data.frame(rownames(x), x[, components, drop = FALSE])
  colnames(x)[1] <- "var"
  plotobj <- melt(x, id = "var")
  colnames(plotobj) <- c("var", "PC", "Loading")
  plotobj <- unique(as.character(plotobj$var))
  return(plotobj)
}
