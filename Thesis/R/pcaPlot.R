pcaPlot <- function(dfnew,Key,title){
  require(ggplot2)
  as.matrix(dfnew)
  pcaANAL <- prcomp(t(dfnew), scale. = T)
  pcaPlotDF <- data.frame(pcaANAL$x[,1:2]) 
  
  ggplot(data = pcaPlotDF, aes(x = PC1, y = PC2, label = row.names(pcaPlotDF),
                               color = Key))+
    geom_point()+
    geom_text(aes(label = row.names(pcaPlotDF)),hjust=0.6, vjust=0, size = 5)+
    labs(title = paste0(title))+
    theme_bw()
}