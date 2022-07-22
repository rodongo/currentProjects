setwd('G:/Projects/NIH/Input/Animal Models/GEM Simulations/Mouse  - Alzheimer/iMAT_Fisher/')
load('iMAT_Fisher.RData')
library(readxl)
library(openxlsx)
library(dplyr)
library(VarfromPDB)
library(RISmed)
library(stringi)

resDir <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/Analyzed/'
saveDir <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/Analyzed/Analysed B/'
# data1 <- 'GSE163857_E3E4_CON.xlsx'
# data1 <- 'GSE163857_E3E4_FAD.xlsx'

files <- list.files(path = paste0(resDir),pattern = '*.xlsx')
ftMetAlt <- list()
ftMetAlt1 <- list()
for(j in 1:length(files)){
  wb <- loadWorkbook(paste0(resDir,files[j]))
  sheetNames <- names(wb)
  samps <- as.matrix(table(readxl::read_xlsx(path = paste0(resDir,files[j]), sheet = sheetNames[2])$Condition))
  df <- readxl::read_xlsx(path = paste0(resDir,files[j]), sheet = sheetNames[1])
  ftMetAlt1[[paste0(gsub('.xlsx','',files[j],''))]] <- df
  dfRes <- data.frame()
  if(samps[[1]]<4 | samps[[2]]<4){
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==100 &
                                    as.numeric(df$`% number of rnxs active in AD`)==0),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==2/samps[2,]*100 &
                                    as.numeric(df$`% number of rnxs active in AD`)==0),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==1/samps[2,]*100 &
                                    as.numeric(df$`% number of rnxs active in AD`)==100),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==0 &
                                    as.numeric(df$`% number of rnxs active in AD`)==100),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==0.0 &
                                    as.numeric(df$`% number of rnxs active in AD`)==80),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==100 &
                                    as.numeric(df$`% number of rnxs active in AD`)==20),])
    dfRes <- rbind(dfRes,df[which(as.numeric(df$`% number of rnxs active in control`)==0 &
                                    as.numeric(df$`% number of rnxs active in AD`)==40),])
  }
  else{
    df$`Odds Ratio`[is.na(df$`Odds Ratio`)] <- 0
    df$`Odds Ratio`[df$`Odds Ratio`=='Inf'] <- 0
    dfRes <- df[order(df$`Fishers P-Value`, decreasing = F),] %>% head(100)
  }
  ftMetAlt[[paste0(gsub('.xlsx','',files[j],''))]] <- dfRes
  # ftMetAlt[[paste0(gsub('.xlsx','',files[j]))]] <- df[order(df$`Fishers P-Value`, decreasing = F),] %>% head(50)
}

# Add the leading genes for each reaction computed from the gene p-value in the MOOMIN weight function
library(plyr)
fn <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/leadingGenes.xlsx'
wb <- loadWorkbook(paste0(fn))
sheetNames <- names(wb)
ftMetAlt2 <- list()
ftMetAlt3 <- list()
for(j in 1:length(ftMetAlt)){
  fn2 <- names(ftMetAlt)[j]
  sheetNamesJ <- sheetNames[which(paste0(fn2)==sheetNames)]
  df <- readxl::read_xlsx(path = paste0(fn), sheet = paste0(sheetNamesJ))
  df2 <- ftMetAlt[[j]]
  df3 <- ftMetAlt1[[j]]
  ftMetAlt2[[paste0(fn2)]] <- join(df2,df,by='Reaction')
  ftMetAlt3[[paste0(fn2)]] <- join(df3,df,by='Reaction')
}

# Separate p-value based from ratio based
ftMetAltRatio <- list()
ftMetAltPVals <- list()
for(j in 1:length(ftMetAlt2)){
  df <- ftMetAlt2[[j]]
  if(nrow(df)!=100){
    ftMetAltRatio[[paste0(names(ftMetAlt2)[j])]] <- df
  }
  else{
    ftMetAltPVals[[paste0(names(ftMetAlt2)[j])]] <- df
  }
}
# Perform Metabolic Pathway Enrichment Based on Fisher's p-values
library(collapse)

sigLevel <- 0.5
MetEnrichRes <- list()
for(k in 1:length(ftMetAltPVals)){
  myDF <- ftMetAltPVals[[k]]
  paths <- unique(ftMetAltPVals[[k]]$Subsystem)
  myPaths <- data.frame(matrix(data = NA, nrow = length(paths), ncol = 5))
  for(i in 1:length(paths)){
    p <- paths[i]
    pathDF <- myDF[myDF$Subsystem %in% p,]
    if(pathDF$Reaction %in% myDF$Reaction){
      lG <- myDF[which(myDF$Reaction %in% pathDF$Reaction),]
      lG <- lG[order(lG$`Reaction Weight`,decreasing = T),] %>% head(1)
      lR <- lG$Reaction
      lG <- lG$LeadingGenes
      lE <- pathDF$`Metabolic Equation`[pathDF$Reaction==lR]
    }
    else{
      lG <- ''
      lR <- ''
      lE <- ''
    }
    ct <- data.frame(rbind(cbind(n11 = length(which(myDF[which(myDF$`Fishers P-Value`<=sigLevel),]$Subsystem %in% pathDF$Subsystem)),
                                 n12 = length(which(myDF[which(myDF$`Fishers P-Value`>sigLevel),]$Subsystem %in% pathDF$Subsystem))), 
                           cbind(n21 = length(which(myDF[which(myDF$`Fishers P-Value`<=sigLevel),]$Subsystem %!in% pathDF$Subsystem)), 
                                 n22 = length(which(myDF[which(myDF$`Fishers P-Value`>sigLevel),]$Subsystem %!in% pathDF$Subsystem)))))
    ct <- rbind(ct, apply(ct,2,sum))
    ct$Total <- apply(ct,1,sum)
    hyperPval <- phyper(ct[1,1], ct[1,3], ct[3,3]-ct[1,3], ct[3,1], lower.tail = TRUE, log.p = FALSE)
    # hyperPval <- mcnemar.test(ct[1,1], ct[1,3], ct[3,3]-ct[1,3], ct[3,1], lower.tail = TRUE, log.p = FALSE)
    myPaths[i,] <- c(paste(p),hyperPval,paste0(lR),paste0(lG),paste0(lE))
  }
  colnames(myPaths) <- c('Pathway','P-Value','Leading Reaction','Leading Gene','Metabolic Reaction Equation')
  MetEnrichRes[[paste0(names(ftMetAltPVals)[k])]] <- myPaths
}

top5 <- lapply(MetEnrichRes, function(x){x <- x[order(x$`P-Value`),] %>% head(10)})

top5 %>% writexl::write_xlsx(path = paste0(resDirB,'Top10_MetPaths.xlsx'))


MetEnrichRes.1 <- list()
for(i in 1:length(MetEnrichRes)){
  df <- MetEnrichRes[[i]][,c(1:2)]
  colnames(df)[2] <- paste0(names(MetEnrichRes)[i],'_PValue')
  MetEnrichRes.1[[paste0(names(MetEnrichRes)[i])]] <- df
}
library(plyr)
MetEnrichRes.1 <- join_all(MetEnrichRes.1, by = 'Pathway', type = 'full',match = 'all')
# Pathway enrichment based on ratios
ftMetAltEn <- list()
for(i in 1:length(ftMetAltRatio)){
  dataI <- ftMetAltRatio[[i]]
  dataI <- table(dataI$Subsystem) %>%
    as.data.frame() 
  dataI <- dataI[order(dataI$Freq, decreasing = T),]
  colnames(dataI) <- c('Subsystem','Number of affected reactions')
  ftMetAltEn[[paste0(names(ftMetAltRatio)[i])]] <- dataI %>% head(10)
}
ftMetAltEn %>% writexl::write_xlsx(path = paste0(resDirB,'Top10_MetPathsRatios.xlsx'))

ftMetAltEn <- list()
for(i in 1:length(ftMetAltPVals)){
  dataI <- ftMetAltPVals[[i]]
  dataI <- table(dataI$Subsystem) %>%
    as.data.frame() 
  dataI <- dataI[order(dataI$Freq, decreasing = T),]
  colnames(dataI) <- c('Subsystem','Number of affected reactions')
  ftMetAltEn[[paste0(names(ftMetAltPVals)[i])]] <- dataI %>% head(10)
}
ftMetAltEn %>% writexl::write_xlsx(path = paste0(resDirB,'Top10_MetPathsPV.xlsx'))
# myPaths <- myPaths[order(as.numeric(myPaths$`P-Value`),decreasing = F),]
# |
#   (`% number of rnxs active in control`==100 & `% number of rnxs active in AD`==11.11111)
wb <- loadWorkbook('MouseiMAT_FishersTest2.xlsx')
sheetNames <- names(wb)
ftMetAlt2 <- list()
for(i in 1:length(sheetNames)){
  df <- readxl::read_xlsx(path = 'MouseiMAT_FishersTest2.xlsx', sheet = sheetNames[i])
  df$OR[is.na(df$OR)] <- 0
  df$OR[df$OR=='Inf'] <- 0
  df <- df[which(df$`P-Value`<0.05),] 
  ftMetAlt2[[paste(sheetNames[i])]] <- df[order(df$`P-Value`, decreasing = F),]
}

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringdist)

plots <- list()
for(p in 1:length(dbList)){
  dataP <- dbList[[p]]
  dataP$Intersect <- dataP$intersection_size
  dataP$term_name <- factor(dataP$term_name, 
                            levels = dataP$term_name[order(dataP$Intersect)])
  plots[[paste(names(dbList)[p])]]<- ggplot(data = dataP, aes(x = term_name, y = Intersect, col = adjusted_p_value))+
    geom_point(aes(size = intersection_size))+
    coord_flip()+
    ylab('Size of Intersection')+
    xlab('Enriched Pathway')+
    theme_classic()
}

visE3E4 <- ftMetAlt

datDF1 <- data.frame()
for(i in 1:length(visE3E4)){
  datDF <- visE3E4[[i]]
  # for(k in 1:length(dbsource)){
  #   dfT <- df %>% filter(df$source == dbsource[k]) %>% head(10)
  #   datDF <- rbind(datDF,dfT)
  # }
  datDF$Group <- paste(names(visE3E4)[i])
  datDF <- datDF[order(datDF$`Fishers P-Value`, decreasing = F),]
  datDF1 <- rbind(datDF1,datDF %>% head(30))
}
colnames(datDF1)[7] <- 'PValue'
# datDF1$OR <- as.numeric(datDF1$OR)
datDF1$Group <- gsub('GSE163857_','APO',datDF1$Group)
datDF1$'-log10(P.Value)' <- -log10(as.numeric(datDF1$PValue))
library(ggplot2)
library(ggpubr)
png(filename = paste('Metabolic-Reactions-APOE3-4-FAD.png'),
    width = 4.1, height = 7, units = 'in', res = 600)
ggplot(data = datDF1, aes(x= Group, y=Reaction, size = `Odds Ratio`, color = `-log10(P.Value)`, group=Group)) + 
  geom_point(alpha = 0.8) + 
  ylab('Reactions')+
  xlab('Genetic Background')+
  ggtitle('GSE163857: APOE3/4 and AD')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1),
        text = element_text(size = 8))
dev.off()
p0Plots <- list()
p1Plots <- list()
for(j in 1:length(visE3E4)){
  df <- visE3E4[[j]]
  df <- df[order(df$`Fishers P-Value`, decreasing = F),] %>% head(30)
  df <- table(df$Subsystem) %>% as.data.frame()
  df$Freq <- round(df$Freq/sum(df$Freq)*100,digits = 1)
  df <- df %>% mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)
  df <- df[order(df$Freq, decreasing = T),] %>% head(10)
  colnames(df) <- c('Subsystem','Frequency','lab.ypos')
  theme_set(theme_pubclean())
  p0 <- ggplot(df, aes(x = "", y = Frequency, fill = Subsystem)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    geom_text(aes(y = lab.ypos, label = Frequency), color = "white")+
    coord_polar("y", start = 0)+
    ggpubr::fill_palette("jco")+
    theme_void()
  p1 <- ggpie(
    df, x = "Frequency", label = "Frequency",
    lab.pos = "out", lab.font = list(color = "white"), 
    fill = "Subsystem", color = "white",
    palette = "jco",
    legend = "right"
  )
  p0Plots[[names(visE3E4)[j]]] <- p0
  p1Plots[[names(visE3E4)[j]]] <- p1
}
for(k in 1:length(p1Plots)){
  png(filename = paste0(names(p1Plots)[k],'.png'),
      width = 6.5, height = 6.5, units = 'in', res = 600)
  print(p1Plots[[k]])
  dev.off()
}

# #####################################################################################################
# #####################################################################################################
# ML
# Load data
library(dplyr)
library(tidyr)
library(tidyverse)
fn <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/Analyzed/Analysed A/'

files <- list.files(path = paste0(fn),pattern = '*.xlsx')
# ftMetAltBR <- list()
# ftMetAltGB <- list()
# ftMetAltMN <- list()
ftMetAltBR_CON <- list()
ftMetAltBR_TRE <- list()
ftMetAltGB_TRE <- list()
ftMetAltGB_CON <- list()
for(i in 1:length(files)){
  # i <- 1
  wb <- loadWorkbook(paste0(fn,files[i]))
  sheetNames <- names(wb)
  df.con <- readxl::read_xlsx(path = paste0(fn,files[i]), sheet = paste0(sheetNames[1])) %>% 
    select(c(Var1_dataTable1,starts_with('data')))
  df.tre <- readxl::read_xlsx(path = paste0(fn,files[i]), sheet = paste0(sheetNames[2]))%>% 
    select(c(Var1_dataTable2,starts_with('data')))
  # df.McNem <- data.frame(matrix(data = NA, nrow = nrow(df.tre), ncol = 2))
  # for(j in 1:nrow(df.tre)){
  #   j <- 17
  #   Control <- as.matrix(df.con[j,-1])# c(length(which(df.con[j,-1]!=0)),length(which(df.con[j,-1]==0)))
  #   AD <- as.matrix(df.tre[j,-1]) #c(length(which(df.tre[j,-1]!=0)),length(which(df.tre[j,-1]==0)))
  #   mcNem <- mcnemar.test(table(Control,AD))
  #   df.McNem[j,] <- c(df.tre[j,1],mcNem$p.value)
  # }
  
  fnames <- sapply(files[i], function(x){strsplit(x,'.',fixed = T)[[1]][1]})[[1]]
  fnames1 <- sapply(files[i], function(x){strsplit(x,'_',fixed = T)[[1]][1]})[[1]]
  fnames2 <- sapply(files[i], function(x){strsplit(x,'_',fixed = T)[[1]][2]})[[1]]
  fnames3 <- sapply(files[i], function(x){strsplit(x,'_',fixed = T)[[1]][3]})[[1]]
  
  fnames2 <- gsub('.xlsx','',fnames2)
  fnames2 <- gsub('.xlsx','',fnames2)
  fnames3 <- gsub('.xlsx','',fnames3)
  
  colnames(df.con)[-1] <- paste0(fnames1,'_',fnames2,'_CON')
  colnames(df.tre)[-1] <- paste0(fnames1,'_',fnames2,'_TRE')
  
  colnames(df.con)[1] <- 'Reaction'
  colnames(df.tre)[1] <- 'Reaction'
  
  ftMetAltBR_TRE[[paste0(fnames)]] <- df.tre
  ftMetAltBR_CON[[paste0(fnames)]] <- df.con
  colnames(df.con)[-1] <- paste0(fnames1,'_',fnames3,'_CON')
  colnames(df.tre)[-1] <- paste0(fnames1,'_',fnames3,'_TRE')
  
  colnames(df.con)[1] <- 'Reaction'
  colnames(df.tre)[1] <- 'Reaction'
  
  ftMetAltGB_TRE[[paste0(fnames)]] <- df.tre
  ftMetAltGB_CON[[paste0(fnames)]] <- df.con
  # ftMetAltMN[[paste0(fnames)]] <- df.McNem
}
# Combine the datasets: 1-Genetic Background
library(plyr)
df1 <- ftMetAltGB_CON[[1]]
df2 <- ftMetAltGB_TRE[[1]]
for(i in 2:length(ftMetAltGB_CON)){
  df1 <- cbind(df1,ftMetAltGB_CON[[i]][,-1])
  df2 <- cbind(df2,ftMetAltGB_TRE[[i]][,-1])
}

df <- cbind(df1,df2[,-1])
df <- df[!duplicated(df$Reaction),]
df <- df[complete.cases(df$Reaction),]
# rxns <- unique(df$Reaction)
dff <- df %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Reaction')
labs <- sapply(colnames(dff), 
               function(x){strsplit(x,'_',fixed = T)[[1]][2]})#[[1]]
labs2 <- sapply(colnames(dff), 
                function(x){strsplit(x,'_',fixed = T)[[1]][3][1]})
labs2 <- gsub('\\.[0-9]','',labs2)
labs3 <- paste(labs,labs2,sep = '_')
labs3 <- gsub('[0-9]','',labs3)

# Combine the datasets: 2-Brain Region

df1.1 <- ftMetAltBR_CON[[1]]
df2.1 <- ftMetAltBR_TRE[[1]]
for(i in 2:length(ftMetAltBR_CON)){
  df1.1 <- cbind(df1.1,ftMetAltBR_CON[[i]][,-1])
  df2.1 <- cbind(df2.1,ftMetAltBR_TRE[[i]][,-1])
}

df.1 <- cbind(df1.1,df2.1[,-1])
df.1 <- df.1[!duplicated(df.1$Reaction),]
df.1 <- df.1[complete.cases(df.1$Reaction),]
# rxns <- unique(df$Reaction)
dff.1 <- df.1 %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Reaction')
labs.1 <- sapply(colnames(dff.1), 
                 function(x){strsplit(x,'_',fixed = T)[[1]][2]})#[[1]]
labs2.1 <- sapply(colnames(dff.1), 
                  function(x){strsplit(x,'_',fixed = T)[[1]][3][1]})
labs2.1 <- gsub('\\.[0-9]','',labs2.1)
labs3.1 <- paste(labs.1,labs2.1,sep = '_')
labs3.1 <- gsub('[0-9]','',labs3.1)
# ftMetAltGB_TRE2 <- join_all(ftMetAltGB_TRE, by = 'Reaction')
# ftMetAltGB_CON2 <- join_all(ftMetAltGB_CON, by = 'Reaction')
# 
# fnames <- sapply(files[1], function(x){strsplit(x,'_',fixed = T)[[1]][1]})[[1]]
# fnames2 <- sapply(files[1], function(x){strsplit(x,'_',fixed = T)[[1]][2]})[[1]]
# fnames3 <- sapply(files[1], function(x){strsplit(x,'_',fixed = T)[[1]][3]})[[1]]
# Load required libraries

library(glmnet)
library(dplyr)
library(doParallel)
library(parallel)
set.seed(5)

lassoFeatExt <- function(dataLoc, datLabs){
  
  # Partition
  inpData <- as.matrix(t(dataLoc))
  
  # Fit the LASSO model (Lasso: Alpha = 1)
  cv.lasso <- cv.glmnet(inpData, datLabs, family = 'binomial', 
                        alpha = 1, parallel = TRUE, 
                        standardize = TRUE, 
                        type.measure = 'class')
  df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
  df_coef <- data.frame(rxnName = row.names(df_coef), Coeff = as.numeric(df_coef[,1]))
  # See all contributing variables
  df_coef <- df_coef[order(df_coef$Coeff, decreasing = T),] 
  df_coef <- rbind(head(df_coef,30), tail(df_coef,30))
  
  a <- seq(0.1, 0.9, 0.05)
  search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(inpData, datLabs, family = 'binomial', nfold = 10, type.measure = 'class', paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  md3 <- glmnet(inpData, datLabs, family = 'binomial', 
                lambda = cv3$lambda.1se, alpha = cv3$alpha)
  
  
  df_coef.en <- coef(md3)
  df_coef.en <- data.frame(rxnName = row.names(df_coef.en), Coeff = as.numeric(df_coef.en[,1]))
  # See all contributing variables
  df_coef.en <- df_coef.en[order(df_coef.en$Coeff, decreasing = T),] 
  df_coef.en <- rbind(head(df_coef.en,30), tail(df_coef.en,30))
  
  res <- list(cv.lasso = cv.lasso,
              coefMat.lasso = df_coef,
              elasnet_Alpha = 1,
              elasnet_Lambda = cv.lasso$lambda.min,
              cv.elasnet = md3,
              coefMat.elasnet = df_coef.en,
              elasnet_Alpha = cv3$alpha,
              elasnet_Lambda=cv3$lambda.1se)
  return(res)
}


labelsClass <- as.integer(factor(labs3))
labelsClass2 <- sapply(labs3, 
                       function(x){strsplit(x,'_',fixed = T)[[1]][2][1]})
labelsClass2 <- as.integer(factor(labelsClass2))
mlRes <- lassoFeatExt(dff,labelsClass2)

# ElsaticNet Coefficients# 
coeffMat <- mlRes$coefMat.elasnet
metIn <- df.con[,c(1,6,7,8)]
colnames(metIn) <- c('Reaction','EC','Equation','Subsystem')

metIn2 <- merge(metIn,coeffMat, by.x = 'Reaction',by.y = 'rxnName')
metIn2 <- metIn2[which(metIn2$Subsystem!='Drug metabolism'),] 
metIn2 <- metIn2[order(metIn2$Coeff,decreasing = T),]



# 


cons <- df1.1[which(df1.1$Reaction==metIn2$Reaction[1]),]
tre <- df2.1[which(df2.1$Reaction==metIn2$Reaction[1]),]
df.1 <- cbind(df1.1,df2.1[,-1])




writexl::write_xlsx(metIn2,'iMAT-Fishers-tr-con_ML_features.xlsx')
# Genetic Background
labelsClass <- as.integer(factor(sapply(labs3, 
                                        function(x){strsplit(x,'_',fixed = T)[[1]][1][1]})))
MouseExpressionLabels <- read_excel('Tx Data/Processed/Expression/MouseExpressionLabels.xlsx')$Var1
Labs <- data.frame(DataSet = colnames(dff),
                   GeneticBackground = labs,
                   GeneticBackgroundLab = as.integer(factor(labs)),
                   BrainRegion = labs.1,
                   BrainRegionLab = as.integer(factor(labs.1))
)
actMat <- dff[rownames(dff) %in% coeffMat$rxnName,]
# Q. Which metabolic pathways/ metabolites are enriched in affected reactions?
# Met Pathways/Reactions
Human_GEM <- read_excel('G:/Thesis Data/Human-GEM-main/model/Human-GEM.xlsx')
paths.actMat <- table(Human_GEM$SUBSYSTEM[which(Human_GEM$ID %in% row.names(actMat))]) %>% as.data.frame()
mets.actMat <- Human_GEM[which(Human_GEM$ID %in% row.names(actMat)),] %>%
  select(ID,EQUATION,SUBSYSTEM)
mets.actMat$Metabolites <- sapply(mets.actMat$EQUATION, 
                                  function(x){strsplit(x,'=',fixed = T)[[1]][1][1]})
mets.actMat$Metabolites <- gsub('<','',mets.actMat$Metabolites)
mets.actMat$Metabolites <- gsub('\\+',',',mets.actMat$Metabolites)
affPaths <- unique(mets.actMat$SUBSYSTEM)
res.Data <- data.frame(matrix(data = NA, ncol = 4,nrow = length(affPaths)))
colnames(res.Data) <- c('Subsystem','Reactions','Metabolites','Enrichment')
for(i in 1:length(affPaths)){
  # i <- 1
  df <- mets.actMat %>% filter(.$SUBSYSTEM==affPaths[i])
  res.Data[i,] <- c(affPaths[i],
                    toString(df$ID),
                    toString(df$Metabolites),
                    nrow(df))
}
# View transport reactions
Human_GEM.compart <- read_excel('E:/Human-GEM.xlsx', 
                                sheet = 'COMPS')
mouseComps <- read_excel('G:/Projects/NIH/Input/Animal Models/Models/Mouse-GEM-main/model/mouseComps.xlsx')
Human_GEM.mets <- read_excel('E:/Human-GEM.xlsx', 
                             sheet = 'METS')

t.Reactions <- res.Data %>% filter(res.Data$Subsystem == 'Transport reactions')
comps <- strsplit(t.Reactions$Metabolites,',')[[1]]
for(j in 1:length(comps)){
  comps[j] <- stringr::str_extract(comps[j],'(?<=\\[).*(?=\\])')
}
comps <- comps[!is.na(comps)]
comps <- table(mouseComps$Names[match(comps,mouseComps$Abbreviations)]) %>%
  as.data.frame()
colnames(comps) <- c('Compartment','Frequency')
# Blood group
b.Reactions <- res.Data %>% filter(res.Data$Subsystem == 'Blood group biosynthesis')
b.comps <- strsplit(b.Reactions$Metabolites,',')[[1]]
for(j in 1:length(b.comps)){
  b.comps[j] <- stringr::str_extract(b.comps[j],'(?<=\\[).*(?=\\])')
}
b.comps <- b.comps[!is.na(b.comps)]
b.comps <- table(mouseComps$Names[match(b.comps,mouseComps$Abbreviations)]) %>%
  as.data.frame()
colnames(b.comps) <- c('Compartment','Frequency')
as_tibble(b.comps)

# Generate Pie Charts
library(scales)
png(filename = paste0('Results/metabolites-transport-compartment.png'),
    width = 6.5, height = 5.5, res = 600, units = 'in')
print(piePlot(comps,'Metabolites in Transport Reactions'))
dev.off()
# Exchange reactions
e.Reactions <- res.Data %>% filter(res.Data$Subsystem == 'Exchange/demand reactions')

comps2 <- strsplit(e.Reactions$Metabolites,',')[[1]]
for(j in 1:length(comps2)){
  comps2[j] <- stringr::str_extract(comps2[j],'(?<=\\[).*(?=\\])')
}
comps2 <- comps2[!is.na(comps2)]
comps2 <- table(mouseComps$Names[match(comps2,mouseComps$Abbreviations)]) %>%
  as.data.frame()
colnames(comps2) <- c('Compartment','Frequency')

metDet <- function(mets){
  
  for(i in 1:length(mets)){
    mets[i] <- stringr::str_remove(mets[i],'(?<=\\[).*(?=\\])')
    mets[i] <- gsub('\\[','',mets[i])
    mets[i] <- gsub('\\]','',mets[i])
  }
  mets <- table(mets) %>%
    as.data.frame()
  colnames(mets) <- c('Metabolite','Frequency')
  mets <- mets[as.character(mets$Metabolite)!=' ',]
  mets <- mets[order(mets$Frequency,decreasing = T),]
  return(mets)
}
comps <- strsplit(t.Reactions$Metabolites,',')[[1]]
comps.mets <- metDet(comps)

comps2 <- strsplit(e.Reactions$Metabolites,',')[[1]]

comps2.mets <- metDet(comps2) %>%
  as.data.frame() 
colnames(comps2.mets) <- c('Metabolite','Frequency')

for(i in 1:nrow(comps2.mets)){
  i <- 1
  metI <- as.character(comps2.mets$Metabolite[i][1])
  st <- contains(metI,ignore.case = T,Human_GEM$EQUATION)
  subSys <- Human_GEM[grepl(paste0(metI),Human_GEM$EQUATION,ignore.case = T)==T,]
}
comps2.mets$MetPathway <- paste0(toString(Human_GEM.mets$NAME))
############# Correlation Analysis################
library(WGCNA)

metExtractor <- function(mets){
  mets <- gsub('<','',mets)
  mets <- gsub('\\+',',',mets)
  res.mets <- data.frame(matrix(data = NA,ncol = 1,nrow = length(mets)))
  for(k in 1:length(mets)){
    mets.k <- mets[k]
    mets.k <- strsplit(mets.k,',')[[1]]
    # for(j in 1:length(mets.k)){
    #   mets.k[j] <- stringr::str_extract(mets.k[j],'(?<=\\[).*(?=\\])')
    # }
    # mets.k <- mets.k[!is.na(mets.k)]
    
    # mets.k <- strsplit(mets.k,',')[[1]]
    for(l in 1:length(mets.k)){
      mets.k[l] <- stringr::str_remove(mets.k[l],'(?<=\\[).*(?=\\])')
      mets.k[l] <- gsub('\\[','',mets.k[l])
      mets.k[l] <- gsub('\\]','',mets.k[l])
    }
    mets.k <- table(mets.k) %>%
      as.data.frame()
    # colnames(mets.k) <- c('Metabolite','Frequency')
    # mets.k <- mets.k[as.character(mets.k$Metabolite)!=' ',]
    res.mets[k,] <- paste0(toString(unique(mets.k[,1])))
  }
  return(res.mets)
}

actMat.cor <- cor(t(actMat),
                  Labs$GeneticBackgroundLab, 
                  use = 'p',
                  method = 'spearman') %>% 
  as.data.frame()
colnames(actMat.cor) <- 'Genetic.Background'
actMat.cor <- actMat.cor %>% mutate(Brain.Region = cor(t(actMat), 
                                                       Labs$BrainRegionLab, use = 'p',
                                                       method = 'spearman')) 
library(tidyverse)
actMat.cor.1 <- actMat.cor %>% mutate(SubSystem = Human_GEM$SUBSYSTEM[match(rownames(actMat.cor), Human_GEM$ID)]) %>%
  rownames_to_column(var = "ID")

metEqs <- sapply(mets.actMat$EQUATION, function(x){strsplit(x,'=',fixed = T)[[1]][1][1]})
metEqs.1 <- data.frame(mets.actMat,Metabolites2 = metExtractor(metEqs))
actMat.cor.2 <- merge(actMat.cor.1,metEqs.1, by = 'ID') [,c(1:7)]

write.xlsx(actMat.cor.2, 'Results/iMAT_Fisher_Overall_ML_Cor.xlsx')
actMat.lm <- data.frame(t(actMat))
actMat.lm$GenBack <- Labs$GeneticBackgroundLab

library(dplyr)
library(tidyr)
actMat.lm.1 <- actMat.lm %>% pivot_longer(cols = starts_with('M'))
# actMat.lm.1 <- actMat.lm[1,]
actMat.lm.1 <- lm(GenBack~ 0 + value,data = actMat.lm.1)
summary(actMat.lm.1)
coef(actMat.lm.1)
# KMeans
library(cluster)    # clustering algorithms
library(factoextra) 
actMat2 <- actMat
colnames(actMat2) <- colnames(df)[-1]
k2 <- kmeans(t(actMat), centers = 8, nstart = 25)
fviz_cluster(k2, data = t(actMat))
# # modelFeat <- readxl::read_excel(path = 'G:/desk2/TUSEB/iBrain/MODEL_21_04_20.xlsx')[,c(1:2)]
# regrFeat2 <- list()
# for(r in 1:length(regrFeat)){
#   dataR <- regrFeat[[r]]
#   rxnNams <- dataR$rxnName
#   # dfff <- gsub('^(n_)','', rxnNams)
#   dataR$Subsystem <- hcData$Subsystem[match(rxnNams, hcData$Reaction)]
#   dataR$Equation <- hcData$`Met Reaction Equation`[match(rxnNams, hcData$Reaction)]
#   dataR$EC <- hcData$`Enzyme EC`[match(rxnNams, hcData$Reaction)]
#   dataR <- dataR[complete.cases(dataR$Subsystem),]
#   df <- subset(dataR, Coeff!=0)
#   df <- cbind(df$rxnName,df$EC,df$Equation,df$Subsystem,df$Coeff) %>%as.data.frame()
#   
#   colnames(df) <- c('Reaction','EC','Equation','Subsystem','Coefficient')
#   regrFeat2[[names(regrFeat)[r]]] <- subset(df, Coefficient!=0) 
# }
# regrFeat2 %>% writexl::write_xlsx(path = 'regreFeatureSelection.xlsx')
# library(ggplot2)
# # # Basic piechart
# # ggplot(data, aes(x="", y=prop, fill=Var1)) +
# #   geom_bar(stat="identity", width=1, color="white") +
# #   coord_polar("y", start=0) +
# #   theme_void() + 
# #   theme(legend.position="none") +
# #   
# #   geom_text(aes(y = ypos, label = Var1), color = "white", size=6) +
# #   scale_fill_brewer(palette="RdYlBu")
# 
# for(g in 1:length(regrFeat2)){
#   g <- 4
#   tyh <- as.data.frame(table(regrFeat2[[g]][3]))
#   data <- tyh %>% 
#     arrange(desc(Var1)) %>%
#     mutate(prop = Freq / sum(tyh$Freq) *100) %>%
#     mutate(ypos = cumsum(prop)- 0.5*prop )
#   png(filename = paste0('G:/ParkinsonsDatasets/Fluxome ML/GSMM Results/',names(regrFeat2)[g], '.png'),
#       width = 6.5, height = 5.5, res = 600, units = 'in')
#   ggplot(data, aes (x="", y = Freq, fill = factor(Var1))) + 
#     geom_col(position = 'stack', width = 1) +
#     geom_text(aes(label = paste(round(Freq / sum(Freq) * 100, 1), "%"), x = 1.3),
#               position = position_stack(vjust = 0.5), size =4.2) +
#     theme_classic() +
#     theme(plot.title = element_text(hjust=0.5),
#           axis.line = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank()) +
#     labs(fill = 'Metabolic Pathway',
#          x = NULL,
#          y = NULL,
#          title = paste(names(regrFeat2)[g])) + 
#     coord_polar("y") +
#     scale_fill_brewer(palette='RdYlBu')
#   dev.off()
# }
# 
# regrFeat2 %>% writexl::write_xlsx(path = 'G:/ParkinsonsDatasets/Fluxome ML/GSMM Results/regreFeatures.xlsx')
save.image('iMAT_Fisher.RData')
