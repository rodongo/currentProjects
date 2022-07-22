setwd('G:/Thesis Data/OMIC Data/Priori/')
load('BrendaDb.RData')
BiocManager::install('brendaDb', dependencies=TRUE)
library(brendaDb)
library(dplyr)
library(readxl)
Human_GEM <- read_excel('G:/Thesis Data/Human-GEM-main/model/Human-GEM.xlsx', sheet = 4)
ec_code_data_hGEM <- read_excel('G:/Thesis Data/OMIC Data/Expression/AD/Processed/ec_code_data_hGEM.xlsx')

# Li et al 2022 - Aging Cell GPMM
acel13595_sup_0002_tables1_s7 <- read_excel('C:/Users/ROdongo/Downloads/acel13595-sup-0002-tables1-s7.xlsx', 
                                            sheet = 'Table S1', skip = 1)
ec_code_data_hGEM$KCat <- acel13595_sup_0002_tables1_s7$Kcat[match(ec_code_data_hGEM$ecData_4,acel13595_sup_0002_tables1_s7$EC_number)]

ec_code_data_hGEM <- ec_code_data_hGEM[!is.na(ec_code_data_hGEM$ecData_4),]

ec_code_data_hGEM.1 <- ec_code_data_hGEM[is.na(ec_code_data_hGEM$KCat),]
brenda.filepath <- DownloadBrenda()
df <- ReadBrenda(brenda.filepath)
saveRDS(df,'BrendaData.RDS')
# df <- df %>% filter(df$field=='SPECIFIC_ACTIVITY')

humEnzKcats <- c()
humEnzEC <- list()
ecs <- c()
for(i in 1338:nrow(ec_code_data_hGEM.1)){
  # i <- 2
  ec.i <- unlist(strsplit(ec_code_data_hGEM.1$ecData_4[i],';')[[1]])
  res <- QueryBrenda(df, EC = paste(ec.i),  organisms = "Homo sapiens", n.core = 3)
  humEnzEC[[paste0(i)]] <- ExtractField(res, field = 'parameters$turnover.number') %>% 
    filter(.$organism=='Homo sapiens')
  # if(nrow(kcat.vals)>0){
  #   kcat.vals <- kcat.vals %>%
  #     mutate(description = as.numeric(.$description))
  #   kcat.vals <- kcat.vals[!is.na(kcat.vals$description),] 
  #   kcat.vals <- kcat.vals %>% 
  #     select(ec,description) %>%
  #     group_by(ec) %>%
  #     summarise(KCal = median(description))
  #   KCal <- kcat.vals$KCal
  #   ec <- kcat.vals$ec
  # }
  # else{
  #   kcat.vals.1 <- c(kcat.vals.1,0)
  #   ec <- ec.i
  # }
  # humEnzKcats <- c(humEnzKcats,paste(KCal))
  # ecs <- c(ecs,paste(ec.i))
}
kcat.vals.1 <- c()
for(j in 1:length(humEnzEC)){
  # j <- 1
  dat <- humEnzEC[[j]]
  if(nrow(dat)>0){
    kcat.vals <- dat %>%
      mutate(description = as.numeric(.$description))
    kcat.vals <- kcat.vals[!is.na(kcat.vals$description),]
    kcat.vals <- kcat.vals %>%
      select(ec,description) %>%
      group_by(ec) %>%
      summarise(KCal = median(description))
    # KCal <- kcat.vals$KCal
    # ec <- kcat.vals$ec
    kcat.vals.1 <- rbind(kcat.vals.1,kcat.vals)
  }
}
ecNums2 <- kcat.vals.1 %>% as.data.frame()


for(l in 1:nrow(ec_code_data_hGEM.1)){
  # l <- 1
  ecs <- ec_code_data_hGEM.1[l,]
  ecs <- strsplit(ecs$ecData_4,split = '\\;')[[1]]
  ecs.kcat <- ecNums2 %>% filter(ec %in% ecs)
  ecs.kcat <- ecs.kcat[!duplicated(ecs.kcat$ec),]
  ec_code_data_hGEM.1$KCat[l] <- paste0(toString(ecs.kcat$KCal))
}
geneLength <- readRDS(file = 'G:/Thesis Data/OMIC Data/Priori/geneLengthDF.RDS')
Human_GEM.Genes <- geneLength$SYMBOL[match(Human_GEM$NAME,geneLength$ENSEMBLID)]  
brenda2data <- ID2Enzyme(brenda = df, ids = Human_GEM.Genes)
brenda2data <- QueryBrenda(df, EC = unique(df$ID), n.core = 0, organisms = "Homo sapiens")
save.image('BrendaDb.RData')
