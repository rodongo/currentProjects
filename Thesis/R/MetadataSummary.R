





# Summarize the samples

Average.Age <- list(Male.Control = gsub('+','',CON_Samples$ageDeath[CON_Samples$sex=='male' & 
                                                                      CON_Samples$individualID %in% CO_SamplesUpdated]),
                    Female.Control = gsub('+','',CON_Samples$ageDeath[CON_Samples$sex=='female' & 
                                                                        CON_Samples$individualID %in% CO_SamplesUpdated]),
                    Male.AD = gsub('+','',AD_Samples$ageDeath[AD_Samples$sex=='male' & 
                                                                AD_Samples$individualID %in% AD_SamplesUpdated]),
                    Female.AD = gsub('+','',AD_Samples$ageDeath[AD_Samples$sex=='female' & 
                                                                  AD_Samples$individualID %in% AD_SamplesUpdated])) 
Average.Age <- lapply(Average.Age, function(x){x = as.numeric(gsub('\\+','',x))})
Average.Age <- lapply(Average.Age, mean)


APOE.Status <- list(Male.AD = gsub('+','',AD_Samples$apoeGenotype[AD_Samples$sex=='male' & AD_Samples$individualID %in% AD_SamplesUpdated]),
                    Female.AD = gsub('+','',AD_Samples$apoeGenotype[AD_Samples$sex=='female' & AD_Samples$individualID %in% AD_SamplesUpdated])) 
APOE.Status <- lapply(APOE.Status, table)
APOE.Status <- lapply(APOE.Status, as.data.frame)


BRAAK.Stage <- list(Male.AD = AD_Samples$Braak[AD_Samples$sex=='male' & AD_Samples$individualID %in% AD_SamplesUpdated],
                    Female.AD = AD_Samples$Braak[AD_Samples$sex=='female' & AD_Samples$individualID %in% AD_SamplesUpdated]) 

BRAAK.Stage <- lapply(BRAAK.Stage, table)
BRAAK.Stage <- lapply(BRAAK.Stage, as.data.frame)
