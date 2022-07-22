rm(list=ls()); cat("\014"); # Matlab'teki clear ve clc kombinasyonu
# Lib_loc <- unique(c("T:/Portable Files/Portable-Apps/R/Library", .libPaths())) # kullandigim k?t?phaler
Func_Folder <- "R" # Kendi yazdigim veya kullandigim kodlarin lokasyonu
#
Loc <- getwd()
# .libPaths(Lib_loc) # calismada kullandigim paket k?t?phanesi fazla kafa yormaya gerek yok
# .Options
options(digits=10, scipen=3)

###############################################################
library(stringr)
library(stringi)
library(limma)
# library(clusterProfiler)
# library(hgug4112a.db)
# library(annotate)
###############################################################
###############################################################
###############################################################
###############################################################
Targets <- readTargets("Target_1.txt", row.names=1, sep="\t") # Target dosyasi ne nedir bilgisini i?eriyor
setwd(paste0(Loc, "/Raw_Data") )
Targets$Files_Name[!file.exists(Targets$Files_Name)]
Targets$Files_Name[file.exists(Targets$Files_Name)]
Targets <- Targets[file.exists(Targets$Files_Name), ]
# 
x <- read.maimages(Targets$Files_Name, source="agilent", green.only=TRUE) # Belirtilen dosyalari sisteme y?kl?yor
#
colnames(x$E) <- Targets$ColNames # kolon isimlerini sisteme belirtim yoksa dosya isimlerini kullaniyor
#
setwd(Loc) # Normal calisma klas?r?ne geri d?nd?m
png(file="Output/Boxplot_Non_Normalize.png",	width = 480, height = 480, units = "px")
boxplot(log2(x$E) )#: Visualize the eigenvalues
dev.off()



y <- backgroundCorrect(x, method="normexp", offset=10)

# y <- normalizeWithinArrays(y, method="loess")
#
# for (i in 1:dim(y$targets)[1]) {
# 	png(file=paste0("Output/plotMD_", as.character(i), ".png"), width=1600, height=1200)
# 	plotMD(y, column=i, status=y$genes$ControlType)
# 	dev.off()
# }
# rm(i)

NormData <- normalizeBetweenArrays(y, method="quantile")
rm(x, y)
# 
Probe_ID <- NormData$genes$ProbeName
SystematicName <- NormData$genes$SystematicName
ControlType <- NormData$genes$ControlType
#
#
GPL21185 <- read.delim("Data/GPL21185-21174.txt", header=TRUE)
#
source(paste0(Func_Folder, "/", "ID_Converter_X.R"))

SYMBOL <- ID_Converter_X(Probe_ID, Columns_Number = 1, ID_List_Org = GPL21185[ , c(1, 6)], Remove_Unavailable_IDs = TRUE)
ENTREZID <- ID_Converter_X(Probe_ID, Columns_Number = 1, ID_List_Org = GPL21185[ , c(1, 5)], Remove_Unavailable_IDs = TRUE)
rm(GPL21185)
#
Norm_Data_Last <- data.frame(Probe_ID, SystematicName, ControlType, ENTREZID, SYMBOL, NormData$E)
colnames(Norm_Data_Last)[4:5] <- c("ENTREZID", "SYMBOL")
k_1 <- ! is.element(as.matrix(ENTREZID), "")
k_2 <- ! is.element(as.matrix(SYMBOL), "")
k <- which( (k_1 | k_2) & ControlType == 0)
rm(k_1, k_2)
#DEG filter

write.table(Norm_Data_Last, file = "Output/Norm_Data_Last.txt", sep="\t", col.names=TRUE, row.names=FALSE, qmethod="double", dec = ".")
###############################################################
write.table(Norm_Data_Last[k, ], file = "Output/Norm_Data_Last_Cleaned.txt", sep="\t", col.names=TRUE, row.names=FALSE, qmethod="double", dec = ".")
###############################################################
k_X <- as.matrix(which(! is.element(as.matrix(SYMBOL), ""))) 
Norm_Data_Last <- Norm_Data_Last[k_X, ]
rownames(Norm_Data_Last) <- NULL
rm(k, k_X)
rm(Probe_ID, SystematicName, ControlType, ENTREZID, SYMBOL, NormData)

###############################################################


library(stringr)
library(stringi)
library(limma)

Dim <- dim(Norm_Data_Last)
Data_Array <- Norm_Data_Last[ , 6:Dim[2]]
# Temp <- colnames(Data_Array) == Data_Pheno$geo_accession
Data_Expriment <- colnames(Norm_Data_Last)[6:Dim[2]]
Data_Expriment <- Targets
Data_GeneInfo <- Norm_Data_Last[ , 1:5]
rm(Targets, Norm_Data_Last)
# Temp <- rownames(Data_Array) == Data_GeneInfo$ID
# X6ohda_12h_rep
# Control_12h_rep
# 
# X6ohda_18h_rep
# Control_18h_rep
#

png(file="Output/Boxplot_Cleaned.png", width=1600, height=1200)
boxplot(Data_Array)#: Visualize the eigenvalues
dev.off()

#
Group_Control_12 <- which(str_detect(Data_Expriment$ColNames, "Control_12h_rep"))
Group_Control_18 <- which(str_detect(Data_Expriment$ColNames, "Control_18h_rep"))
#
Group_Disease_12 <- which(str_detect(Data_Expriment$ColNames, "n6_ohda_12h_rep"))
Group_Disease_18 <- which(str_detect(Data_Expriment$ColNames, "n6_ohda_18h_rep"))
#
Data_Array <- Data_Array[ , c(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18)]
Data_Expriment <- Data_Expriment[c(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18), ]
#
#
Group_Control_12 <- which(str_detect(Data_Expriment$ColNames, "Control_12h_rep"))
Group_Control_18 <- which(str_detect(Data_Expriment$ColNames, "Control_18h_rep"))
#
Group_Disease_12 <- which(str_detect(Data_Expriment$ColNames, "n6_ohda_12h_rep"))
Group_Disease_18 <- which(str_detect(Data_Expriment$ColNames, "n6_ohda_18h_rep"))
# 
Level <- factor(c(rep(0, each=length(c(Group_Control_12, Group_Control_18))), rep(1, each=length(c(Group_Disease_12, Group_Disease_18)))))
design <- cbind(Intercept=1, Grp2vs1=Level)
fit <- lmFit(as.matrix(Data_Array[ , c(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18)]), design)
cont.matrix <- makeContrasts(Grp2vs1, levels=design)
fit <- contrasts.fit(fit, cont.matrix)
#
fit2 <- eBayes(fit)
tT <- topTable(fit2, adjust="BH", number=dim(Data_Array)[1], sort.by="none")
rm(Level, design, fit, cont.matrix, fit2)
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

###############################################################
###############################################################
setwd(Loc)
save.image("BackUp.Rdata")
# rm(list=ls()); cat("\014")
# load(file="BackUp.Rdata")
###############################################################
###############################################################

Gene_Result <- data.frame(Data_GeneInfo$SYMBOL, tT[c("P.Value", "adj.P.Val", "AveExpr", "logFC")])

Gene_Result$tTest_P <- matrix(apply(Data_Array, 1, function(x){t.test(as.numeric(x[c(Group_Control_12, Group_Control_18)]), as.numeric(x[c(Group_Disease_12, Group_Disease_18)]))$p.value}), ncol=1)
Gene_Result$Avg_Control <- matrix(apply(Data_Array, 1, function(x){ mean( x[c(Group_Control_12, Group_Control_18)], na.rm = TRUE) }), ncol=1)
Gene_Result$Avg_Disease <- matrix(apply(Data_Array, 1, function(x){ mean( x[c(Group_Disease_12, Group_Disease_18)], na.rm = TRUE) }), ncol=1)
Gene_Result$FC  <- 2^Gene_Result$logFC

rm(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18)
#
# rm(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18)
rm(tT)

# kk <- which(!is.element(Data_GeneInfo$SYMBOL, "")) 
colnames(Gene_Result) <- c("Gene_Symbol", "Limma_P.Value", "Limma_adjP.Value", "AveExpr", "logFC", "tTest_P.Value", "Avg_Control", "Avg_Disease", "FC")

rownames(Gene_Result) <- NULL
#
write.table(Gene_Result, file="Output/Gene_Result.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
#
##############################################################################################################################
PPI_Network <- read.delim("Data/PPI_Network.txt", header=TRUE)#[ , c(1, 3)]
k <- which(!(as.character(PPI_Network[ , 1]) == as.character(PPI_Network[ , 3])))
PPI_Network <- PPI_Network[k, ]
rm(k)

write.table(PPI_Network, file="KPM/PPI_Network.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)

###############################################################

X <- data.frame(Gene_Result[ , c("Gene_Symbol")], ((Gene_Result[ , c("Limma_adjP.Value")] <= 0.05) & (Gene_Result[ , c("FC")] >= 1.20) )*1 )
colnames(X) <- c("Gene_Symbol", "Sig_Genes")
X <- aggregate(X$Sig_Genes, by = list(X$Gene_Symbol), FUN = sum)
colnames(X) <- c("Gene_Symbol", "Sig_Genes")
X$Sig_Genes <- (X$Sig_Genes>0)*1
write.table(X, file="KPM/Data/KPM_Input_BH_005_Up_Unique.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)



Y <- data.frame(Gene_Result[ , c("Gene_Symbol")], ((Gene_Result[ , c("Limma_adjP.Value")] <= 0.05) & (Gene_Result[ , c("FC")] <= (1/1.20)) )*1 )
colnames(Y) <- c("Gene_Symbol", "Sig_Genes")
Y <- aggregate(Y$Sig_Genes, by = list(Y$Gene_Symbol), FUN = sum)
colnames(Y) <- c("Gene_Symbol", "Sig_Genes")
Y$Sig_Genes <- (Y$Sig_Genes>0)*1
write.table(Y, file="KPM/Data/KPM_Input_BH_005_Down_Unique.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)

rm(X, Y)

setwd(paste0(Loc, "/KPM") )
###############################################################
###############################################################
CMD <- "java -jar -Xmx1G KPM-4.0.jar -strategy=INES -algo=ACO -K=5 -L1=0 -graphFile=PPI_Network.txt -matrix1=Data/KPM_Input_BH_005_Up_Unique.txt"

system(CMD, intern = TRUE, wait = TRUE)
#
dir.create("Outputs/BH_005_Up_Unique", showWarnings = FALSE)
file.copy("gene_stats-.txt", "Outputs/BH_005_Up_Unique/gene_stats-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.copy("summary-.txt", "Outputs/BH_005_Up_Unique/summary-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.copy("pathways_stats-.txt", "Outputs/BH_005_Up_Unique/pathways_stats-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.remove( c("gene_stats-.txt", "summary-.txt", "pathways_stats-.txt") )
#
file.copy(paste0("results/", list.files("results/")), paste0("Outputs/BH_005_Up_Unique/", list.files("results/")), overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.remove( paste0("results/", list.files("results/")) )
###############################################################
###############################################################
CMD <- "java -jar -Xmx1G KPM-4.0.jar -strategy=INES -algo=ACO -K=5 -L1=0 -graphFile=PPI_Network.txt -matrix1=Data/KPM_Input_BH_005_Down_Unique.txt"
system(CMD, intern = TRUE, wait = TRUE)
#
dir.create("Outputs/BH_005_Down_Unique", showWarnings = FALSE)
file.copy("gene_stats-.txt", "Outputs/BH_005_Down_Unique/gene_stats-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.copy("summary-.txt", "Outputs/BH_005_Down_Unique/summary-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.copy("pathways_stats-.txt", "Outputs/BH_005_Down_Unique/pathways_stats-.txt", overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.remove( c("gene_stats-.txt", "summary-.txt", "pathways_stats-.txt") )
#
file.copy(paste0("results/", list.files("results/")), paste0("Outputs/BH_005_Down_Unique/", list.files("results/")), overwrite = TRUE, recursive = FALSE,
		  copy.mode = TRUE, copy.date = FALSE)
file.remove( paste0("results/", list.files("results/")) )
rm(CMD)
setwd(Loc)
###############################################################
save.image("BackUp2.Rdata")
# rm(list=ls()); cat("\014")
# load(file="BackUp2.Rdata")
###############################################################


CNS_DrugTarget_List <- read.delim("Data/DrugTarget_List_CNS_Drugs.txt", header=TRUE)
All_DrugTarget_List <- read.delim("Data/DrugTarget_List_Human.txt", header=TRUE)
k <- which(is.element(CNS_DrugTarget_List[ , 2], Data_GeneInfo$SYMBOL))
CNS_DrugTarget_List <- CNS_DrugTarget_List[k, c(2, 4)]
#
k <- which(is.element(All_DrugTarget_List[ , 2], Data_GeneInfo$SYMBOL))
All_DrugTarget_List <- All_DrugTarget_List[k, c(2, 4)]
rm(k)



Gene_Result$CNS_DrugTarget <- ""
Gene_Result$All_DrugTarget <- ""
#
Uni_X <- unique(as.character(All_DrugTarget_List[ ,1]))
for (i in 1:length(Uni_X)) {
	Temp_A <- which(is.element(as.character(Gene_Result[ , 1]), Uni_X[i]))
	Temp_B <- which(is.element(as.character(CNS_DrugTarget_List[ ,1]), Uni_X[i]))
	Temp_C <- which(is.element(as.character(All_DrugTarget_List[ ,1]), Uni_X[i]))

	Gene_Result[Temp_A , "CNS_DrugTarget"] <- paste(CNS_DrugTarget_List[Temp_B ,2], sep = "", collapse = " ; ")
	Gene_Result[Temp_A , "All_DrugTarget"] <- paste(All_DrugTarget_List[Temp_C ,2], sep = "", collapse = " ; ")
	
	rm(Temp_A, Temp_B, Temp_C)
}
rm(i, Uni_X, All_DrugTarget_List, CNS_DrugTarget_List)
###############################################################
###############################################################
library(igraph)
library(stringr)

Files_List <- list.files("KPM/Outputs/BH_005_Up_Unique/")

Files_List <- Files_List[which(str_detect(Files_List , "-NODES-.txt"))]

Nodes_List_Up <- read.delim(paste0("KPM/Outputs/BH_005_Up_Unique/", Files_List[1]), header=FALSE)

for (i in 1:length(Files_List)) {
	Nodes_List_Up <- rbind.data.frame(Nodes_List_Up, read.delim(paste0("KPM/Outputs/BH_005_Up_Unique/", Files_List[i]), header=FALSE))
}
rm(i, Files_List)
Nodes_List_Up <- unique.data.frame(Nodes_List_Up)
####
Files_List <- list.files("KPM/Outputs/BH_005_Down_Unique/")

Files_List <- Files_List[which(str_detect(Files_List , "-NODES-.txt"))]

Nodes_List_Down <- read.delim(paste0("KPM/Outputs/BH_005_Down_Unique/", Files_List[1]), header=FALSE)

for (i in 1:length(Files_List)) {
	Nodes_List_Down <- rbind.data.frame(Nodes_List_Down, read.delim(paste0("KPM/Outputs/BH_005_Down_Unique/", Files_List[i]), header=FALSE))
}
rm(i, Files_List)
Nodes_List_Down <- unique.data.frame(Nodes_List_Down)
###############################################################

PPI_Network_Up <- PPI_Network[which((is.element(PPI_Network[ , 1], as.matrix(Nodes_List_Up) )==TRUE) & (is.element(PPI_Network[ , 3], as.matrix(Nodes_List_Up) )==TRUE)), ]
PPI_Network_Down <- PPI_Network[which((is.element(PPI_Network[ , 1], as.matrix(Nodes_List_Down) )==TRUE) & (is.element(PPI_Network[ , 3], as.matrix(Nodes_List_Down) )==TRUE)), ]

write.table(PPI_Network_Up, file="KPM/Final/KPM_BH_005_Up_Edges.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
write.table(PPI_Network_Down, file="KPM/Final/KPM_BH_005_Down_Edges.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
rm(PPI_Network_Up, PPI_Network_Down)

###############################################################
Gene_Result_Up <- Gene_Result[is.element(Gene_Result[ , 1], as.matrix(Nodes_List_Up) ), ]
Gene_Result_Down <- Gene_Result[is.element(Gene_Result[ , 1], as.matrix(Nodes_List_Down) ), ]
write.table(Gene_Result_Up, file="KPM/Final/KPM_BH_005_Up_Nodes.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
write.table(Gene_Result_Down, file="KPM/Final/KPM_BH_005_Down_Nodes.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)

###########
Gene_Result_Up_Edit <- data.frame(matrix(NA, ncol=dim(Gene_Result_Up)[2]+1, nrow=0))
k_Up <- c()

for (i in 1:length(Nodes_List_Up[, 1]) ) {
	k <- which(is.element(Gene_Result_Up$Gene_Symbol, Nodes_List_Up[i, ] ))
	if (length(k) == 0) {
		rm(k)
		next
	} else if (length(k) == 1) {
		k_Up <- c(k_Up, k)
	} else if (length(k) > 1) {
		Sig_Genes <- which(Gene_Result_Up$Limma_adjP.Value[k] <= 0.05 & (Gene_Result_Up$FC[k] <= (1/1.2) | Gene_Result_Up$FC[k] >= 1.2))
		if ( length(Sig_Genes) == 0 ) {
			k_Up <- c(k_Up, k[which.min(Gene_Result_Up$Limma_adjP.Value[k])])
		}	else {
			k_Up <- c(k_Up, k[ Gene_Result_Up$Limma_adjP.Value[k]==min(Gene_Result_Up$Limma_adjP.Value[k[Sig_Genes]])])
		}
		rm(Sig_Genes)	
	}
	rm(k)	
}
rm(i)

Gene_Result_Up <- Gene_Result_Up[k_Up, ]
###########
Gene_Result_Down_Edit <- data.frame(matrix(NA, ncol=dim(Gene_Result_Down)[2]+1, nrow=0))
k_Down <- c()

for (i in 1:length(Nodes_List_Down[, 1]) ) {
	k <- which(is.element(Gene_Result_Down$Gene_Symbol, Nodes_List_Down[i, ] ))
	if (length(k) == 0) {
		rm(k)
		next
	} else if (length(k) == 1) {
		k_Down <- c(k_Down, k)
	} else if (length(k) > 1) {
		Sig_Genes <- which(Gene_Result_Down$Limma_adjP.Value[k] <= 0.05 & (Gene_Result_Down$FC[k] <= (1/1.2) | Gene_Result_Down$FC[k] >= 1.2))
		if ( length(Sig_Genes) == 0 ) {
			k_Down <- c(k_Down, k[which.min(Gene_Result_Down$Limma_adjP.Value[k])])
		}	else {
			k_Down <- c(k_Down, k[ Gene_Result_Down$Limma_adjP.Value[k]==min(Gene_Result_Down$Limma_adjP.Value[k[Sig_Genes]])])
		}
		rm(Sig_Genes)	
	}
	rm(k)	
}
rm(i)

Gene_Result_Down <- Gene_Result_Down[k_Down, ]

write.table(Gene_Result_Up, file="KPM/Final/KPM_BH_005_Up_Nodes_OneProbe.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
write.table(Gene_Result_Down, file="KPM/Final/KPM_BH_005_Down_Nodes_OneProbe.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
###########



Temp_N1 <- read.delim("KPM/Outputs/BH_005_Up_Unique/Pathway-01-NODES-.txt", header=FALSE)
write.table(Gene_Result_Up[is.element(Gene_Result_Up$Gene_Symbol, Temp_N1$V1) , ], file="KPM/Final/KPM_BH_005_Up_Nodes_OneProbe_N1.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)

Temp_N1 <- read.delim("KPM/Outputs/BH_005_Down_Unique/Pathway-01-NODES-.txt", header=FALSE)
write.table(Gene_Result_Down[is.element(Gene_Result_Down$Gene_Symbol, Temp_N1$V1) , ], file="KPM/Final/KPM_BH_005_Down_Nodes_OneProbe_N1.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".", quote = FALSE)
rm(Temp_N1)
# rm(Gene_Result_Up, Gene_Result_Down)
###############################################################

setwd(Loc)


PreFilter_pValues <- 0.05
Sig_Genes <- as.vector(Gene_Result$Gene_Symbol[which(Gene_Result$tTest_P.Value < PreFilter_pValues)])
k <- as.matrix(which(Gene_Result$tTest_P.Value < PreFilter_pValues)) 




Group_Control_12 <- which(str_detect(colnames(Data_Array), "Control_12h_rep"))
Group_Control_18 <- which(str_detect(colnames(Data_Array), "Control_18h_rep"))
#
Group_Disease_12 <- which(str_detect(colnames(Data_Array), "n6_ohda_12h_rep"))
Group_Disease_18 <- which(str_detect(colnames(Data_Array), "n6_ohda_18h_rep"))
Group_Control <- c(Group_Control_12, Group_Control_18)
Group_Disease <- c(Group_Disease_12, Group_Disease_18)

# Group_Control <- c(Group_Control_18)
# Group_Disease <- c(Group_Disease_18)





library(factoextra)
#
GSE.PCA <- prcomp(t(Data_Array[ , c(Group_Control, Group_Disease)]), scale = TRUE)
###############################################################
PCA_Dim <- get_eigenvalue(GSE.PCA)#: Extract the eigenvalues/variances of principal components
write.table(PCA_Dim, file="Output/PCA_Dim.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".")
rm(PCA_Dim)
#
png(file="Output/eigenvalue.png", width=1600, height=1200)
fviz_eig(GSE.PCA)#: Visualize the eigenvalues
dev.off()
#
Temp_Name <- factor(c(rep("Group_Control_12", each=length(Group_Control_12)), c(rep("Group_Control_18", each=length(Group_Control_18)), 
					rep("Group_Disease_12", each=length(Group_Disease_12)), rep("Group_Disease_18", each=length(Group_Disease_18)))))
Temp_Index <- as.character(c(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18))#
png(file="Output/PCA_2D_Matrix.png", width=1600, height=1200)
Data <- GSE.PCA$x[ , 1:3]
species_labels <- as.factor(Temp_Name)
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(length(unique(species_labels))))[as.numeric(species_labels)]
# Plot a SPLOM:
pairs(Data, col = species_col, 
	  lower.panel = NULL, 
	  cex.labels=2, pch=19, cex = 4)
# Add a legend
par(xpd = TRUE)
legend(x = 0.1, y = 0.3, cex = 3,
	   legend = as.character(levels(species_labels)),
	   fill = unique(species_col))
par(xpd = NA)
dev.off()
rm(Data, species_labels, species_col, Temp_Name, Temp_Index)
###############################################################
GSE.PCA <- prcomp(t(Data_Array[k, c(Group_Control, Group_Disease)]), scale = TRUE)
###############################################################
PCA_Dim <- get_eigenvalue(GSE.PCA)#: Extract the eigenvalues/variances of principal components
write.table(PCA_Dim, file="Output/PCA_Dim_k.txt", sep="\t", col.names=TRUE, row.names = FALSE, qmethod="double", dec = ".")
rm(PCA_Dim)
#
png(file="Output/eigenvalue_k.png", width=1600, height=1200)
fviz_eig(GSE.PCA)#: Visualize the eigenvalues
dev.off()
#
Temp_Name <- factor(c(rep("Group_Control_12", each=length(Group_Control_12)), c(rep("Group_Control_18", each=length(Group_Control_18)), 
																				rep("Group_Disease_12", each=length(Group_Disease_12)), rep("Group_Disease_18", each=length(Group_Disease_18)))))
Temp_Index <- as.character(c(Group_Control_12, Group_Control_18, Group_Disease_12, Group_Disease_18))#
png(file="Output/PCA_2D_Matrix_k.png", width=1600, height=1200)
Data <- GSE.PCA$x[ , 1:3]
species_labels <- as.factor(Temp_Name)
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(length(unique(species_labels))))[as.numeric(species_labels)]
# Plot a SPLOM:
pairs(Data, col = species_col, 
	  lower.panel = NULL, 
	  cex.labels=2, pch=19, cex = 4)
# Add a legend
par(xpd = TRUE)
legend(x = 0.1, y = 0.3, cex = 3,
	   legend = as.character(levels(species_labels)),
	   fill = unique(species_col))
par(xpd = NA)
dev.off()
rm(Data, species_labels, species_col, Temp_Name, Temp_Index)
###############################################################
###############################################################
rm(GSE.PCA)
###############################################################
###############################################################
###############################################################
require(graphics); require(grDevices)
#
Selected_Gene_Index <- 1:length(k)

Data_HeatMap <- Data_Array[k[Selected_Gene_Index], c(Group_Control, Group_Disease)]
Col_Name <- factor(c(rep("Group_Control_12", each=length(Group_Control_12)), c(rep("Group_Control_18", each=length(Group_Control_18)),
					rep("Group_Disease_12", each=length(Group_Disease_12)), rep("Group_Disease_18", each=length(Group_Disease_18)))))
#
#
Row_Name <- Data_GeneInfo[k[Selected_Gene_Index], 5]
#
#
png(file="Output/PD_HeatMap_1.png", width=1600, height=1200)
heatmap(as.matrix(Data_HeatMap), labRow = Row_Name, labCol = Col_Name, scale = "row", margins = c(15, 3), main = "PD_HeatMap", xlab = "Samples", ylab =  "Genes")
dev.off()
# #
Col_Name <- factor(Data_Expriment[c(Group_Control, Group_Disease), 3])
#
png(file="Output/PD_HeatMap_2.png", width=1600, height=1200)
heatmap(as.matrix(Data_HeatMap), labRow = Row_Name, labCol = Col_Name, scale = "row", margins = c(15, 3), main = "PD_HeatMap", xlab = "Samples", ylab =  "Genes")
dev.off()


rm(Selected_Gene_Index, Data_HeatMap, Col_Name, Row_Name)
#
rm(PreFilter_pValues, k, Dim)
rm(Group_Control, Group_Control_12, Group_Control_18, Group_Disease, Group_Disease_12, Group_Disease_18)

