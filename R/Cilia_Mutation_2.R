#      Created 05/08/2024             #
#                                     #
#      Author: Nicholas Hadas         #
#             [Redacted]              #


####-------Install Packages-------####
install.packages("ggraph")
install.packages("igraph")
install.packages('BiocManager')
install.packages("clustree")
BiocManager::install('multtest')
install.packages(c('dplyr', 'Seurat', 'patchwork', 'devtools', 'tidyverse', 
                   'gridExtra', 'harmony', 'metap', 'pals'))
devtools::install_github('satijalab/seurat-data')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
install.packages('hdf5r')
BiocManager::install("zellkonverter")

library(DoubletFinder)
library(Seurat)
library(SeuratData)
library(patchwork)
library(devtools)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)
library(sctransform)
library(BiocManager)
library(multtest)
library(metap)
library(ggraph)
library(clustree)
library(pals)
library(hdf5r)
library(SeuratDisk)
library(zellkonverter)
library(dplyr)


####-------Set a Working Directory-------####
setwd('/Users/labuser/Documents/SR004328') # Work Mac
setwd('/home/haider/Documents/Nico_Data/SR004328') # Supercomputer

####-------Load Raw Data-------####

# PCD171 #
PCD171 <- Read10X("Raw_Data/SR004328_PCD171/outs/filtered_feature_bc_matrix")

# PCD171M #
PCD171M <- Read10X("Raw_Data/SR004328_PCD171M/outs/filtered_feature_bc_matrix")

# PCD172 #
PCD172 <- Read10X("Raw_Data/SR004328_PCD172/outs/filtered_feature_bc_matrix")

# PCD177 #
PCD177 <- Read10X("Raw_Data/SR004328_PCD177/outs/filtered_feature_bc_matrix")


####-------Convert to Seurat Object-------####

PCD171 <- CreateSeuratObject(counts = PCD171, project = "SR004328", min.cells = 3, min.features = 200)
PCD171M <- CreateSeuratObject(counts = PCD171M, project = "SR004328", min.cells = 3, min.features = 200)
PCD172 <- CreateSeuratObject(counts = PCD172, project = "SR004328", min.cells = 3, min.features = 200)
PCD177 <- CreateSeuratObject(counts = PCD177, project = "SR004328", min.cells = 3, min.features = 200)


####-------Add Metadata-------####

# PCD #
PCD171@meta.data$stim <- "PCD"
PCD172@meta.data$stim <- "PCD"
PCD177@meta.data$stim <- "PCD"

# Mother #
PCD171M@meta.data$stim <- "Mother"

# Sample ID #
PCD171@meta.data$sampleid <- "PCD171"
PCD172@meta.data$sampleid <- "PCD172"
PCD177@meta.data$sampleid <- "PCD177"
PCD171M@meta.data$sampleid <- "PCD171M"



####-------QC & Filtering-------####

# PCD171 #
PCD171[["percent.mt"]] <- PercentageFeatureSet(PCD171, pattern = "^MT-")
PCD171[["percent.rp"]] <- PercentageFeatureSet(PCD171, pattern = "^RP[SL]")
VlnPlot(PCD171, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
        ncol = 3)
plot1 <- FeatureScatter(PCD171, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(PCD171, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
PCD171 <- subset(PCD171, subset = nFeature_RNA < 11000 & percent.mt < 30 & percent.rp < 25)

# PCD171M #
PCD171M[["percent.mt"]] <- PercentageFeatureSet(PCD171M, pattern = "^MT-")
PCD171M[["percent.rp"]] <- PercentageFeatureSet(PCD171M, pattern = "^RP[SL]")
VlnPlot(PCD171M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
        ncol = 3)
plot1 <- FeatureScatter(PCD171M, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(PCD171M, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
PCD171M <- subset(PCD171M, subset = nFeature_RNA < 10000 & percent.mt < 25 & percent.rp < 11)

# PCD172 #
PCD172[["percent.mt"]] <- PercentageFeatureSet(PCD172, pattern = "^MT-")
PCD172[["percent.rp"]] <- PercentageFeatureSet(PCD172, pattern = "^RP[SL]")
VlnPlot(PCD172, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
        ncol = 3)
plot1 <- FeatureScatter(PCD172, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(PCD172, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
PCD172 <- subset(PCD172, subset = nFeature_RNA < 10000 & percent.mt < 35 & percent.rp < 17)

# PCD177 #
PCD177[["percent.mt"]] <- PercentageFeatureSet(PCD177, pattern = "^MT-")
PCD177[["percent.rp"]] <- PercentageFeatureSet(PCD177, pattern = "^RP[SL]")
VlnPlot(PCD177, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
        ncol = 3)
plot1 <- FeatureScatter(PCD177, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(PCD177, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
PCD177 <- subset(PCD177, subset = nFeature_RNA < 12500 & percent.mt < 50 & percent.rp < 20)


####-------Doublet Finder-------####

# Preprocessing #
a <- c(PCD171, PCD171M, PCD172, PCD177)
for(x in 1:length(a)) {
  a[[x]] <- NormalizeData(a[[x]])
  a[[x]] <- FindVariableFeatures(a[[x]])
  a[[x]] <- ScaleData(a[[x]])
  a[[x]] <- RunPCA(a[[x]])
  print(ElbowPlot(a[[x]]))
}

# Create List of PCs #
b <- c(18, 17, 18, 20)

# Run Doublet Finder #
for (m in 1:length(a)) {
  sweep.res.list_nsclc <- paramSweep(a[[m]], PCs = 1:b[[m]], sct = FALSE)
  sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
  
  pK <- bcmvn_nsclc %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- a[[m]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(a[[m]]@meta.data)) # Assuming 7.5% doublet formation rate
  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  a[[m]] <- doubletFinder(a[[m]],
                          PCs = 1:b[[m]],
                          pN = 0.25,
                          pK = pK,
                          nExp = nExp_poi.adj,
                          reuse.pANN = FALSE, sct = FALSE)
}

# Visualize Doublets #
for (l in 1:length(a)) {
  colnames(a[[l]]@meta.data)[8] <- "DoubletFinder" 
  a[[l]]@meta.data <- a[[l]]@meta.data[, -7]
  print(table(a[[l]]@meta.data$DoubletFinder))
}

# Remove Doublets #
PCD171 <- subset(x = a[[1]], subset = DoubletFinder == "Singlet")
PCD171M <- subset(x = a[[2]], subset = DoubletFinder == "Singlet")
PCD172 <- subset(x = a[[3]], subset = DoubletFinder == "Singlet")
PCD177 <- subset(x = a[[4]], subset = DoubletFinder == "Singlet")

# Save QCd and Doublet Removed Samples #
Doublet_List_SR004328 <- c(PCD171, PCD171M, PCD172, PCD177)
save(Doublet_List_SR004328, file = "Doublet_List_SR004328.rds")


####-------Merge-------####

# Merge Original Object with Published Data #
merged_object <- merge(PCD171, y = c(PCD171M, PCD172, PCD177), 
                       add.cell.ids = c("PCD171", "PCD171M", "PCD172", "PCD177"), project = "SR004328")

# Save #
saveRDS(merged_object, "merged_SR004328_object.rds")


####-------SCTransform-------####

# Load Object #
merged_object <- readRDS("Objects/merged_SR004328_object.rds")

# Join Layers #
merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])

# Split Object by Individual Donors #
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f = merged_object$sampleid)

# Install glmGamPoi for More Efficient Implementation #
# BiocManager::install("glmGamPoi")
# library(glmGamPoi)

# Run SCTransform
merged_object_SCT <- SCTransform(merged_object)

# Run PCA #
merged_object_SCT <- RunPCA(merged_object_SCT)

# Run UMAP #
merged_object_SCT <- RunUMAP(merged_object_SCT, dims = 1:30)

# Visualize Unintegrated UMAP #
DimPlot(merged_object_SCT, reduction = "umap", group.by = c("sampleid"))


####-------CCA Integration-------####

# Perform CCA Integration #
merged_object_SCT_Integrated <- IntegrateLayers(object = merged_object_SCT, method = CCAIntegration, normalization.method = "SCT", verbose = F)

# Run Find Neighbours #
merged_object_SCT_Integrated <- FindNeighbors(merged_object_SCT_Integrated, reduction = "integrated.dr", dims = 1:30)

# Find Clusters #
merged_object_SCT_Integrated <- FindClusters(merged_object_SCT_Integrated, resolution = 1.0)

# Run UMAP #
merged_object_SCT_Integrated <- RunUMAP(merged_object_SCT_Integrated, dims = 1:30, reduction = "integrated.dr")

# Save #
saveRDS(merged_object_SCT_Integrated, "merged_object_SCT_Integrated_SR004328.rds")

# Visualize UMAP #
DimPlot(merged_object_SCT_Integrated, reduction = "umap", group.by = c("seurat_clusters"))
DimPlot(merged_object_SCT_Integrated, reduction = "umap", label = TRUE)

# Join Layers #
merged_object_SCT_Integrated[["RNA"]] <- JoinLayers(merged_object_SCT_Integrated[["RNA"]])


####-------Cluster Annotation-------####

# Load Object #
merged_object_SCT_Integrated <- readRDS("Objects/merged_object_SCT_Integrated_SR004328.rds")

# Change Default Assay #
DefaultAssay(merged_object_SCT_Integrated) <- 'RNA'

# Run FindAllMarkers #
All.Markers <- FindAllMarkers(merged_object_SCT_Integrated, only.pos = TRUE, max.cells.per.ident = 1000)

# Save as CSV #
write.csv(All.Markers, "SR004328_pub_markers.csv")

# Pertinent Markers #
FeaturePlot(merged_object_SCT_Integrated,"JUNB") # Stress/Ribosomal Contamination
FeaturePlot(merged_object_SCT_Integrated, "nFeature_RNA") # RNA Levels
FeaturePlot(merged_object_SCT_Integrated, "DNAH5") # Ciliated
FeaturePlot(merged_object_SCT_Integrated, "BANK1") # B Cells
FeaturePlot(merged_object_SCT_Integrated, "JCHAIN") # Germinal B cells
FeaturePlot(merged_object_SCT_Integrated, "C1QA") # IM
FeaturePlot(merged_object_SCT_Integrated, "CD300E") # ncMONO
FeaturePlot(merged_object_SCT_Integrated, "THEMIS") # T Cells
FeaturePlot(merged_object_SCT_Integrated, "EAR2") # AM
FeaturePlot(merged_object_SCT_Integrated, "TOP2A") # Div
FeaturePlot(merged_object_SCT_Integrated, "CCR2") # cMono
FeaturePlot(merged_object_SCT_Integrated, "KRT5") # Basal
FeaturePlot(merged_object_SCT_Integrated, "CYP2F1", label = TRUE) # Secretory

# Set Assay to RNA #
DefaultAssay(merged_object_SCT_Integrated) <- "RNA"

# Rename Idents #
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "0" = "PMN") 
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "1" = "NK")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "2" = "AM")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "3" = "Secretory 6")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "4" = "Secretory 8")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "5" = "Basal 1")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "6" = "Secretory 9")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "7" = "Secretory 10")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "8" = "Secretory 5")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "9" = "T Cells")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "10" = "Secretory 3")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "11" = "Secretory 4")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "12" = "Secretory 9") # Subcluster
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "13" = "IM+ncMono+cMono") # Subcluster
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "14" = "Basal 2")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "15" = "Secretory 11")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "16" = "MT")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "17" = "Ciliated 2")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "18" = "PMN")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "19" = "Ciliated 1")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "20" = "B+GB") # Subcluster
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "21" = "Doublets")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "22" = "Doublets")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "23" = "Secretory 7")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "24" = "Div") # Subcluster
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "25" = "Hemoglobin")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "26" = "Secretory 1") # Subcluster
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "27" = "Doublets")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "28" = "Secretory 2")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "29" = "Ionocyte")

# Visualize #
DimPlot(merged_object_SCT_Integrated, reduction = "umap", label = TRUE) + NoLegend()

# Remove Low Quality Clusters #
merged_object_SCT_Integrated <- subset(merged_object_SCT_Integrated, idents = c("Doublets", "MT", "Hemoglobin"), invert = TRUE)

# Make New Column #
merged_object_SCT_Integrated@meta.data$Cluster.Names <- Idents(merged_object_SCT_Integrated)


##---Secretory 9 Subcluster---##

# Subset #
Secretory_9 <- subset(merged_object_SCT_Integrated, subset=Cluster.Names %in% "Secretory 9")

# Process #
DefaultAssay(Secretory_9) <- 'SCT'
Secretory_9 <- FindNeighbors(Secretory_9, reduction = "integrated.dr", dims = 1:30)
Secretory_9 <- FindClusters(Secretory_9, resolution = 1)

# Visualize #
DimPlot(Secretory_9, label = TRUE)

# Rename Clusters #
DefaultAssay(Secretory_9) <- 'RNA'
Secretory_9 <- RenameIdents(Secretory_9, '0' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '1' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '2' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '3' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '4' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '5' = "Secretory 9")
Secretory_9 <- RenameIdents(Secretory_9, '6' = "Remove")
Secretory_9 <- RenameIdents(Secretory_9, '7' = "Secretory 9")

# Visualize #
DimPlot(Secretory_9, label = TRUE) + NoLegend()

# Create New Column of Idents #
Secretory_9$final_cluster <- Secretory_9@active.ident


##---IM+ncMono+cMono Subcluster---##

# Subset #
IM_Mono <- subset(merged_object_SCT_Integrated, subset=Cluster.Names %in% "IM+ncMono+cMono")

# Process #
DefaultAssay(IM_Mono) <- 'SCT'
IM_Mono <- FindNeighbors(IM_Mono, reduction = "integrated.dr", dims = 1:30)
IM_Mono <- FindClusters(IM_Mono, resolution = 1)

# Visualize #
DimPlot(IM_Mono, label = TRUE)

# Rename Clusters #
DefaultAssay(IM_Mono) <- 'RNA'
IM_Mono <- RenameIdents(IM_Mono, '0' = "IM")
IM_Mono <- RenameIdents(IM_Mono, '1' = "ncMono")
IM_Mono <- RenameIdents(IM_Mono, '2' = "Remove")
IM_Mono <- RenameIdents(IM_Mono, '3' = "cMono")
IM_Mono <- RenameIdents(IM_Mono, '4' = "ncMono")
IM_Mono <- RenameIdents(IM_Mono, '5' = "Remove")
IM_Mono <- RenameIdents(IM_Mono, '6' = "ncMono")
IM_Mono <- RenameIdents(IM_Mono, '7' = "ncMono")
IM_Mono <- RenameIdents(IM_Mono, '8' = "cMono")
IM_Mono <- RenameIdents(IM_Mono, '9' = "IM")

# Visualize #
DimPlot(IM_Mono, label = TRUE) + NoLegend()

# Create New Column of Idents #
IM_Mono$final_cluster <- IM_Mono@active.ident


##---Secretory 1 Subcluster---##

# Subset #
Secretory_1 <- subset(merged_object_SCT_Integrated, subset=Cluster.Names %in% "Secretory 1")

# Process #
DefaultAssay(Secretory_1) <- 'SCT'
Secretory_1 <- FindNeighbors(Secretory_1, reduction = "integrated.dr", dims = 1:30)
Secretory_1 <- FindClusters(Secretory_1, resolution = 1)

# Visualize #
DimPlot(Secretory_1, label = TRUE)

# Rename Clusters #
DefaultAssay(Secretory_1) <- 'RNA'
Secretory_1 <- RenameIdents(Secretory_1, '0' = "Secretory 1")
Secretory_1 <- RenameIdents(Secretory_1, '1' = "Ciliated 3")
Secretory_1 <- RenameIdents(Secretory_1, '2' = "Secretory 1")
Secretory_1 <- RenameIdents(Secretory_1, '3' = "Secretory 1")

# Visualize #
DimPlot(Secretory_1, label = TRUE) + NoLegend()

# Create New Column of Idents #
Secretory_1$final_cluster <- Secretory_1@active.ident


##---B+GB Subcluster---##

# Subset #
B_GB <- subset(merged_object_SCT_Integrated, subset=Cluster.Names %in% "B+GB")

# Process #
DefaultAssay(B_GB) <- 'SCT'
B_GB <- FindNeighbors(B_GB, reduction = "integrated.dr", dims = 1:30)
B_GB <- FindClusters(B_GB, resolution = 1)

# Visualize #
DimPlot(B_GB, label = TRUE)

# Rename Clusters #
DefaultAssay(B_GB) <- 'RNA'
B_GB <- RenameIdents(B_GB, '0' = "B Cells")
B_GB <- RenameIdents(B_GB, '1' = "B Cells")
B_GB <- RenameIdents(B_GB, '2' = "B Cells")
B_GB <- RenameIdents(B_GB, '3' = "B Cells")
B_GB <- RenameIdents(B_GB, '4' = "B Cells")
B_GB <- RenameIdents(B_GB, '5' = "B Cells")
B_GB <- RenameIdents(B_GB, '6' = "GB Cells")

# Visualize #
DimPlot(B_GB, label = TRUE) + NoLegend()

# Create New Column of Idents #
B_GB$final_cluster <- B_GB@active.ident


##---Div Subcluster---##

# Subset #
Div <- subset(merged_object_SCT_Integrated, subset=Cluster.Names %in% "Div")

# Process #
DefaultAssay(Div) <- 'SCT'
Div <- FindNeighbors(Div, reduction = "integrated.dr", dims = 1:30)
Div <- FindClusters(Div, resolution = 1)

# Visualize #
DimPlot(Div, label = TRUE)

# Rename Clusters #
DefaultAssay(Div) <- 'RNA'
Div <- RenameIdents(Div, '0' = "Div")
Div <- RenameIdents(Div, '1' = "Div")
Div <- RenameIdents(Div, '2' = "Div")
Div <- RenameIdents(Div, '3' = "Div")
Div <- RenameIdents(Div, '4' = "Remove")
Div <- RenameIdents(Div, '5' = "Div")


# Visualize #
DimPlot(Div, label = TRUE) + NoLegend()

# Create New Column of Idents #
Div$final_cluster <- Div@active.ident


##---Annotating Back---##

# Create column with Idents #
merged_object_SCT_Integrated$final_cluster <- merged_object_SCT_Integrated@active.ident
merged_object_SCT_Integrated$final_cluster <- as.character(merged_object_SCT_Integrated$final_cluster)

# Paste Idents of subclusters into final_cluster column of merged_object_SCT_Integrated #
merged_object_SCT_Integrated$final_cluster[colnames(merged_object_SCT_Integrated)[colnames(merged_object_SCT_Integrated) %in% colnames(Secretory_9)]] <- paste0(Secretory_9@active.ident)
merged_object_SCT_Integrated$final_cluster[colnames(merged_object_SCT_Integrated)[colnames(merged_object_SCT_Integrated) %in% colnames(IM_Mono)]] <- paste0(IM_Mono@active.ident)
merged_object_SCT_Integrated$final_cluster[colnames(merged_object_SCT_Integrated)[colnames(merged_object_SCT_Integrated) %in% colnames(Secretory_1)]] <- paste0(Secretory_1@active.ident)
merged_object_SCT_Integrated$final_cluster[colnames(merged_object_SCT_Integrated)[colnames(merged_object_SCT_Integrated) %in% colnames(B_GB)]] <- paste0(B_GB@active.ident)
merged_object_SCT_Integrated$final_cluster[colnames(merged_object_SCT_Integrated)[colnames(merged_object_SCT_Integrated) %in% colnames(Div)]] <- paste0(Div@active.ident)

# Set New Idents #
merged_object_SCT_Integrated$final_cluster <- as.factor(merged_object_SCT_Integrated$final_cluster)
merged_object_SCT_Integrated@active.ident <- merged_object_SCT_Integrated$final_cluster

# Visualize #
DimPlot(merged_object_SCT_Integrated, label=TRUE) 

# Remove Unwanted Clusters #
merged_object_SCT_Integrated <- subset(merged_object_SCT_Integrated, idents = "Remove", invert = TRUE)

# Visualize #
DimPlot(merged_object_SCT_Integrated, label=TRUE, label.size = 5) + NoLegend()

# Save #
saveRDS(merged_object_SCT_Integrated, "merged_object_SCT_Integrated_SR004328_numbered.rds")


##---Lump Airway Cell Cluster Names---##

# Rename Idents #
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 1" = "Secretory") 
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 2" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 3" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 4" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 5" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 6" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 7" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 8" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 9" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 10" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory 11" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Ciliated 1" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Ciliated 2" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Ciliated 3" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Basal 1" = "Basal")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Basal 2" = "Basal")

# Visualize #
DimPlot(merged_object_SCT_Integrated, label=TRUE, label.size = 5) + NoLegend()

# Visualize Integration by Sample ID #
DimPlot(merged_object_SCT_Integrated, label.size = 5, group.by = "sampleid")

merged_object_SCT_Integrated$final_cluster <- Idents(merged_object_SCT_Integrated)

# Save #
saveRDS(merged_object_SCT_Integrated, "merged_object_SCT_Integrated_SR004328_annotated.rds")


####-------Integration Analysis-------####

# Number of Cells per Cluster per Identity #
table <- as.data.frame.matrix(table(merged_object_SCT_Integrated@meta.data$final_cluster, merged_object_SCT_Integrated@meta.data$sampleid))

# Export as CSV #
write.csv(table, "cells_per_ident.csv")


####-------Cell Cluster Analysis-------####

# Install dittoSeq #
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")
library("dittoSeq")

# Run Bar Plot #
dittoBarPlot(merged_object_SCT_Integrated, "final_cluster", scale = "percent", main = "Percent Cell Counts", group.by = "sampleid") + 
  theme(title = element_text(face="bold", size=20), plot.title = element_text(hjust=0.5))


####-------DEG Analysis-------####

# Load Object #
integrated_object <- readRDS("Objects/merged_object_SCT_Integrated_SR004328_annotated.rds")

# Join Layers #
integrated_object[["RNA"]] <- JoinLayers(integrated_object[["RNA"]])

# Create New Column for PCD vs Mother #
integrated_object$stim2 <- integrated_object$sampleid

# Reformat Cluster Names #
integrated_object$stim2 <- as.character(integrated_object$stim2)

# Change Metadata Elements #
integrated_object$stim2[integrated_object$stim2 %in% c("PCD171", "PCD172")] <- paste0("PCD")
integrated_object$stim2[integrated_object$stim2 %in% 'PCD171M'] <- paste0('Mother')
integrated_object$stim2[integrated_object$stim2 %in% 'PCD177'] <- paste0('Possible_Control')

# Reformat Cluster Names #
integrated_object$stim2 <- as.factor(integrated_object$stim2)

# Add Column Discerning Control v PCD #
integrated_object$stim3 <- paste(integrated_object$final_cluster, integrated_object$stim2, sep='_')

# Set Idents #
Idents(integrated_object) <- integrated_object$stim3

# Create List of Cell Types #
pcd_list <- levels(integrated_object)[1:14]
pcd_list <- pcd_list[-c(11, 13)] # Remove GB Cells & Ionocytes 
pcd_list <- pcd_list[order(names(setNames(pcd_list, pcd_list)))]

mother_list <- levels(integrated_object)[15:28]
mother_list <- mother_list[-c(11, 14)] # Remove GB Cells & Ionocytes 
mother_list <- mother_list[order(names(setNames(mother_list, mother_list)))]

possible_control_list <- levels(integrated_object)[29:42]
possible_control_list <- possible_control[-c(11, 14)] # Remove GB Cells & Ionocytes 
possible_control_list <- possible_control[order(names(setNames(possible_control, possible_control)))]

# Run PrepSCTFindMarkers #
integrated_object <- PrepSCTFindMarkers(integrated_object)


####-------DE Analysis (PCD vs Mother) - Wilcox - SCT -------####

# Run FindMarkers #
for (i in 1:12) {
  name = FindMarkers(integrated_object, ident.1 = pcd_list[[i]], ident.2 = mother_list[[i]], verbose = TRUE)
  file = paste0(unlist((strsplit(mother_list[[i]], split = "_")))[1], "_wilcox_", "[Redacted]_", "SCT",  ".csv")
  write.csv(name, file)
}


####-------DE Analysis (PCD vs Mother) - MAST - SCT-------####

# Run FindMarkers #
for (i in 1:12) {
  name = FindMarkers(integrated_object, ident.1 = pcd_list[[i]], ident.2 = mother_list[[i]], verbose = TRUE, test.use = "MAST")
  file = paste0(unlist((strsplit(mother_list[[i]], split = "_")))[1], "_MAST_", "[Redacted]_", "SCT", ".csv")
  write.csv(name, file)
}


####-------DE Analysis (PCD vs Potential Control) - Wilcox - SCT -------####

# Run FindMarkers #
for (i in 1:12) {
  name = FindMarkers(integrated_object, ident.1 = pcd_list[[i]], ident.2 = possible_control_list[[i]], verbose = TRUE)
  file = paste0(unlist((strsplit(mother_list[[i]], split = "_")))[1], "_wilcox_", "[Redacted]_", "possible_control_", "SCT",  ".csv")
  write.csv(name, file)
}


####-------DE Analysis (PCD vs Potential Control) - MAST - SCT-------####

# Run FindMarkers #
for (i in 1:12) {
  name = FindMarkers(integrated_object, ident.1 = pcd_list[[i]], ident.2 = possible_control_list[[i]], verbose = TRUE, test.use = "MAST")
  file = paste0(unlist((strsplit(mother_list[[i]], split = "_")))[1], "_MAST_", "[Redacted]_", "possible_control_", "SCT", ".csv")
  write.csv(name, file)
}


####-------Pathway Analysis - Set-Up-------####

# Install Packages #
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("rWikiPathways")
library(clusterProfiler)
library("org.Hs.eg.db")
library(dplyr)
library(tidyverse)
library(rWikiPathways)

# Download GMT File #
wp.hs.gmt <- rWikiPathways::readPathwayGMT("WP/wikipathways-20240510-gmt-Homo_sapiens.gmt")

# Create Necessary DataFrames for Enricher #
wpid2gene <- dplyr::select(wp.hs.gmt, wpid, gene) #TERM2GENE
wpid2name <- dplyr::select(wp.hs.gmt, wpid, name) #TERM2NAME

##---Alveolar Macrophages - geneList Example---##

# Read CSV #
AM.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/AM_wilcox_[Redacted]_SCT.csv", header = TRUE)

# Filter LogFC #
AM.Markers <- AM.Markers[AM.Markers$avg_log2FC > 1.5,]

# Filter P-values #
AM.Markers <- AM.Markers[AM.Markers$p_val < 0.01,]

# Rename Column 1 #
colnames(AM.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- AM.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
AM.Markers <- inner_join(genes_to_test, AM.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# Create geneList for clusterProfiler
geneList <- AM.Markers$avg_log2FC 
names(geneList)<- as.character(AM.Markers$ENTREZID)
geneList <-geneList[!duplicated(names(geneList))] # Remove Duplicates
geneList<- sort(geneList,decreasing=TRUE) # Sort Decreasing


####-------Pathway Analysis-------####

##---Alveolar Macrophages - WikiPathways---##

# Read CSV #
AM.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/AM_wilcox_[Redacted]_SCT.csv", header = TRUE)

# Rename Column 1 #
colnames(AM.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- AM.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
AM.Markers <- inner_join(genes_to_test, AM.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
AM.WP.Up.Genes <- AM.Markers[AM.Markers$avg_log2FC > 1 & AM.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
AM.WP.Bkgd.Genes <- AM.Markers[["ENTREZID"]]

# Run Enrich #
AM.WP.UP <- clusterProfiler::enricher(AM.WP.Up.Genes, universe = AM.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(AM.WP.UP, showCategory = 10) +
  ggtitle("AM PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Alveolar Macrophages - CompareCluster---##

# Read CSV #
AM.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/AM_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(AM.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- AM.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
AM.Markers.Mother <- inner_join(genes_to_test.Mother, AM.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
AM.CC.Up.Genes.Mother <- AM.Markers.Mother[AM.Markers.Mother$avg_log2FC > 1 & AM.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
AM.CC.Down.Genes.Mother <- AM.Markers.Mother[AM.Markers.Mother$avg_log2FC < -1 & AM.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
AM.CC.List <- list(PCDvMOM.UP = AM.CC.Up.Genes.Mother, PCDvMOM.DOWN = AM.CC.Down.Genes.Mother)

# Run Compare Cluster #
AM.CC <- compareCluster(geneCluster = AM.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(AM.CC, showCategory = 5) +
  ggtitle("Alveolar Macrophage") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---B Cells - WikiPathways---##

# Read CSV #
BCells.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/B_Cells_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(BCells.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- BCells.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
BCells.Markers <- inner_join(genes_to_test, BCells.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
BCells.WP.Up.Genes <- BCells.Markers[BCells.Markers$avg_log2FC > 1 & BCells.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
BCells.WP.Bkgd.Genes <- BCells.Markers[["ENTREZID"]]

# Run Enrich #
BCells.WP.UP <- clusterProfiler::enricher(BCells.WP.Up.Genes, universe = BCells.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(BCells.WP.UP, showCategory = 10) +
  ggtitle("B Cells PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---B Cells - CompareCluster---##

# Read CSV #
BCells.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/B_Cells_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(BCells.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- BCells.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
BCells.Markers.Mother <- inner_join(genes_to_test.Mother, BCells.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
BCells.CC.Up.Genes.Mother <- BCells.Markers.Mother[BCells.Markers.Mother$avg_log2FC > 1 & BCells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
BCells.CC.Down.Genes.Mother <- BCells.Markers.Mother[BCells.Markers.Mother$avg_log2FC < -1 & BCells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
BCells.CC.List <- list(PCDvMOM.UP = BCells.CC.Up.Genes.Mother, PCDvMOM.DOWN = BCells.CC.Down.Genes.Mother)

# Run Compare Cluster #
BCells.CC <- compareCluster(geneCluster = BCells.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(BCells.CC, showCategory = 5) +
  ggtitle("B Cells") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Basal Cells - WikiPathways---##

# Read CSV #
BasalCells.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Basal_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(BasalCells.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- BasalCells.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
BasalCells.Markers <- inner_join(genes_to_test, BasalCells.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
BasalCells.WP.Up.Genes <- BasalCells.Markers[BasalCells.Markers$avg_log2FC > 1 & BasalCells.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
BasalCells.WP.Bkgd.Genes <- BasalCells.Markers[["ENTREZID"]]

# Run Enrich #
BasalCells.WP.UP <- clusterProfiler::enricher(BasalCells.WP.Up.Genes, universe = BasalCells.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(BasalCells.WP.UP, showCategory = 10) +
  ggtitle("Basal Cells PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Basal Cells - CompareCluster---##

# Read CSV #
BasalCells.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Basal_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(BasalCells.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- BasalCells.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
BasalCells.Markers.Mother <- inner_join(genes_to_test.Mother, BasalCells.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
BasalCells.CC.Up.Genes.Mother <- BasalCells.Markers.Mother[BasalCells.Markers.Mother$avg_log2FC > 1 & BasalCells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
BasalCells.CC.Down.Genes.Mother <- BasalCells.Markers.Mother[BasalCells.Markers.Mother$avg_log2FC < -1 & BasalCells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
BasalCells.CC.List <- list(PCDvMOM.UP = BasalCells.CC.Up.Genes.Mother, PCDvMOM.DOWN = BasalCells.CC.Down.Genes.Mother)

# Run Compare Cluster #
BasalCells.CC <- compareCluster(geneCluster = BasalCells.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(BasalCells.CC, showCategory = 5) +
  ggtitle("Basal Cells") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Ciliated - WikiPathways---##

# Read CSV #
Ciliated.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Ciliated_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Ciliated.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Ciliated.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Ciliated.Markers <- inner_join(genes_to_test, Ciliated.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Ciliated.WP.Up.Genes <- Ciliated.Markers[Ciliated.Markers$avg_log2FC > 1 & Ciliated.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Ciliated.WP.Bkgd.Genes <- Ciliated.Markers[["ENTREZID"]]

# Run Enrich #
Ciliated.WP.UP <- clusterProfiler::enricher(Ciliated.WP.Up.Genes, universe = Ciliated.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Ciliated.WP.UP, showCategory = 10) +
  ggtitle("Ciliated PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Ciliated - CompareCluster---##

# Read CSV #
Ciliated.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Ciliated_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Ciliated.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Ciliated.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
Ciliated.Markers.Mother <- inner_join(genes_to_test.Mother, Ciliated.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
Ciliated.CC.Up.Genes.Mother <- Ciliated.Markers.Mother[Ciliated.Markers.Mother$avg_log2FC > 1 & Ciliated.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Ciliated.CC.Down.Genes.Mother <- Ciliated.Markers.Mother[Ciliated.Markers.Mother$avg_log2FC < -1 & Ciliated.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
Ciliated.CC.List <- list(PCDvMOM.UP = Ciliated.CC.Up.Genes.Mother, PCDvMOM.DOWN = Ciliated.CC.Down.Genes.Mother)

# Run Compare Cluster #
Ciliated.CC <- compareCluster(geneCluster = Ciliated.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Ciliated.CC, showCategory = 5) +
  ggtitle("Ciliated") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---cMono - WikiPathways---##

# Read CSV #
cMono.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/cMono_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(cMono.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- cMono.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
cMono.Markers <- inner_join(genes_to_test, cMono.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
cMono.WP.Up.Genes <- cMono.Markers[cMono.Markers$avg_log2FC > 1 & cMono.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
cMono.WP.Bkgd.Genes <- cMono.Markers[["ENTREZID"]]

# Run Enrich #
cMono.WP.UP <- clusterProfiler::enricher(cMono.WP.Up.Genes, universe = cMono.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(cMono.WP.UP, showCategory = 10) +
  ggtitle("cMono PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---cMono - CompareCluster---##

# Read CSV #
cMono.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/cMono_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(cMono.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- cMono.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
cMono.Markers.Mother <- inner_join(genes_to_test.Mother, cMono.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
cMono.CC.Up.Genes.Mother <- cMono.Markers.Mother[cMono.Markers.Mother$avg_log2FC > 1 & cMono.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
cMono.CC.Down.Genes.Mother <- cMono.Markers.Mother[cMono.Markers.Mother$avg_log2FC < -1 & cMono.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
cMono.CC.List <- list(PCDvMOM.UP = cMono.CC.Up.Genes.Mother, PCDvMOM.DOWN = cMono.CC.Down.Genes.Mother)

# Run Compare Cluster #
cMono.CC <- compareCluster(geneCluster = cMono.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(cMono.CC, showCategory = 5) +
  ggtitle("cMono") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Div - WikiPathways---##

# Read CSV #
Div.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Div_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Div.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Div.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Div.Markers <- inner_join(genes_to_test, Div.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Div.WP.Up.Genes <- Div.Markers[Div.Markers$avg_log2FC > 1 & Div.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Div.WP.Bkgd.Genes <- Div.Markers[["ENTREZID"]]

# Run Enrich #
Div.WP.UP <- clusterProfiler::enricher(Div.WP.Up.Genes, universe = Div.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Div.WP.UP, showCategory = 10) +
  ggtitle("Div PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Div - CompareCluster---##

# Read CSV #
Div.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Div_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Div.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Div.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
Div.Markers.Mother <- inner_join(genes_to_test.Mother, Div.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
Div.CC.Up.Genes.Mother <- Div.Markers.Mother[Div.Markers.Mother$avg_log2FC > 1 & Div.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Div.CC.Down.Genes.Mother <- Div.Markers.Mother[Div.Markers.Mother$avg_log2FC < -1 & Div.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
Div.CC.List <- list(PCDvMOM.UP = Div.CC.Up.Genes.Mother, PCDvMOM.DOWN = Div.CC.Down.Genes.Mother)

# Run Compare Cluster #
Div.CC <- compareCluster(geneCluster = Div.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Div.CC, showCategory = 5) +
  ggtitle("Divider") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---IM - WikiPathways---##

# Read CSV #
IM.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/IM_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(IM.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- IM.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
IM.Markers <- inner_join(genes_to_test, IM.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
IM.WP.Up.Genes <- IM.Markers[IM.Markers$avg_log2FC > 1 & IM.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
IM.WP.Bkgd.Genes <- IM.Markers[["ENTREZID"]]

# Run Enrich #
IM.WP.UP <- clusterProfiler::enricher(IM.WP.Up.Genes, universe = IM.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(IM.WP.UP, showCategory = 10) +
  ggtitle("IM PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---IM - CompareCluster---##

# Read CSV #
IM.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/IM_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(IM.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- IM.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
IM.Markers.Mother <- inner_join(genes_to_test.Mother, IM.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
IM.CC.Up.Genes.Mother <- IM.Markers.Mother[IM.Markers.Mother$avg_log2FC > 1 & IM.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
IM.CC.Down.Genes.Mother <- IM.Markers.Mother[IM.Markers.Mother$avg_log2FC < -1 & IM.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
IM.CC.List <- list(PCDvMOM.UP = IM.CC.Up.Genes.Mother, PCDvMOM.DOWN = IM.CC.Down.Genes.Mother)

# Run Compare Cluster #
IM.CC <- compareCluster(geneCluster = IM.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(IM.CC, showCategory = 5) +
  ggtitle("IM") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---ncMono - WikiPathways---##

# Read CSV #
ncMono.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/ncMono_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(ncMono.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- ncMono.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
ncMono.Markers <- inner_join(genes_to_test, ncMono.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
ncMono.WP.Up.Genes <- ncMono.Markers[ncMono.Markers$avg_log2FC > 1 & ncMono.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
ncMono.WP.Bkgd.Genes <- ncMono.Markers[["ENTREZID"]]

# Run Enrich #
ncMono.WP.UP <- clusterProfiler::enricher(ncMono.WP.Up.Genes, universe = ncMono.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(ncMono.WP.UP, showCategory = 10) +
  ggtitle("ncMono PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---ncMono - CompareCluster---##

# Read CSV #
ncMono.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/ncMono_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(ncMono.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- ncMono.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
ncMono.Markers.Mother <- inner_join(genes_to_test.Mother, ncMono.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
ncMono.CC.Up.Genes.Mother <- ncMono.Markers.Mother[ncMono.Markers.Mother$avg_log2FC > 1 & ncMono.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
ncMono.CC.Down.Genes.Mother <- ncMono.Markers.Mother[ncMono.Markers.Mother$avg_log2FC < -1 & ncMono.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
ncMono.CC.List <- list(PCDvMOM.UP = ncMono.CC.Up.Genes.Mother, PCDvMOM.DOWN = ncMono.CC.Down.Genes.Mother)

# Run Compare Cluster #
ncMono.CC <- compareCluster(geneCluster = ncMono.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(ncMono.CC, showCategory = 5) +
  ggtitle("ncMono") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---NK - WikiPathways---##

# Read CSV #
NK.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/NK_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(NK.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- NK.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
NK.Markers <- inner_join(genes_to_test, NK.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
NK.WP.Up.Genes <- NK.Markers[NK.Markers$avg_log2FC > 1 & NK.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
NK.WP.Bkgd.Genes <- NK.Markers[["ENTREZID"]]

# Run Enrich #
NK.WP.UP <- clusterProfiler::enricher(NK.WP.Up.Genes, universe = NK.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(NK.WP.UP, showCategory = 10) +
  ggtitle("NK PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---NK - CompareCluster---##

# Read CSV #
NK.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/NK_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(NK.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- NK.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
NK.Markers.Mother <- inner_join(genes_to_test.Mother, NK.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
NK.CC.Up.Genes.Mother <- NK.Markers.Mother[NK.Markers.Mother$avg_log2FC > 1 & NK.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
NK.CC.Down.Genes.Mother <- NK.Markers.Mother[NK.Markers.Mother$avg_log2FC < -1 & NK.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
NK.CC.List <- list(PCDvMOM.UP = NK.CC.Up.Genes.Mother, PCDvMOM.DOWN = NK.CC.Down.Genes.Mother)

# Run Compare Cluster #
NK.CC <- compareCluster(geneCluster = NK.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(NK.CC, showCategory = 5) +
  ggtitle("NK") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---PMN - WikiPathways---##

# Read CSV #
PMN.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/PMN_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(PMN.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- PMN.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
PMN.Markers <- inner_join(genes_to_test, PMN.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
PMN.WP.Up.Genes <- PMN.Markers[PMN.Markers$avg_log2FC > 1 & PMN.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
PMN.WP.Bkgd.Genes <- PMN.Markers[["ENTREZID"]]

# Run Enrich #
PMN.WP.UP <- clusterProfiler::enricher(PMN.WP.Up.Genes, universe = PMN.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(PMN.WP.UP, showCategory = 10) +
  ggtitle("PMN PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---PMN - CompareCluster---##

# Read CSV #
PMN.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/PMN_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(PMN.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- PMN.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
PMN.Markers.Mother <- inner_join(genes_to_test.Mother, PMN.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
PMN.CC.Up.Genes.Mother <- PMN.Markers.Mother[PMN.Markers.Mother$avg_log2FC > 1 & PMN.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
PMN.CC.Down.Genes.Mother <- PMN.Markers.Mother[PMN.Markers.Mother$avg_log2FC < -1 & PMN.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
PMN.CC.List <- list(PCDvMOM.UP = PMN.CC.Up.Genes.Mother, PCDvMOM.DOWN = PMN.CC.Down.Genes.Mother)

# Run Compare Cluster #
PMN.CC <- compareCluster(geneCluster = PMN.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(PMN.CC, showCategory = 5) +
  ggtitle("PMN") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Secretory - WikiPathways---##

# Read CSV #
Secretory.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Secretory_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Secretory.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Secretory.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Secretory.Markers <- inner_join(genes_to_test, Secretory.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Secretory.WP.Up.Genes <- Secretory.Markers[Secretory.Markers$avg_log2FC > 1 & Secretory.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Secretory.WP.Bkgd.Genes <- Secretory.Markers[["ENTREZID"]]

# Run Enrich #
Secretory.WP.UP <- clusterProfiler::enricher(Secretory.WP.Up.Genes, universe = Secretory.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Secretory.WP.UP, showCategory = 10) +
  ggtitle("Secretory PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Secretory - CompareCluster---##

# Read CSV #
Secretory.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Secretory_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Secretory.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Secretory.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
Secretory.Markers.Mother <- inner_join(genes_to_test.Mother, Secretory.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
Secretory.CC.Up.Genes.Mother <- Secretory.Markers.Mother[Secretory.Markers.Mother$avg_log2FC > 1 & Secretory.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Secretory.CC.Down.Genes.Mother <- Secretory.Markers.Mother[Secretory.Markers.Mother$avg_log2FC < -1 & Secretory.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
Secretory.CC.List <- list(PCDvMOM.UP = Secretory.CC.Up.Genes.Mother, PCDvMOM.DOWN = Secretory.CC.Down.Genes.Mother)

# Run Compare Cluster #
Secretory.CC <- compareCluster(geneCluster = Secretory.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Secretory.CC, showCategory = 5) +
  ggtitle("Secretory") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---T Cells - WikiPathways---##

# Read CSV #
T_Cells.Markers <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/T_Cells_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(T_Cells.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- T_Cells.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
T_Cells.Markers <- inner_join(genes_to_test, T_Cells.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
T_Cells.WP.Up.Genes <- T_Cells.Markers[T_Cells.Markers$avg_log2FC > 1 & T_Cells.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
T_Cells.WP.Bkgd.Genes <- T_Cells.Markers[["ENTREZID"]]

# Run Enrich #
T_Cells.WP.UP <- clusterProfiler::enricher(T_Cells.WP.Up.Genes, universe = T_Cells.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(T_Cells.WP.UP, showCategory = 10) +
  ggtitle("T Cells PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---T Cells - CompareCluster---##

# Read CSV #
T_Cells.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/T_Cells_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(T_Cells.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- T_Cells.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
T_Cells.Markers.Mother <- inner_join(genes_to_test.Mother, T_Cells.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
T_Cells.CC.Up.Genes.Mother <- T_Cells.Markers.Mother[T_Cells.Markers.Mother$avg_log2FC > 1 & T_Cells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
T_Cells.CC.Down.Genes.Mother <- T_Cells.Markers.Mother[T_Cells.Markers.Mother$avg_log2FC < -1 & T_Cells.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
T_Cells.CC.List <- list(PCDvMOM.UP = T_Cells.CC.Up.Genes.Mother, PCDvMOM.DOWN = T_Cells.CC.Down.Genes.Mother)

# Run Compare Cluster #
T_Cells.CC <- compareCluster(geneCluster = T_Cells.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(T_Cells.CC, showCategory = 5) +
  ggtitle("T Cells") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


####-------DEG Analysis - Immune Cells-------####

# Load Object #
integrated_object <- readRDS("Objects/merged_object_SCT_Integrated_SR004328_annotated.rds")

# Join Layers #
integrated_object[["RNA"]] <- JoinLayers(integrated_object[["RNA"]])

# Create New Column for PCD vs Mother #
integrated_object$stim5 <- integrated_object$sampleid

# Reformat Cluster Names #
integrated_object$stim5 <- as.character(integrated_object$stim5)

# Change Metadata Elements #
integrated_object$stim5[integrated_object$stim5 %in% c("PCD171", "PCD172")] <- paste0("PCD")
integrated_object$stim5[integrated_object$stim5 %in% 'PCD171M'] <- paste0('Mother')
integrated_object$stim5[integrated_object$stim5 %in% 'PCD177'] <- paste0('Possible_Control')

# Reformat Cluster Names #
integrated_object$stim5 <- as.factor(integrated_object$stim5)

# Group Immune CLuster Names #
integrated_object$inames <- integrated_object$final_cluster
integrated_object$inames <- as.character(integrated_object$inames)
integrated_object$inames[integrated_object$inames %in% 'AM'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'B Cells'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'cMono'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'GB Cells'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'IM'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'ncMono'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'NK'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'PMN'] <- paste0('Immune')
integrated_object$inames[integrated_object$inames %in% 'T Cells'] <- paste0('Immune')
integrated_object$inames <- as.factor(integrated_object$inames)

# Add Column Discerning Control v PCD #
integrated_object$stim3 <- paste(integrated_object$inames, integrated_object$stim5, sep='_')

# Set Idents #
Idents(integrated_object) <- integrated_object$stim3

# Run PrepSCTFindMarkers #
integrated_object <- PrepSCTFindMarkers(integrated_object)


####-------DE Analysis - Immune (PCD vs Mother) - Wilcox - SCT -------####

# Run FindMarkers #
name = FindMarkers(integrated_object, ident.1 = "Immune_PCD", ident.2 = "Immune_Mother", verbose = TRUE)
file = paste0("PCD", "_wilcox_", "[Redacted]_", "SCT", "_immune",  ".csv")
write.csv(name, file)


####-------DE Analysis - Immune (PCD vs Mother) - MAST - SCT-------####

# Run FindMarkers #
name = FindMarkers(integrated_object, ident.1 = "Immune_PCD", ident.2 = "Immune_Mother", verbose = TRUE, test.use = "MAST")
file = paste0("PCD", "_MAST_", "[Redacted]_", "SCT", "_immune",  ".csv")
write.csv(name, file)


####-------DE Analysis - Immune (PCD vs Potential Control) - Wilcox - SCT -------####

# Run FindMarkers #
name = FindMarkers(integrated_object, ident.1 = "Immune_PCD", ident.2 = "Immune_Possible_Control", verbose = TRUE)
file = paste0("PCD", "_wilcox_", "[Redacted]_", "possible_control_", "SCT", "_immune",  ".csv")
write.csv(name, file)


####-------DE Analysis - Immune (PCD vs Potential Control) - MAST - SCT-------####

# Run FindMarkers #
name = FindMarkers(integrated_object, ident.1 = "Immune_PCD", ident.2 = "Immune_Possible_Control", verbose = TRUE, test.use = "MAST")
file = paste0("PCD", "_MAST_", "[Redacted]_", "possible_control_", "SCT", "_immune",  ".csv")
write.csv(name, file)


####-------Pathway Analysis - Set-Up - Immune-------####

# Install Packages #
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("rWikiPathways")
library(clusterProfiler)
library("org.Hs.eg.db")
library(dplyr)
library(tidyverse)
library(rWikiPathways)

# Download GMT File #
wp.hs.gmt <- rWikiPathways::readPathwayGMT("WP/wikipathways-20240510-gmt-Homo_sapiens.gmt")

# Create Necessary DataFrames for Enricher #
wpid2gene <- dplyr::select(wp.hs.gmt, wpid, gene) #TERM2GENE
wpid2name <- dplyr::select(wp.hs.gmt, wpid, name) #TERM2NAME



####-------Pathway Analysis - Immune-------####

##---Immune - WikiPathways---##

# Read CSV #
Immune.Markers <- read.csv("Differential_Expression/Immune/Mother/Wilcoxon/PCD_wilcox_[Redacted]_SCT_immune.csv", header = TRUE)

# Rename Column 1 #
colnames(Immune.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Immune.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Immune.Markers <- inner_join(genes_to_test, Immune.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Immune.WP.Up.Genes <- Immune.Markers[Immune.Markers$avg_log2FC > 1 & Immune.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Immune.WP.Bkgd.Genes <- Immune.Markers[["ENTREZID"]]

# Run Enrich #
Immune.WP.UP <- clusterProfiler::enricher(Immune.WP.Up.Genes, universe = Immune.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Immune.WP.UP, showCategory = 10) +
  ggtitle("Immune PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Immune - CompareCluster---##

# Read CSV #
Immune.Markers.Mother <- read.csv("Differential_Expression/Immune/Mother/Wilcoxon/PCD_wilcox_[Redacted]_SCT_immune.csv", header = TRUE) # Mother

# Rename Column 1 #
colnames(Immune.Markers.Mother)[1] <- "GeneNames" # Mother

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Immune.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

# Join Both DataFrames by GeneNames #
Immune.Markers.Mother <- inner_join(genes_to_test.Mother, Immune.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother

# List Upregulated Genes #
Immune.CC.Up.Genes.Mother <- Immune.Markers.Mother[Immune.Markers.Mother$avg_log2FC > 1 & Immune.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Immune.CC.Down.Genes.Mother <- Immune.Markers.Mother[Immune.Markers.Mother$avg_log2FC < -1 & Immune.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

# Create Named List #
Immune.CC.List <- list(PCDvMOM.UP = Immune.CC.Up.Genes.Mother, PCDvMOM.DOWN = Immune.CC.Down.Genes.Mother)

# Run Compare Cluster #
Immune.CC <- compareCluster(geneCluster = Immune.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Immune.CC, showCategory = 5) +
  ggtitle("Immune") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


####-------DEG Analysis - Airway Cells-------####

# Load Object #
integrated_object <- readRDS("Objects/merged_object_SCT_Integrated_SR004328_annotated.rds")

# Join Layers #
integrated_object[["RNA"]] <- JoinLayers(integrated_object[["RNA"]])

# Create New Column for PCD vs Mother #
integrated_object$stim5 <- integrated_object$sampleid

# Reformat Cluster Names #
integrated_object$stim5 <- as.character(integrated_object$stim5)

# Change Metadata Elements #
integrated_object$stim5[integrated_object$stim5 %in% c("PCD171", "PCD172")] <- paste0("PCD")
integrated_object$stim5[integrated_object$stim5 %in% 'PCD171M'] <- paste0('Mother')
integrated_object$stim5[integrated_object$stim5 %in% 'PCD177'] <- paste0('Possible_Control')

# Reformat Cluster Names #
integrated_object$stim5 <- as.factor(integrated_object$stim5)

# Group Immune CLuster Names #
integrated_object$inames <- integrated_object$final_cluster
integrated_object$inames <- as.character(integrated_object$inames)
integrated_object$inames[integrated_object$inames %in% 'Ciliated'] <- paste0('Airway')
integrated_object$inames[integrated_object$inames %in% 'Basal'] <- paste0('Airway')
integrated_object$inames[integrated_object$inames %in% 'Secretory'] <- paste0('Airway')
integrated_object$inames <- as.factor(integrated_object$inames)

# Add Column Discerning Control v PCD #
integrated_object$stim3 <- paste(integrated_object$inames, integrated_object$stim5, sep='_')

# Set Idents #
Idents(integrated_object) <- integrated_object$stim3

# Run PrepSCTFindMarkers #
integrated_object <- PrepSCTFindMarkers(integrated_object)


####-------DE Analysis - Immune (PCD vs Mother) - Wilcox - SCT -------####

# Run FindMarkers #
name = FindMarkers(integrated_object, ident.1 = "Airway_PCD", ident.2 = "Airway_Mother", verbose = TRUE)
file = paste0("PCD", "_wilcox_", "[Redacted]_", "SCT", "_airway",  ".csv")
write.csv(name, file)


####-------DNAH5 - Prebuilt - Set-Up-------####

# Load Raw Data #
object.control <- readRDS("Objects/Rpca_13_SINGLETS_final.rds")

# Subset out PCD and extra controls #
object.control <- subset(object.control, subset = stim2 %in% c("HNEC102", "HNEC103"))

# Remove Integrated Assay #
Assays(object.control)
object.control[['integrated']] <- NULL

# Change Column Name #
object.control$sampleid <- object.control$stim2

####-------[Redacted] Object Set-Up-------####

# Load Object #
object.pcd <- readRDS("Objects/merged_object_SCT_Integrated_SR004328_annotated.rds")

# Join Layers #
object.pcd[["RNA"]] <- JoinLayers(object.pcd[["RNA"]])

# Subset out PCD and extra controls #
object.pcd <- subset(object.pcd, subset = sampleid %in% c("PCD171", "PCD172"))

# Subset out Immune Cells #
object.pcd <- subset(object.pcd, subset = final_cluster %in% c("AM", "T Cells", "B Cells", "PMN", "cMono", "GB Cells", "IM", "ncMono", "NK"), invert = TRUE)

# Remove Integrated Assay #
Assays(object.pcd)
object.pcd[['SCT']] <- NULL


####-------Merge-------####

# Merge Objects#
merged_object <- merge(object.control, y = c(object.pcd), 
                       add.cell.ids = c("control", "pcd"), project = "[Redacted]")

# Save #
saveRDS(merged_object, "merged_object_pcd_HNEC102_HNEC103.rds") 


####-------Update Meta Data-------####

# Create New Column for PCD vs Mother #
merged_object$stim6 <- integrated_object$stim2

# Reformat Cluster Names #
integrated_object$stim5 <- as.character(integrated_object$stim5)

# Change Metadata Elements #
integrated_object$stim5[integrated_object$stim5 %in% c("PCD171", "PCD172")] <- paste0("PCD")
integrated_object$stim5[integrated_object$stim5 %in% 'PCD171M'] <- paste0('Mother')
integrated_object$stim5[integrated_object$stim5 %in% 'PCD177'] <- paste0('Possible_Control')

# Reformat Cluster Names #
integrated_object$stim5 <- as.factor(integrated_object$stim5)


####-------SCTransform-------####

# Load Object #
merged_object <- readRDS("Objects/merged_object_pcd_HNEC102_HNEC103.rds")

# Join Layers #
merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])

# Split Object by Individual Donors #
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f = merged_object$sampleid)

# Install glmGamPoi for More Efficient Implementation #
# BiocManager::install("glmGamPoi")
# library(glmGamPoi)

# Run SCTransform
merged_object_SCT <- SCTransform(merged_object)

# Run PCA #
merged_object_SCT <- RunPCA(merged_object_SCT)

# Run UMAP #
merged_object_SCT <- RunUMAP(merged_object_SCT, dims = 1:30)

# Visualize Unintegrated UMAP #
DimPlot(merged_object_SCT, reduction = "umap", group.by = c("sampleid"))


####-------CCA Integration-------####

# Perform CCA Integration #
merged_object_SCT_Integrated <- IntegrateLayers(object = merged_object_SCT, method = CCAIntegration, normalization.method = "SCT", verbose = F)

# Run Find Neighbours #
merged_object_SCT_Integrated <- FindNeighbors(merged_object_SCT_Integrated, reduction = "integrated.dr", dims = 1:30)

# Find Clusters #
merged_object_SCT_Integrated <- FindClusters(merged_object_SCT_Integrated, resolution = 1.0)

# Run UMAP #
merged_object_SCT_Integrated <- RunUMAP(merged_object_SCT_Integrated, dims = 1:30, reduction = "integrated.dr")

# Save #
saveRDS(merged_object_SCT_Integrated, "merged_object_SCT_Integrated_[Redacted]_DNAH5.rds")

# Visualize UMAP #
DimPlot(merged_object_SCT_Integrated, reduction = "umap", group.by = c("seurat_clusters"))
DimPlot(merged_object_SCT_Integrated, reduction = "umap", label = TRUE)

merged_object_SCT_Integrated[["RNA"]] <- JoinLayers(merged_object_SCT_Integrated[["RNA"]])

####-------Cluster Annotation-------####

# Load Object #
merged_object_SCT_Integrated <- readRDS("Objects/merged_object_SCT_Integrated_[Redacted]_DNAH5.rds")

# Change Default Assay #
DefaultAssay(merged_object_SCT_Integrated) <- 'RNA'

# Run FindAllMarkers #
All.Markers <- FindAllMarkers(merged_object_SCT_Integrated, only.pos = TRUE, max.cells.per.ident = 1000)

# Save as CSV #
write.csv(All.Markers, "[Redacted]_DNAH5_markers.csv")

# Pertinent Markers #
FeaturePlot(merged_object_SCT_Integrated,"JUNB") # Stress/Ribosomal Contamination
FeaturePlot(merged_object_SCT_Integrated, "nFeature_RNA") # RNA Levels
FeaturePlot(merged_object_SCT_Integrated, "DNAH5", label = TRUE) # Ciliated
FeaturePlot(merged_object_SCT_Integrated, "TOP2A", label = TRUE) # Div
FeaturePlot(merged_object_SCT_Integrated, "KRT5", label = TRUE) # Basal
FeaturePlot(merged_object_SCT_Integrated, "CYP2F1", label = TRUE) # Secretory


# Set Assay to RNA #
DefaultAssay(merged_object_SCT_Integrated) <- "RNA"

# Rename Idents #
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "0" = "Basal") 
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "1" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "2" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "3" = "Basal")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "4" = "Basal")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "5" = "Basal")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "6" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "7" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "8" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "9" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "10" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "11" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "12" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "13" = "Secretory") 
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "14" = "Basal")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "15" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "16" = "Div")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "17" = "Ciliated")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "18" = "Secretory")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "19" = "Secretory")

# Create final_cluster column #
merged_object_SCT_Integrated$final_cluster <- Idents(merged_object_SCT_Integrated)

# Visualize
DimPlot(merged_object_SCT_Integrated, label = TRUE)


####-------DE Analysis - [Redacted]/DNAH5 (PCD vs Control) - Set-Up-------####

# Add Column Discerning Control v PCD #
merged_object_SCT_Integrated$stim10 <- paste(merged_object_SCT_Integrated$final_cluster, merged_object_SCT_Integrated$stim, sep='_')

# Set Idents #
Idents(merged_object_SCT_Integrated) <- merged_object_SCT_Integrated$stim10

# Create List of Cell Types #
pcd_list <- levels(merged_object_SCT_Integrated)[5:8]
pcd_list <- pcd_list[order(names(setNames(pcd_list, pcd_list)))]

control_list <- levels(merged_object_SCT_Integrated)[1:4]
control_list <- mother_list[order(names(setNames(mother_list, mother_list)))]


# Run PrepSCTFindMarkers #
merged_object_SCT_Integrated <- PrepSCTFindMarkers(merged_object_SCT_Integrated)


####-------DE Analysis - DNAH5 (PCD vs Control) - Wilcox - SCT -------####

# Run FindMarkers #
for (i in 1:4) {
  name = FindMarkers(merged_object_SCT_Integrated, ident.1 = pcd_list[[i]], ident.2 = control_list[[i]], verbose = TRUE)
  file = paste0(unlist((strsplit(control_list[[i]], split = "_")))[1], "_wilcox_", "[Redacted]_", "SCT", "DNAH5_control", ".csv")
  write.csv(name, file)
}


####-------DE Analysis - DNAH5 (PCD vs Control) - MAST - SCT-------####

# Run FindMarkers #
for (i in 1:4) {
  name = FindMarkers(merged_object_SCT_Integrated, ident.1 = pcd_list[[i]], ident.2 = control_list[[i]], verbose = TRUE, test.use = "MAST")
  file = paste0(unlist((strsplit(control_list[[i]], split = "_")))[1], "_MAST_", "[Redacted]_", "SCT", "DNAH5_control", ".csv")
  write.csv(name, file)
}


####-------Pathway Analysis - Set-Up-------####

# Install Packages #
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("rWikiPathways")
library(clusterProfiler)
library("org.Hs.eg.db")
library(dplyr)
library(tidyverse)
library(rWikiPathways)

# Download GMT File #
wp.hs.gmt <- rWikiPathways::readPathwayGMT("WP/wikipathways-20240510-gmt-Homo_sapiens.gmt")

# Create Necessary DataFrames for Enricher #
wpid2gene <- dplyr::select(wp.hs.gmt, wpid, gene) #TERM2GENE
wpid2name <- dplyr::select(wp.hs.gmt, wpid, name) #TERM2NAME


####-------Pathway Analysis - DNAH5 Control-------####

##---Basal - WikiPathways---##

# Read CSV #
Basal.Markers <- read.csv("Differential_Expression/Control/Wilcoxon/Basal_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE)

# Rename Column 1 #
colnames(Basal.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Basal.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Basal.Markers <- inner_join(genes_to_test, Basal.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Basal.WP.Up.Genes <- Basal.Markers[Basal.Markers$avg_log2FC > 1 & Basal.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Basal.WP.Bkgd.Genes <- Basal.Markers[["ENTREZID"]]

# Run Enrich #
Basal.WP.UP <- clusterProfiler::enricher(Basal.WP.Up.Genes, universe = Basal.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Basal.WP.UP, showCategory = 10) +
  ggtitle("Basal PCD vs Control WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Basal - CompareCluster---##

# Read CSV #
Basal.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Basal_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother
Basal.Markers.Control <- read.csv("Differential_Expression/Control/Wilcoxon/Basal_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE) # Possible Control

# Rename Column 1 #
colnames(Basal.Markers.Mother)[1] <- "GeneNames" # Mother
colnames(Basal.Markers.Control)[1] <- "GeneNames" # Possible Control

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Basal.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

genes_to_test.Control <- Basal.Markers.Control$GeneNames # Possible Control
genes_to_test.Control <- bitr(genes_to_test.Control, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Possible Control
colnames(genes_to_test.Control)[1]<-"GeneNames" # Possible Control

# Join Both DataFrames by GeneNames #
Basal.Markers.Mother <- inner_join(genes_to_test.Mother, Basal.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother
Basal.Markers.Control <- inner_join(genes_to_test.Control, Basal.Markers.Control, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Possible Control

# List Upregulated Genes #
Basal.CC.Up.Genes.Mother <- Basal.Markers.Mother[Basal.Markers.Mother$avg_log2FC > 1 & Basal.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Basal.CC.Down.Genes.Mother <- Basal.Markers.Mother[Basal.Markers.Mother$avg_log2FC < -1 & Basal.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

Basal.CC.Up.Genes.Control <- Basal.Markers.Control[Basal.Markers.Control$avg_log2FC > 1 & Basal.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control
Basal.CC.Down.Genes.Control <- Basal.Markers.Control[Basal.Markers.Control$avg_log2FC < -1 & Basal.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control

# Create Named List #
Basal.CC.List <- list(PCDvMOM.UP = Basal.CC.Up.Genes.Mother, PCDvMOM.DOWN = Basal.CC.Down.Genes.Mother, PCDvCTRL.UP = Basal.CC.Up.Genes.Control , PCDvCTRL.DOWN = Basal.CC.Down.Genes.Control)

# Run Compare Cluster #
Basal.CC <- compareCluster(geneCluster = Basal.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Basal.CC, showCategory = 5) +
  ggtitle("Basal") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Ciliated - WikiPathways---##

# Read CSV #
Ciliated.Markers <- read.csv("Differential_Expression/Control/Wilcoxon/Ciliated_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE)

# Rename Column 1 #
colnames(Ciliated.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Ciliated.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Ciliated.Markers <- inner_join(genes_to_test, Ciliated.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Ciliated.WP.Up.Genes <- Ciliated.Markers[Ciliated.Markers$avg_log2FC > 1 & Ciliated.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Ciliated.WP.Bkgd.Genes <- Ciliated.Markers[["ENTREZID"]]

# Run Enrich #
Ciliated.WP.UP <- clusterProfiler::enricher(Ciliated.WP.Up.Genes, universe = Ciliated.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Ciliated.WP.UP, showCategory = 10) +
  ggtitle("Ciliated PCD vs Control WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Ciliated - CompareCluster---##

# Read CSV #
Ciliated.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Ciliated_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother
Ciliated.Markers.Control <- read.csv("Differential_Expression/Control/Wilcoxon/Ciliated_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE) # Possible Control

# Rename Column 1 #
colnames(Ciliated.Markers.Mother)[1] <- "GeneNames" # Mother
colnames(Ciliated.Markers.Control)[1] <- "GeneNames" # Possible Control

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Ciliated.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

genes_to_test.Control <- Ciliated.Markers.Control$GeneNames # Possible Control
genes_to_test.Control <- bitr(genes_to_test.Control, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Possible Control
colnames(genes_to_test.Control)[1]<-"GeneNames" # Possible Control

# Join Both DataFrames by GeneNames #
Ciliated.Markers.Mother <- inner_join(genes_to_test.Mother, Ciliated.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother
Ciliated.Markers.Control <- inner_join(genes_to_test.Control, Ciliated.Markers.Control, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Possible Control

# List Upregulated Genes #
Ciliated.CC.Up.Genes.Mother <- Ciliated.Markers.Mother[Ciliated.Markers.Mother$avg_log2FC > 1 & Ciliated.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Ciliated.CC.Down.Genes.Mother <- Ciliated.Markers.Mother[Ciliated.Markers.Mother$avg_log2FC < -1 & Ciliated.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

Ciliated.CC.Up.Genes.Control <- Ciliated.Markers.Control[Ciliated.Markers.Control$avg_log2FC > 1 & Ciliated.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control
Ciliated.CC.Down.Genes.Control <- Ciliated.Markers.Control[Ciliated.Markers.Control$avg_log2FC < -1 & Ciliated.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control

# Create Named List #
Ciliated.CC.List <- list(PCDvMOM.UP = Ciliated.CC.Up.Genes.Mother, PCDvMOM.DOWN = Ciliated.CC.Down.Genes.Mother, PCDvCTRL.UP = Ciliated.CC.Up.Genes.Control , PCDvCTRL.DOWN = Ciliated.CC.Down.Genes.Control)

# Run Compare Cluster #
Ciliated.CC <- compareCluster(geneCluster = Ciliated.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Ciliated.CC, showCategory = 5) +
  ggtitle("Ciliated") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Secretory - WikiPathways---##

# Read CSV #
Secretory.Markers <- read.csv("Differential_Expression/Control/Wilcoxon/Secretory_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE)

# Rename Column 1 #
colnames(Secretory.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Secretory.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Secretory.Markers <- inner_join(genes_to_test, Secretory.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Secretory.WP.Up.Genes <- Secretory.Markers[Secretory.Markers$avg_log2FC > 1 & Secretory.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Secretory.WP.Bkgd.Genes <- Secretory.Markers[["ENTREZID"]]

# Run Enrich #
Secretory.WP.UP <- clusterProfiler::enricher(Secretory.WP.Up.Genes, universe = Secretory.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Secretory.WP.UP, showCategory = 10) +
  ggtitle("Secretory PCD vs Control WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Secretory - CompareCluster---##

# Read CSV #
Secretory.Markers.Mother <- read.csv("Differential_Expression/Mother/Wilcoxon_SCT/Secretory_wilcox_[Redacted]_SCT.csv", header = TRUE) # Mother
Secretory.Markers.Control <- read.csv("Differential_Expression/Control/Wilcoxon/Secretory_wilcox_[Redacted]_SCTDNAH5_control.csv", header = TRUE) # Possible Control

# Rename Column 1 #
colnames(Secretory.Markers.Mother)[1] <- "GeneNames" # Mother
colnames(Secretory.Markers.Control)[1] <- "GeneNames" # Possible Control

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Secretory.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

genes_to_test.Control <- Secretory.Markers.Control$GeneNames # Possible Control
genes_to_test.Control <- bitr(genes_to_test.Control, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Possible Control
colnames(genes_to_test.Control)[1]<-"GeneNames" # Possible Control

# Join Both DataFrames by GeneNames #
Secretory.Markers.Mother <- inner_join(genes_to_test.Mother, Secretory.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother
Secretory.Markers.Control <- inner_join(genes_to_test.Control, Secretory.Markers.Control, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Possible Control

# List Upregulated Genes #
Secretory.CC.Up.Genes.Mother <- Secretory.Markers.Mother[Secretory.Markers.Mother$avg_log2FC > 1 & Secretory.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Secretory.CC.Down.Genes.Mother <- Secretory.Markers.Mother[Secretory.Markers.Mother$avg_log2FC < -1 & Secretory.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

Secretory.CC.Up.Genes.Control <- Secretory.Markers.Control[Secretory.Markers.Control$avg_log2FC > 1 & Secretory.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control
Secretory.CC.Down.Genes.Control <- Secretory.Markers.Control[Secretory.Markers.Control$avg_log2FC < -1 & Secretory.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control

# Create Named List #
Secretory.CC.List <- list(PCDvMOM.UP = Secretory.CC.Up.Genes.Mother, PCDvMOM.DOWN = Secretory.CC.Down.Genes.Mother, PCDvCTRL.UP = Secretory.CC.Up.Genes.Control , PCDvCTRL.DOWN = Secretory.CC.Down.Genes.Control)

# Run Compare Cluster #
Secretory.CC <- compareCluster(geneCluster = Secretory.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Secretory.CC, showCategory = 5) +
  ggtitle("Secretory") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


####------- DE Airway Analysis - DNAH5 -------####

# Rename Idents #
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Basal" = "Airway") 
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Ciliated" = "Airway")
merged_object_SCT_Integrated <- RenameIdents(merged_object_SCT_Integrated, "Secretory" = "Airway")

# Create final_cluster Column #
merged_object_SCT_Integrated$final_cluster <- Idents(merged_object_SCT_Integrated)

# Visualize
DimPlot(merged_object_SCT_Integrated, label = TRUE)

####-------DE Analysis - Airway - [Redacted]/DNAH5 (PCD vs Control) - Set-Up-------####

# Add Column Discerning Control v PCD #
merged_object_SCT_Integrated$stim10 <- paste(merged_object_SCT_Integrated$final_cluster, merged_object_SCT_Integrated$stim, sep='_')

# Set Idents #
Idents(merged_object_SCT_Integrated) <- merged_object_SCT_Integrated$stim10

# Run PrepSCTFindMarkers #
merged_object_SCT_Integrated <- PrepSCTFindMarkers(merged_object_SCT_Integrated)


####-------DE Analysis - DNAH5 (PCD vs Control) - Wilcox - SCT -------####

# Run FindMarkers #
name = FindMarkers(merged_object_SCT_Integrated, ident.1 = "Airway_PCD", ident.2 = "Airway_CTRL", verbose = TRUE)
file = paste0("PCD", "_wilcox_", "[Redacted]_", "SCT", "DNAH5_control", "_Airway", ".csv")
write.csv(name, file)


##---Airway - WikiPathways - Mother---##

# Read CSV #
Airway.Markers <- read.csv("Differential_Expression/Airway/Mother/Wilcoxon/PCD_wilcox_[Redacted]_SCT_airway.csv", header = TRUE)

# Rename Column 1 #
colnames(Airway.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Airway.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Airway.Markers <- inner_join(genes_to_test, Airway.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Airway.WP.Up.Genes <- Airway.Markers[Airway.Markers$avg_log2FC > 1 & Airway.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Airway.WP.Bkgd.Genes <- Airway.Markers[["ENTREZID"]]

# Run Enrich #
Airway.WP.UP <- clusterProfiler::enricher(Airway.WP.Up.Genes, universe = Airway.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Airway.WP.UP, showCategory = 10) +
  ggtitle("Airway PCD vs Mother WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Airway - WikiPathways - Control---##

# Read CSV #
Airway.Markers <- read.csv("Differential_Expression/Airway/Control/Wilcoxon/PCD_wilcox_[Redacted]_SCTDNAH5_control_Airway.csv", header = TRUE)

# Rename Column 1 #
colnames(Airway.Markers)[1] <- "GeneNames"

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test <- Airway.Markers$GeneNames
genes_to_test <- bitr(genes_to_test, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
colnames(genes_to_test)[1]<-"GeneNames"

# Join Both DataFrames by GeneNames #
Airway.Markers <- inner_join(genes_to_test, Airway.Markers, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything())

# List Upregulated Genes #
Airway.WP.Up.Genes <- Airway.Markers[Airway.Markers$avg_log2FC > 1 & Airway.Markers$p_val < 0.05,][["ENTREZID"]]

# List Background Genes #
Airway.WP.Bkgd.Genes <- Airway.Markers[["ENTREZID"]]

# Run Enrich #
Airway.WP.UP <- clusterProfiler::enricher(Airway.WP.Up.Genes, universe = Airway.WP.Bkgd.Genes, pAdjustMethod = "fdr", TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Visualize #
barplot(Airway.WP.UP, showCategory = 6) +
  ggtitle("Airway PCD vs Control WikiPathway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))


##---Airway - CompareCluster---##

# Read CSV #
Airway.Markers.Mother <- read.csv("Differential_Expression/Airway/Mother/Wilcoxon/PCD_wilcox_[Redacted]_SCT_airway.csv", header = TRUE) # Mother
Airway.Markers.Control <- read.csv("Differential_Expression/Airway/Control/Wilcoxon/PCD_wilcox_[Redacted]_SCTDNAH5_control_Airway.csv", header = TRUE) # Possible Control

# Rename Column 1 #
colnames(Airway.Markers.Mother)[1] <- "GeneNames" # Mother
colnames(Airway.Markers.Control)[1] <- "GeneNames" # Possible Control

# Create New DataFrame of SYMBOLS and Convert to ENTREZIDs #
genes_to_test.Mother <- Airway.Markers.Mother$GeneNames # Mother
genes_to_test.Mother <- bitr(genes_to_test.Mother, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Mother
colnames(genes_to_test.Mother)[1]<-"GeneNames" # Mother

genes_to_test.Control <- Airway.Markers.Control$GeneNames # Possible Control
genes_to_test.Control <- bitr(genes_to_test.Control, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Possible Control
colnames(genes_to_test.Control)[1]<-"GeneNames" # Possible Control

# Join Both DataFrames by GeneNames #
Airway.Markers.Mother <- inner_join(genes_to_test.Mother, Airway.Markers.Mother, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Mother
Airway.Markers.Control <- inner_join(genes_to_test.Control, Airway.Markers.Control, by='GeneNames') %>%  dplyr::select(ENTREZID, avg_log2FC, everything()) # Possible Control

# List Upregulated Genes #
Airway.CC.Up.Genes.Mother <- Airway.Markers.Mother[Airway.Markers.Mother$avg_log2FC > 1 & Airway.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother
Airway.CC.Down.Genes.Mother <- Airway.Markers.Mother[Airway.Markers.Mother$avg_log2FC < -1 & Airway.Markers.Mother$p_val < 0.05,][["ENTREZID"]] # Mother

Airway.CC.Up.Genes.Control <- Airway.Markers.Control[Airway.Markers.Control$avg_log2FC > 1 & Airway.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control
Airway.CC.Down.Genes.Control <- Airway.Markers.Control[Airway.Markers.Control$avg_log2FC < -1 & Airway.Markers.Control$p_val < 0.05,][["ENTREZID"]] # Possible Control

# Create Named List #
Airway.CC.List <- list(PCDvMOM.UP = Airway.CC.Up.Genes.Mother, PCDvMOM.DOWN = Airway.CC.Down.Genes.Mother, PCDvCTRL.UP = Airway.CC.Up.Genes.Control , PCDvCTRL.DOWN = Airway.CC.Down.Genes.Control)

# Run Compare Cluster #
Airway.CC <- compareCluster(geneCluster = Airway.CC.List, fun = enrichWP, organism = "Homo sapien")

# Run DotPlot #
dotplot(Airway.CC, showCategory = 5) +
  ggtitle("Airway") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face = "bold"))

