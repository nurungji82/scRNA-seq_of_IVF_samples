###Leukocytes in human follicular aspirates obtaineid from 4 IVF patients

  
library(devtools)
library(htmltools) 
library(Seurat)
library(SeuratData)
library(patchwork)
library(Matrix)
library(xlsx)
library(umap)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(dittoSeq)
library(tidyverse)
library(devtools)
library(scater)
library(patchwork)
library(SeuratDisk)
library(readr)


## Setup the Seurat Object
### Setup the Seurat Object-IVF1
IVF1.data <- Read10X(data.dir = "D:/scRNA-seq data/Processed data/IVF1/outs/filtered_feature_bc_matrix/")
IVF1.atleastone <- apply(IVF1.data, 2, function(x) sum(x>0))
hist(IVF1.atleastone, breaks = 100,
     main = "Distribution of detected genes - IVF1",
     xlab = "Genes with at least one tag")
IVF1.tmp <- apply(IVF1.data, 1, function(x) sum (x>0))
table(IVF1.tmp>=3)
IVF1.keep <- IVF1.tmp>=3
IVF1.tmp <- IVF1.data[IVF1.keep,]
IVFI.atleastone <- apply(IVF1.tmp, 2, function(x) sum(x>0))
summary(IVF1.atleastone)
rownames(x = IVF1.data) <- gsub(pattern = '_', replacement = '-', x= rownames(x = IVF1.data))
IVF1 <- CreateSeuratObject(counts = IVF1.data, project = "1IVF", min.cell = 3, min.features = 200)
IVF1$patient <- "IVF#1"
IVF1 <- subset(x = IVF1, subset = nFeature_RNA > 500)
IVF1 <- NormalizeData(object = IVF1, verbose = FALSE)
IVF1 <- FindVariableFeatures(object = IVF1, selection.method = "vst", nfeatures = 2000)

### Setup the Seurat Object-IVF2
IVF2.data <- Read10X(data.dir = "D:/scRNA-seq data/Processed data/IVF2/outs/filtered_feature_bc_matrix/")
IVF2.atleastone <- apply(IVF2.data, 2, function(x) sum(x>0))
hist(IVF2.atleastone, breaks = 100,
     main = "Distribution of detected genes - IVF2",
     xlab = "Genes with at least one tag")
IVF2.tmp <- apply(IVF2.data, 1, function(x) sum (x>0))
table(IVF2.tmp>=3)
IVF2.keep <- IVF2.tmp>=3
IVF2.tmp <- IVF2.data[IVF2.keep,]
IVF2.atleastone <- apply(IVF2.tmp, 2, function(x) sum(x>0))
summary(IVF2.atleastone)
rownames(x = IVF2.data) <- gsub(pattern = '_', replacement = '-', x= rownames(x = IVF2.data))
IVF2 <- CreateSeuratObject(counts = IVF2.data, project = "2IVF", min.cell = 3, min.features = 200)
IVF2$patient <- "IVF#2"
IVF2 <- subset(x = IVF2, subset = nFeature_RNA > 500)
IVF2 <- NormalizeData(object = IVF2, verbose = FALSE)
IVF2 <- FindVariableFeatures(object = IVF2, selection.method = "vst", nfeatures = 2000)

### Setup the Seurat Object-IVF3
IVF3.data <- Read10X(data.dir = "D:/scRNA-seq data/Processed data/IVF3/outs/filtered_feature_bc_matrix/")
IVF3.atleastone <- apply(IVF3.data, 2, function(x) sum(x>0))
hist(IVF3.atleastone, breaks = 100,
     main = "Distribution of detected genes - IVF3",
     xlab = "Genes with at least one tag")
IVF3.tmp <- apply(IVF3.data, 1, function(x) sum (x>0))
table(IVF3.tmp>=3)
IVF3.keep <- IVF3.tmp>=3
IVF3.tmp <- IVF3.data[IVF3.keep,]
IVF3.atleastone <- apply(IVF3.tmp, 2, function(x) sum(x>0))
summary(IVF3.atleastone)
rownames(x = IVF3.data) <- gsub(pattern = '_', replacement = '-', x= rownames(x = IVF3.data))
IVF3 <- CreateSeuratObject(counts = IVF3.data, project = "3IVF", min.cell = 3, min.features = 200)
IVF3$patient <- "IVF#3"
IVF3 <- subset(x = IVF3, subset = nFeature_RNA > 500)
IVF3 <- NormalizeData(object = IVF3, verbose = FALSE)
IVF3 <- FindVariableFeatures(object = IVF3, selection.method = "vst", nfeatures = 2000)

### Setup the Seurat Object-IVF4
IVF4.data <- Read10X(data.dir = "D:/scRNA-seq data/Processed data/IVF4/outs/filtered_feature_bc_matrix/")
IVF4.atleastone <- apply(IVF4.data, 2, function(x) sum(x>0))
hist(IVF4.atleastone, breaks = 100,
     main = "Distribution of detected genes - IVF4",
     xlab = "Genes with at least one tag")
IVF4.tmp <- apply(IVF4.data, 1, function(x) sum (x>0))
table(IVF4.tmp>=3)
IVF4.keep <- IVF4.tmp>=3
IVF4.tmp <- IVF4.data[IVF4.keep,]
IVF4.atleastone <- apply(IVF4.tmp, 2, function(x) sum(x>0))
summary(IVF4.atleastone)
rownames(x = IVF4.data) <- gsub(pattern = '_', replacement = '-', x= rownames(x = IVF4.data))
IVF4 <- CreateSeuratObject(counts = IVF4.data, project = "4IVF", min.cell = 3, min.features = 200)
IVF4$patient <- "IVF#4"
IVF4 <- subset(x = IVF4, subset = nFeature_RNA > 500)
IVF4 <- NormalizeData(object = IVF4, verbose = FALSE)
IVF4 <- FindVariableFeatures(object = IVF4, selection.method = "vst", nfeatures = 2000)

## Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(IVF1, IVF2, IVF3, IVF4), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"

## visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
DefaultAssay(object = combined) <- "RNA"
saveRDS(combined, file="integrated_IVFs.rds")


readRDS("integrated_IVFs.rds")

## Plotting
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "patient")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1)
plot_grid(p2)
DimPlot(combined, reduction = "umap", pt.size = 1, label = TRUE, label.size = 9) + NoLegend()

## UMAPs showing 4 individual IVF patients
DimPlot(object = combined, reduction = "umap", split.by = "patient")


## Cell-specific markers for re-clustering
###Follicular cells vs. leukocytes vs. RBCs
FeaturePlot(combined, features = c("CYP11A1", "PTPRC", "HBA1"))

###Follicular cells
VlnPlot(combined, features = c("CYP11A1", "HSD3B2"), slot = "counts", log = TRUE)

###Granulosa cells
VlnPlot(combined, features = c("CYP19A1", "INHA"), slot = "counts", log = TRUE)

###Theca cells
VlnPlot(combined, features = c("COL1A1", "COL1A2", "COL3A1"), slot = "counts", log = TRUE)

###Leukocytes, RBCs, & endothelial cells
VlnPlot(combined, features = c("PTPRC", "HBA1", "TIE1", "VWF"), slot = "counts", log = TRUE)


###Subclustering leukocytes
####Cytotoxic & Helper T cells
VlnPlot(combined, features = c("CD3E", "CD8A", "CD4"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD3E", "CD8A", "CD4"))

####NK cells
VlnPlot(combined, features = c("NKG7", "NCR1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("NKG7", "NCR1"))

####B cells
VlnPlot(combined, features = c("CD19", "MS4A1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD19", "MS4A1"))

####RBCs & Platelets
VlnPlot(combined, features = c("HBA1", "PF4", "PPBP"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("PF4", "PPBP"))

####Neutrophils
VlnPlot(combined, features = c("FCER1A", "CRLF2"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("FCER1A", "CRLF2"))

####baso/eosinophils
VlnPlot(combined, features = c("CCR3", "PTGDR2"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CCR3", "PTGDR2"))

####Macrophages
VlnPlot(combined, features = c("CD14", "CD68"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD14", "CD68"))
#####M1
VlnPlot(combined, features = c("ITGAX", "HLA-DRA", "HLA-DRB1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("ITGAX", "HLA-DRA", "HLA-DRB1"))
#####M2
VlnPlot(combined, features = c("CD163", "MSR1", "MRC1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD163", "MSR1", "MRC1"))



##Re-cluster by localization of cell-specific markers
new.cluster.ids <- c("0- GC_3", 
                     "1- T-helper", 
                     "2- GC_4", 
                     "3- M1-Macrophage", 
                     "4- Cytotoxic T", 
                     "5- M2-Macrophage", 
                     "6- NK", 
                     "7- NKT", 
                     "8- Neutrophil", 
                     "9- GC_1", 
                     "10- Baso/eosinophil", 
                     "11- GC_5", 
                     "12- GC_2", 
                     "13- M2-Macrophage", 
                     "14- M1-Macrophage", 
                     "15- B", 
                     "16- RBC/Platelet", 
                     "17- Dendritic", 
                     "18- TC")
names(new.cluster.ids) <- levels(combined)
combined1 <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined1, reduction = "umap", label = TRUE, label.size = 9, pt.size = 1) + NoLegend()

##Re-order the clusters by manually

reordered <- c("9- GC_1",
               "12- GC_2",
               "0- GC_3",
               "2- GC_4",
               "11- GC_5",
               "18- TC",
               "6- NK", 
               "7- NKT",
               "3- M1-Macrophage",
               "14- M1-Macrophage",
               "5- M2-Macrophage",
               "13- M2-Macrophage",
               "1- T-helper", 
               "4- Cytotoxic T",
               "8- Neutrophil",
               "10- Baso/eosinophil",  
               "15- B",
               "17- Dendritic",
               "16- RBC/Platelet" 
)


combined1@active.ident <- factor(combined1@active.ident, levels = reordered)
##UMAP showing reordered clusters
DimPlot(combined1, reduction = "umap", label = TRUE, label.size = 9, pt.size = 1)

##Re-cluster by localization of cell-specific markers
new.cluster.ids1<- c("GC_1",
                     "GC_2",
                     "GC_3",
                     "GC_4",
                     "GC_5",
                     "TC",
                     "NK", 
                     "NKT",
                     "M1-Macrophage",
                     "M1-Macrophage",
                     "M2-Macrophage",
                     "M2-Macrophage",
                     "T-helper", 
                     "Cytotoxic T",
                     "Neutrophil",
                     "Baso/eosinophil",  
                     "B",
                     "Dendritic",
                     "RBC/Platelet" 
)
names(new.cluster.ids1) <- levels(combined1)
combined1 <- RenameIdents(combined1, new.cluster.ids1)
DimPlot(combined1, reduction = "umap", pt.size = 1) + NoLegend()

saveRDS(combined1, file="combined1.rds")
readRDS("combined1.rds")



## UMAPs showing 4 individual IVF patients after re-clustering
DimPlot(object = combined1, reduction = "umap", split.by = "patient")

##Frequencies of clusters or cell types
###Barplot representing percentages of cell subpopulations
dittoBarPlot(combined1, "ident", group.by = "patient")
###Barplot representing counts of cell subpopulations
dittoBarPlot(combined1, "ident", group.by = "patient", scale = "count")


##sub-clustering

###all FC
FC <- subset(x = combined1, idents = c("GC_1", "GC_2", "GC_3","GC_4", "GC_5", "TC"))
FC
FC <- ScaleData(object = FC, verbose = FALSE)
FC <- RunPCA(object = FC, npcs = 30, verbose = FALSE)
FC <- RunUMAP(object = FC, reduction = "pca", dims = 1:25)
FC <- FindNeighbors(object = FC, reduction = "pca", dims = 1:25)
FC <- FindClusters(FC, resolution = 0.7)
DimPlot(object = FC, reduction = "umap", label.size = 5)
DimPlot(object = FC, reduction = "umap", label.size = 5, split.by = "patient")



###all GC
all_GC <- subset(x = combined1, idents = c("GC_1", "GC_2", "GC_3", "GC_4", "GC_5"))
all_GC
all_GC <- ScaleData(object = all_GC, verbose = FALSE)
all_GC <- RunPCA(object = all_GC, npcs = 30, verbose = FALSE)
all_GC <- RunUMAP(object = all_GC, reduction = "pca", dims = 1:25)
all_GC <- FindNeighbors(object = all_GC, reduction = "pca", dims = 1:25)
all_GC <- FindClusters(all_GC, resolution = 0.7)
DefaultAssay(object = all_GC) <- "RNA"
DimPlot(object = all_GC, reduction = "umap", label.size = 5)
DimPlot(object = all_GC, reduction = "umap", label.size = 5, split.by = "patient")

###GC-1
GC_1 <- subset(x = combined1, idents = c("GC_1"))
GC_1
GC_1 <- ScaleData(object = GC_1, verbose = FALSE)
GC_1 <- RunPCA(object = GC_1, npcs = 30, verbose = FALSE)
GC_1 <- RunUMAP(object = GC_1, reduction = "pca", dims = 1:25)
GC_1 <- FindNeighbors(object = GC_1, reduction = "pca", dims = 1:25)
GC_1 <- FindClusters(GC_1, resolution = 0.7)
DimPlot(object = GC_1, reduction = "umap", label.size = 5)

###GC
GC_2 <- subset(x = combined1, idents = c("GC_2"))
GC_2
GC_2 <- ScaleData(object = GC_2, verbose = FALSE)
GC_2 <- RunPCA(object = GC_2, npcs = 30, verbose = FALSE)
GC_2 <- RunUMAP(object = GC_2, reduction = "pca", dims = 1:25)
GC_2 <- FindNeighbors(object = GC_2, reduction = "pca", dims = 1:25)
GC_2 <- FindClusters(GC_2, resolution = 0.7)
DimPlot(object = GC_2, reduction = "umap", label.size = 5)

###GC
GC_3 <- subset(x = combined1, idents = c("GC_3"))
GC_3
GC_3 <- ScaleData(object = GC_3, verbose = FALSE)
GC_3 <- RunPCA(object = GC_3, npcs = 30, verbose = FALSE)
GC_3 <- RunUMAP(object = GC_3, reduction = "pca", dims = 1:25)
GC_3 <- FindNeighbors(object = GC_3, reduction = "pca", dims = 1:25)
GC_3 <- FindClusters(GC_3, resolution = 0.7)
DimPlot(object = GC_3, reduction = "umap", label.size = 5)

###GC
GC_4 <- subset(x = combined1, idents = c("GC_4"))
GC_4
GC_4 <- ScaleData(object = GC_4, verbose = FALSE)
GC_4 <- RunPCA(object = GC_4, npcs = 30, verbose = FALSE)
GC_4 <- RunUMAP(object = GC_4, reduction = "pca", dims = 1:25)
GC_4 <- FindNeighbors(object = GC_4, reduction = "pca", dims = 1:25)
GC_4 <- FindClusters(GC_4, resolution = 0.7)
DimPlot(object = GC_4, reduction = "umap", label.size = 5)

###GC_5
GC_5 <- subset(x = combined1, idents = c("GC_5"))
GC_5
GC_5 <- ScaleData(object = GC_2, verbose = FALSE)
GC_5 <- RunPCA(object = GC_5, npcs = 30, verbose = FALSE)
GC_5 <- RunUMAP(object = GC_5, reduction = "pca", dims = 1:25)
GC_5 <- FindNeighbors(object = GC_5, reduction = "pca", dims = 1:25)
GC_5 <- FindClusters(GC_5, resolution = 0.7)
DimPlot(object = GC_5, reduction = "umap", label.size = 5)


###TC
TC <- subset(x = combined1, idents = c("TC"))
TC
TC <- ScaleData(object = TC, verbose = FALSE)
TC <- RunPCA(object = TC, npcs = 30, verbose = FALSE)
TC <- RunUMAP(object = TC, reduction = "pca", dims = 1:25)
TC <- FindNeighbors(object = TC, reduction = "pca", dims = 1:25)
TC <- FindClusters(TC, resolution = 0.7)
DimPlot(object = TC, reduction = "umap", label.size = 5)

###leukocytes
leukocyte <- subset(x = combined1, idents = c("T-helper",
                                              "M1-Macrophage",
                                              "M2-Macrophage",
                                              "Cytotoxic T",
                                              "NK",
                                              "NKT",
                                              "Neutrophil",
                                              "Baso/eosinophil",
                                              "B",
                                              "Dendritic"))
leukocyte
leukocyte <- ScaleData(object = leukocyte, verbose = FALSE)
leukocyte <- RunPCA(object = leukocyte, npcs = 30, verbose = FALSE)
leukocyte <- RunUMAP(object = leukocyte, reduction = "pca", dims = 1:25)
leukocyte <- FindNeighbors(object = leukocyte, reduction = "pca", dims = 1:25)
leukocyte <- FindClusters(leukocyte, resolution = 0.7)
DimPlot(object = leukocyte, reduction = "umap", label.size = 5)
DimPlot(object = leukocyte, reduction = "umap", label.size = 5, split.by = "patient")

##Macrophages
Mac <- subset(x = combined, idents = c("3- M1-Macrophage", "5- M2-Macrophage", "13- M2-Macrophage", "14- M1-Macrophage"))
Mac
Mac <- ScaleData(object = Mac, verbose = FALSE)
Mac <- RunPCA(object = Mac, npcs = 30, verbose = FALSE)
Mac <- RunUMAP(object = Mac, reduction = "pca", dims = 1:25)
Mac <- FindNeighbors(object = Mac, reduction = "pca", dims = 1:25)
Mac <- FindClusters(Mac, resolution = 0.7)
DimPlot(object = Mac, reduction = "umap", label.size = 5)

### M1-Macrophages
M1 <- subset(x = combined1, idents = c("M1-Macrophage"))
M1
M1 <- ScaleData(object = M1, verbose = FALSE)
M1 <- RunPCA(object = M1, npcs = 30, verbose = FALSE)
M1 <- RunUMAP(object = M1, reduction = "pca", dims = 1:25)
M1 <- FindNeighbors(object = M1, reduction = "pca", dims = 1:25)
M1 <- FindClusters(M1, resolution = 0.7)
DimPlot(object = M1, reduction = "umap", label.size = 5)

### M2-Macrophages
M2 <- subset(x = combined1, idents = c("M2-Macrophage"))
M2
M2 <- ScaleData(object = M2, verbose = FALSE)
M2 <- RunPCA(object = M2, npcs = 30, verbose = FALSE)
M2 <- RunUMAP(object = M2, reduction = "pca", dims = 1:25)
M2 <- FindNeighbors(object = M2, reduction = "pca", dims = 1:25)
M2 <- FindClusters(M2, resolution = 0.7)
DimPlot(object = M2, reduction = "umap", label.size = 5)

### T-helper cells
Thelper <- subset(x = combined1, idents = c("T-helper"))
Thelper
Thelper <- ScaleData(object = Thelper, verbose = FALSE)
Thelper <- RunPCA(object = Thelper, npcs = 30, verbose = FALSE)
Thelper <- RunUMAP(object = Thelper, reduction = "pca", dims = 1:25)
Thelper <- FindNeighbors(object = Thelper, reduction = "pca", dims = 1:25)
Thelper <- FindClusters(Thelper, resolution = 0.7)
DimPlot(object = Thelper, reduction = "umap", label.size = 5)

### Cytotoxic T cells
CytoT <- subset(x = combined1, idents = c("Cytotoxic T"))
CytoT
CytoT <- ScaleData(object = CytoT, verbose = FALSE)
CytoT <- RunPCA(object = CytoT, npcs = 30, verbose = FALSE)
CytoT <- RunUMAP(object = CytoT, reduction = "pca", dims = 1:25)
CytoT <- FindNeighbors(object = CytoT, reduction = "pca", dims = 1:25)
CytoT <- FindClusters(CytoT, resolution = 0.7)
DimPlot(object = CytoT, reduction = "umap", label.size = 5)


### NK cells
NK <- subset(x = combined1, idents = c("NK"))
NK
NK <- ScaleData(object = NK, verbose = FALSE)
NK <- RunPCA(object = NK, npcs = 30, verbose = FALSE)
NK <- RunUMAP(object = NK, reduction = "pca", dims = 1:25)
NK <- FindNeighbors(object = NK, reduction = "pca", dims = 1:25)
NK <- FindClusters(NK, resolution = 0.7)
DimPlot(object = NK, reduction = "umap", label.size = 5)

### NKT cells
NKT <- subset(x = combined1, idents = c("NKT"))
NKT
NKT <- ScaleData(object = NKT, verbose = FALSE)
NKT <- RunPCA(object = NKT, npcs = 30, verbose = FALSE)
NKT <- RunUMAP(object = NKT, reduction = "pca", dims = 1:25)
NKT <- FindNeighbors(object = NKT, reduction = "pca", dims = 1:25)
NKT <- FindClusters(NKT, resolution = 0.7)
DimPlot(object = NKT, reduction = "umap", label.size = 5)

### Neutrophils
Neutro <- subset(x = combined1, idents = c("Neutrophil"))
Neutro
Neutro <- ScaleData(object = Neutro, verbose = FALSE)
Neutro <- RunPCA(object = Neutro, npcs = 30, verbose = FALSE)
Neutro <- RunUMAP(object = Neutro, reduction = "pca", dims = 1:25)
Neutro <- FindNeighbors(object = Neutro, reduction = "pca", dims = 1:25)
Neutro <- FindClusters(Neutro, resolution = 0.7)
DimPlot(object = Neutro, reduction = "umap", label.size = 5)

### Baso/eosinophils
BasoEosino <- subset(x = combined1, idents = c("Baso/eosinophil"))
BasoEosino
BasoEosino <- ScaleData(object = BasoEosino, verbose = FALSE)
BasoEosino <- RunPCA(object = BasoEosino, npcs = 30, verbose = FALSE)
BasoEosino <- RunUMAP(object = BasoEosino, reduction = "pca", dims = 1:25)
BasoEosino <- FindNeighbors(object = BasoEosino, reduction = "pca", dims = 1:25)
BasoEosino <- FindClusters(BasoEosino, resolution = 0.7)
DimPlot(object = BasoEosino, reduction = "umap", label.size = 5)

### B-cells
B <- subset(x = combined1, idents = c("B"))
B
B <- ScaleData(object = B, verbose = FALSE)
B <- RunPCA(object = B, npcs = 30, verbose = FALSE)
B <- RunUMAP(object = B, reduction = "pca", dims = 1:25)
B <- FindNeighbors(object = B, reduction = "pca", dims = 1:25)
B <- FindClusters(B, resolution = 0.7)
DimPlot(object = B, reduction = "umap", label.size = 5)

### Dendritic cells
Dendritic <- subset(x = combined1, idents = c("Dendritic"))
Dendritic
Dendritic <- ScaleData(object = Dendritic, verbose = FALSE)
Dendritic <- RunPCA(object = Dendritic, npcs = 30, verbose = FALSE)
Dendritic <- RunUMAP(object = Dendritic, reduction = "pca", dims = 1:25)
Dendritic <- FindNeighbors(object = Dendritic, reduction = "pca", dims = 1:25)
Dendritic <- FindClusters(Dendritic, resolution = 0.7)
DimPlot(object = Dendritic, reduction = "umap", label.size = 5)

### RBCs and Platelets
RBCPlatelet <- subset(x = combined1, idents = c("RBC/Platelet"))
RBCPlatelet
RBCPlatelet <- ScaleData(object = RBCPlatelet, verbose = FALSE)
RBCPlatelet <- RunPCA(object = RBCPlatelet, npcs = 30, verbose = FALSE)
RBCPlatelet <- RunUMAP(object = RBCPlatelet, reduction = "pca", dims = 1:25)
RBCPlatelet <- FindNeighbors(object = RBCPlatelet, reduction = "pca", dims = 1:25)
RBCPlatelet <- FindClusters(RBCPlatelet, resolution = 0.7)
DimPlot(object = RBCPlatelet, reduction = "umap", label.size = 5)


### Proportions of FC and leukocytes
FC_Leuko <- c("FC",
              "FC",
              "FC",
              "FC",
              "FC",
              "FC",
              "Leukocyte", 
              "Leukocyte",
              "Leukocyte",
              "Leukocyte",
              "Leukocyte",
              "Leukocyte",
              "Leukocyte", 
              "Leukocyte",
              "Leukocyte",
              "Leukocyte",  
              "Leukocyte",
              "Leukocyte",
              "RBC/Platelet")
names(FC_Leuko) <- levels(combined1)
combined3 <- RenameIdents(combined1, FC_Leuko)
DimPlot(combined3, reduction = "umap", label = TRUE, label.size = 9, pt.size = 0.5) + NoLegend()


##Frequencies of clusters or cell types
###Barplot representing percentages of cell subpopulations
dittoBarPlot(combined3, "ident", group.by = "patient", data.out = TRUE)



###Barplot representing counts of cell subpopulations
dittoBarPlot(combined3, "ident", group.by = "patient", scale = "count")




##add a 'cell_type' variable to GC clusters
all_GC$cell_type <- as.character(Idents(all_GC))
combined1$cell_type <- as.character(Idents(combined1))
leukocyte$cell_type <- as.character(Idents(leukocyte))


saveRDS("all_GC.rds")

readRDS("combined1.RDs")





#create looms
library(SeuratDisk)
library(SeuratData)



# convert the all_GC obeect to loom
GC_all5.loom <- as.loom(all_GC, filename = "GC_all5.loom", verbose = FALSE)
GC_all5.loom

# close looms
GC_all5.loom$close_all()




##boxplot by RNA velocity

###GC125
GC125 <- subset(x = combined1, idents = c("GC_1", "GC_2", "GC_5"))
GC125
GC125 <- ScaleData(object = GC125, verbose = FALSE)
GC125 <- RunPCA(object = GC125, npcs = 30, verbose = FALSE)
GC125 <- RunUMAP(object = GC125, reduction = "pca", dims = 1:25)
GC125 <- FindNeighbors(object = GC125, reduction = "pca", dims = 1:25)
GC125 <- FindClusters(GC125, resolution = 0.7)
DefaultAssay(object = GC125) <- "RNA"



###GC1345
GC1345 <- subset(x = combined1, idents = c("GC_1", "GC_3", "GC_4", "GC_5"))
GC1345
GC1345 <- ScaleData(object = GC1345, verbose = FALSE)
GC1345 <- RunPCA(object = GC1345, npcs = 30, verbose = FALSE)
GC1345 <- RunUMAP(object = GC1345, reduction = "pca", dims = 1:25)
GC1345 <- FindNeighbors(object = GC1345, reduction = "pca", dims = 1:25)
GC1345 <- FindClusters(GC1345, resolution = 0.7)
DefaultAssay(object = GC1345) <- "RNA"

dittoHeatmap(all_GC, 
             c("AMH", "FSHR", "PGR", "INHBA", "INHBB", "BMP3", "RUNX2", "ADAMTS1", "HSD11B1", "EDN2", "SLCO2A1", "CD24"),
             group.by = "idents")

DoHeatmap(all_GC, 
          features = c("AMH", "FSHR", "PGR", "INHBA", "INHBB", "BMP3", "DHCR24", "IL6ST", "HSD11B1", "SORBS2", "TSHZ2", "DUXAR8", "HSPG2", "FLNB", "RUNX2", "ADAMTS1", "HSD11B1", "EDN2", "SLCO2A1", "MT2A", "STAR", "CMAS"))




##save objects in rds files
saveRDS(combined1, "combined1")
saveRDS(all_GC, "all_GC.rds")
saveRDS(leukocyte, "leukocyte.rds")

readRDS("all_GC.rds")


##create a dot plot showing expression levels and the percentage of cells in a cluster


###Cytokines and receptors
dittoDotPlot(combined1, c("IFNGR2", "IFNGR1", "IFNG", "CD9", "IL16", "IL1RAP", "IL1R1", "IL1B"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("IL1B", "IL1R1", "IL1RAP"))
FeaturePlot(combined1, features = c("IL16", "CD9"))
FeaturePlot(combined1, features = c("IFNG", "IFNGR1", "IFNGR2"))

###Death ligands and receptors
dittoDotPlot(combined1, c("TNFRSF10C", "TNFRSF10B", "TNFRSF10A", "TNFSF10", "LTBR", "LTB", "LTA", "TNFRSF1B", "TNFRSF1A", "TNF"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("TNF", "TNFRSF1A", "TNFRSF1B"))
FeaturePlot(combined1, features = c("LTA", "LTB", "LTBR"))
FeaturePlot(combined1, features = c("TNFSF10", "TNFRSF10A", "TNFRSF10B"))


###PLAU & PLAUR
dittoDotPlot(combined1, c("F2R", "SERPINE1", "PLAUR", "PLAU"),
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined, features = c("PLAU", "PLAUR", "SERPINE1", "F2R"))


###NGR1 & its receptors
dittoDotPlot(combined1, c("ERBB4", "ERBB3", "ERBB2", "NRG1"),
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined, features = c("NRG1", "ERBB2", "ERBB3", "ERBB4"))


dittoDotPlot(combined1, c("PLAU", "PLAT"),
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined, features = c("PLAT", "PLAU"))

###MMPs
dittoDotPlot(combined1, c("MMP25", "MMP23B", "MMP19", "MMP17", "MMP14", "MMP10", "MMP9", "MMP2"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("MMP2", "MMP9"))
FeaturePlot(combined1, features = c("MMP10", "MMP14"))
FeaturePlot(combined1, features = c("MMP17", "MMP19"))
FeaturePlot(combined1, features = c("MMP23B", "MMP25"))

###ADAMs
dittoDotPlot(combined1, c("ADAM28", "ADAM19", "ADAM15", "ADAM10", "ADAM9", "ADAM8"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("ADAM8", "ADAM9"))
FeaturePlot(combined1, features = c("ADAM10", "ADAM15"))
FeaturePlot(combined1, features = c("ADAM19", "ADAM28"))

###ADAMTSs
dittoDotPlot(combined1, c("ADAMTS17", "ADAMTS10", "ADAMTS6", "ADAMTS5", "ADAMTS4", "ADAMTS1"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("ADAMTS1", "ADAMTS4"))
FeaturePlot(combined1, features = c("ADAMTS5", "ADAMTS6"))
FeaturePlot(combined1, features = c("ADAMTS10", "ADAMTS17"))

###TIMPs
dittoDotPlot(combined1, c("TIMP4", "TIMP3", "TIMP2", "TIMP1"), 
             group.by = "ident") + scale_x_discrete(position = "bottom") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined1, features = c("TIMP1", "TIMP2", "TIMP3", "TIMP4"))

###VEGFs
dittoDotPlot(combined1, c("THBS1", "IGF1", "FGF2", "PGF", "VEGFC", "VEGFB", "VEGFA"),
             group.by = "ident") + scale_x_discrete(position = "bottom") + coord_flip()
FeaturePlot(combined, features = c("VEGFA", "VEGFB", "VEGFC", "PGF", "FGF2", "IGF1", "THBS1"))
