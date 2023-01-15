library(dplyr)
library(Seurat)
library(Matrix)
library(DoubletFinder)
library(ggplot2)


##load datasets

syn5844h.data <- read.delim("P5844_Hippo_syn_transcript.tsv", row.names=1)
syn5844h <- CreateSeuratObject(counts = syn5844h.data, project = "syn5844h", min.cells = 3, min.features = 100)
##QC
syn5844h[["percent.mt"]] <- PercentageFeatureSet(syn5844h, pattern = "^MT-")
syn5844h[["UMI_count"]] <-syn5844h[["nCount_RNA"]] 
syn5844h[["Gene_count"]] <-syn5844h[["nFeature_RNA"]] 
VlnPlot(syn5844h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
plot1 <- FeatureScatter(syn5844h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(syn5844h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

syn5844h <- subset(syn5844h, subset = nCount_RNA<quantile(a, 0.99))
syn5844h[["patient_ID"]]="P5844"

syn5844h <- SCTransform(syn5844h, method = "glmGamPoi", verbose = TRUE)

syn5844h <- RunPCA(syn5844h)
ElbowPlot(syn5844h,ndims=50)
syn5844h <- RunUMAP(syn5844h, dims = 1:10)
DimPlot(syn5844h,reduction = "umap")
p1 <- DimPlot(syn5844h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn5844h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot(syn5844h, reduction="umap",features = c("FNDC1", "SEMA5A","GAD1","AQP4","SPOCK1"), pt.size=0.5)

#pK_selection
sweep.res.list_syn5844h <- paramSweep_v3(syn5844h, PCs = 1:10, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn5844h <- summarizeSweep(sweep.res.list_syn5844h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn5844h <- find.pK(sweep.stats_syn5844h)
mpK<-as.numeric(as.vector(bcmvn_syn5844h$pK[which.max(bcmvn_syn5844h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn5844h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn5844h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn5844h <- doubletFinder_v3(syn5844h, PCs = 1:10, 
                             pN = 0.25, pK = 0.09, 
                             nExp = 132, 
                             reuse.pANN = FALSE ,
                             sct = T)
table(syn5844h@meta.data[["DF.classifications_0.25_0.29_132"]])
syn5844h[["doublet"]] <-syn5844h[["DF.classifications_0.25_0.09_132"]] 



syn5818h.data <- read.delim("P5818_Hippo_syn_transcript.tsv", row.names=1)
syn5818h <- CreateSeuratObject(counts = syn5818h.data, project = "syn5818h", min.cells = 3, min.features = 100)
##QC
syn5818h[["percent.mt"]] <- PercentageFeatureSet(syn5818h, pattern = "^MT-")
syn5818h[["UMI_count"]] <-syn5818h[["nCount_RNA"]] 
syn5818h[["Gene_count"]] <-syn5818h[["nFeature_RNA"]] 
VlnPlot(syn5818h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
plot1 <- FeatureScatter(syn5818h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(syn5818h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
a=syn5818h[["UMI_count"]]$UMI_count
quantile(a, 0.99) 
syn5818h <- subset(syn5818h, subset = nCount_RNA<quantile(a, 0.99))
syn5818h[["patient_ID"]]="P5818"


syn5818h <- SCTransform(syn5818h, method = "glmGamPoi", verbose = TRUE)

syn5818h <- RunPCA(syn5818h)
ElbowPlot(syn5818h,ndims=50)

syn5818h <- RunUMAP(syn5818h, dims = 1:10,spread=3,min.dist=0.1)
DimPlot(syn5818h,reduction = "umap")
p1 <- DimPlot(syn5818h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn5818h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot
#pK_selection
sweep.res.list_syn5818h <- paramSweep_v3(syn5818h, PCs = 1:10, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn5818h <- summarizeSweep(sweep.res.list_syn5818h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn5818h <- find.pK(sweep.stats_syn5818h)
mpK<-as.numeric(as.vector(bcmvn_syn5818h$pK[which.max(bcmvn_syn5818h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn5818h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn5818h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn5818h <- doubletFinder_v3(syn5818h, PCs = 1:10, 
                             pN = 0.25, pK = mpK, 
                             nExp = nExp_poi.adj, 
                             reuse.pANN = FALSE ,
                             sct = T)
table(syn5818h@meta.data[["DF.classifications_0.25_0.09_241"]])
syn5818h[["doublet"]] <-syn5818h[["DF.classifications_0.25_0.05_241"]]



##Merge datasets
synh.combined <- merge(syn5844h, y = syn5818h, add.cell.ids = c("syn5844h", "syn5818h"), project = "synh")


# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(synh.combined, split.by = "patient_ID")


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x<- SCTransform(x, method = "glmGamPoi", verbose = TRUE)
})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

###Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined.sct) <- "integrated"
immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 50, verbose = FALSE)
ElbowPlot(immune.combined.sct,ndims = 50)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:14, 
                               spread=2,min.dist=0.15,metric="cosine")


# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "patient_ID",shuffle=T)
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "doublet")
p3 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 +p3



