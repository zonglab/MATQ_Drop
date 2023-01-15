library(dplyr)
library(Seurat)
library(Matrix)
library(DoubletFinder)
library(ggplot2)


syn417h.data <- read.delim("Combined_417syn_transcript_counts_rRNA_MT_filtered_gene_symbol_noDup.tsv", row.names=1)
syn417h <- CreateSeuratObject(counts = syn417h.data, project = "syn417h", min.cells = 3, min.features = 100)
##QC
syn417h[["percent.mt"]] <- PercentageFeatureSet(syn417h, pattern = "^MT-")
syn417h[["UMI_count"]] <-syn417h[["nCount_RNA"]] 
syn417h[["Gene_count"]] <-syn417h[["nFeature_RNA"]] 
a=syn417h[["UMI_count"]]$UMI_count
quantile(a, 0.99) 
syn417h <- subset(syn417h, subset = nCount_RNA<quantile(a, 0.99))
syn417h[["patient_ID"]]="P417"

syn417h <- SCTransform(syn417h, method = "glmGamPoi", verbose = TRUE,ncells = 20000)

syn417h <- RunPCA(syn417h)
ElbowPlot(syn417h,ndims=50)

syn417h <- RunUMAP(syn417h, dims = 1:15,spread=3,min.dist=0.1)
DimPlot(syn417h,reduction = "umap")
p1 <- DimPlot(syn417h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn417h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot(syn417h,features=c("UMI_count"),max.cutoff = 2000)
FeaturePlot(syn417h,features=c("Csmd1"))
#pK_selection
sweep.res.list_syn417h <- paramSweep_v3(syn417h, PCs = 1:15, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn417h <- summarizeSweep(sweep.res.list_syn417h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn417h <- find.pK(sweep.stats_syn417h)
mpK<-as.numeric(as.vector(bcmvn_syn417h$pK[which.max(bcmvn_syn417h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn417h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn417h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn417h <- doubletFinder_v3(syn417h, PCs = 1:15, 
                            pN = 0.25, pK = mpK, 
                            nExp = nExp_poi.adj, 
                            reuse.pANN = FALSE ,
                            sct = T)
table(syn417h@meta.data[["DF.classifications_0.25_0.01_352"]])
syn417h[["doublet"]] <-syn417h[["DF.classifications_0.25_0.01_352"]]


########
#########
syn418h.data <- read.delim("Combined_418syn_transcript_counts_rRNA_MT_filtered_gene_symbol_noDup.tsv", row.names=1)
syn418h <- CreateSeuratObject(counts = syn418h.data, project = "syn418h", min.cells = 3, min.features = 100)
##QC
syn418h[["percent.mt"]] <- PercentageFeatureSet(syn418h, pattern = "^MT-")
syn418h[["UMI_count"]] <-syn418h[["nCount_RNA"]] 
syn418h[["Gene_count"]] <-syn418h[["nFeature_RNA"]] 
a=syn418h[["UMI_count"]]$UMI_count
quantile(a, 0.99) 
syn418h <- subset(syn418h, subset = nCount_RNA<quantile(a, 0.99))
syn418h[["patient_ID"]]="P418"

syn418h <- SCTransform(syn418h, method = "glmGamPoi", verbose = TRUE,ncells = 20000)

syn418h <- RunPCA(syn418h)
ElbowPlot(syn418h,ndims=50)

syn418h <- RunUMAP(syn418h, dims = 1:15,spread=3,min.dist=0.1)
DimPlot(syn418h,reduction = "umap")
p1 <- DimPlot(syn418h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn418h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot(syn418h,features=c("UMI_count"),max.cutoff = 150)
FeaturePlot(syn418h,features=c("Zbtb20"))
#pK_selection
sweep.res.list_syn418h <- paramSweep_v3(syn418h, PCs = 1:15, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn418h <- summarizeSweep(sweep.res.list_syn418h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn418h <- find.pK(sweep.stats_syn418h)
mpK<-as.numeric(as.vector(bcmvn_syn418h$pK[which.max(bcmvn_syn418h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn418h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn418h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn418h <- doubletFinder_v3(syn418h, PCs = 1:15, 
                            pN = 0.25, pK = mpK, 
                            nExp = nExp_poi.adj, 
                            reuse.pANN = FALSE ,
                            sct = T)
table(syn418h@meta.data[["DF.classifications_0.25_0.01_324"]])
syn418h[["doublet"]] <-syn418h[["DF.classifications_0.25_0.01_324"]]


######
######
syn480h.data <- read.delim("Combined_480syn_transcript_counts_rRNA_MT_filtered_gene_symbol_noDup.tsv", row.names=1)
syn480h <- CreateSeuratObject(counts = syn480h.data, project = "syn480h", min.cells = 3, min.features = 100)
##QC
syn480h[["percent.mt"]] <- PercentageFeatureSet(syn480h, pattern = "^MT-")
syn480h[["UMI_count"]] <-syn480h[["nCount_RNA"]] 
syn480h[["Gene_count"]] <-syn480h[["nFeature_RNA"]] 
a=syn480h[["UMI_count"]]$UMI_count
quantile(a, 0.99) 
syn480h <- subset(syn480h, subset = nCount_RNA<quantile(a, 0.99))
syn480h[["patient_ID"]]="P480"

syn480h <- SCTransform(syn480h, method = "glmGamPoi", verbose = TRUE,ncells = 20000)

syn480h <- RunPCA(syn480h)
ElbowPlot(syn480h,ndims=50)

syn480h <- RunUMAP(syn480h, dims = 1:15,spread=3,min.dist=0.1)
DimPlot(syn480h,reduction = "umap")
p1 <- DimPlot(syn480h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn480h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot(syn480h,features=c("UMI_count"),max.cutoff = 2000)
FeaturePlot(syn480h,features=c("Hpca"))
#pK_selection
sweep.res.list_syn480h <- paramSweep_v3(syn480h, PCs = 1:15, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn480h <- summarizeSweep(sweep.res.list_syn480h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn480h <- find.pK(sweep.stats_syn480h)
mpK<-as.numeric(as.vector(bcmvn_syn480h$pK[which.max(bcmvn_syn480h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn480h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn480h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn480h <- doubletFinder_v3(syn480h, PCs = 1:15, 
                            pN = 0.25, pK = mpK, 
                            nExp = nExp_poi.adj, 
                            reuse.pANN = FALSE ,
                            sct = T)
table(syn480h@meta.data[["DF.classifications_0.25_0.06_199"]])
syn480h[["doublet"]] <-syn480h[["DF.classifications_0.25_0.06_199"]]


#########
#########
syn514h.data <- read.delim("Combined_514syn_transcript_counts_rRNA_MT_filtered_gene_symbol_noDup.tsv", row.names=1)
syn514h <- CreateSeuratObject(counts = syn514h.data, project = "syn514h", min.cells = 3, min.features = 100)
##QC
syn514h[["percent.mt"]] <- PercentageFeatureSet(syn514h, pattern = "^MT-")
syn514h[["UMI_count"]] <-syn514h[["nCount_RNA"]] 
syn514h[["Gene_count"]] <-syn514h[["nFeature_RNA"]] 
a=syn514h[["UMI_count"]]$UMI_count
quantile(a, 0.99) 
syn514h <- subset(syn514h, subset = nCount_RNA<quantile(a, 0.99))
syn514h[["patient_ID"]]="P514"

syn514h <- SCTransform(syn514h, method = "glmGamPoi", verbose = TRUE,ncells = 20000)

syn514h <- RunPCA(syn514h)
ElbowPlot(syn514h,ndims=50)

syn514h <- RunUMAP(syn514h, dims = 1:15,spread=3,min.dist=0.1)
DimPlot(syn514h,reduction = "umap")
p1 <- DimPlot(syn514h, reduction = "umap", group.by = "patient_ID")
p2 <- DimPlot(syn514h, reduction = "umap", group.by = "doublet")
p1 + p2
FeaturePlot(syn514h,features=c("UMI_count"),max.cutoff = 2000)
FeaturePlot(syn514h,features=c("Hpca"))
#pK_selection
sweep.res.list_syn514h <- paramSweep_v3(syn514h, PCs = 1:15, sct = T)
#head(sweep.res.list_syn5844h)
sweep.stats_syn514h <- summarizeSweep(sweep.res.list_syn514h, GT = FALSE)
#head(sweep.stats_syn5844h)
bcmvn_syn514h <- find.pK(sweep.stats_syn514h)
mpK<-as.numeric(as.vector(bcmvn_syn514h$pK[which.max(bcmvn_syn514h$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- syn514h@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(syn514h@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

syn514h <- doubletFinder_v3(syn514h, PCs = 1:15, 
                            pN = 0.25, pK = mpK, 
                            nExp = nExp_poi.adj, 
                            reuse.pANN = FALSE ,
                            sct = T)
table(syn514h@meta.data[["DF.classifications_0.25_0.08_201"]])
syn514h[["doublet"]] <-syn514h[["DF.classifications_0.25_0.08_201"]]


#########
#########
##Merge datasets
syn.combined <- merge(syn417h, y = c(syn418h,syn480h, syn514h), project = "syn.combined")

#syn.combined <- merge(syn418h, y = c(syn480h,syn514h), project = "syn.combined")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(syn.combined, split.by = "patient_ID")

# normalize and identify variable features for each dataset independently
#ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 2000)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
#})
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x<- SCTransform(x, method = "glmGamPoi", verbose = TRUE,ncells=20000)
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
immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 50, verbose = T)
ElbowPlot(immune.combined.sct,ndims = 50)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:16, 
                               spread=3,min.dist=0.1,metric="cosine")


# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "patient_ID",shuffle=T)
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "doublet")
p3 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 +p3



