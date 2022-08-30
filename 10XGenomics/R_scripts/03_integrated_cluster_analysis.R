## Integration of tBM and utBM Brugia malayi mf datasets based on Seurat integration tutorial

install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
install.packages('pals')
install.packages('pheatmap')
BiocManager::install("dittoSeq")
install.packages("aplot")
---------------------------
  
library(multtest)
library(metap)
library(pals)
library(ggtext)
library(pheat)
library(dittoSeq)
library(aplot)


## Perform the integration starting with the utBM_subset and tBM_subset processed Seurat objects
utBM_subset <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BMSinglecell-ms/10XGenomics/scmulti_utBM.rds")
tBM_subset <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BMSinglecell-ms/10XGenomics/scmulti_tBM.rds")

# Identify anchors in both datasets with the FindIntegrationAnchors function 
# Completes canonical correlation analysis (CCA) by using PCA to find the greaest sources of shared variation across the groups and then identifies the anchors across the datasets
# The anchors or mutual nearest neighbors (MNNs) are cellular 'best buddies' --> looking for a cells closest neighbor in the other condition based on gene expression and if the reciprocal cell also identifies the original cell as a neighbor, the cells are marked as anchors
# incorrect anchors are removed by comparing the similarity between anchor pairs by the overlap in the local neighborhoods - do the adjacent cell shave "best buddies" that are adjacent to each other?

#merge the two seurat objects
combined <- merge(utBM_subset, y = tBM_subset)

# create split object list
obj.list <- SplitObject(combined, split.by = "orig.ident")

# log normalize and selecte variable features from each dataset
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = TRUE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
})

# select the features
features <- SelectIntegrationFeatures(object.list = obj.list)

# run PCA on each object 
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# determine the "best buddy" features for integration based on the features identified above
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, verbose = FALSE)

# integrate the datasets based on anchor features
new_combined <- IntegrateData(anchorset = anchors, verbose = FALSE)


DefaultAssay(new_combined) <- "integrated"
new_combined <- ScaleData(new_combined, verbose = FALSE)
new_combined <- RunPCA(new_combined, features = VariableFeatures(new_combined),verbose = FALSE)
new_combined <- RunUMAP(new_combined, reduction = "pca", dims = 1:30, verbose=FALSE)

DimPlot(new_combined, group.by = "orig.ident")


# cluster the cells
new_combined <- FindNeighbors(new_combined, reduction = "pca", dims = 1:50, verbose = FALSE)
new_combined <- FindClusters(new_combined, resolution = 0.5, verbose = FALSE)
p1 <- DimPlot(new_combined, pt.size = 0.3, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(new_combined, pt.size = 0.3, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


#saveRDS(new_combined, "~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BmSinglecell-ms/10XGenomics/scmulti_integrated.RDS")

new_combined <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BmSinglecell-ms/10XGenomics/scmulti_integrated.RDS")

# how many cells are in each cluster and what proportion of those cells is represented by tBM/utBM
table(new_combined@active.ident, group_by = new_combined@meta.data$orig.ident)




################################################
#Differential Gene Expression Exploration
# DEG tests should be run on the unintegrated data since integration inherently introduces dependencies between data points and violates assumptions of the statistical tests used for DEG.

# Find differential expressed genes (markers) for identity classes. Uses Wilcoxon Rank Sum test with bonferroni p-value correction to account for false positives. Find conserved markers will find markers that are conserved between tBM and utBM but are differential expressed compared to all other clusters. 
BM_marker_0 <- FindConservedMarkers(new_combined, ident.1 = 0, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_1 <- FindConservedMarkers(new_combined, ident.1 = 1, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_2 <- FindConservedMarkers(new_combined, ident.1 = 2, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_3 <- FindConservedMarkers(new_combined, ident.1 = 3, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_4 <- FindConservedMarkers(new_combined, ident.1 = 4, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_5 <- FindConservedMarkers(new_combined, ident.1 = 5, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_6 <- FindConservedMarkers(new_combined, ident.1 = 6, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_7 <- FindConservedMarkers(new_combined, ident.1 = 7, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_8 <- FindConservedMarkers(new_combined, ident.1 = 8, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_9 <- FindConservedMarkers(new_combined, ident.1 = 9, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_10 <- FindConservedMarkers(new_combined, ident.1 = 10, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_11 <- FindConservedMarkers(new_combined, ident.1 = 11, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_12 <- FindConservedMarkers(new_combined, ident.1 = 12, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_13 <- FindConservedMarkers(new_combined, ident.1 = 13, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_14 <- FindConservedMarkers(new_combined, ident.1 = 14, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_15 <- FindConservedMarkers(new_combined, ident.1 = 15, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_16 <- FindConservedMarkers(new_combined, ident.1 = 16, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_17 <- FindConservedMarkers(new_combined, ident.1 = 17, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_18 <- FindConservedMarkers(new_combined, ident.1 = 18, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_19 <- FindConservedMarkers(new_combined, ident.1 = 19, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_20 <- FindConservedMarkers(new_combined, ident.1 = 20, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_21 <- FindConservedMarkers(new_combined, ident.1 = 21, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_22 <- FindConservedMarkers(new_combined, ident.1 = 22, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_23 <- FindConservedMarkers(new_combined, ident.1 = 23, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_24 <- FindConservedMarkers(new_combined, ident.1 = 24, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_25 <- FindConservedMarkers(new_combined, ident.1 = 25, grouping.var = "orig.ident", verbose = TRUE)
BM_marker_26 <- FindConservedMarkers(new_combined, ident.1 = 26, grouping.var = "orig.ident", verbose = TRUE)


BM_marker_0<- rownames_to_column(BM_marker_0, var = "gene_id")
BM_marker_0['Cluster'] <- "1"

BM_marker_1<- rownames_to_column(BM_marker_1, var = "gene_id")
BM_marker_1['Cluster'] <- "2"

BM_marker_2<- rownames_to_column(BM_marker_2, var = "gene_id")
BM_marker_2['Cluster'] <- "3"

BM_marker_3<- rownames_to_column(BM_marker_3, var = "gene_id")
BM_marker_3['Cluster'] <- "4"

BM_marker_4<- rownames_to_column(BM_marker_4, var = "gene_id")
BM_marker_4['Cluster'] <- "5"

BM_marker_5<- rownames_to_column(BM_marker_5, var = "gene_id")
BM_marker_5['Cluster'] <- "6"

BM_marker_6<- rownames_to_column(BM_marker_6, var = "gene_id")
BM_marker_6['Cluster'] <- "7"

BM_marker_7<- rownames_to_column(BM_marker_7, var = "gene_id")
BM_marker_7['Cluster'] <- "8"

BM_marker_8<- rownames_to_column(BM_marker_8, var = "gene_id")
BM_marker_8['Cluster'] <- "9"

BM_marker_9<- rownames_to_column(BM_marker_9, var = "gene_id")
BM_marker_9['Cluster'] <- "10"

BM_marker_10<- rownames_to_column(BM_marker_10, var = "gene_id")
BM_marker_10['Cluster'] <- "11"

BM_marker_11<- rownames_to_column(BM_marker_11, var = "gene_id")
BM_marker_11['Cluster'] <- "12"

BM_marker_12<- rownames_to_column(BM_marker_12, var = "gene_id")
BM_marker_12['Cluster'] <- "13"

BM_marker_13<- rownames_to_column(BM_marker_13, var = "gene_id")
BM_marker_13['Cluster'] <- "14"

BM_marker_14<- rownames_to_column(BM_marker_14, var = "gene_id")
BM_marker_14['Cluster'] <- "15"

BM_marker_15<- rownames_to_column(BM_marker_15, var = "gene_id")
BM_marker_15['Cluster'] <- "16"

BM_marker_16<- rownames_to_column(BM_marker_16, var = "gene_id")
BM_marker_16['Cluster'] <- "17"

BM_marker_17<- rownames_to_column(BM_marker_17, var = "gene_id")
BM_marker_17['Cluster'] <- "18"

BM_marker_18<- rownames_to_column(BM_marker_18, var = "gene_id")
BM_marker_18['Cluster'] <- "19"

BM_marker_19<- rownames_to_column(BM_marker_19, var = "gene_id")
BM_marker_19['Cluster'] <- "20"

BM_marker_20<- rownames_to_column(BM_marker_20, var = "gene_id")
BM_marker_20['Cluster'] <- "21"

BM_marker_21<- rownames_to_column(BM_marker_21, var = "gene_id")
BM_marker_21['Cluster'] <- "22"

BM_marker_22<- rownames_to_column(BM_marker_22, var = "gene_id")
BM_marker_22['Cluster'] <- "23"

BM_marker_23<- rownames_to_column(BM_marker_23, var = "gene_id")
BM_marker_23['Cluster'] <- "24"

BM_marker_24<- rownames_to_column(BM_marker_24, var = "gene_id")
BM_marker_24['Cluster'] <- "25"

BM_marker_25<- rownames_to_column(BM_marker_25, var = "gene_id")
BM_marker_25['Cluster'] <- "26"

BM_marker_26<- rownames_to_column(BM_marker_26, var = "gene_id")
BM_marker_26['Cluster'] <- "27"



#combine all identified markers into single data frame
conserved_markers <- rbind(BM_marker_0, BM_marker_1, BM_marker_2, BM_marker_3, BM_marker_4, BM_marker_5, BM_marker_6, BM_marker_7, BM_marker_8, BM_marker_9, BM_marker_10, BM_marker_11, BM_marker_12, BM_marker_13, BM_marker_14, BM_marker_15, BM_marker_16, BM_marker_17, BM_marker_18, BM_marker_19, BM_marker_20, BM_marker_21, BM_marker_22, BM_marker_23, BM_marker_24, BM_marker_25, BM_marker_26)


# export data table as csv
write.csv(conserved_markers, "~/Desktop/BM_combined_conservedmarkers.csv")

orthos <- select(orthos, c(-bma_genome_project, -cel_gene_name))
colnames(orthos) <- c("gene_id", "cel_gene_id")











