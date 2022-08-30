# sub-setting the untreated Bma sc UMAP data to contain only the genes that are orthologs to Cel. for input into Scanorama

# install packages
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
}
#remotes::install_github("mojaveazure/seurat-disk")
#devtools::install_github('satijalab/seurat-data')

# load required libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# read in untreated bma (utBM) object
utBM_subset <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BmSinglecell-ms/10XGenomics/scmulti_utBM_subset.rds")


#import csv of cel and bma genes
orthos <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/sc_R_scripts/cel_orthos.csv")


#remove rows that are not 1:1 orthologus in cel
orthos_tmp <- orthos %>% 
  na_if("") %>% 
  na.omit %>% 
  select(-bma_genome_project, -cel_gene_name)%>%
  orthos_tmp[!duplicated(orthos_tmp$bma_gene_ID),]%>%
  select(orthos_tmp, bma_gene_ID)


#create vector to subset seurat object
bma_IDs <- bma_IDs$bma_gene_ID


#subset seurat object based on 1:1 orthologs using the bma gene IDs
utBM_tmp <-subset(x = utBM_subset, features = bma_IDs)



# replace Bma WBGene IDs in the seurat object with the C. elegans orthologs WBGene IDs

X <- data.frame(utBM_tmp@assays$RNA@counts@Dimnames[[1]])

colnames(X) <- "bma_gene_ID"

tmp <- X %>% 
  left_join(orthos_tmp, by = "bma_gene_ID") %>% 
  mutate(X2 = ifelse(is.na(cel_gene_ID), bma_gene_ID, cel_gene_ID))

tmp <- tmp$X2

utBM_tmp@assays$RNA@counts@Dimnames[[1]] <- tmp





# export new seurat object with Cel WBGene IDs and convert to AnnData structure (h5ad)

SaveH5Seurat(utBM_tmp, filename = "~/Desktop/utBM_tmp.h5Seurat")
Convert("~/Desktop/utBM_tmp.h5Seurat", dest = "h5ad")



