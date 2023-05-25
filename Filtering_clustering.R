###Merge Bb all weeks -------------
.libPaths(c("/home/mw232/~R/-4.0.1", "/n/app/R/4.1.1/lib64/R/library"))
setwd("/home/mw232/Bb_scrnaseq/")
.libPaths()
library(Seurat)
library(tidyverse)
library(DropletUtils)
set.seed(661851) # needed for empty drops reproducibility

#Function to convert matrix files to empty droplets, filter out empty drops and load Seurat
droplet_to_seurat <- function(data.dir, FDR = 0.01, project = "Seurat", min.cells = 3, min.features = 0) {
  X_mat <- Read10X(data.dir = data.dir)
  emptydrops <- emptyDrops(X_mat)
  filtered_cells <- rownames(emptydrops)[emptydrops$FDR < FDR]
  filtered_cells <- na.omit(filtered_cells)
  filtered_matrix <- X_mat[,filtered_cells]
  Seurat <- CreateSeuratObject(filtered_matrix, project = project,
                               min.cells = min.cells, min.features = min.features)
  return(Seurat)
}

wk0pos <- droplet_to_seurat(data.dir = "cd45_pos_0wk/raw_feature_bc_matrix/", project = "wk0pos")
wk0neg <- droplet_to_seurat(data.dir = "cd45_neg_0wk/raw_feature_bc_matrix/", project = "wk0neg")
wk2pos <- droplet_to_seurat(data.dir = "cd45_pos_2wk/raw_feature_bc_matrix/", project = "wk2pos")
wk2neg <- droplet_to_seurat(data.dir = "cd45_neg_2wk/raw_feature_bc_matrix/", project = "wk2neg")
wk4pos <- droplet_to_seurat(data.dir = "cd45_pos_4wk/raw_feature_bc_matrix/", project = "wk4pos")
wk4neg <- droplet_to_seurat(data.dir = "cd45_neg_4wk/raw_feature_bc_matrix/", project = "wk4neg")
wk6pos <- droplet_to_seurat(data.dir = "cd45_pos_6wk/raw_feature_bc_matrix/", project = "wk6pos")
wk6neg <- droplet_to_seurat(data.dir = "cd45_neg_6wk/raw_feature_bc_matrix/", project = "wk6neg")
wk8pos <- droplet_to_seurat(data.dir = "cd45_pos_8wk/raw_feature_bc_matrix/", project = "wk8pos")
wk8neg <- droplet_to_seurat(data.dir = "cd45_neg_8wk/raw_feature_bc_matrix/", project = "wk8neg")

combined <- merge(x = wk0pos,
                  y = c(wk0neg, wk2pos, wk2neg, wk4pos, wk4neg, wk6pos, wk6neg, wk8pos, wk8neg),
                  add.cell.ids = c("wk0_CD45+", "wk0_CD45-", "wk2_CD45+", "wk2_CD45-", "wk4_CD45+", "wk4_CD45-", "wk6_CD45+", "wk6_CD45-", "wk8_CD45+", "wk8_CD45-"),
                  project = "combined")


saveRDS(combined, file = "all_weeks_combined_MW.rds")

print("Saved RDS file following merging of dropletUtils-filtered matrices")

# Filtering ---------------------------------------

#mitochodrial genes and total gene filtering
paste("Pre-filtering = ", paste0(dim(combined), collapse = " x ")) #20642 x 167725
combined$perc_mito <- PercentageFeatureSet(combined, "^mt-")
mito_cutoff <- mean(combined$perc_mito) + 2 * sd(combined$perc_mito)
paste("Mitochondrial cutoff = ", mito_cutoff) #12.812607928373
combined <- combined[,combined$perc_mito <= mito_cutoff] #remove cells above cutoff (12.81)
paste("Post mitochondrial filtering = ", paste0(dim(combined), collapse = " x ")) #20642 x 163351
combined <- combined[,combined$nFeature_RNA > 500] #remove cells that have fewer than 500 genes
paste("Post mitochondrial and gene filtering = ", paste0(dim(combined), collapse = " x ")) #20642 x 119086

combined_list <- SplitObject(combined, split.by = "orig.ident")
combined_list <- lapply(X = combined_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
print("Seurat object split and normalized")

features <- SelectIntegrationFeatures(object.list = combined_list)
combined_list <- lapply(X = combined_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = combined_list, reference = c(1, 2), reduction = "rpca", dims = 1:50)

# Integration -------------------------------------------------------
combined_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
print("Integration complete")

saveRDS(combined_integrated, file = "all_weeks_integrated_MW.rds")
remove(list = setdiff(ls(), "combined_integrated"))
gc()
print("Saved integrated file")

combined_integrated <- ScaleData(combined_integrated, verbose = FALSE)
print("Integrated Data Scaled")
combined_integrated <- RunPCA(combined_integrated, verbose = FALSE)
print("PCA run")
combined_integrated <- RunUMAP(combined_integrated, dims = 1:50)
print("UMAP run")

saveRDS(combined_integrated, file = "all_weeks_integrated_MW.rds")
print("File saved")

combined_integrated <- FindNeighbors(combined_integrated, dims = 1:50)
combined_integrated <- FindClusters(combined_integrated)
print("Find Neighbors and Clusters")

DefaultAssay(combined_integrated) <- "RNA"

markers <- FindAllMarkers(combined_integrated, only.pos = TRUE)
combined_integrated@misc$clusters_df <- markers
saveRDS(combined_integrated, "all_weeks_integrated_MW.rds")
print("Clusters found and file saved")

# add week metadata
combined_integrated@meta.data$week <- str_extract(rownames(combined_integrated@meta.data),"wk\\d")

# Merge clusters ----------------------------------------------------------
merged_clusters <- combined_integrated
merged_clusters <- RenameIdents(merged_clusters, '0' = "Fibroblast", "2" = "Fibroblast", "3" = "Fibroblast", "4" = "Fibroblast",
                                "5" = "Fibroblast", "8" = "Fibroblast", "9" = "Fibroblast",
                                "11" = "Fibroblast", "18" = "Fibroblast", "23" = "Fibroblast",
                                "25" = "Fibroblast", "32" = "Fibroblast", "34" = "Fibroblast", "35" = "Fibroblast")
merged_clusters <- RenameIdents(merged_clusters, "1" = "Myeloid", "6" = "Myeloid", "7" = "Myeloid",
                                "12" = "Myeloid", "13" = "Myeloid", "19" = "Myeloid", 
                                "22" = "Myeloid", "28" = "Myeloid", "29" = "Myeloid",
                                "36" = "Myeloid", "37" = "Myeloid")
merged_clusters <- RenameIdents(merged_clusters, "10" = "Lymphoid", "14" = "Lymphoid",
                                "15" = "Lymphoid", "17" = "Lymphoid", "20" = "Lymphoid",
                                "26" = "Lymphoid", "27" = "Lymphoid", "30" = "Lymphoid",
                                "31" = "Lymphoid", "33" = "Lymphoid")
merged_clusters <- RenameIdents(merged_clusters, "16" = "Vascular",
                                "21" = "Cycling cells", "24" = "Vascular", "38" = "Vascular")

#Reorder based on size of population
levels(merged_clusters@active.ident)
my_levels <- c("Fibroblast", "Myeloid", "Lymphoid", "Vascular", "Cycling cells")
merged_clusters@active.ident <- factor(x = merged_clusters@active.ident, levels = my_levels)

#Rename clusters
merged_clusters$finalcluster <- Idents(merged_clusters)

# Subset vascular cells -----------------------------------------------
vascular <- subset(merged_clusters, idents = "Vascular")
Idents(vascular) <- "integrated_snn_res.0.8"
vascular <- RenameIdents(vascular, "16" = "Adipocytes", "24" = "Endothelial cells", "38" = "Red blood cells")
vascular$vascular_final <- Idents(vascular)
my_levels <- c("Adipocytes", "Endothelial cells", "Red blood cells")
vascular@active.ident <- factor(x = vascular@active.ident, levels = my_levels)
vascular$vascular_final <- Idents(vascular)

# Subclustering myeloid cells ---------------------------------------------

myeloid <- subset(merged_clusters, idents = "Myeloid")
DefaultAssay(myeloid) <- "RNA"
myeloid <- FindNeighbors(myeloid, dims = 1:30)
DefaultAssay(myeloid) <- "integrated"
myeloid <- FindClusters(myeloid, resolution = 0.15)
myeloid <- RenameIdents(myeloid, "0" = "Macrophages/Monocytes", "1" = "Macrophages/Monocytes",
                        "2" = "Macrophages/Monocytes", "3" = "Neutrophils",
                        "4" = "Neutrophils", "5" = "Dendritic cells", 
                        "6" = "Macrophages/Monocytes", "7" = "Neutrophils", 
                        "8" = "Neutrophils", "9" = "Macrophages/Monocytes")
myeloid$myeloid_final <- Idents(myeloid)

## Mac/Mono initial subclustering (doublet removal) --------------------------
mac_mono <- subset(myeloid, ident = "Macrophages/Monocytes")
DefaultAssay(mac_mono) <- "RNA"
mac_mono <- FindVariableFeatures(mac_mono, selection.method = "vst", nfeatures = 2000)
DefaultAssay(mac_mono) <- "integrated"
mac_mono <- ScaleData(mac_mono, vars.to.regress = c("nFeature_RNA", "perc_mito"))
mac_mono <- RunPCA(mac_mono)
mac_mono <- FindNeighbors(mac_mono, dims = 1:30)
mac_mono <- FindClusters(mac_mono, resolution = 0.2)

mac_mono <- RenameIdents(mac_mono, "0" = "Macrophage/Monocyte",
                         "1" = "Macrophage/Monocyte",
                         "2" = "Macrophage/Monocyte", "3" = "Macrophage/Monocyte", 
                         "4" = "Doublets", "5" = "Doublets",
                         "6" = "Macrophage/Monocyte",
                         "7" = "Doublets", "8" = "Macrophage/Monocyte")

#### Mac/Mono reclustering post-doublets -----------------------------------

mac_nodoublets <- subset(mac_mono, ident = "Macrophage/Monocyte")
DefaultAssay(mac_nodoublets) <- "RNA"
mac_nodoublets <- FindVariableFeatures(mac_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(mac_nodoublets) <- "integrated"
mac_nodoublets <- ScaleData(mac_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
mac_nodoublets <- RunPCA(mac_nodoublets)
mac_nodoublets <- FindNeighbors(mac_nodoublets, dims = 1:30)
mac_nodoublets <- FindClusters(mac_nodoublets, resolution = 0.2)
mac_nodoublets <- RunUMAP(mac_nodoublets, dims = 1:30)
#merge very similar clusters
mac_nodoublets <- RenameIdents(mac_nodoublets, "0" = "0", "1" = "1", "2" = "2", "3" = "2", "4" = "3", "5" = "1", "6" = "4", "7" = "5", "8" = "4")
# Final naming
mac_nodoublets <- RenameIdents(mac_nodoublets, "2" = "MacMono.1", 
                               "0" = "MacMono.2", "1" = "MacMono.3", "3" = "MacMono.4",
                               "4" = "MacMono.5", "5" = "MacMono.6")
mac_nodoublets$mac_nd_clusters <- Idents(mac_nodoublets)
DefaultAssay(mac_nodoublets) <- "RNA"
# find new markers
macrophage_markers <- FindAllMarkers(mac_nodoublets, only.pos = TRUE)
mac_nodoublets@misc$macrophage_markers <- macrophage_markers # add to seurat object to save

saveRDS(mac_nodoublets, file = "mac_nodoublets.rds")


## Neutrophil initial subclustering (doublet removal) ----------------------------------------------
neut <- subset(myeloid, ident = "Neutrophils")
DefaultAssay(neut) <- "RNA"
neut <- FindVariableFeatures(neut, selection.method = "vst", nfeatures = 2000)
DefaultAssay(neut) <- "integrated"
neut <- ScaleData(neut, vars.to.regress = c("nFeature_RNA", "perc_mito"))
neut <- RunPCA(neut)
neut <- FindNeighbors(neut, dims = 1:30)
neut <- FindClusters(neut, resolution = 0.3)
neut <- RenameIdents(neut, "0" = "Neutrophils", "1" = "Neutrophils", "2" = "Doublets",
                     "3" = "Doublets", "4" = "Neutrophils", "5" = "Doublets", "6" = "Neutrophils", "7" = "Neutrophils", "8" = "Neutrophils", "9" = "Doublets")


### neutrophil reclustering post-doublets -----------------------------------
neut_nodoublets <- subset(neut, ident = "Neutrophils")
DefaultAssay(neut_nodoublets) <- "RNA"
neut_nodoublets <- FindVariableFeatures(neut_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(neut_nodoublets) <- "integrated"
neut_nodoublets <- ScaleData(neut_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
neut_nodoublets <- RunPCA(neut_nodoublets)
neut_nodoublets <- FindNeighbors(neut_nodoublets, dims = 1:30)
neut_nodoublets <- FindClusters(neut_nodoublets, resolution = 0.2)
neut_nodoublets <- RunUMAP(neut_nodoublets, dims = 1:30)

neut_nodoublets <- RenameIdents(neut_nodoublets, "0" = "Neut.1",
                                "1" = "Neut.2", "2" = "Neut.3",
                                "3" = "Neut.4", "4" = "Neut.5", "5" = "Neut.6")
neut_nodoublets$neut_nd_clusters <- Idents(neut_nodoublets)
# find new markers
DefaultAssay(neut_nodoublets) <- "RNA"
Neut_markers <- FindAllMarkers(neut_nodoublets, only.pos = TRUE)
neut_nodoublets@misc$neut_nd_clusters <- Neut_markers # add to seurat object to save
saveRDS(neut_nodoublets, file = "neut_nodoublets.rds")

neutrophil_markers <- FindAllMarkers(neut_nodoublets, only.pos = TRUE)



# Subclustering lymphoid cells ---------------------------------------------
lymphoid <- subset(merged_clusters, idents = "Lymphoid")
Idents(lymphoid) <- "integrated_snn_res.0.8"
lymphoid <- RenameIdents(lymphoid, "10" = "B cells", "14" = "B cells",
                         "15" = "B cells", "17" = "T cells",
                         "20" = "B cells", "26" = "NK cells",
                         "27" = "B cells", "30" = "Fibroblasts", # upon closer inspection
                         "31" = "T cells", "33" = "ILCs")
lymphoid$lymphoid_final <- Idents(lymphoid)

## T cell initial subclustering (doublet removal) --------------------------
Tcells <- subset(lymphoid, ident = "T cells")
DefaultAssay(Tcells) <- "RNA"
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Tcells) <- "integrated"
Tcells <- ScaleData(Tcells, vars.to.regress = c("nFeature_RNA", "perc_mito"))
Tcells <- RunPCA(Tcells)
Tcells <- FindNeighbors(Tcells, dims = 1:30)
Tcells <- FindClusters(Tcells, resolution = 0.3)
Tcells <- RenameIdents(Tcells, "0" = "T cells", "1" = "T cells", "2" = "T cells",
                       "3" = "Doublets", "4" = "Doublets", "5" = "Doublets",
                       "6" = "T cells", "7" = "T cells", "8" = "Doublets")


### T cell reclustering post-doublets -----------------------------------
set.seed(952850)
Tcell_nodoublets <- subset(Tcells, ident = "T cells")
DefaultAssay(Tcell_nodoublets) <- "RNA"
Tcell_nodoublets <- FindVariableFeatures(Tcell_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Tcell_nodoublets) <- "integrated"
Tcell_nodoublets <- ScaleData(Tcell_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
Tcell_nodoublets <- RunPCA(Tcell_nodoublets)
Tcell_nodoublets <- FindNeighbors(Tcell_nodoublets, dims = 1:30)
Tcell_nodoublets <- FindClusters(Tcell_nodoublets, resolution = 0.2)
Tcell_nodoublets <- RunUMAP(Tcell_nodoublets, dims = 1:30)
Tcell_nodoublets <- RenameIdents(Tcell_nodoublets, "0" = "Tcell.1",
                                 "1" = "Tcell.2", "2" = "Tcell.3", "3" = "Tcell.4")
Tcell_nodoublets$Tcell_nd_clusters <- Idents(Tcell_nodoublets)
DefaultAssay(Tcell_nodoublets) <- "RNA"
#find new markers
Tcell_markers <- FindAllMarkers(Tcell_nodoublets, only.pos = TRUE)
Tcell_nodoublets@misc$Tcell_markers <- Tcell_markers # save to seurat object
saveRDS(Tcell_nodoublets, file = "Tcell_nodoublets.rds")

## B cell initial subclustering (doublet removal) --------------------------
Bcells <- subset(lymphoid, ident = "B cells")
DefaultAssay(Bcells) <- "RNA"
Bcells <- FindVariableFeatures(Bcells, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Bcells) <- "integrated"
Bcells <- ScaleData(Bcells, vars.to.regress = c("nFeature_RNA", "perc_mito"))
Bcells <- RunPCA(Bcells)
Bcells <- FindNeighbors(Bcells, dims = 1:30)
Bcells <- FindClusters(Bcells, resolution = 0.2)

Bcells <- RenameIdents(Bcells, "0" = "B cells", "1" = "B cells", 
                       "2" = "B cells", "3" = "B cells",
                       "4" = "B cells", "5" = "Doublets",
                       "6" = "Doublets", "7" = "Doublets")

# B cell reclustering post-doublets -----------------------------------
Bcell_nodoublets <- subset(Bcells, ident = "B cells")
DefaultAssay(Bcell_nodoublets) <- "RNA"
Bcell_nodoublets <- FindVariableFeatures(Bcell_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Bcell_nodoublets) <- "integrated"
Bcell_nodoublets <- ScaleData(Bcell_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
Bcell_nodoublets <- RunPCA(Bcell_nodoublets)
Bcell_nodoublets <- FindNeighbors(Bcell_nodoublets, dims = 1:30)
Bcell_nodoublets <- FindClusters(Bcell_nodoublets, resolution = 0.3)
Bcell_nodoublets <- RunUMAP(Bcell_nodoublets, dims = 1:30)

#merge very similar clusters
Bcell_nodoublets <- RenameIdents(Bcell_nodoublets, "0" = "0", "1" = "1", 
                                 "2" = "2", "3" = "3", "4" = "4", "5" = "5", "6" = "1", "7" = "2")
# Final naming
Bcell_nodoublets <- RenameIdents(Bcell_nodoublets, "0" = "Bcell.1", "1" = "Bcell.2",
                                 "2" = "Bcell.3", "3" = "Bcell.4",
                                 "4" = "Bcell.5", "5" = "Bcell.6")

Bcell_nodoublets$Bcell_nd_clusters <- Idents(Bcell_nodoublets)
DefaultAssay(Bcell_nodoublets) <- "RNA"
#find new markers
Bcell_markers <- FindAllMarkers(Bcell_nodoublets, only.pos = TRUE)
Bcell_nodoublets@misc$Bcell_markers <- Bcell_markers # save to seurat object

saveRDS(Bcell_nodoublets, file = "Bcell_nodoublets.rds")


# Subclustering Non-immune/fibroblast-like cells  --------------------------

fibroblasts <- subset(merged_clusters, idents = "Fibroblast")
DefaultAssay(fibroblasts) <- "RNA"
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:30)
DefaultAssay(fibroblasts) <- "integrated"
fibroblasts <- FindClusters(fibroblasts, resolution = 0.2)

fibroblasts <- RenameIdents(fibroblasts, "0" = "Fibroblasts", "1" = "Fibroblasts",
                            "2" = "Fibroblasts", 
                            "3" = "Fibroblasts",
                            "4" = "Synoviocytes",
                            "5" = "Muscle cells", "6" = "Osteoblasts",
                            "7" = "Fibroblasts", "8" = "Chondrocytes")
fibroblasts$fibroblast_final <- Idents(fibroblasts)


## Fibroblast initial subclustering (doublet removal) ----------------------------------
fib_only <- subset(fibroblasts, ident = "Fibroblasts")
DefaultAssay(fib_only) <- "RNA"
fib_only <- FindVariableFeatures(fib_only, selection.method = "vst", nfeatures = 2000)
DefaultAssay(fib_only) <- "integrated"
fib_only <- ScaleData(fib_only, vars.to.regress = c("nFeature_RNA", "perc_mito"))
fib_only <- RunPCA(fib_only)
fib_only <- FindNeighbors(fib_only, dims = 1:30)
fib_only <- FindClusters(fib_only, resolution = c0.2)
fib_only <- RenameIdents(fib_only, "0" = "Fibroblast", "1" = "Fibroblast",
                         "2" = "Doublets", "3" = "Fibroblast", 
                         "4" = "Fibroblast", "5" = "Fibroblast",
                         "6" = "Doublets", "7" = "Doublets", "8" = "Doublets", "9" = "Doublets")

### fibroblast reclustering post-doublets -----------------------------------
fib_nodoublets <- subset(fib_only, ident = "Fibroblast")
DefaultAssay(fib_nodoublets) <- "RNA"
fib_nodoublets <- FindVariableFeatures(fib_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(fib_nodoublets) <- "integrated"
fib_nodoublets <- ScaleData(fib_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
fib_nodoublets <- RunPCA(fib_nodoublets)
fib_nodoublets <- FindNeighbors(fib_nodoublets, dims = 1:30)
fib_nodoublets <- FindClusters(fib_nodoublets, resolution = 0.2)

fib_nodoublets <- RunUMAP(fib_nodoublets, dims = 1:30)

fib_nodoublets <- RenameIdents(fib_nodoublets, "0" = "Fib.1", 
                               "1" = "Fib.2", "2" = "Fib.3",
                               "3" = "Fib.4", "4" = "Fib.5", 
                               "5" = "Fib.6", "6" = "Fib.7", "7" = "Fib.8")

fib_nodoublets$fib_nd_clusters <- Idents(fib_nodoublets)
#find new markers
DefaultAssay(fib_nodoublets) <- "RNA"
Fib_markers <- FindAllMarkers(fib_nodoublets, only.pos = T)
fib_nodoublets@misc$Fib_markers <- Fib_markers # add markers to seurat object

saveRDS(fib_nodoublets, file = "fib_nodoublets.rds")

# synoviocyte initial reclustering (doublet removal) ----------------------------------
synov <- subset(fibroblasts, ident = "Synoviocytes")
DefaultAssay(synov) <- "RNA"
synov <- FindVariableFeatures(synov, selection.method = "vst", nfeatures = 2000)
DefaultAssay(synov) <- "integrated"
synov <- ScaleData(synov, vars.to.regress = c("nFeature_RNA", "perc_mito"))
synov <- RunPCA(synov)
synov <- FindNeighbors(synov, dims = 1:30)
synov <- FindClusters(synov, resolution = 0.25)

synov <- RenameIdents(synov, "0" = "Synoviocytes",
                      "1" = "Synoviocytes", "2" = "Synoviocytes",
                      "3" = "Doublets", "4" = "Doublets",
                      "5" = "Doublets", "6" = "Doublets", "7" = "Doublets")

# Synoviocyte reclustering post-doublets -----------------------------------
synov_nodoublets <- subset(synov, ident = "Synoviocytes")
DefaultAssay(synov_nodoublets) <- "RNA"
synov_nodoublets <- FindVariableFeatures(synov_nodoublets, selection.method = "vst", nfeatures = 2000)
DefaultAssay(synov_nodoublets) <- "integrated"
synov_nodoublets <- ScaleData(synov_nodoublets, vars.to.regress = c("nFeature_RNA", "perc_mito"))
synov_nodoublets <- RunPCA(synov_nodoublets)
synov_nodoublets <- FindNeighbors(synov_nodoublets, dims = 1:30)
synov_nodoublets <- FindClusters(synov_nodoublets, resolution = 0.2)
synov_nodoublets <- RunUMAP(synov_nodoublets, dims = 1:30)

synov_nodoublets <- RenameIdents(synov_nodoublets, "0" = "Synoviocyte.1", "1" = "Synoviocyte.2", "2" = "Synoviocyte.3", "3" = "Synoviocyte.4", "4" = "Synoviocyte.5")
synov_nodoublets$synov_nd_clusters <- Idents(synov_nodoublets)

DefaultAssay(synov_nodoublets) <- "RNA"
Syn_markers <- FindAllMarkers(synov_nodoublets, only.pos = T)
synov_nodoublets@misc$Syn_markers <- Syn_markers

saveRDS(synov_nodoublets, file = "synov_nodoublets.rds")


# Merge cell types back into combined seurat object -----------------------

f_clusters <- data.frame(no_doublets = fibroblasts@meta.data$fibroblast_final, 
                         cell = rownames(fibroblasts@meta.data))

l_clusters <- data.frame(no_doublets = lymphoid@meta.data$lymphoid_final, 
                         cell = rownames(lymphoid@meta.data))

m_clusters <- data.frame(no_doublets = myeloid@meta.data$myeloid_final, 
                         cell = rownames(myeloid@meta.data))

v_clusters <- data.frame(no_doublets = vascular@meta.data$vascular_final, 
                         cell = rownames(vascular@meta.data))

c_clusters <- data.frame(no_doublets = rep("Cycling cells",
                                           length(rownames(merged_clusters@meta.data)[merged_clusters@meta.data$finalcluster == "Cycling cells"])),
                         cell = rownames(merged_clusters@meta.data)[merged_clusters@meta.data$finalcluster == "Cycling cells"])


no_doublets <- rbind(f_clusters, l_clusters, m_clusters, v_clusters, c_clusters)



my_levels <- c("Fibroblasts", "Macrophages/Monocytes", "B cells",
               "Neutrophils", "Synoviocytes", "T cells",
               "Dendritic cells", "Adipocytes", "Cycling cells",
               "Muscle cells", "Endothelial cells", 
               "NK cells", "Osteoblasts", "ILCs", "Chondrocytes", "Red blood cells")

no_doublets$no_doublets <- factor(no_doublets,levels = my_levels)

if(nrow(no_doublets) == ncol(merged_clusters)){
  merged_clusters$no_doublets <- no_doublets$no_doublets[match(colnames(merged_clusters), no_doublets$cell)]
  
}

Idents(merged_clusters) <- "no_doublets"
saveRDS(merged_clusters, file  = "merged_clusters.rds")
