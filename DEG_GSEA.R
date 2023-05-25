library(tidyverse)
library(Seurat)
library(DESeq2)
library(xlsx)
library(fgsea)

setwd("~/Documents/Hu lab/R/Hu lab joint single cell/all weeks analysis")

# Function to make the pseudobulk
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype==ct]
    pseudobulk <- rowSums(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

#pathways ------------------------
pathwaysG <- gmtPathways("/Users/Jennifer/Downloads/c5.go.bp.v7.5.1.symbols.gmt")
pathwaysG <- lapply(pathwaysG, str_to_title)
pathwaysG <- lapply(pathwaysG, function(x) gsub(pattern = "Mt-", replacement = "mt-", x))

pathwaysH <- gmtPathways("/Users/Jennifer/Downloads/h.all.v7.5.1.symbols.gmt")
pathwaysH <- lapply(pathwaysH, str_to_title)
pathwaysH <- lapply(pathwaysH, function(x) gsub(pattern = "Mt-", replacement = "mt-", x))


#Macrophages Pseudobulk  ---------------
data <- readRDS("mac_nodoublets.rds")
mat <- data$RNA@counts
data@meta.data$samp_clust <- paste(stringr::str_extract(rownames(data@meta.data), "[^_]*_[^_]*"),
                                   data@meta.data$mac_nd_clusters, sep = "_") #change cluster name
data@meta.data$samp_clust <- factor(data@meta.data$samp_clust)
samp_clust_list <- data@meta.data$samp_clust
names(samp_clust_list) <- colnames(mat)

#filter to remove samples with few cells
howmanycells = table(samp_clust_list)
remove = names(howmanycells)[howmanycells == 1]
remove_indx = rownames(data@meta.data)[data@meta.data$samp_clust %in% remove]

new_mat = mat[,!colnames(mat) %in% remove_indx]
new_list = samp_clust_list[!names(samp_clust_list) %in% remove_indx]
new_list <- droplevels(new_list)

pseudo_mtx <- getPseudobulk(new_mat, new_list)

meta_data <- data.frame(row.names = colnames(pseudo_mtx))
meta_data$sample <- stringr::str_extract(rownames(meta_data), "[^_]*_[^_]*")
meta_data$cluster <- sub("^[^_]*_[^_]*_", "", rownames(meta_data))
meta_data$tmpt <- str_extract(rownames(meta_data), "[^_]*")
meta_data$cell_number <- as.numeric(table(new_list))
meta_data$time <- as.numeric(str_extract(meta_data$tmpt, "\\d")) / 2

## Mac for loop dataframe -------------------------------
# for loop to make a master fgsea dataframe to save as a csv
# DEGs (wk0 vs wk2) will be saved as an xlsx file
fgsea_mac <- data.frame()
for(cl in names(table(meta_data$cluster))){
  # Lets make a temporary dds file to make a graph
  meta_temp <- meta_data[meta_data$cluster == cl,]
  count_temp <- pseudo_mtx[,rownames(meta_temp)]
  dds <- DESeqDataSetFromMatrix(countData = count_temp,
                                colData = meta_temp,
                                design = ~ cell_number + tmpt) #adding in cell_number to formula 
  dds <- DESeq(dds,
               test = "LRT",
               reduced = ~cell_number) #reduced formula: accounting for cell_number differences
  
  res <- results(dds, contrast = c("tmpt", "wk2", "wk0"))
  
  # Create a dataframe of DEGs from DESeq2
  res_df <- res %>% as.data.frame() %>% arrange(padj) %>% 
    .[1:3000,]
  write.xlsx(res_df, file = "mac_LRT_wk2_cell_number.xlsx", sheetName = cl, append = TRUE)
  
  res <- as.data.frame(res) %>%
    dplyr::filter(!is.na(padj)) %>%
    arrange(pvalue)
  
  ranks <- res$log2FoldChange #rank the genes by log2foldchange
  ranks <- as.numeric(ranks)
  names(ranks) <- rownames(res)
  ranks <- ranks[!is.na(ranks)]
  
  set.seed(85849484) # seed for reproducibility of fgsea
  fgsea <- fgseaMultilevel(pathwaysG, ranks, minSize=5, maxSize = 500, eps = 0) 
  fgsea$cluster <- cl
  fgsea_mac <- rbind(fgsea_mac, fgsea)
  print(paste("Finished", cl))
}

fgsea_mac <- fgsea_mac %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = "|"))

write.csv(fgsea_mac, file = "fgsea/fgsea_go_macmono_week0_week2.csv",
          row.names = F)


# T cell pseudobulk GSEA -------------------------------
data <- readRDS("Tcell_nodoublets.rds")
mat <- data$RNA@counts
data@meta.data$samp_clust <- paste(stringr::str_extract(rownames(data@meta.data), "[^_]*_[^_]*"),
                                   data@meta.data$Tcell_nd_clusters, sep = "_") #change cluster name
data@meta.data$samp_clust <- factor(data@meta.data$samp_clust)
samp_clust_list <- data@meta.data$samp_clust
names(samp_clust_list) <- colnames(mat)

#filter to remove samples with few cells
howmanycells = table(samp_clust_list)
remove = names(howmanycells)[howmanycells == 1]
remove_indx = rownames(data@meta.data)[data@meta.data$samp_clust %in% remove]

new_mat = mat[,!colnames(mat) %in% remove_indx]
new_list = samp_clust_list[!names(samp_clust_list) %in% remove_indx]
new_list <- droplevels(new_list)

pseudo_mtx <- getPseudobulk(new_mat, new_list)

meta_data <- data.frame(row.names = colnames(pseudo_mtx))
meta_data$sample <- stringr::str_extract(rownames(meta_data), "[^_]*_[^_]*")
meta_data$cluster <- sub("^[^_]*_[^_]*_", "", rownames(meta_data))
meta_data$tmpt <- str_extract(rownames(meta_data), "[^_]*")
meta_data$cell_number <- as.numeric(table(new_list))
meta_data$time <- as.numeric(str_extract(meta_data$tmpt, "\\d")) / 2

## T cell for loop dataframe -------------------------------
# for loop to make a master fgsea dataframe to save as a csv
# DEGs (wk0 vs wk2) will be saved as an xlsx file
fgsea_T <- data.frame()
# Only Tcell.2 and Tcell.4 have enough cells to perform pseudobulk DEG analysis
for(cl in c("Tcell.2", "Tcell.4")){
  # Lets make a temporary dds file to make a graph
  meta_temp <- meta_data[meta_data$cluster == cl,]
  count_temp <- pseudo_mtx[,rownames(meta_temp)]
  dds <- DESeqDataSetFromMatrix(countData = count_temp,
                                colData = meta_temp,
                                design = ~ cell_number + tmpt) #adding in cell_number to formula 
  dds <- DESeq(dds,
               test = "LRT",
               reduced = ~cell_number) #reduced formula: accounting for cell_number differences
  
  res <- results(dds, contrast = c("tmpt", "wk2", "wk0"))
  
  # Create a dataframe of DEGs from DESeq2
  res_df <- res %>% as.data.frame() %>% arrange(padj) %>% 
    .[1:3000,]
  write.xlsx(res_df, file = "Tcell_LRT_wk2_cell_number.xlsx", sheetName = cl, append = TRUE)
  
  #remove unwanted genes
  res <- as.data.frame(res) %>%
    .[-grep("^mt-|Gm[0-9]|Rp[ls]|Tr[abgd][vjc]|Tcrg|-ps", rownames(res)),]
  
  
  ranks <- res$log2FoldChange #rank the genes by log2foldchange
  ranks <- as.numeric(ranks)
  names(ranks) <- rownames(res)
  ranks <- ranks[!is.na(ranks)]
  
  set.seed(85849484) # set seed for reproducibility
  fgsea <- fgseaMultilevel(pathwaysH, ranks, minSize=5, maxSize = 500, eps = 0) 
  fgsea$cluster <- cl
  fgsea_T <- rbind(fgsea_T, fgsea)
  print(paste("Finished", cl))
}

fgsea_T <- fgsea_T %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = "|"))

write.csv(fgsea_T, file = "fgsea/fgsea_Hallmark_Tcell_week0_week2.csv",
          row.names = F)


# Fibroblast pseudobulk GSEA -------------------------------
data <- readRDS("fib_nodoublets.rds")
mat <- data$RNA@counts
data@meta.data$samp_clust <- paste(stringr::str_extract(rownames(data@meta.data), "[^_]*_[^_]*"),
                                   data@meta.data$fib_nd_clusters, sep = "_") #change cluster name
data@meta.data$samp_clust <- factor(data@meta.data$samp_clust)
samp_clust_list <- data@meta.data$samp_clust
names(samp_clust_list) <- colnames(mat)

#filter to remove samples with few cells
howmanycells = table(samp_clust_list)
remove = names(howmanycells)[howmanycells == 1]
remove_indx = rownames(data@meta.data)[data@meta.data$samp_clust %in% remove]

new_mat = mat[,!colnames(mat) %in% remove_indx]
new_list = samp_clust_list[!names(samp_clust_list) %in% remove_indx]
new_list <- droplevels(new_list)

pseudo_mtx <- getPseudobulk(new_mat, new_list)

meta_data <- data.frame(row.names = colnames(pseudo_mtx))
meta_data$sample <- stringr::str_extract(rownames(meta_data), "[^_]*_[^_]*")
meta_data$cluster <- sub("^[^_]*_[^_]*_", "", rownames(meta_data))
meta_data$tmpt <- str_extract(rownames(meta_data), "[^_]*")
meta_data$cell_number <- as.numeric(table(new_list))
meta_data$time <- as.numeric(str_extract(meta_data$tmpt, "\\d")) / 2


## fibroblast for loop dataframe -------------------------------
# for loop to make a master fgsea dataframe to save as a csv
# DEGs (wk0 vs wk2) will be saved as an xlsx file
fgsea_fib <- data.frame()
for(cl in names(table(meta_data$cluster))){
  # Lets make a temporary dds file to make a graph
  meta_temp <- meta_data[meta_data$cluster == cl,]
  count_temp <- pseudo_mtx[,rownames(meta_temp)]
  dds <- DESeqDataSetFromMatrix(countData = count_temp,
                                colData = meta_temp,
                                design = ~ cell_number + tmpt) #adding in cell_number to formula 
  dds <- DESeq(dds,
               test = "LRT",
               reduced = ~cell_number) #reduced formula: accounting for cell_number differences
  
  res <- results(dds, contrast = c("tmpt", "wk2", "wk0"))
  
  # Create a dataframe of DEGs from DESeq2
  res_df <- res %>% as.data.frame() %>% arrange(padj) %>% 
    .[1:3000,]
  write.xlsx(res_df, file = "Fib_LRT_wk2_cell_number.xlsx", sheetName = cl, append = TRUE)
  
  #remove unwanted genes
  res <- as.data.frame(res) %>%
    .[-grep("^mt-|Gm[0-9]|Rp[ls]|Tr[abgd][vjc]|Tcrg|-ps", rownames(res)),]
  
  ranks <- res$log2FoldChange #rank the genes by log2foldchange
  ranks <- as.numeric(ranks)
  names(ranks) <- rownames(res)
  ranks <- ranks[!is.na(ranks)]
  
  set.seed(85849484) # set seed for reproducibility
  fgsea <- fgseaMultilevel(pathwaysH, ranks, minSize=5, maxSize = 500, eps = 0) 
  fgsea$cluster <- cl
  fgsea_fib <- rbind(fgsea_fib, fgsea)
  print(paste("Finished", cl))
}

fgsea_fib <- fgsea_fib %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = "|"))

write.csv(fgsea_fib, file = "fgsea/fgsea_Hallmark_fib_week0_week2.csv",
          row.names = F)

#Pseudobulk synoviocytes ----------------------------
data <- readRDS("synov_nodoublets.rds")
mat <- data$RNA@counts
data@meta.data$samp_clust <- paste(stringr::str_extract(rownames(data@meta.data), "[^_]*_[^_]*"),
                                   data@meta.data$synov_nd_clusters, sep = "_") #change cluster name
data@meta.data$samp_clust <- factor(data@meta.data$samp_clust)
samp_clust_list <- data@meta.data$samp_clust
names(samp_clust_list) <- colnames(mat)

#filter to remove samples with few cells
howmanycells = table(samp_clust_list)
remove = names(howmanycells)[howmanycells == 1]
remove_indx = rownames(data@meta.data)[data@meta.data$samp_clust %in% remove]

new_mat = mat[,!colnames(mat) %in% remove_indx]
new_list = samp_clust_list[!names(samp_clust_list) %in% remove_indx]
new_list <- droplevels(new_list)

pseudo_mtx <- getPseudobulk(new_mat, new_list)

meta_data <- data.frame(row.names = colnames(pseudo_mtx))
meta_data$sample <- stringr::str_extract(rownames(meta_data), "[^_]*_[^_]*")
meta_data$cluster <- sub("^[^_]*_[^_]*_", "", rownames(meta_data))
meta_data$tmpt <- str_extract(rownames(meta_data), "[^_]*")
meta_data$cell_number <- as.numeric(table(new_list))
meta_data$time <- as.numeric(str_extract(meta_data$tmpt, "\\d")) / 2

## Synoviocyte for loop dataframe -------------------------------
# for loop to make a master fgsea dataframe to save as a csv
# DEGs (wk0 vs wk2) will be saved as an xlsx file
fgsea_syn <- data.frame()
for(cl in names(table(meta_data$cluster))){
  # Lets make a temporary dds file to make a graph
  meta_temp <- meta_data[meta_data$cluster == cl,]
  count_temp <- pseudo_mtx[,rownames(meta_temp)]
  dds <- DESeqDataSetFromMatrix(countData = count_temp,
                                colData = meta_temp,
                                design = ~ cell_number + tmpt) #adding in cell_number to formula 
  dds <- DESeq(dds,
               test = "LRT",
               reduced = ~cell_number) #reduced formula: accounting for cell_number differences
  
  res <- results(dds, contrast = c("tmpt", "wk2", "wk0"))
  
  # Create a dataframe of DEGs from DESeq2
  res_df <- res %>% as.data.frame() %>% arrange(padj) %>% 
    .[1:3000,]
  write.xlsx(res_df, file = "Syn_LRT_wk2_cell_number.xlsx", sheetName = cl, append = TRUE)
  
  #remove unwanted genes
  res <- as.data.frame(res) %>%
    .[-grep("^mt-|Gm[0-9]|Rp[ls]|Tr[abgd][vjc]|Tcrg|-ps", rownames(res)),]
  
  ranks <- res$log2FoldChange #rank the genes by log2foldchange
  ranks <- as.numeric(ranks)
  names(ranks) <- rownames(res)
  ranks <- ranks[!is.na(ranks)]
  
  set.seed(85849484) # set seed for reproducibility
  fgsea <- fgseaMultilevel(pathwaysH, ranks, minSize=5, maxSize = 500, eps = 0) 
  fgsea$cluster <- cl
  fgsea_syn <- rbind(fgsea_syn, fgsea)
  print(paste("Finished", cl))
}

fgsea_syn <- fgsea_syn %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = "|"))

write.csv(fgsea_syn, file = "fgsea/fgsea_Hallmark_Syn_week0_week2.csv",
          row.names = F)
