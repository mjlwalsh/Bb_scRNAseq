library(OmnipathR)
library(Matrix)
library(stringr)
library(data.table)
library(magrittr)
library(dplyr)
library(glue)
library(rlang)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

omni_full <- import_ligrecextra_interactions(organism = 10090)

# Lets get a list of all unique interations
unique_rec_lig_interactions <- omni_full %>%
  dplyr::select(source_genesymbol, target_genesymbol) %>%
  distinct()

# Read in the DEG info
de_files <- list.files("/projects/home/nealpsmith/projects/hu/data/week_comps/degs", full.names = TRUE)
de_files <- de_files[grep("padj_0.05_logfc_1.csv", de_files)]

all_de <- data.frame()
for (f in de_files){
  df <- read.csv(f)
  all_de <- rbind(all_de, df)
}

de_counts <- all_de %>%
  mutate("direction" = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(cluster, direction) %>%
  summarise(n_genes = n()) %>%
  ungroup() %>%
  complete(cluster, direction) %>%
  replace(is.na(.), 0) %>%
  group_by(cluster) %>%
  mutate(n_total = sum(n_genes)) %>%
  arrange(desc(n_total)) %>%
  mutate(n_genes = ifelse(direction == "down", -n_genes, n_genes))
de_counts$cluster <- factor(de_counts$cluster, levels = rev(unique(de_counts$cluster)))
ggplot(de_counts, aes(x = cluster, y = n_genes, group = direction, fill = direction)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("#FF8000", "#40007F")) +
    scale_y_continuous(labels = abs) +
    ggtitle(glue("Number of DE genes")) +
    ylab("# of DE genes") + xlab("") +
    # geom_text(data =de[plot_df$value != 0,],
    #           aes(label = abs(value)),
    #           position = position_dodge(width = 0.9),
    #           vjust = ifelse(plot_df[plot_df$value != 0,]$value >= 0, -1, 1),
    #           hjust = ifelse(plot_df[plot_df$value != 0,]$value >= 0, -0.15, 1.15)) +
    theme_classic(base_size = 20)
ggsave("/projects/home/nealpsmith/projects/hu/figures/cell_cell_degs/n_total_degs_padj_0.05_week_0_vs_2_logfc_1.pdf",
       width = 8, height = 8,
       useDingbats = FALSE)


de_counts <- lapply(seq_along(de_paths), function(x){
    print(de_paths[x])
    print(str_extract(de_paths[x], "cluster_(\\d)+"))
    print("_________________________________")
    data <- read.csv(de_paths[x])
    up_aa <- nrow(data[data$padj < 0.1 & data$log2FoldChange > 0,])
    up_ana <- nrow(data[data$padj < 0.1 & data$log2FoldChange < 0,])
    info_df <- data.frame("AA" = up_aa,
                          "ANA" = up_ana, row.names = str_extract(de_paths[x], "cluster_(\\d)+"))
  }) %>%
    do.call(rbind, .)
# Filter to genes of interest
all_rec_lig_genes <- c(unique_rec_lig_interactions$target_genesymbol,
                       unique_rec_lig_interactions$source_genesymbol) %>%
  unique()


de_filt <- all_de %>%
  dplyr::filter(gene %in% all_rec_lig_genes, log2FoldChange > 0) # For now, just get the ones higher at week 2

# Now lets get pairs
sig_interactions <- data.frame()
for (i in 1:nrow(unique_rec_lig_interactions)){
  intr <- unique_rec_lig_interactions[i,]
  genes <- c(intr$source_genesymbol, intr$target_genesymbol)
  d <- de_filt %>%
    dplyr::filter(gene %in% genes)
  if (all(genes %in% d$gene)){
    source_clusts <- d %>%
      dplyr::filter(gene == intr$source_genesymbol) %>%
      .$cluster
    target_clusts <- d %>%
      dplyr::filter(gene == intr$target_genesymbol) %>%
      .$cluster
    all_combs <- expand.grid(source_clusts, target_clusts) %>%
      `colnames<-`(c("source_cluster", "target_cluster")) %>%
      mutate(source = intr$source_genesymbol, target = intr$target_genesymbol)
    sig_interactions <- rbind(sig_interactions, all_combs)
  }
}

# Now add back the stats (prob could have done this better)

source_stats <- de_filt %>%
  dplyr::select(gene, cluster, log2FoldChange, pvalue, padj) %>%
  `colnames<-`(c("source", "source_cluster", "source_log2FoldChange", "source_pvalue", "source_padj"))
head(sig_interactions)
sig_interactions %<>%
  dplyr::left_join(source_stats, by = c("source", "source_cluster"))

target_stats <- de_filt %>%
  dplyr::select(gene, cluster, log2FoldChange, pvalue, padj) %>%
  `colnames<-`(c("target", "target_cluster", "target_log2FoldChange", "target_pvalue", "target_padj"))
sig_interactions %<>%
  dplyr::left_join(target_stats, by = c("target", "target_cluster"))
write.csv(sig_interactions, "/projects/home/nealpsmith/projects/hu/data/rec_lig_interactions/sig_rec_lig_interactions_wk_0_vs_2_pseudobulk.csv",
          row.names = FALSE)



# Make a heatmap to see which cluster pairs have the most predicted interactions
unique_clusts <- sort(unique(c(sig_interactions$source_cluster, sig_interactions$target_cluster)))

heatmap_df <- data.frame(matrix(ncol = length(unique_clusts), nrow = length(unique_clusts)),
                         row.names = unique_clusts) %>%
  `colnames<-`(unique_clusts)

for (cl1 in unique_clusts){
  for (cl2 in unique_clusts){
    # Get all unique rows in either orientation
    s1 <- sig_interactions %>%
      dplyr::filter(target_cluster == cl1, source_cluster == cl2) %>%
      nrow()
    if (cl1 != cl2){
     s2 <- sig_interactions %>%
      dplyr::filter(target_cluster == cl2, source_cluster == cl1) %>%
       nrow()
    } else {
      s2 <- 0
    }
    n_inter <- s1 + s2
    heatmap_df[cl1, cl2] <- n_inter
    heatmap_df[cl2, cl1] <- n_inter
  }
}
# heatmap_df[upper.tri(heatmap_df, diag = FALSE)] = NA
sum(heatmap_df[!is.na(heatmap_df)]) == nrow(sig_interactions)

head(heatmap_df)
heatmap_col_fun = colorRamp2(c(min(heatmap_df[!is.na(heatmap_df)]), max(heatmap_df[!is.na(heatmap_df)])),
                             c("#fafafa", '#02024a'))
hmap <- Heatmap(heatmap_df, name = "# interactions", col = heatmap_col_fun, na_col = "black",
                    cluster_rows = FALSE, cluster_columns = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill) {
                   if(i <= j) {
                     grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white"))
                   }
                })

pdf("/projects/home/nealpsmith/projects/hu/figures/cell_cell_degs/n_interaction_heatmap.pdf")
draw(hmap)
dev.off()

