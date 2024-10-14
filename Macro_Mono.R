setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/Macro_Mono")
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(clustree)
library(dplyr)
library(patchwork)
library(ggforce)

merged_filtered_2 <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/第一次注释_2.Data2_Annotated.rds")
Macro_Mono <- subset(merged_filtered_2, subset = Celltype %in% c("Macrophage", "Monocyte"))
Macro_Mono$Sample <- as.character(Macro_Mono$Sample)
Macro_Mono %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            # model.use = c("linear")) %>%
            model.use = c("negbinom")) %>%
  RunPCA() -> Macro_Mono

Macro_Mono <- RunUMAP(Macro_Mono,
              reduction = "pca",
              dims = 1:30, n.neighbors = 50,
              min.dist = 0.6, n.epochs = 1000)
DimPlot(Macro_Mono, group.by = "Celltype")

Macro_Mono %>%
  FindNeighbors(dims = 1:30, reduction = "pca",
                k.param = 50) -> Macro_Mono
Macro_Mono %>%
  FindClusters(resolution = seq(from = 0.1, by = 0.1, to = 2.5),
               group.singletons = TRUE) -> Macro_Mono

DimPlot(Macro_Mono)

pdf("./1.Macro_Mono_UMAP_Sample.pdf", height = 4.5, width = 5.5)
DimPlot(Macro_Mono, group.by = "Sample")
dev.off()

pdf("./2.Macro_Mono_UMAP_DataSets.pdf", height = 4.5, width = 5.5)
DimPlot(Macro_Mono, group.by = "Batch")
dev.off()

pdf("./3.Macro_Mono_UMAP_Clusters.pdf")
for (i in seq(from = 0.1, by = 0.1, to = 2.5)) {
  p <- DimPlot(Macro_Mono, label = T, group.by = paste0("RNA_snn_res.",i))
  print(p)
}
dev.off()

DimPlot(Macro_Mono, group.by = "Celltype")

clustree(Macro_Mono@meta.data, prefix = "RNA_snn_res.")


Macro_Mono %>%
  FindClusters(resolution = 0.9,
               group.singletons = TRUE) -> Macro_Mono

pdf("./4.Macro_Mono_UMAP_Clusters_resolution_0.9.pdf",
    height = 4.5, width = 5.4)
DimPlot(Macro_Mono, label = T)
dev.off()


pdf("./4.Macro_Mono_UMAP_Clusters_resolution_0.9_splited.pdf",
    height = 7, width = 14)
DimPlot(Macro_Mono, label = T,
        split.by = "seurat_clusters", ncol = 4)
dev.off()

exp <- as.matrix(Macro_Mono@assays$RNA@data)
exp <- t(exp)
cluster_mean <- aggregate.data.frame(exp,
                                     by = list(as.character(Macro_Mono@meta.data[rownames(exp), "seurat_clusters"])),
                                     FUN = mean)
cluster_mean$Group.1 <- paste0("Seurat_cluster_", cluster_mean$Group.1)
rownames(cluster_mean) <- cluster_mean$Group.1
cluster_mean <- cluster_mean[,-1]
cluster_mean <- t(cluster_mean)
corr <- cor(cluster_mean, method = "pearson")
pdf("./5.Macro_Mono_seurat_cluster_cor_heatmap_all_genes.pdf",
    width = nrow(corr) * 1 + 1, height = nrow(corr) * 1)
p <- pheatmap::pheatmap(corr, angle_col = 90, treeheight_row = 45,
                        fontsize = 14, treeheight_col = 45, #cutree_rows = 11, cutree_cols = 11,
                        cellwidth = 35, cellheight = 35, display_numbers = TRUE)
dev.off()

Idents(Macro_Mono) <- Macro_Mono$RNA_snn_res.0.9
Cluster_diff <- FindAllMarkers(Macro_Mono, only.pos = T,
                               logfc.threshold = 0.25)
Cluster_diff <- Cluster_diff[order(Cluster_diff$cluster,
                                   -Cluster_diff$avg_log2FC),]
Cluster_diff_list <- split.data.frame(Cluster_diff, f = list(Cluster_diff$cluster))
openxlsx::write.xlsx(Cluster_diff_list,
                     "./6.Cluster_DEG_resolution_0.9_log2FC_0.25.xlsx")

library(clusterProfiler)
library(org.Mm.eg.db)
DEG_Enrich_GO <- function(temp_reclustered_markers, logfc_thr = 0.25, Org = org.Mm.eg.db, prefix = " - Mus musculus (house mouse)") {
  genes_list <- split.data.frame(temp_reclustered_markers, f = list(temp_reclustered_markers$cluster))
  for (i in names(genes_list)) {
    temp <- genes_list[[i]]
    temp <- temp[which(abs(temp$avg_log2FC) > logfc_thr),]
    temp <- temp$gene
    temp <- bitr(temp, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Org)
    genes_list[[i]] <- temp$ENTREZID
  }
  GO_result_clusters <- c()
  for (i in names(genes_list)) {
    print(i)
    genes <- genes_list[[i]]
    temp_go <- enrichGO(genes, keyType = "ENTREZID",
                        ont = "BP", pvalueCutoff = 1,
                        qvalueCutoff = 1, OrgDb = Org,
                        readable = T)
    temp_go_result <- temp_go@result
    temp_go_result$Description <- Hmisc::capitalize(temp_go_result$Description)
    temp_go_result <- temp_go_result[temp_go_result$pvalue < 0.05,]
    temp_go_result$Celltype <- i
    temp_go_result$Description <- unlist(lapply(temp_go_result$Description,
                                                function(x){
                                                  unlist(strsplit(x, prefix, fixed = T))[1]
                                                }))
    temp_go@result <- temp_go_result
    GO_result_clusters <- c(GO_result_clusters,
                            list(temp_go))
  }
  names(GO_result_clusters) <- names(genes_list)
  return(GO_result_clusters)
}

Cluster_GO <- DEG_Enrich_GO(temp_reclustered_markers = Cluster_diff,
                            logfc_thr = 0.25, Org = org.Mm.eg.db, 
                            prefix = " - Mus musculus (house mouse)")
openxlsx::write.xlsx(Cluster_GO, "./7.Macro_Mono_Cluster_GO_resolution_0.9.xlsx")

saveRDS(Macro_Mono, "Macro_Mono_subclusters.rds")

Macro_Mono <- readRDS("Macro_Mono_subclusters.rds")

DimPlot(Macro_Mono, label = T)
DimPlot(Macro_Mono, group.by = "TCR")
DimPlot(Macro_Mono, group.by = "Cluster_Celltype")

Macro_Mono <- RenameIdents(Macro_Mono,
                      "0" = "Mac-C1qb",
                      "1" = "Mono-Chil3",
                      "2" = "Mac-Mif",
                      "3" = "Mac-prolif",
                      "4" = "Mac-Plin2",
                      "5" = "Mac-mt",
                      "6" = "Unidentified",
                      "7" = "Mono-Ace")
Macro_Mono$SubCelltype <- as.character(Idents(Macro_Mono))
sort(unique(Macro_Mono$SubCelltype))

Idents(Macro_Mono) <- factor(Macro_Mono$SubCelltype,
                     levels = sort(unique(Macro_Mono$SubCelltype)))
colors <- ArchR::ArchRPalettes$circus
names(colors) <- NULL
colors <- c(colors[1:7], "gray")
pdf("8.Macro_Mono_Subtypes.pdf", width = 5.2, height = 4)
DimPlot(Macro_Mono, label = T, cols = colors, pt.size = 0.7)
dev.off()

saveRDS(Macro_Mono, "Macro_Mono_annotated.rds")

features <- list("Mac-C1qb" = c("C1qb"),
                 "Mac-Mif" = c("Mif"),
                 "Mac-mt" = c("mt-Co1", "mt-Co2", "mt-Atp8", "mt-Co3"),
                 "Mac-Plin2" = c("Cd14", "Plin2"),
                 "Mac-prolif" = c("Stmn1", "Top2a", "Pclaf"),
                 "Mono-Ace" = c("Ace", "Csf1r"),
                 "Mono-Chil3" = c("Lyz2", "Ccl2", "Chil3", "Mgst1"))

pdf("9.Macro_Mono_Subtypes_DotPlot.pdf", height = 3.5, width = 10)
DotPlot(Macro_Mono, features = features, cols = c("gray", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))
dev.off()




