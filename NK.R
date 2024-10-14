setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/NK")
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(clustree)
library(dplyr)
library(patchwork)
library(ggforce)

merged_filtered_2 <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/第一次注释_2.Data2_Annotated.rds")
NK <- subset(merged_filtered_2, subset = Celltype == "NK")
NK$Sample <- as.character(NK$Sample)
NK %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> NK

NK <- RunUMAP(NK,
                 reduction = "pca",
                 dims = 1:30, n.neighbors = 50,
                 min.dist = 0.25, n.epochs = 1000)
NK %>%
  FindNeighbors(dims = 1:30, reduction = "pca",
                k.param = 50) -> NK
NK %>%
  FindClusters(resolution = seq(from = 0.1, by = 0.1, to = 2.5),
               group.singletons = TRUE) -> NK

DimPlot(NK)

pdf("./1.NK_UMAP_Sample.pdf", height = 4.5, width = 5.5)
DimPlot(NK, group.by = "Sample")
dev.off()

pdf("./2.NK_UMAP_DataSets.pdf", height = 4.5, width = 5.5)
DimPlot(NK, group.by = "Batch")
dev.off()

pdf("./3.NK_UMAP_Clusters.pdf")
for (i in seq(from = 0.1, by = 0.1, to = 2.5)) {
  p <- DimPlot(NK, label = T, group.by = paste0("RNA_snn_res.",i))
  print(p)
}
dev.off()

NK %>%
  FindClusters(resolution = 0.7,
               group.singletons = TRUE) -> NK

pdf("./4.NK_UMAP_Clusters_resolution_0.7.pdf",
    height = 4.5, width = 5.4)
DimPlot(NK, label = T)
dev.off()


pdf("./4.NK_UMAP_Clusters_resolution_0.7_splited.pdf",
    height = 7, width = 10)
DimPlot(NK, label = T,
        split.by = "seurat_clusters", ncol = 4)
dev.off()

exp <- as.matrix(NK@assays$RNA@data)
exp <- t(exp)
cluster_mean <- aggregate.data.frame(exp,
                                     by = list(as.character(NK@meta.data[rownames(exp), "seurat_clusters"])),
                                     FUN = mean)
cluster_mean$Group.1 <- paste0("Seurat_cluster_", cluster_mean$Group.1)
rownames(cluster_mean) <- cluster_mean$Group.1
cluster_mean <- cluster_mean[,-1]
cluster_mean <- t(cluster_mean)
corr <- cor(cluster_mean, method = "pearson")
pdf("./5.NK_seurat_cluster_cor_heatmap_all_genes.pdf",
    width = nrow(corr) * 1 + 1, height = nrow(corr) * 1)
p <- pheatmap::pheatmap(corr, angle_col = 90, treeheight_row = 45,
                        fontsize = 14, treeheight_col = 45, #cutree_rows = 11, cutree_cols = 11,
                        cellwidth = 35, cellheight = 35, display_numbers = TRUE)
dev.off()

Idents(NK) <- NK$RNA_snn_res.0.7
Cluster_diff <- FindAllMarkers(NK, only.pos = T,
                               logfc.threshold = 0.25)
Cluster_diff <- Cluster_diff[order(Cluster_diff$cluster,
                                   -Cluster_diff$avg_log2FC),]
Cluster_diff_list <- split.data.frame(Cluster_diff, f = list(Cluster_diff$cluster))
openxlsx::write.xlsx(Cluster_diff_list,
                     "./6.Cluster_DEG_resolution_0.7_log2FC_0.25.xlsx")

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
openxlsx::write.xlsx(Cluster_GO, "./7.NK_Cluster_GO_resolution_0.7.xlsx")

saveRDS(NK, "NK_subclusters.rds")

NK <- readRDS("NK_subclusters.rds")

DimPlot(NK, label = T)
DimPlot(NK, group.by = "TCR")
DimPlot(NK, group.by = "Cluster_Celltype")

NK <- RenameIdents(NK,
                   "0" = "NK-1",
                   "1" = "NK-4",
                   "2" = "NK-activation",
                   "3" = "NK-3",
                   "4" = "NK-3",
                   "5" = "NK-6",
                   "6" = "NK-2",
                   "7" = "NKT",
                   "8" = "NK-antigen presentation")

NK$SubCelltype <- as.character(Idents(NK))
sort(unique(NK$SubCelltype))

Idents(NK) <- factor(NK$SubCelltype,
                     levels = sort(unique(NK$SubCelltype)))
colors <- ArchR::ArchRPalettes$kelly
names(colors) <- NULL

pdf("8.NK_Subtypes.pdf", width = 6.2, height = 4)
DimPlot(NK, label = T, cols = colors, pt.size = 0.7,
        repel = T)
dev.off()

saveRDS(NK, "NK_annotated.rds")

features <- list("NK-1" = c("Plac8", "Ctla2a", "Fosl2", "Vegfa", "Sult2b1", "Tigit", "Ier5l"),
                 "NK-2" = c("Ltb", "Ly6e", "Emb", "Ly6c2"),
                 "NK-3" = c("Hsp90ab1", "Rps17", "Rps20", "Rpl12", "Rps2"),
                 "NK-4" = c("Ccl4", "Ccl3", "Nr4a1", "Icam1", "Ifng", "Nfkbid", "Nfkbiz", "Nfkbia"),
                 "NK-5" = c("Ccl5", "S100a6", "Lgals1", "Klrg1", "Gzma"),
                 "NK-antigen presentation" = c("H2-Aa", "H2-Eb1", "H2-Ab1"),
                 "NKT" = c("Il7r", "Ifi27l2a", "Ly6a"))

pdf("9.NK_Subtypes_DotPlot.pdf", height = 3.5, width = 14)
DotPlot(NK, features = features, cols = c("gray", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))
dev.off()


