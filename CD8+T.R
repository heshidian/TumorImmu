setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/CD8+")
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(clustree)
library(dplyr)
library(patchwork)
library(ggforce)
load("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Data.RData")

merged_filtered_2 <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/第一次注释_2.Data2_Annotated.rds")
CD8_T <- subset(merged_filtered_2, subset = Celltype == "CD8+ T")
CD8_T$Sample <- as.character(CD8_T$Sample)
CD8_T %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> CD8_T
CD8_T <- RunUMAP(CD8_T,
                 reduction = "pca",
                 dims = 1:30, n.neighbors = 50,
                 min.dist = 0.25, n.epochs = 1000)
CD8_T %>%
  FindNeighbors(dims = 1:30, reduction = "pca",
                k.param = 50) -> CD8_T

CD8_T %>%
  FindClusters(resolution = seq(from = 0.1, by = 0.1, to = 2.5),
               group.singletons = TRUE) -> CD8_T
DimPlot(CD8_T)

clustree(CD8_T@meta.data, prefix = "RNA_snn_res.")

CD8_T %>%
  FindClusters(resolution = 1,
               group.singletons = TRUE) -> CD8_T

pdf("./1.CD8_T_UMAP_Sample.pdf", height = 4.5, width = 5.5)
DimPlot(CD8_T, group.by = "Sample", cols = Sample_color)
dev.off()

pdf("./2.CD8_T_UMAP_DataSets.pdf", height = 4.5, width = 5.5)
DimPlot(CD8_T, group.by = "Batch")
dev.off()

pdf("./3.CD8_T_UMAP_Clusters.pdf")
for (i in seq(from = 0.1, by = 0.1, to = 2.5)) {
  p <- DimPlot(CD8_T, label = T, group.by = paste0("RNA_snn_res.",i))
  print(p)
}
dev.off()

CD8_T %>%
  FindClusters(resolution = 0.4,
               group.singletons = TRUE) -> CD8_T

pdf("./4.CD8_T_UMAP_Clusters_resolution_0.4.pdf",
    height = 4.5, width = 5.4)
DimPlot(CD8_T, label = T)
dev.off()

pdf("./4.CD8_T_UMAP_Clusters_resolution_0.4_splited.pdf",
    height = 7, width = 3.5 * 4)
DimPlot(CD8_T, label = T,
        split.by = "seurat_clusters", ncol = 4)
dev.off()

exp <- as.matrix(CD8_T@assays$RNA@data)
exp <- t(exp)
cluster_mean <- aggregate.data.frame(exp,
                                     by = list(as.character(CD8_T@meta.data[rownames(exp), "seurat_clusters"])),
                                     FUN = mean)
cluster_mean$Group.1 <- paste0("Seurat_cluster_", cluster_mean$Group.1)
rownames(cluster_mean) <- cluster_mean$Group.1
cluster_mean <- cluster_mean[,-1]
cluster_mean <- t(cluster_mean)
corr <- cor(cluster_mean, method = "pearson")
pdf("./5.CD8_T_seurat_cluster_cor_heatmap_all_genes.pdf",
    width = nrow(corr) * 1 + 1, height = nrow(corr) * 1)
p <- pheatmap::pheatmap(corr, angle_col = 90, treeheight_row = 45,
                        fontsize = 14, treeheight_col = 45, #cutree_rows = 11, cutree_cols = 11,
                        cellwidth = 35, cellheight = 35, display_numbers = TRUE)
dev.off()

Idents(CD8_T) <- CD8_T$RNA_snn_res.0.5
Cluster_diff <- FindAllMarkers(CD8_T, only.pos = T,
                               logfc.threshold = 0.25)
Cluster_diff <- Cluster_diff[order(Cluster_diff$cluster,
                                   -Cluster_diff$avg_log2FC),]
Cluster_diff_list <- split.data.frame(Cluster_diff, f = list(Cluster_diff$cluster))
openxlsx::write.xlsx(Cluster_diff_list,
                     "./6.Cluster_DEG_resolution_0.5_log2FC_0.25.xlsx")

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
openxlsx::write.xlsx(Cluster_GO, "./7.CD8_T_Cluster_GO_resolution_0.5.xlsx")

saveRDS(CD8_T, "CD8_T_subclusters.rds")

CD8_T <- readRDS("CD8_T_subclusters.rds")
Idents(CD8_T) <- CD8_T$RNA_snn_res.0.5

DimPlot(CD8_T, label = T)

CD8_T$TCR
DimPlot(CD8_T, group.by = "TCR")

CD8_T <- RenameIdents(CD8_T,
                      "0" = "CD8 Teff",
                      "1" = "CD8 Tex",
                      "2" = "CD8 IFN Response",
                      "3" = "CD8 Proliferating-Mki67",
                      "4" = "CD8 Tem",
                      "5" = "CD8 Proliferating-Cdc20",
                      "6" = "CD8 Naive/Tcm",
                      "7" = "CD8 Tpex",
                      "8" = "CD8 Texpanding"
                      )

CD8_T$SubCelltype <- as.character(Idents(CD8_T))
Idents(CD8_T) <- factor(CD8_T$SubCelltype,
                        levels = sort(unique(CD8_T$SubCelltype)))
colors <- ArchR::ArchRPalettes$circus
colors <- colors[1:length(unique(CD8_T$SubCelltype))]
names(colors) <- levels(Idents(CD8_T))
pdf("8.CD8_T_Subtypes.pdf",
    height = 4.5, width = 7)
DimPlot(CD8_T, label = T, repel = F,
        cols = colors)
dev.off()

features <- list("CD8 IFN Response" = c("Ifit3", "Ifit1", "Isg15", "Ifit3b",
                                        "Cd69", "Ifngr1", "Bst2", "Rtp4"),
                 "CD8 Naive/Tcm" = c("Lef1", "Sell", "Tcf7", "Ccr7", "Il7r", "Dapl1"),
                 "CD8 Proliferating-Cdc20" = c("Cdc20", "Ccnb2", "Ccnb1", "Hmgb2",
                                               "Dut", "Plk1"),
                 "CD8 Proliferating-Mki67" = c("Mki67", "Top2a", "Tubb5", "Hist1h1b",
                                               "Hist1h2ap", "Tuba1b"),
                 "CD8 Teff" = c("Xcl1", "Ifng", "Prf1", "Gzmb", "Gzmf", "Ccl3"),
                 "CD8 Tem" = c("Gzmk", "Cd8b1", "Cd7", "Ly6c2", "Gzma"),
                 "CD8 Tex" = c("Tox", "Bhlhe40", "Rbpj", "Nr4a2", "Nfil3", "Spp1",
                               "Lag3", "Tigit", "Pdcd1", "Tnfrsf9", "Havcr2", "Entpd1"),
                 "CD8 Texpanding" = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"),
                 "CD8 Tpex" = c("Stmn1", "Hist1h1e", "Birc5"))

pdf("9.CD8_T_Subtypes_DotPlot.pdf", height = 3.5, width = 19)
DotPlot(CD8_T, features = features, cols = c("gray", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))
dev.off()

saveRDS(CD8_T, "CD8_T_annotated.rds")

##################


sub_type_markers <- list("Naive-like" = c("Lef1", "Hnf1a", "Mal",
                                          "Ifit2", "Oasl1", "Oasl2", "Ifit3", "Isg15"),
                         "TEM" = c("Tcf7", "Cd27", "Cd7",
                                   "Cd52", "Ccl5", "Gzmk"),
                         "Proliferating T cell" = c("Mki67", "Top2a", "H1f5", "H3c2",
                                                    "Tyms", "Pcna", "Kifc1"),
                         "Teff" = c("Gzma", "Gzmb", "Nkg7", "Xcl1", "Prf1", "Gnly", "Fox", "Ltb"),
                         "Tpex" = c("Il7r", "Il2ra", "Myb", "Sell", "Slamf6", "Id3"),
                         "Tex" = c("Tigit", "Lag3", "Ctla4", "Cxcr6",
                                   "Krt18", "Entpd1", "Tox", "Havcr2",
                                   "Krt86", "Klrb1", "Pdcd1", "Cxcl13", "Layn"),
                         "Tact" = c("Hspa1b"),
                         "Tn" = c("Ccr7", "Sesn3", "S1pr1"))

DotPlot(CD8_T, features = sub_type_markers) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

CD8_T <- RenameIdents(CD8_T,
                      "0" = "",
                      "1" = "",
                      "2" = "",
                      "3" = "",
                      "4" = "",
                      "5" = "",
                      "6" = "",
                      "7" = "",
                      "8" = ""
                      )

Cluster8_cells <- colnames(CD8_T)[CD8_T$RNA_snn_res.0.5 == "8"]
  
merged_filtered_2$CD8_cluster <- "Others"
merged_filtered_2@meta.data[colnames(CD8_T), "CD8_cluster"] <- CD8_T$RNA_snn_res.0.5

DimPlot(merged_filtered_2, group.by = "CD8_cluster")

table(merged_filtered_2$Cluster_Celltype, merged_filtered_2$CD8_cluster)
