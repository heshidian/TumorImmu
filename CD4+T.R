setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/CD4+")
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
CD4_T <- subset(merged_filtered_2, subset = Celltype == "CD4+ T")
CD4_T$Sample <- as.character(CD4_T$Sample)
CD4_T %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> CD4_T
CD4_T <- RunUMAP(CD4_T,
                 reduction = "pca",
                 dims = 1:30, n.neighbors = 50,
                 min.dist = 0.25, n.epochs = 1000)
CD4_T %>%
  FindNeighbors(dims = 1:30, reduction = "pca",
                k.param = 50) -> CD4_T

CD4_T %>%
  FindClusters(resolution = seq(from = 0.5, by = 0.1, to = 2.5),
               group.singletons = TRUE) -> CD4_T
DimPlot(CD4_T)

clustree(CD4_T@meta.data, prefix = "RNA_snn_res.")

CD4_T %>%
  FindClusters(resolution = 1,
               group.singletons = TRUE) -> CD4_T

pdf("./1.CD4_T_UMAP_Sample.pdf", height = 4.5, width = 5.5)
DimPlot(CD4_T, group.by = "Sample", cols = Sample_color)
dev.off()

pdf("./2.CD4_T_UMAP_DataSets.pdf", height = 4.5, width = 5.5)
DimPlot(CD4_T, group.by = "Batch")
dev.off()

pdf("./3.CD4_T_UMAP_Clusters.pdf")
for (i in seq(from = 0.5, by = 0.1, to = 2.5)) {
  p <- DimPlot(CD4_T, label = T, group.by = paste0("RNA_snn_res.",i))
  print(p)
}
dev.off()

CD4_T %>%
  FindClusters(resolution = 0.5,
               group.singletons = TRUE) -> CD4_T

pdf("./4.CD4_T_UMAP_Clusters_resolution_0.5.pdf",
    height = 4.5, width = 5.4)
DimPlot(CD4_T, label = T)
dev.off()

pdf("./4.CD4_T_UMAP_Clusters_resolution_0.5_splited.pdf",
    height = 7, width = 3.5 * 4)
DimPlot(CD4_T, label = T,
        split.by = "seurat_clusters", ncol = 4)
dev.off()

exp <- as.matrix(CD4_T@assays$RNA@data)
exp <- t(exp)
cluster_mean <- aggregate.data.frame(exp,
                                     by = list(as.character(CD4_T@meta.data[rownames(exp), "seurat_clusters"])),
                                     FUN = mean)
cluster_mean$Group.1 <- paste0("Seurat_cluster_", cluster_mean$Group.1)
rownames(cluster_mean) <- cluster_mean$Group.1
cluster_mean <- cluster_mean[,-1]
cluster_mean <- t(cluster_mean)
corr <- cor(cluster_mean, method = "pearson")
pdf("./5.CD4_T_seurat_cluster_cor_heatmap_all_genes.pdf",
    width = nrow(corr) * 0.6 + 1, height = nrow(corr) * 0.6)
p <- pheatmap::pheatmap(corr, angle_col = 90, treeheight_row = 45,
                        fontsize = 14, treeheight_col = 45, #cutree_rows = 11, cutree_cols = 11,
                        cellwidth = 35, cellheight = 35, display_numbers = TRUE)
dev.off()

Cluster_diff <- FindAllMarkers(CD4_T, only.pos = T,
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
openxlsx::write.xlsx(Cluster_GO, "./7.CD4_T_Cluster_GO_resolution_0.5.xlsx")

saveRDS(CD4_T, "CD4_T_subclusters.rds")

CD4_T <- readRDS("CD4_T_subclusters.rds")
Idents(CD4_T) <- CD4_T$RNA_snn_res.0.5

DotPlot(CD4_T, features = c("Isg15", "Ifit1", # CD4 IFN response
                            "Mki67", "Top2a", "Pcna", # CD4 proliferating
                            "Foxp3", "Ctla4", "Ccr8", "Tnfrsf9", # CD4 Treg
                            "Ccr7", "Lef1", "Foxp1", "Tcf7", "Il7r", # CD4 Naive
                            "Cxcr5", "Cd27", # CD4 TEM
                            "Pdcd1", "Prdm1", "Havcr2", # CD4 Exh
                            "Ifitm3", # CD4 effector/Th
                            "Ifng", "Ifngr1", "Cxcr3", "Tbx21", "Stat4", "Stat1", # Th1
                            "Gata3", "Ccr3", "Ccr4", # Th2
                            "Spi1", "Il9", "Il10", # Th9
                            "Il23a", "Ccr6", "Il1r1", "Il1r2", "Klrb1c", "Klrb1a", "Klrb1b", "Stat3", # Th17
                            "Ccr10", # Th22
                            "Il21", "Il4" # Tfh
                            )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

FeaturePlot(CD4_T, features = c("Foxp3", "Ctla4", "Ccr8", "Tnfrsf9", "Prdm1","Pdcd1"), order = T)


CD4_T <- readRDS("CD4_T_subclusters.rds")
DimPlot(CD4_T)
CD4_T <- RenameIdents(CD4_T,
                      "0" = "CD4 Naive/Memory",
                      "1" = "CD4 Tpex",
                      "2" = "CD4 Treg",
                      "3" = "CD4 Texpanding",
                      "4" = "CD4 Tex",
                      "5" = "CD4 IFN Response",
                      "6" = "CD4 Tem",
                      "7" = "CD4 Teff"
                      )
CD4_T$SubCelltype <- as.character(Idents(CD4_T))
Idents(CD4_T) <- factor(CD4_T$SubCelltype,
                        levels = sort(unique(CD4_T$SubCelltype)))
colors <- ArchR::ArchRPalettes$bear
colors <- colors[1:length(unique(CD4_T$SubCelltype))]
names(colors) <- levels(Idents(CD4_T))
pdf("8.CD4_T_Subtypes.pdf",
    height = 4.5, width = 7)
DimPlot(CD4_T, label = T, repel = T,
        cols = colors)
dev.off()

saveRDS(CD4_T, "CD4_T_annotated.rds")

features <- list("CD4 IFN Response" = c("Ifit3", "Ifit3b", "Isg15", "Ifit1"),
                 "CD4 Naive/Memory" = c("Lef1", "Sell", "Tcf7", "Ccr7", "Slamf6", "Il7r"),
                 "CD4 Teff" = c("Gzmb", "Prf1", "Xcl1", "Ifng"),
                 "CD4 Tem" = c("Tnfsf11", "Gata3", "Itgb1", "Cd40lg", "Junb", "Fos", "Atf3"),
                 "CD4 Tex" = c("Tox", "Fyn", "Havcr2", "Gzma"),
                 "CD4 Texpanding" = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Ifitm2", "Ifitm3"),
                 "CD4 Tpex" = c("Pdcd1", "Cd7", "Bcl2", "Ccl5"),
                 "CD4 Treg" = c("Foxp3", "Ctla4", "Il2ra", "Tnfrsf4", "Ikzf2",
                                "Lag3", "Tnfrsf9", "Pglyrp1", "Tigit"))

pdf("9.CD4_T_Subtypes_DotPlot.pdf", height = 3.5, width = 17)
DotPlot(CD4_T, features = features, cols = c("gray", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))
dev.off()


##### Monocle2
library(monocle)
data <- as.matrix(CD4_T@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CD4_T@meta.data)
fData <- data.frame(gene_short_name = row.names(data),
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4)
mycds <- setOrderingFilter(mycds, VariableFeatures(CD4_T))
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2,
                         reduction_method = 'DDRTree',
                         residualModelFormulaStr = "~Sample")
mycds <- orderCells(mycds)

plot1 <- plot_cell_trajectory(mycds, cell_size = 0.5,
                              color_by = "SubCelltype") +
  scale_color_manual(values = colors) +
  facet_grid(".~SubCelltype")
pdf("10.Monocle2_CD4_T_Subtypes.pdf",
    height = 4, width = 17)
plot1
dev.off()

# Monocle3
library(monocle3)
##创建CDS对象并预处理数据
data <- GetAssayData(CD4_T, assay = 'RNA', slot = 'counts')
cell_metadata <- CD4_T@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 30)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP",
                 color_cells_by="SubCelltype") +
  ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CD4_T, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP",
                 color_cells_by="SubCelltype") +
  ggtitle('int.umap')
cds <- cluster_cells(cds, resolution = 0.001)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) +
  ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition",
                 show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p <- patchwork::wrap_plots(p1, p2)
cds <- learn_graph(cds)
cds <- order_cells(cds)
pdf("5.monocle3.pdf", height = 4.5, width = 5.5)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()



# ########
# ########
# ########
# library(ProjecTILs)
# library(Seurat)
# 
# ref <- readRDS("../ref_TILAtlas_mouse_v1.rds")
# refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
#              "#87f6a5", "#e812dd")
# DimPlot(ref, label = T, cols = refCols)
# 
# query.projected <- Run.ProjecTILs(CD4_T, ref = ref, filter.cells = FALSE)
# 
# 
# table(query.projected@meta.data[["functional.cluster"]])
# 
# plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5)
# table(project_meta$functional.cluster)
# project_meta <- query.projected@meta.data
# quantile(project_meta$functional.cluster.conf)
# ggplot(data = project_meta, aes(x = functional.cluster,
#                                 y = functional.cluster.conf)) +
#   geom_boxplot()
# 
# genes4radar <- c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Sell", "Gzmb", "Gzmk", "Pdcd1",
#                  "Havcr2", "Tox", "Mki67")
# 
# plot.states.radar(ref, query = query.projected, genes4radar = genes4radar, min.cells = 20)
# 
# project_meta_filter <- project_meta[project_meta$functional.cluster.conf > 0.8,]
# table(project_meta_filter$functional.cluster)
# 
