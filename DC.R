setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/DC")
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(clustree)
library(dplyr)
library(patchwork)
library(ggforce)

percent_bar <- function(sce, Ident, Group,
                        fill_color = NULL, fill_label = "",
                        flow = FALSE) {
  if (flow == TRUE) {
    temp <- data.frame(Ident = sce@meta.data[,Ident],
                       Group = sce@meta.data[,Group])
    temp <- as.data.frame.array(table(temp$Ident, temp$Group))
    for (i in 1:ncol(temp)) {
      temp[,i] <- temp[,i] / sum(temp[,i])
    }
    temp <- data.frame(Ident = rownames(temp),
                       temp)
    dat <- reshape::melt(temp, id = "Ident")
    dat$Ident <- factor(dat$Ident, levels = sort(unique(sce@meta.data[,Ident])))
    
    p <- ggplot(dat, aes(x = variable, y = value, fill = Ident, 
                         stratum = Ident, alluvium = Ident)) +
      ggalluvial::geom_stratum(width = 0.55) +  #代替 geom_col() 绘制堆叠柱形图
      ggalluvial::geom_flow(width = 0.55, alpha = 1) +  #绘制同类别之间的连接线\
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(x = "", y = "Percentage", fill = fill_label) +
      theme(axis.text = element_text(size = 13, color = "black"),
            axis.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 13, color = "black"),
            legend.title = element_text(size = 13, color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
  } else {
    temp <- data.frame(Ident = sce@meta.data[,Ident],
                       Group = sce@meta.data[,Group])
    temp$Ident <- factor(temp$Ident,
                         levels = sort(unique((sce@meta.data[,Ident]))))
    temp$Group <- factor(temp$Group,
                         levels = sort(unique((sce@meta.data[,Group]))))
    p <- ggplot(data = temp, aes(x = Group, fill = Ident)) +
      geom_bar(stat = "count", position = "fill", width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(x = "", y = "Percentage", fill = fill_label) +
      theme(axis.text = element_text(size = 13, color = "black"),
            axis.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 13, color = "black"),
            legend.title = element_text(size = 13, color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  if (!is.null(fill_color)) {
    p <- p + scale_fill_manual(values = fill_color)
  }
  p
}

merged_filtered_2 <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/第一次注释_2.Data2_Annotated.rds")
table(merged_filtered_2$Celltype)
DC <- subset(merged_filtered_2, subset = Celltype %in% c("DC1", "DC2", "mregDC", "pDC"))
DC$Sample <- as.character(DC$Sample)
DC %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> DC

DC <- RunUMAP(DC,
                 reduction = "pca",
                 dims = 1:30, n.neighbors = 50,
                 min.dist = 1.5, n.epochs = 1000)

DimPlot(DC, group.by = "Celltype")

DC %>%
  FindNeighbors(dims = 1:30, reduction = "pca",
                k.param = 50) -> DC
DC %>%
  FindClusters(resolution = seq(from = 0.1, by = 0.1, to = 2.5),
               group.singletons = TRUE) -> DC

DimPlot(DC)

pdf("./1.DC_UMAP_Sample.pdf", height = 4.5, width = 5.5)
DimPlot(DC, group.by = "Sample")
dev.off()

pdf("./2.DC_UMAP_DataSets.pdf", height = 4.5, width = 5.5)
DimPlot(DC, group.by = "Batch")
dev.off()

pdf("./3.DC_UMAP_Clusters.pdf")
for (i in seq(from = 0.1, by = 0.1, to = 2.5)) {
  p <- DimPlot(DC, label = T, group.by = paste0("RNA_snn_res.",i))
  print(p)
}
dev.off()

DC %>%
  FindClusters(resolution = 1.5,
               group.singletons = TRUE) -> DC

pdf("./4.DC_UMAP_Clusters_resolution_1.5.pdf",
    height = 4.5, width = 5.4)
DimPlot(DC, label = T)
dev.off()


pdf("./4.DC_UMAP_Clusters_resolution_1.5_splited.pdf",
    height = 7, width = 10)
DimPlot(DC, label = T,
        split.by = "seurat_clusters", ncol = 4)
dev.off()

exp <- as.matrix(DC@assays$RNA@data)
exp <- t(exp)
cluster_mean <- aggregate.data.frame(exp,
                                     by = list(as.character(DC@meta.data[rownames(exp), "seurat_clusters"])),
                                     FUN = mean)
cluster_mean$Group.1 <- paste0("Seurat_cluster_", cluster_mean$Group.1)
rownames(cluster_mean) <- cluster_mean$Group.1
cluster_mean <- cluster_mean[,-1]
cluster_mean <- t(cluster_mean)
corr <- cor(cluster_mean, method = "pearson")
pdf("./5.DC_seurat_cluster_cor_heatmap_all_genes.pdf",
    width = nrow(corr) * 1 + 1, height = nrow(corr) * 1)
p <- pheatmap::pheatmap(corr, angle_col = 90, treeheight_row = 45,
                        fontsize = 14, treeheight_col = 45, #cutree_rows = 11, cutree_cols = 11,
                        cellwidth = 35, cellheight = 35, display_numbers = TRUE)
dev.off()

Idents(DC) <- DC$RNA_snn_res.1.5
Cluster_diff <- FindAllMarkers(DC, only.pos = T,
                               logfc.threshold = 0.25)
Cluster_diff <- Cluster_diff[order(Cluster_diff$cluster,
                                   -Cluster_diff$avg_log2FC),]
Cluster_diff_list <- split.data.frame(Cluster_diff, f = list(Cluster_diff$cluster))
openxlsx::write.xlsx(Cluster_diff_list,
                     "./6.Cluster_DEG_resolution_1.5_log2FC_0.25.xlsx")

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
openxlsx::write.xlsx(Cluster_GO, "./7.DC_Cluster_GO_resolution_1.5.xlsx")

saveRDS(DC, "DC_subclusters.rds")

DC <- readRDS("DC_subclusters.rds")
DimPlot(DC, label = T)
DimPlot(DC, group.by = "Cluster_Celltype")
DimPlot(DC, group.by = "TCR")

DC <- RenameIdents(DC,
                   "0" = "cDC1",
                   "1" = "mregDC1",
                   "2" = "cDC1",
                   "3" = "cDC2",
                   "4" = "DC-prolif",
                   "5" = "cDC1",
                   "6" = "Unidentified",
                   "7" = "mregDC2",
                   "8" = "pDC",
                   "9" = "mo-DC"
                   )
DC$SubCelltype <- as.character(Idents(DC))
sort(unique(DC$SubCelltype))

Idents(DC) <- factor(DC$SubCelltype,
                     levels = sort(unique(DC$SubCelltype)))
# colors <- ArchR::ArchRPalettes$kelly
# names(colors) <- NULL
# colors <- c(colors[1:7], "gray")
colors <- scales::hue_pal()(7)
colors <- c(colors[1:7], "gray")
names(colors) <- sort(unique(DC$SubCelltype))
pdf("8.DC_Subtypes.pdf", width = 5.2, height = 4)
DimPlot(DC, label = T, cols = colors, pt.size = 0.7)
dev.off()

saveRDS(DC, "DC_annotated.rds")

features <- list("DC" = c("H2-Ab1", "H2-Aa", "Cd74", "H2-DMb1", "H2-DMa"),
                 "cDC1" = c("Clec9a", "Xcr1"),
                 "cDC2" = c("Cd209a", "Itgam", "Ifitm1"),
                 "DC−prolif" = c("Stmn1", "Top2a"),
                 "mregDC" = c("Fscn1", "Ccr7", "Ccl22"),
                 "pDC" = c("Cox6a2", "Siglech", "Ccr9", "Cd300c"))

pdf("9.DC_Subtypes_DotPlot.pdf", height = 3.5, width = 11)
DotPlot(DC, features = features, cols = c("gray", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))
dev.off()

DC <- readRDS("DC_annotated.rds")
DimPlot(DC)

table(DC$SubCelltype)
DC_sub <- subset(DC, subset = SubCelltype %in% c("cDC1", "cDC2",
                                                 "DC-prolif", "mregDC1", "mregDC2"))
DC_sub %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> DC_sub
DC_sub <- RunUMAP(DC_sub,
                  reduction = "pca",
                  dims = 1:30, n.neighbors = 50,
                  min.dist = 1.5, n.epochs = 1000)
DimPlot(DC_sub, group.by = "Sample")
Sample_rename <- data.frame(original = c("Vector-T_2", "XCL1-T_2", "Fx-T_2", "CD40-T_2"),
                            rename = c("OT-I",
                                       "OT-I/XCL1",
                                       "OT-I/FX",
                                       "OT-I/FX + CD40 mAb"))
rownames(Sample_rename) <- Sample_rename$original
DC_sub$Sample_rename <- Sample_rename[DC_sub$Sample, "rename"]
sample_color <- c("#1E699D", "#8F94C0", "#D87177", "#B42225")
names(sample_color) <- c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb")

DC_sub$Sample_rename <- factor(DC_sub$Sample_rename,
                               levels = c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb"))
pdf("补充.1.经典DC_sample_umap.pdf",
    width = 6.2, height = 4.5)
DimPlot(DC_sub, cols = sample_color, group.by = "Sample_rename") +
  ggtitle("Sample")
dev.off()

pdf("补充.2.经典DC_celltype_umap.pdf",
    width = 5.5, height = 4.5)
DimPlot(DC_sub, cols = colors,
        group.by = "SubCelltype") +
  ggtitle("Celltype")
dev.off()

DC_sub_2 <- subset(DC_sub, subset = Sample_rename != "OT-I/FX + CD40 mAb")
Idents(DC_sub_2) <- DC_sub_2$Sample_rename
pdf("补充_V2.7.经典DC_dotplot_geneset_1.pdf",
    height = 4, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Lif", "Il1b", "Il12b", "Il15",
                         "Tnf", "Cxcl9", "Cxcl10",
                         "Cxcl16", "Ccl17", "Ccl22"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_2.pdf",
    height = 3, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Relb", "Cd86", "Cd83", "Cd80", "Cd40"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_3.pdf",
    height = 4, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Myo1g", "Marcks11", "Marcks", "Icam1", "Fscn1", "Cxcl16", "Ccr7", "Adam8"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_4.pdf",
    height = 4.5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Tlr9", "Tlr8", "Tlr7",
                         "Tlr6", "Tlr4", "Tlr3", "Tlr2", "Tlr1",
                         "Myd88", "Mavs"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_5.pdf",
    height = 5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Wdfy4", "Tapbp1", "Tapbp", "Tap2",
                         "Tap1", "Rab7", "Psme3", "Mfge8", "Lamp2",
                         "Lamp1", "Ctss", "Ciita", "Cd74"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_6.pdf",
    height = 6, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Tnfaip3", "Slfn2", "Pycard", "Pstpip1",
                         "Nlrp3", "Nlrc4", "Il1b", "Il18bp",
                         "Il18", "Ctsl", "Ctsb", "Cd14", "Casp4",
                         "Casp1", "Aim2"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_maturation.pdf",
    height = 3, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Cd40", "Cd80", "Cd86", "Cd83", "Relb"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_regulatory.pdf",
    height = 5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Cd274", "Pdcd1lg2", "Cd200", "Fas",
                         "Aldh1a2", "Socs1", "Socs2", "Ido1", "Il4i1",
                         "Cd63", "Cd39"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_TLR.pdf",
    height = 5.5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Myd88", "Mavs", "Tlr1", "Tlr2", "Tlr3",
                         "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr8",
                         "Tlr9", "Tlr10", "Tlr11"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充_V2.3.经典DC_dotplot_geneset_migration.pdf",
    height = 4.8, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Ccr7", "Myo1g", "Cxcl16", "Adam8",
                         "Icam1", "Fscn1", "Marcks", "Marcksl1",
                         "Il15", "Cxcl9", "Cxcl10", "Il12b"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_Th2.pdf",
    height = 4.3, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("Il4ra", "Il4i1", "Ccl17", "Ccl22", "Tnfrsf4",
                         "Stat6", "Bcl2l1", "Irf4"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_MHC.pdf",
    height = 4.3, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("H2-D1", "H2-K1", "B2m", "H2-Aa",
                         "H2-Ab1", "H2-Eb1", "H2-Eb2", "Cd74"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_MHC_I.pdf",
    height = 2, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("H2-D1", "H2-K1", "B2m"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

pdf("补充.3.经典DC_dotplot_geneset_MHC_II.pdf",
    height = 3.5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("H2-Aa", "H2-Ab1",
                         "H2-Eb1", "H2-Eb2", "Cd74"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()



DC_sub_2$Sample_rename <- factor(as.character(DC_sub_2$Sample_rename),
                                 levels = c("OT-I", "OT-I/XCL1", "OT-I/FX"))
pdf("补充.4.经典DC_celltype_percent_cDC.pdf", width = 4, height = 3.5)
percent_bar(DC_sub_2, Ident = "SubCelltype",
            Group = "Sample_rename",
            fill_color = colors, fill_label = "",
            flow = FALSE) +
  labs(y = "Percentage of total cDC")
dev.off()

DC$Sample_rename <- Sample_rename[DC$Sample, "rename"]
DC$Sample_rename <- factor(DC$Sample_rename,
                           levels = c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb"))
temp <- as.data.frame.array(table(DC$SubCelltype, DC$Sample_rename))
temp <- t(temp) / colSums(temp)
temp <- reshape2::melt(temp)
temp <- temp[temp$Var2 %in% c("cDC1", "cDC2",
                              "DC-prolif", "mregDC1", "mregDC2"),]
temp <- temp[temp$Var1 != "OT-I/FX + CD40 mAb",]
pdf("补充.4.经典DC_celltype_percent_all_DC.pdf", width = 4, height = 3.5)
ggplot(data = temp, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  theme_classic() +
  labs(x = "", y = "Percentage of total DC", fill = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

library(AUCell)
AUCell_score <- AUCell_run(exprMat = DC_sub_2@assays$RNA@data,
           geneSets = list("mregDC signature" = c("Cd40", "Cd80", "Cd86", "Cd83", "Relb",
                                                  "Cd274", "Pdcd1lg2", "Cd200", "Fas",
                                                  "Aldh1a2", "Socs1", "Socs2", "Ido1", "Il4i1",
                                                  "Cd63", "Cd39"),
                           "MHC class I signature" = c("H2-D1", "H2-K1", "B2m"),
                           "MHC class II signature" = c("H2-Aa", "H2-Ab1", "H2-Eb1",
                                                        "H2-Eb2", "Cd74"),
                           "MHC class I&II signature" = c("H2-D1", "H2-K1", "B2m", "H2-Aa",
                                                          "H2-Ab1", "H2-Eb1", "H2-Eb2", "Cd74")))
AUCell_score <- AUCell_score@assays@data@listData[["AUC"]]
AUCell_score <- data.frame(t(AUCell_score), check.rows = F, check.names = F)
AUCell_score$Sample <- DC_sub_2@meta.data[rownames(AUCell_score),
                                          "Sample_rename"]
pdf("补充.5.经典DC_ECDF_mregDC signature.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `mregDC signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("mregDC signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15))
dev.off()

pdf("补充.5.经典DC_ECDF_MHC class I signature.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `MHC class I signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("MHC class I signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15))
dev.off()

pdf("补充.5.经典DC_ECDF_MHC class II signature.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `MHC class II signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("MHC class II signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15))
dev.off()

pdf("补充.5.经典DC_ECDF_MHC class I&II signature.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `MHC class I&II signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("MHC class I&II signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15))
dev.off()


DC2 <- subset(DC_sub_2, subset = SubCelltype == "cDC2")
DC2 <- NormalizeData(DC2)
DC2 <- ScaleData(DC2, features = rownames(DC2))
pdf("补充_V2.2.cDC2_heatmap.pdf", height = 5, width = 5)
DoHeatmap(DC2, features = c("Ffar2", "Ltb", "H2-DMb2", "Ifngr1", "H2-Eb1",
                            "Cd72", "Dab2", "H2-Oa", "Klrd1", "Cfp", "Il1r2",
                            "Cd300c2", "Clec4a3", "Ifitm1", "Mgl2", "Ifit1bl1",
                            "Ifit3b", "Ifit3", "Ifit1", "Ifi204", "Cxcl10",
                            "Ms4a4c", "Ifit2", "Isg15", "Fcgr1", "Rsad2", "Phf11d",
                            "Isg20", "Usp18", "Ly6a", "Axl", "Cd5"), group.by = "Sample_rename",
          group.colors = sample_color,
          disp.min = -1.8,
          disp.max = 1.8,
          size = 4)
dev.off()


library(AUCell)
AUCell_score <- AUCell_run(exprMat = DC_sub_2@assays$RNA@data,
                           geneSets = list("mregDC signature" = c("Cd40", "Cd80", "Cd86", "Cd83", "Relb",
                                                                  "Cd274", "Pdcd1lg2", "Cd200", "Fas",
                                                                  "Aldh1a2", "Socs1", "Socs2", "Ido1", "Il4i1",
                                                                  "Cd63", "Cd39"),
                                           "MHC class I signature" = c("H2-Dma", "H2-M2", "H2-T23", "Tap1", "Tap2",
                                                                       "Tapbpl", "Ctsl", "Ctss", "Ctsb", "Ciita",
                                                                       "Lamp1", "Lamp2", "Cybb", "Wdfy4", "Cd74",
                                                                       "Apol7c", "Mfge8", "Psme3"),
                                           "Inflammasome" = c("Aim2", "Nlrp3", "Ctsl", "Ctsb", "Casp1", "Casp4",
                                                              "Il1b", "Il18", "Il18bp", "Tnfaip3", "Slfn2", "Pycard",
                                                              "Pstpip1"))
                           )
AUCell_score <- AUCell_score@assays@data@listData[["AUC"]]
AUCell_score <- data.frame(t(AUCell_score), check.rows = F, check.names = F)
AUCell_score$Sample <- DC_sub_2@meta.data[rownames(AUCell_score),
                                          "Sample_rename"]
AUCell_score$SubCelltype <- DC_sub_2@meta.data[rownames(AUCell_score),
                                               "SubCelltype"]
temp <- AUCell_score[AUCell_score$SubCelltype %in% c("mregDC1", "mregDC2"),]
leveneTest(`mregDC signature` ~ Sample, data = temp, center = mean)
tapply(temp$`mregDC signature`, temp$Sample, shapiro.test)
anova_result <- aov(`mregDC signature` ~ Sample, data = temp)
summary(anova_result)
# library(goftest)
# XCL1_OT_I <- ks.test(temp$`mregDC signature`[temp$Sample == "OT-I/XCL1"],
#                          temp$`mregDC signature`[temp$Sample == "OT-I"])
# XCL1_OT_I$p.value
# FX_OT_I <- ks.test(temp$`mregDC signature`[temp$Sample == "OT-I/FX"],
#                        temp$`mregDC signature`[temp$Sample == "OT-I"])
# FX_OT_I$p.value
pdf("补充_V3.3.mregDC1_2_ECDF_mregDC signature.pdf",
    height = 4, width = 5.4)
ggplot(temp, aes(x = `mregDC signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("mregDC signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  annotate("text", x = 0.3, y = 0.12, label = expression(italic("P")["ANOVA test"]~" = 0.0387"),
           color = "black", size = 4, hjust = 0, vjust = 0)
dev.off()


# XCL1_OT_I <- wilcox.test(AUCell_score$`Inflammasome`[AUCell_score$Sample == "OT-I/XCL1"],
#                          AUCell_score$`Inflammasome`[AUCell_score$Sample == "OT-I"])
# XCL1_OT_I$pvalue
# FX_OT_I <- StatComp18087::cvm.test(AUCell_score$`Inflammasome`[AUCell_score$Sample == "OT-I/FX"],
#                                    AUCell_score$`Inflammasome`[AUCell_score$Sample == "OT-I"])
# FX_OT_I$pvalue
anova_result <- aov(`Inflammasome` ~ Sample, data = AUCell_score)
summary(anova_result)
pdf("补充_V3.4.经典DC_ECDF_Inflammasome.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `Inflammasome`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("Inflammasome") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  annotate("text", x = 0.2, y = 0.12, label = expression(italic("P")["ANOVA test"]~" = 0.0199"),
           color = "black", size = 4, hjust = 0, vjust = 0)
dev.off()

anova_result <- aov(`MHC class I signature` ~ Sample, data = AUCell_score)
summary(anova_result)
pdf("补充_V3.5.经典DC_ECDF_MHC class I signature.pdf",
    height = 4, width = 5.4)
ggplot(AUCell_score, aes(x = `MHC class I signature`, colour = Sample)) +
  stat_ecdf(geom = "step", size = 0.7) +
  scale_color_manual(values = sample_color) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "AUCell score", y = "ECDF") +
  ggtitle("MHC class I signature") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  annotate("text", x = 0.25, y = 0.12, label = expression(italic("P")["ANOVA test"]~" = 0.215"),
           color = "black", size = 4, hjust = 0, vjust = 0)
dev.off()


pdf("补充_V2.6.经典DC_dotplot_geneset_MHC_I.pdf",
    height = 6.5, width = 3.5)
DotPlot(DC_sub_2,
        features = rev(c("H2-Dma", "H2-M2", "H2-T23", "Tap1", "Tap2",
                         "Tapbpl", "Ctsl", "Ctss", "Ctsb", "Ciita",
                         "Lamp1", "Lamp2", "Cybb", "Wdfy4", "Cd74",
                         "Apol7c", "Mfge8", "Psme3"))) +
  theme_bw() +
  coord_flip() +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  
  labs(x = "", y = "")
dev.off()

