{
  library(Seurat)
  library(ggplot2)
  library(dplyr)
}
All_seu <- readRDS("第一次注释_3.Data2_Annotated_with_CD45_classification.rds")
unique(All_seu$Celltype)
Idents(All_seu) <- All_seu$Celltype
All_seu <- RenameIdents(All_seu,
                        "CD8+ T" = "CD8",
                        "NK" = "NK",
                        "Macrophage" = "Mono/Macs",
                        "CD4+ T" = "CD4",
                        "Neutrophil" = "Neutrophil",
                        "DPT" = "DPT",
                        "Monocyte" = "Mono/Macs",
                        "pDC" = "pDC",
                        "B cell" = "B",
                        "Cancer cell" = "Melanocyte",
                        "DC2" = "DC",
                        "Mast cell" = "Mast cell",
                        "mregDC" = "DC",
                        "DC1" = "DC",
                        "Erythroblast" = "Erythroblast"
                        )
All_seu$Celltype_rename <- as.character(Idents(All_seu))
All_seu$Celltype_rename <- factor(All_seu$Celltype_rename,
                                  levels = c("B", "CD4", "CD8", "DPT", "NK",
                                             "DC", "pDC", "Mono/Macs", "Neutrophil",
                                             "Erythroblast", "Mast cell", "Melanocyte"))
Idents(All_seu) <- All_seu$Celltype_rename
All_seu <- RunTSNE(All_seu,
                   reduction = "pca",
                   dims = 1:30)
set.seed(136)
Celltype_color <- randomcoloR::distinctColorPalette(k = 12)
names(Celltype_color) <- levels(All_seu$Celltype_rename)
pdf("./绘图/1.所有细胞_UMAP.pdf",
    height = 4.2, width = 5.7)
DimPlot(All_seu, cols = Celltype_color) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(size = 10 , colour = "black"))
dev.off()

pdf("./绘图/1.所有细胞_TSNE.pdf",
    height = 4.2, width = 5.7)
DimPlot(All_seu, cols = Celltype_color, reduction = "tsne") +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(size = 10 , colour = "black"))
dev.off()


All_seu$Cluster_Celltype <- paste0(as.character(All_seu$seurat_clusters),
                                   " ", as.character(All_seu$Celltype_rename))
# All_seu$Cluster_Celltype <- factor(All_seu$Cluster_Celltype,
#                                    levels = c("16 B",
#                                               "2 CD4", "7 CD4", "10 CD4", "13 CD4", "20 CD4",
#                                               "5 CD8", "6 CD8", "12 CD8", "15 CD8", "19 CD8", "26 CD8",
#                                               "27 DPT",
#                                               "1 NK", "3 NK", "8 NK", "18 NK", "21 NK",
#                                               "14 DC", "24 DC", "28 DC",
#                                               "30 pDC",
#                                               "0 Mono/Macs", "4 Mono/Macs", "9 Mono/Macs", "11 Mono/Macs", "17 Mono/Macs", "22 Mono/Macs",
#                                               "25 Neutrophil",
#                                               "29 Erythroblast",
#                                               "31 Mast cell",
#                                               "23 Melanocyte"
#                                               ))

All_seu$Cluster_Celltype <- factor(All_seu$Cluster_Celltype,
                                   levels = c("0 Mono/Macs", "1 NK", "2 CD4", "3 NK",
                                              "4 Mono/Macs","5 CD8", "6 CD8", "7 CD4",
                                              "8 NK", "9 Mono/Macs", "10 CD4", "11 Mono/Macs",
                                              "12 CD8", "13 CD4", "14 DC", "15 CD8", "16 B",
                                              "17 Mono/Macs", "18 NK", "19 CD8", "20 CD4",
                                              "21 NK", "22 Mono/Macs", "23 Melanocyte",
                                              "24 DC", "25 Neutrophil", "26 CD8", "27 DPT",
                                              "28 DC", "29 Erythroblast", "30 pDC", "31 Mast cell"
                                              ))

set.seed(126)
cluster_color <- randomcoloR::distinctColorPalette(k = 32)
UMAP_meta <- Embeddings(All_seu, reduction = "umap")
UMAP_meta <- as.data.frame(UMAP_meta)
UMAP_meta$Cluster <- as.character(All_seu@meta.data[rownames(UMAP_meta), "seurat_clusters"])
UMAP_meta <- aggregate.data.frame(UMAP_meta[,1:2],
                                  by = list(UMAP_meta$Cluster),
                                  FUN = mean)
colnames(UMAP_meta)[1] <- "Clusters"

TSNE_meta <- Embeddings(All_seu, reduction = "tsne")
TSNE_meta <- as.data.frame(TSNE_meta)
TSNE_meta$Cluster <- as.character(All_seu@meta.data[rownames(TSNE_meta), "seurat_clusters"])
TSNE_meta <- aggregate.data.frame(TSNE_meta[,1:2],
                                  by = list(TSNE_meta$Cluster),
                                  FUN = mean)
colnames(TSNE_meta)[1] <- "Clusters"

pdf("./绘图/2.所有细胞_cluster_TSNE.pdf",
    width = 6.5, height = 4.2)
DimPlot(All_seu, group.by = "Cluster_Celltype",
        reduction = "tsne", cols = cluster_color) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(size = 10 , colour = "black")) +
  ggtitle("") +
  geom_point(data = TSNE_meta, aes(x = tSNE_1, y = tSNE_2),
             size = 5, alpha = 0.7, color = "#EEE9E9") +
  geom_text(data = TSNE_meta, aes(x = tSNE_1, y = tSNE_2, label = Clusters),
             size = 3.3, color = "black")
dev.off()

pdf("./绘图/2.所有细胞_cluster_UMAP.pdf",
    width = 6.5, height = 4.2)
DimPlot(All_seu, group.by = "Cluster_Celltype",
        reduction = "umap", cols = cluster_color) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(size = 10 , colour = "black")) +
  ggtitle("") +
  geom_point(data = UMAP_meta, aes(x = UMAP_1, y = UMAP_2),
             size = 5, alpha = 0.7, color = "#EEE9E9") +
  geom_text(data = UMAP_meta, aes(x = UMAP_1, y = UMAP_2, label = Clusters),
            size = 3.3, color = "black")
dev.off()


Celltype_markers_list <- list("B" = c("Ly6d", "Cd79a", "Cd79b", "Ms4a1"),
                              "T" = c("Cd3e", "Cd3d", "Cd3g"),
                              "CD4" = c("Cd4", "Tcf7", "Bcl11b"),
                              "CD8" = c("Cd8a", "Cd8b1", "Tox"),
                              "NK" = c("Klre1", "Ncr1", "Prf1", "Gzma", "Klrb1c"),
                              "DC" = c("H2-Ab1", "H2-DMa", "H2-Aa", "H2-DMb1", "Cd74"),
                              "pDC" = c("Cox6a2", "Siglech", "Bst2", "Ccr9", "Cd300c"),
                              "Mono/Macs" = c("Cd68", "Ms4a7", "Apoe", "C1qa",
                                              "C1qb", "Ccl8", "C1qc", "Lyz2"),
                              "Neutrophil" = c("S100a9", "S100a8", "G0s2", "Csf3r"),
                              "Erythroblast" = c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2"),
                              "Mast cell" = c("Gata2", "Cpa3"),
                              "Melanocyte" = c("Pmel", "Mlana", "Ptgds", "Dct")
                              )

Idents(All_seu) <- All_seu$Celltype_rename
pdf("./绘图/3所有细胞_celltype_marker_dotplot.pdf",
    height = 5.5, width = 15)
DotPlot(All_seu,
        features = Celltype_markers_list) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10, colour = "black"),
        strip.background = element_rect(color = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "right") +
  scale_color_gradient2(low = "#1E90FF", high = "#EE2C2C") +
  labs(y = "", x = "")
dev.off()

Idents(All_seu) <- All_seu$Sample
All_seu <- RenameIdents(All_seu,
                        "CD40-T_2" = "OT-I/FX + CD40 mAb",
                        "Fx-T_2" = "OT-I/FX",
                        "XCL1-T_2" = "OT-I/XCL1",
                        "Vector-T_2" = "OT-I"
                        )
All_seu$Sample_rename <- factor(as.character(Idents(All_seu)),
                                levels = c("OT-I",
                                           "OT-I/XCL1",
                                           "OT-I/FX",
                                           "OT-I/FX + CD40 mAb"))
Idents(All_seu) <- All_seu$Celltype_rename
pdf("./绘图/4.所有细胞_UMAP_split_by_Sample.pdf",
    height = 4, width = 14)
DimPlot(All_seu, cols = Celltype_color, split.by = "Sample_rename") +
  theme_classic() +
  ggtitle("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        strip.text.x.top = element_text(size = 10, colour = "black", face = "bold"),
        strip.background = element_rect(color = "white", fill = "white"),
        panel.grid = element_blank())
dev.off()

pdf("./绘图/4.所有细胞_TSNE_split_by_Sample.pdf",
    height = 4, width = 14)
DimPlot(All_seu, cols = Celltype_color,
        split.by = "Sample_rename", reduction = "tsne") +
  theme_classic() +
  ggtitle("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        strip.text.x.top = element_text(size = 10, colour = "black", face = "bold"),
        strip.background = element_rect(color = "white", fill = "white"),
        panel.grid = element_blank())
dev.off()

percent_bar <- function(sce, Ident, Group, fill_color = NULL, fill_label = "",
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
pdf("./绘图/5.所有细胞_细胞类型比例.pdf",
    width = 5, height = 5)
percent_bar(sce = All_seu, Ident = "Celltype_rename",
            Group = "Sample_rename",
            fill_color = Celltype_color)
dev.off()

genes <- c("Cd79a", "Cd3d", "Cd8a", "Cd4", "Foxp3",
           "Ncr1", "Cd74", "H2-Ab1", "Clec9a", "Cd209a",
           "Fscn1", "Siglech", "Cd68", "Lyz2", "G0s2", "Hbb-bs",
           "Cpa3", "Mlana")
feature_plot_list <- lapply(genes, function(x){
  FeaturePlot(All_seu, features = c(x), order = T,
              min.cutoff = 1) +
    scale_color_gradientn(colours = c("gray", "#124A9B")) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          plot.title = element_text(size = 14, colour = "black",
                                    face = "bold", hjust = 0.5),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank())
})
pdf("./绘图/6.所有细胞_细胞marker_featurePlot.pdf",
    width = 10.5, height = 18)
patchwork::wrap_plots(feature_plot_list,
                      byrow = F, nrow = 6, ncol = 3)
dev.off()

Idents(All_seu) <- All_seu$Celltype_rename
Celltype_DEGs <- FindAllMarkers(All_seu, only.pos = T, logfc.threshold = 1)

library(dplyr)
Celltype_DEGs %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
top20 <- as.data.frame(top20)
top20 <- top20[-which(top20$cluster == "DPT"),]

Celltype_DEGs %>%
  group_by(cluster) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3
top3 <- as.data.frame(top3)
top3 <- top3[-which(top3$cluster == "DPT"),]
top3 <- top3[!duplicated(top3$gene),]

library(ComplexHeatmap)
library(RColorBrewer)
Cell_group_heatmap <- function(seu,
                               group.by = "Celltype_rename",
                               group.color = Celltype_color,
                               group.label = "Celltype",
                               genes = unique(top5$gene),
                               label_genes = unique(top5$gene)[unique(top5$gene) %in% unlist(lapply(Celltype_markers_list, FUN = as.character))],
                               limit = c(-2, 2)){
  exp <- seu@assays$RNA@data[genes,]
  exp <- t(as.matrix(exp))
  exp <- aggregate.data.frame(exp,
                              by = list(seu@meta.data[,group.by]),
                              FUN = mean)
  rownames(exp) <- exp$Group.1
  exp <- exp[,-1]
  for (i in 1:ncol(exp)) {
    exp[,i] <- (exp[,i] - mean(exp[,i])) / sd(exp[,i])
  }
  exp <- as.matrix(t(exp))
  # exp <- exp[genes_order,]
  exp[exp > limit[2]] <- limit[2]
  exp[exp < limit[1]] <- limit[1]
  neg_length <- 0 - min(exp)
  pos_length <- max(exp) - 0
  sum <- neg_length + pos_length
  neg_length <- 100 * (neg_length / sum)
  neg_length <- round(neg_length)
  pos_length <- 100 * (pos_length / sum)
  pos_length <- round(pos_length)
  neg_color <- colorRampPalette(c("#174C9E", "white"))(neg_length)
  pos_color <- colorRampPalette(c("white", "#E8402B"))(pos_length)
  color_bar <- c(neg_color, pos_color)
  color_map <- circlize::colorRamp2(seq(from = min(exp), to = max(exp), length = length(color_bar)),
                          color_bar)
  group <- apply(exp, 1, function(x){which.max(x)})
  exp_list <- split.data.frame(exp, f = group, drop = F)
  for (i in 1:length(exp_list)) {
    exp_list[[i]] <- exp_list[[i]][order(exp_list[[i]][,i],
                                         decreasing = T),,drop = F]
  }
  exp <- yulab.utils::rbindlist(exp_list)
  col_anno <- HeatmapAnnotation(df = data.frame(Group = colnames(exp)),
                                col = list(Group = group.color),
                                annotation_label = group.label,
                                show_annotation_name = FALSE)
  if (!is.null(label_genes)) {
    row_label <- rowAnnotation(foo = anno_mark(at = match(label_genes, rownames(exp)),
                                               labels = label_genes,
                                               labels_gp = gpar(fontsize = 10),
                                               padding = unit(0, "mm")))
    
    Heatmap(exp, col = color_map, color_space = "RGB",
            cluster_rows = F, cluster_columns = F,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            show_row_names = F,
            right_annotation = row_label, top_annotation = col_anno,
            heatmap_legend_param = list(title = "Scaled\nExpression"#,
                                        # at = c(round(min(exp),1), 0, max(exp)),
                                        # labels = c(round(min(exp),1), 0, max(exp)))
                                        ))
  } else {
    Heatmap(exp, col = color_map, color_space = "RGB",
            cluster_rows = F, cluster_columns = F,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            top_annotation = col_anno,
            heatmap_legend_param = list(title = "Scaled\nExpression"
                                        # at = c(round(min(exp),1), 0, max(exp)),
                                        # labels = c(round(min(exp),1), 0, max(exp))))
            ))
  }
  
}
pdf("./绘图/7.所有细胞_Celltype_top20_差异基因_heatmap.pdf", height = 8, width = 6)
Cell_group_heatmap(seu = All_seu,
                   group.by = "Celltype_rename",
                   group.color = Celltype_color,
                   group.label = "Celltype",
                   genes = unique(top20$gene),
                   label_genes = unique(top3$gene),
                   limit = c(-2, 2))
dev.off()

library(clusterProfiler)
library(org.Mm.eg.db)
unique(All_seu$Celltype_rename)
select_cells <- c("CD4", "CD8", "Mono/Macs", "DC", "Neutrophil")
diff_group <- c("OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb")
Celltype_GO <- c()
for (c in select_cells) {
  print(c)
  temp <- subset(All_seu, subset = Celltype_rename == c)
  Idents(temp) <- temp$Sample_rename
  temp_celltype <- c()
  diff_group_2 <- names(table(temp$Sample_rename))[table(temp$Sample_rename) > 5]
  diff_group_2 <- diff_group[diff_group %in% diff_group_2]
  for (g in diff_group_2) {
    temp_diff <- FindMarkers(temp, ident.1 = g, ident.2 = "OT-I",
                             only.pos = T)
    temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
    temp_gene <- bitr(rownames(temp_diff), fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    temp_gene <- temp_gene$ENTREZID
    temp_go <- enrichGO(temp_gene, keyType = "ENTREZID",
                        ont = "BP", pvalueCutoff = 1,
                        qvalueCutoff = 1, OrgDb = "org.Mm.eg.db")
    temp_go_result <- temp_go@result
    temp_go_result$Description <- Hmisc::capitalize(temp_go_result$Description)
    temp_go_result <- temp_go_result[temp_go_result$pvalue < 0.05,]
    temp_go_result$Group <- g
    temp_celltype <- as.data.frame(rbind(temp_celltype,
                                         temp_go_result))
  }
  Celltype_GO <- c(Celltype_GO,
                   list(temp_celltype))
}
names(Celltype_GO) <- select_cells
length((unique(Celltype_GO$CD4$Group)))
length((unique(Celltype_GO$CD8$Group)))
length((unique(Celltype_GO$`Mono/Macs`$Group)))
length((unique(Celltype_GO$DC$Group)))
length((unique(Celltype_GO$Neutrophil$Group)))

Common_pathway <- c()
for (i in names(Celltype_GO)) {
  temp <- as.data.frame.array(table(Celltype_GO[[i]]$Description, Celltype_GO[[i]]$Group))
  pathways <- rownames(temp)[rowSums(temp) >= 2]
  temp <- Celltype_GO[[i]][Celltype_GO[[i]]$Description %in% pathways,]
  temp$log10P <- -log10(temp$pvalue)
  average <- aggregate.data.frame(temp[,c("pvalue", "Count", "log10P")],
                                  by = list(temp$Description),
                                  FUN = mean)
  average <- average[order(average$log10P, decreasing = T),]
  colnames(average) <- c("Pathway", "pvalue", "Count", "log10P")
  Common_pathway <- c(Common_pathway,
                      list(average))
}
names(Common_pathway) <- names(Celltype_GO)
names(Common_pathway) <- c("CD4", "CD8", "Mono_Macs", "DC", "Neutrophil")
openxlsx::write.xlsx(Common_pathway, "./绘图/8.Common_pathway.xlsx")


Celltype_percent <- as.data.frame.array(table(All_seu$Celltype_rename,
                                              All_seu$Sample_rename))
for (i in 1:ncol(Celltype_percent)) {
  Celltype_percent[,i] <- Celltype_percent[,i] / sum(Celltype_percent[,i])
}
colSums(Celltype_percent)
Celltype_percent <- data.frame(Celltype = rownames(Celltype_percent),
                               Celltype_percent,
                               check.rows = F, check.names = F)
Celltype_percent <- reshape::melt.data.frame(Celltype_percent,
                                             id.vars = "Celltype")
temp <- split.data.frame(Celltype_percent, f = list(Celltype_percent$Celltype))

sample_color <- c("#1E699D", "#8F94C0", "#D87177", "#B42225")
names(sample_color) <- c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb")
percent_list <- lapply(temp, function(X){
  max(X$value)
  ggplot(X,
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_classic() +
    scale_y_continuous(labels = scales::percent,
                       # limits = c(0, range(X$value, na.rm = TRUE)[2]),
                       # expand = c(0, 0),
                       # breaks = seq(from = 0,
                       #              to = range(X$value, na.rm = TRUE)[2],
                       #              length = 5)
                       ) +
    facet_wrap(".~Celltype", scales = "free") +
    labs(y = "Percentage", x = "") +
    scale_fill_manual(values = sample_color) +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, colour = "black"),
          strip.text.x.top = element_text(size = 10, colour = "black", face = "bold"),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank(),
          legend.position = "none")
})
pdf("./绘图/9.所有细胞_细胞比例.pdf",
    width = 10, height = 10)
patchwork::wrap_plots(percent_list,
                      byrow = F, nrow = 3, ncol = 4)
dev.off()

CD4_pathway <- openxlsx::read.xlsx("./绘图/CD4.xlsx", colNames = F)
CD4_pathway <- Common_pathway$CD4[Common_pathway$CD4$Pathway %in% CD4_pathway[,1],]
CD4_pathway$Celltype <- "CD4"
CD8_pathway <- openxlsx::read.xlsx("./绘图/CD8.xlsx", colNames = F)
CD8_pathway <- Common_pathway$CD8[Common_pathway$CD8$Pathway %in% CD8_pathway[,1],]
CD8_pathway$Celltype <- "CD8"
DC_pathway <- openxlsx::read.xlsx("./绘图/DC.xlsx", colNames = F)
DC_pathway <- Common_pathway$DC[Common_pathway$DC$Pathway %in% DC_pathway[,1],]
DC_pathway$Celltype <- "DC"
Macro_pathway <- openxlsx::read.xlsx("./绘图/Mono_Macs.xlsx", colNames = F)
Macro_pathway <- Common_pathway$Mono_Macs[Common_pathway$Mono_Macs$Pathway %in% Macro_pathway[,1],]
Macro_pathway$Celltype <- "Mono/Macs"
Neutro_pathway <- openxlsx::read.xlsx("./绘图/Neutro.xlsx", colNames = F)
Neutro_pathway <- Common_pathway$Neutrophil[Common_pathway$Neutrophil$Pathway %in% Neutro_pathway[,1],]
Neutro_pathway$Celltype <- "Neutrophil"

Celltype_pathway <- as.data.frame(rbind(CD4_pathway, CD8_pathway,
                                        DC_pathway, Macro_pathway, Neutro_pathway))
Celltype_pathway$Pathway <- paste0(Celltype_pathway$Celltype, "_",
                                   Celltype_pathway$Pathway)
# Celltype_pathway <- Celltype_pathway[!duplicated(Celltype_pathway$Pathway),]
Celltype_pathway$Pathway <- factor(Celltype_pathway$Pathway,
                                   levels = rev(Celltype_pathway$Pathway))
Celltype_pathway_list <- split.data.frame(Celltype_pathway, f = Celltype_pathway$Celltype)
Celltype_pathway_list_plot <- lapply(Celltype_pathway_list, function(X){
  ggplot(data = X, aes(x = log10P, y = Pathway, fill = Celltype)) +
    geom_bar(stat = "identity", width = 0.7) +
    ggtitle(unique(X$Celltype)) +
    labs(x = "Average -log10(P-value)", y = "") +
    scale_fill_manual(values = Celltype_color) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          plot.title = element_text(size = 12, colour = "black",
                                    face = "bold", hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank())
})

pdf("./绘图/10所有细胞_细胞类型通路.pdf",
    width = 11.5, height = 10)
patchwork::wrap_plots(Celltype_pathway_list_plot,
                      byrow = F, nrow = 5, ncol = 1,
                      heights = table(Celltype_pathway$Celltype) * 0.3)
dev.off()

saveRDS(All_seu, "第一次注释_4.Data2_Annotated_with_CD45_classification_renamed.rds")

table(All_seu$Celltype_rename)

All_seu$SubCelltype <- as.character(All_seu$Celltype_rename)

CD4_T <- readRDS("./注释/CD4+/CD4_T_annotated.rds")
table(CD4_T$SubCelltype)
CD8_T <- readRDS("./注释/CD8+/CD8_T_annotated.rds")
table(CD8_T$SubCelltype)

All_seu@meta.data[colnames(CD4_T), "SubCelltype"] <- CD4_T$SubCelltype
All_seu@meta.data[colnames(CD8_T), "SubCelltype"] <- CD8_T$SubCelltype

table(All_seu$SubCelltype)

Anti_tumor_cells <- c("NK", "CD4 IFN Response", "CD4 Naive/Memory",
                      "CD4 Teff", "CD4 Tem", "CD4 Texpanding",
                      "CD4 Tpex", "DPT", "CD8 IFN Response", "CD8 Naive/Tcm",
                      "CD8 Proliferating-Cdc20", "CD8 Proliferating-Mki67",
                      "CD8 Teff", "CD8 Tem", "CD8 Texpanding", "CD8 Tpex")
Pro_tumor_cells <- c("CD4 Treg", "CD4 Tex", "CD8 Tex")

All_seu_meta <- All_seu@meta.data
All_seu_meta <- All_seu_meta[which(!c(All_seu_meta$Celltype_rename %in% c("Melanocyte"))),]
All_seu_meta$Anti_Pro <- "Not"
All_seu_meta$Anti_Pro[All_seu_meta$SubCelltype %in% Anti_tumor_cells] <- "Anti"
All_seu_meta$Anti_Pro[All_seu_meta$SubCelltype %in% Pro_tumor_cells] <- "Pro"
temp <- as.data.frame.array(table(All_seu_meta$Sample_rename, All_seu_meta$Anti_Pro))
temp <- data.frame(t(temp), check.rows = F, check.names = F)
for (i in 1:ncol(temp)) {
  temp[,i] <- temp[,i] / sum(temp[,i])
}
temp2 <- data.frame(Sample = colnames(temp),
                    Ratio = 0)
for (i in 1:nrow(temp2)) {
  temp2[i,2] <- temp[1,temp2[i,1]] / temp[3,temp2[i,1]]
}
temp2$Sample <- factor(temp2$Sample,
                       levels = c("OT-I", "OT-I/XCL1",
                                  "OT-I/FX", "OT-I/FX + CD40 mAb"))
pdf("./绘图/11.Anti_pro_ratio.pdf", width = 3.5, height = 6)
ggplot(data = temp2, aes(x = Sample, y = Ratio, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_color) +
  theme_classic() +
  labs(x = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

temp <- as.data.frame.array(table(All_seu_meta$Sample_rename,
                                  All_seu_meta$Anti_Pro))
temp <- data.frame(t(temp), check.rows = F, check.names = F)
for (i in 1:ncol(temp)) {
  temp[,i] <- temp[,i] / sum(temp[,i])
}
temp <- as.data.frame(t(temp))
temp <- data.frame(Sample = rownames(temp),
                   temp)
temp$Sample <- factor(temp$Sample,
                      levels = temp$Sample)
pdf("./绘图/11.Anti.pdf", width = 3.5, height = 6)
ggplot(data = temp, aes(x = Sample, y = Anti, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_color) +
  theme_classic() +
  labs(x = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
pdf("./绘图/11.Pro.pdf", width = 3.5, height = 6)
ggplot(data = temp, aes(x = Sample, y = Pro, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_color) +
  theme_classic() +
  labs(x = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# CD8_T
{
  CD8_T <- readRDS("./注释/CD8+/CD8_T_annotated.rds")
  table(CD8_T$SubCelltype)
  CD8_T$SubCelltype <- factor(CD8_T$SubCelltype,
                              levels = c("CD8 Tex", "CD8 Proliferating-Mki67",
                                         "CD8 Proliferating-Cdc20", "CD8 Tpex",
                                         "CD8 IFN Response", "CD8 Tem", "CD8 Teff",
                                         "CD8 Texpanding", "CD8 Naive/Tcm"))
  CD8_T_colors <- ArchR::ArchRPalettes$circus
  CD8_T_colors <- CD8_T_colors[1:length(unique(CD8_T$SubCelltype))]
  names(CD8_T_colors) <- levels(CD8_T$SubCelltype)
  pdf("./绘图/12.CD8_T_亚群UMAP.pdf",
      height = 4.5, width = 6.7)
  DimPlot(CD8_T, group.by = "SubCelltype", cols = CD8_T_colors) +
    ggtitle("")
  dev.off()
  
  SubCelltype <- data.frame(Sub = CD8_T$SubCelltype,
                            Cluster = CD8_T$RNA_snn_res.0.5)
  SubCelltype <- SubCelltype[!duplicated(SubCelltype),]
  rownames(SubCelltype) <- SubCelltype$Sub
  SubCelltype <- SubCelltype[names(CD8_T_colors),]
  cluster_color <- CD8_T_colors
  names(cluster_color) <- SubCelltype$Cluster
  pdf("./绘图/12.CD8_T_亚群UMAP_cluster.pdf",
      height = 4.5, width = 4.8)
  DimPlot(CD8_T, group.by = "RNA_snn_res.0.5", label = T,
          cols = cluster_color) +
    ggtitle("")
  dev.off()
  
  Idents(CD8_T) <- CD8_T$SubCelltype
  pdf("./绘图/13.CD8_T_dotplot.pdf", height = 3, width = 11)
  DotPlot(CD8_T, features = c("Cd8a", "Cd8b1", "Lef1", "Sell", "Ccr7", "Tcf7",
                              "Xcl1", "Prf1", "Ifng", "Gzma", "Gzmb", "Gzmk",
                              "Itga4", "Cd7", "Ifit3", "Ifit3b", "Stmn1", "Tubb5",
                              "Top2a", "Mki67", "Cdc20", "Ccnb1", "Hist1h1b",
                              "Tuba1b", "Tox", "Lag3", "Pdcd1", "Havcr2", "Cd74", "H2-Ab1")) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10, colour = "black"),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank(),
          legend.position = "right") +
    scale_color_gradient2(low = "#1E90FF", high = "#EE2C2C") +
    labs(y = "", x = "")
  dev.off()
  
  
  genes <- c("Cd4", "Cd74", "Ifng", "Gzmk", "Ifit1", "Tubb5",
             "Cdc20", "Mki67", "Tox", "Tcf7", "Pdcd1", "Havcr2")
  feature_plot_list <- lapply(genes, function(x){
    FeaturePlot(CD8_T, features = c(x), order = T,
                min.cutoff = 0) +
      scale_color_gradientn(colours = c("gray", "#124A9B")) +
      theme_classic() +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            plot.title = element_text(size = 14, colour = "black",
                                      face = "bold", hjust = 0.5),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid = element_blank())
  })
  pdf("./绘图/14.CD8_T_细胞marker_featurePlot.pdf",
      width = 10, height = 18)
  patchwork::wrap_plots(feature_plot_list,
                        byrow = F, nrow = 6, ncol = 3)
  dev.off()
  
  CD8_T$Sample_rename <- All_seu@meta.data[colnames(CD8_T), "Sample_rename"]
  pdf("./绘图/15.CD8_T_细胞类型比例.pdf",
      width = 6, height = 5)
  percent_bar(sce = CD8_T, Ident = "SubCelltype",
              Group = "Sample_rename",
              fill_color = CD8_T_colors)
  dev.off()
  temp <- as.data.frame.array(table(CD8_T$SubCelltype, CD8_T$Sample_rename))
  for (i in 1:ncol(temp)) {
    temp[,i] <- temp[,i] / colSums(temp)[i]
  }
  temp <- data.frame(Celltype = rownames(temp),
                     temp)
  openxlsx::write.xlsx(temp,
                       "./绘图/15.CD8_T_细胞类型比例.xlsx")
  
  library(monocle)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  {
    data <- as.matrix(CD8_T@assays$RNA@counts)
    data <- as(data, 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = CD8_T@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)
    mycds <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily = negbinomial.size())
    
    mycds <- estimateSizeFactors(mycds)
    mycds <- estimateDispersions(mycds, cores = 4)
    disp_table <- dispersionTable(mycds)
    disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    mycds <- setOrderingFilter(mycds, VariableFeatures(CD8_T))
    plot_ordering_genes(mycds)
    mycds <- reduceDimension(mycds, max_components = 2,
                             reduction_method = "DDRTree")
    mycds <- orderCells(mycds)
    pdf("./绘图/16.CD8_T_monocle轨迹.pdf",
        width = 11, height = 9)
    plot_cell_trajectory(mycds, color_by = "SubCelltype",
                         cell_size = 1, show_branch_points = FALSE) +
      scale_color_manual(values = CD8_T_colors) +
      theme_classic() +
      facet_wrap(".~SubCelltype") +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            strip.text.x.top = element_text(size = 10, colour = "black"),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid = element_blank(),
            legend.position = "right") +
      guides(color = guide_legend(override.aes = list(size = 2)))
    dev.off()
  }
  
  
  CD8_T_DEGs <- FindAllMarkers(CD8_T, only.pos = T)
  CD8_T_DEGs %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  top20 <- as.data.frame(top20)
  top20 <- top20[!duplicated(top20$gene),]
  
  CD8_T_DEGs %>%
    group_by(cluster) %>%
    slice_head(n = 3) %>%
    ungroup() -> top3
  top3 <- as.data.frame(top3)
  top3 <- top3[!duplicated(top3$gene),]
  library(ComplexHeatmap)
  pdf("./绘图/17.CD8_T_SubCelltype_top20_差异基因_heatmap.pdf", height = 8, width = 6)
  Cell_group_heatmap(seu = CD8_T,
                     group.by = "SubCelltype",
                     group.color = CD8_T_colors,
                     group.label = "SubCelltype",
                     genes = unique(top20$gene),
                     limit = c(-2, 2),
                     label_genes = unique(top3$gene))
  dev.off()
  
  Idents(CD8_T) <- CD8_T$Sample_rename
  CD8_T_Group_DEGs <- FindAllMarkers(CD8_T, only.pos = T)
  CD8_T_Group_DEGs <- CD8_T_Group_DEGs[order(CD8_T_Group_DEGs$cluster,
                                             -CD8_T_Group_DEGs$avg_log2FC),]
  CD8_T_Group_DEGs_list <- split.data.frame(CD8_T_Group_DEGs,
                                            f = list(CD8_T_Group_DEGs$cluster))
  names(CD8_T_Group_DEGs_list) <- gsub(pattern = "/", replacement = "_",
                                       fixed = T, x = names(CD8_T_Group_DEGs_list))
  openxlsx::write.xlsx(CD8_T_Group_DEGs_list,
                       "./绘图/18.CD8_T_所有亚型_4个组别的高表达基因.xlsx")
  CD8_T_Group_DEGs <- c()
  for (i in c("OT-I/FX + CD40 mAb",
              "OT-I/FX", "OT-I/XCL1")) {
    temp_diff <- FindMarkers(CD8_T, only.pos = T,
                             ident.1 = i,
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[order(temp_diff$avg_log2FC,
                                 decreasing = T),]
    temp_diff$cluster <- i
    temp_diff$gene <- rownames(temp_diff)
    CD8_T_Group_DEGs <- c(CD8_T_Group_DEGs,
                          list(temp_diff))
  }
  names(CD8_T_Group_DEGs) <- c("OT-I_FX + CD40 mAb",
                               "OT-I_FX", "OT-I_XCL1")
  openxlsx::write.xlsx(CD8_T_Group_DEGs,
                       "./绘图/18.CD8_T_所有亚型_3个处理组_分别vs_Vector_高表达基因.xlsx")
  
  
  
  Cell_group_heatmap2 <- function(seu,
                                  group.by = "Celltype_rename",
                                  group.color = Celltype_color,
                                  group.label = "Celltype",
                                  genes = unique(top5$gene),
                                  label_genes = unique(top5$gene)[unique(top5$gene) %in% unlist(lapply(Celltype_markers_list, FUN = as.character))],
                                  limit = c(-2, 2),
                                  flip = F,
                                  show_row_names = T,
                                  show_column_names = T,
                                  cluster_columns = T,
                                  cluster_rows = T,
                                  row_names_side = "left",
                                  column_names_side = "bottom",
                                  column_names_rot = 90,
                                  row_names_rot = 0,
                                  body_width = unit(8, "cm"),
                                  body_height = unit(8, "cm"),
                                  cols = c("#174C9E", "#E8402B"),
                                  genes_order = NULL
  ){
    exp <- seu@assays$RNA@data[genes,]
    exp <- t(as.matrix(exp))
    exp <- aggregate.data.frame(exp,
                                by = list(seu@meta.data[,group.by]),
                                FUN = mean)
    rownames(exp) <- exp$Group.1
    exp <- exp[,-1]
    exp <- exp[,colSums(exp) != 0]
    for (i in 1:ncol(exp)) {
      exp[,i] <- (exp[,i] - mean(exp[,i])) / sd(exp[,i])
    }
    exp <- as.matrix(t(exp))
    if (!is.null(genes_order)) {
      exp <- exp[genes_order,]
    }
    exp[exp > limit[2]] <- limit[2]
    exp[exp < limit[1]] <- limit[1]
    neg_length <- 0 - min(exp)
    pos_length <- max(exp) - 0
    sum <- neg_length + pos_length
    neg_length <- 100 * (neg_length / sum)
    neg_length <- round(neg_length)
    pos_length <- 100 * (pos_length / sum)
    pos_length <- round(pos_length)
    neg_color <- colorRampPalette(c(cols[1], "white"))(neg_length)
    pos_color <- colorRampPalette(c("white", cols[2]))(pos_length)
    color_bar <- c(neg_color, pos_color)
    color_map <- circlize::colorRamp2(seq(from = min(exp), to = max(exp), length = length(color_bar)),
                                      color_bar)
    if (flip == T) {
      exp <- as.matrix(data.frame(t(exp),
                                  check.rows = F,
                                  check.names = F))
    }
    if (!is.null(genes_order)) {
      cluster_rows = FALSE
      cluster_columns = FALSE
    }
    Heatmap(exp, col = color_map, color_space = "RGB",
            cluster_rows = cluster_rows, cluster_columns = cluster_columns,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            # top_annotation = col_anno,
            show_row_names = show_row_names,
            show_column_names = show_column_names,
            column_dend_height = unit(0, "mm"),
            row_dend_width = unit(0, "mm"),
            row_names_side = row_names_side,
            column_names_side = column_names_side,
            show_row_dend = F, show_column_dend = F,
            column_names_rot = column_names_rot,
            row_names_rot = row_names_rot,
            width = body_width, height = body_height,
            heatmap_legend_param = list(title = "Scaled\nExpression"
                                        # at = c(round(min(exp),1), 0, max(exp)),
                                        # labels = c(round(min(exp),1), 0, max(exp))))
            ))
    
  }
  pdf("./绘图/18.CD8_T_co-stimu_co-inhibi_热图.pdf",
      height = 3.5, width = 8)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      # genes = c("Tnfrsf10b", "Cd244a", "Tnfsf10", "Tnfrsf18", "Havcr2", "Ctla4", "Pdcd1", "Tnfrsf9", "Icos",
                      #           "Cd160", "Tnfrsf4", "Lag3", "Cd27", "Cd28", "Tigit"),
                      genes = c("Cd160", "Icos", "Cd27", "Cd28", "Fas", "Tnfrsf4",
                                "Tnfrsf9", "Tnfrsf10b", "Tnfrsf18",
                                "Pdcd1", "Ctla4", "Lag3", "Havcr2", "Cd244a", "Tigit"),
                      limit = c(-2, 2),
                      flip = T,
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = T,
                      show_column_names = T,
                      cluster_columns = T,
                      cluster_rows = F,
                      column_names_rot = 45,
                      row_names_rot = 0,
                      body_width = unit(15 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"),
                      cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/18.CD8_T_TCR_signaling_热图.pdf",
      height = 3.5, width = 8)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      # genes = c("Nr4a1", "Nfatc1", "Nfkbia", "Nfatc2",
                      #           "Nfkbib", "Nfkb1", "Nfatc3", "Jun", "Fos"),
                      genes = c("Ifng", "Nr4a1", "Nr4a2", "Nfatc1", "Nfatc2",
                                "Nfkbib", "Nfkbia", "Nfatc4", "Nfat5", "Nfkb1",
                                "Nfatc3", "Jun", "Fos"),
                      genes_order = c("Nfatc2", "Nr4a1",
                                      "Nfat5", "Nfkbia", "Nr4a2", "Ifng", "Nfkbib",
                                      "Nfatc1", "Jun", "Fos", "Nfatc3", "Nfkb1"),
                      limit = c(-2, 2),
                      flip = T,
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = T,
                      show_column_names = T,
                      cluster_columns = T,
                      cluster_rows = F,
                      column_names_rot = 45,
                      row_names_rot = 0,
                      body_width = unit(13 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"),
                      cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  pdf("./绘图/18.CD8_T_Activation_Differentiation_热图.pdf",
      height = 3.5, width = 13)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      # genes = c("Cd44", "Tbx21", "Cd69", "Klrg1", "Runx3", "Cd101",
                      #           "Entpd1", "Il2rb", "Il2ra", "Klrc1", "Klrd1", "Mki67",
                      #           "Id2", "Cd38", "Eomes", "Lef1", "Bcl2", "Bach2", "Tcf7",
                      #           "Id3", "Tox", "Slamf6", "Sell", "Il7r", "Bcl6", "Runx1"),
                      genes = c("Id2", "Klrg1", "Klrc1", "Klrd1", "Nkg7", "Klf2",
                                "Foxp1", "Tbx21", "Batf", "Tox", "Il2ra", "Il2rb", "Il21r",
                                "Stat5a", "Stat5b", "Entpd1", "Cd101", "Cd38", "Eomes",
                                "Mki67", "Sell", "Cd69", "Bcl6", "Id3", "Il7r", "Lef1", "Bcl2",
                                "Bach2", "Slamf6", "Runx1", "Runx3"),
                      genes_order = c("Tbx21", "Klrc1",
                                      "Klrd1", "Tox",
                                      "Stat5a", "Il2ra", "Entpd1", "Cd101", "Cd38", "Il21r", "Il2rb",
                                      "Eomes", "Batf", "Mki67", "Nkg7", "Sell", "Cd69", "Bcl6", "Id3",
                                      "Klf2", "Slamf6", "Runx3", "Il7r", "Bcl2", "Bach2", "Foxp1",
                                      "Lef1", "Runx1", "Stat5b", "Id2", "Klrg1"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(31 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/18.CD8_T_IFN_signaling_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      # genes = c("Ifitm3", "Ifitm1", "Ifitm2", "Ifitm5", "Stat1"),
                      genes = c("Ifng", "Jak1", "Jak2", "Stat1", "Stat3", "Ifitm1",
                                "Isg15", "Ifitm3", "Ifitm2", "Ifitm5"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(10 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/18.CD8_T_Effector_molecules_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      genes = c("Ifng", "Gzmb", "Prf1", "Tnf", "Gzmd",
                                "Gzme", "Gzmk", "Il2", "Gzma", "Lamp1"),
                      genes_order = c("Il2", "Gzma", "Tnf", "Gzmb", "Ifng",
                                      "Gzmk", "Prf1", "Gzmd", "Lamp1", "Gzme"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(10 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/18.CD8_T_Chemokine_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD8_T,
                      group.by = "Sample_rename",
                      group.color = CD8_T_colors,
                      group.label = "Sample",
                      # genes = c("Cx3cr1", "Ccr5", "Ccr2", "Cxcr3",
                      #           "Cxcr6", "Cxcr5", "Ccr7"),
                      genes = c("Ccl3", "Xcl1", "Ccr7",
                                "Ccr5", "Ccr2", "Cxcr3", "Cxcr6", "Cxcr5",
                                "Cx3cr1"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(11 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  library(scRNAtoolVis)
  CD8_T_sub_diff <- c()
  table(CD8_T$SubCelltype, CD8_T$Sample_rename)
  CD8_SubCelltype <- unique(CD8_T$SubCelltype)
  for (c in CD8_SubCelltype) {
    temp <- subset(CD8_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
    temp_diff$cluster <- c
    temp_diff$gene <- rownames(temp_diff)
    CD8_T_sub_diff <- as.data.frame(rbind(CD8_T_sub_diff,
                                          temp_diff))
  }
  jjVolcano <- function(diffData = NULL, myMarkers = NULL, order.by = c("avg_log2FC"), 
                         log2FC.cutoff = 0.25, pvalue.cutoff = 0.05, adjustP.cutoff = 0.01, 
                         topGeneN = 5, col.type = "updown", back.col = "grey93", pSize = 0.75, 
                         aesCol = c("#0099CC", "#CC3333"), legend.position = c(0.7, 0.9), base_size = 10,
                         tile.col = jjAnno::useMyCol("paired", n = 9), 
                         cluster.order = NULL, polar = FALSE, expand = c(-1, 1), flip = FALSE, ...) 
  {
    diff.marker <- diffData %>% dplyr::filter(abs(avg_log2FC) >= 
                                                log2FC.cutoff & p_val < pvalue.cutoff)
    diff.marker <- diff.marker %>% dplyr::mutate(type = ifelse(avg_log2FC >= 
                                                                 log2FC.cutoff, "sigUp", "sigDown")) %>% dplyr::mutate(type2 = ifelse(p_val_adj < 
                                                                                                                                        adjustP.cutoff, paste("adjust Pvalue < ", adjustP.cutoff, 
                                                                                                                                                              sep = ""), paste("adjust Pvalue >= ", adjustP.cutoff, 
                                                                                                                                                                               sep = "")))
    if (!is.null(cluster.order)) {
      diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
    }
    back.data <- purrr::map_df(unique(diff.marker$cluster), function(x) {
      tmp <- diff.marker %>% dplyr::filter(cluster == x)
      new.tmp <- data.frame(cluster = x, min = min(tmp$avg_log2FC) - 
                              0.2, max = max(tmp$avg_log2FC) + 0.2)
      return(new.tmp)
    })
    top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
    top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = topGeneN, 
                                                          order_by = get(order.by))
    top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = topGeneN, 
                                                          order_by = get(order.by))
    top.marker <- rbind(top.marker.max, top.marker.min)
    if (!is.null(myMarkers)) {
      top.marker <- diff.marker %>% dplyr::filter(gene %in% 
                                                    myMarkers)
    }
    else {
      top.marker <- top.marker
    }
    p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, 
                                                    y = avg_log2FC)) + ggplot2::geom_col(data = back.data, 
                                                                                         ggplot2::aes(x = cluster, y = min), fill = back.col) + 
      ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, 
                                                       y = max), fill = back.col)
    if (col.type == "updown") {
      p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), 
                                      size = pSize) + ggplot2::scale_color_manual(values = c(sigDown = aesCol[1], 
                                                                                             sigUp = aesCol[2]))
    }
    else if (col.type == "adjustP") {
      p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2), 
                                      size = pSize) + ggplot2::scale_color_manual(values = c(aesCol[2], 
                                                                                             aesCol[1]))
    }
    p3 <- p2 + ggplot2::scale_y_continuous(n.breaks = 6) + ggplot2::theme_classic(base_size = base_size) + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                     legend.position = legend.position, legend.title = ggplot2::element_blank(), 
                     legend.background = ggplot2::element_blank()) + ggplot2::xlab("Clusters") + 
      ggplot2::ylab("Average log2FoldChange") + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
    p4 <- p3 + ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, 
                                               fill = cluster), color = "black", height = log2FC.cutoff * 
                                    2, alpha = 0.3, show.legend = FALSE) + ggplot2::scale_fill_manual(values = tile.col) + 
      ggrepel::geom_text_repel(data = top.marker, ggplot2::aes(x = cluster, 
                                                               y = avg_log2FC, label = gene), max.overlaps = 50, 
                               ...)
    if (polar == TRUE) {
      p5 <- p4 + geomtextpath::geom_textpath(ggplot2::aes(x = cluster, 
                                                          y = 0, label = cluster)) + ggplot2::scale_y_continuous(n.breaks = 6, 
                                                                                                                 expand = ggplot2::expansion(mult = expand)) + ggplot2::theme_void(base_size = base_size) + 
        ggplot2::theme(legend.position = legend.position, 
                       legend.title = ggplot2::element_blank()) + ggplot2::coord_polar(clip = "off", 
                                                                                       theta = "x")
    }
    else {
      if (flip == TRUE) {
        p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
          ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                           y = 0, label = cluster), color = "black", size = 3.3) + ggplot2::theme(axis.line.y = ggplot2::element_blank(), 
                                                                                                                axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
          ggplot2::coord_flip()
      }
      else {
        p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
          ggplot2::geom_text(ggplot2::aes(x = cluster, 
                                          y = 0, label = cluster), color = "black", size = 3.3) + ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                                                                                                               axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
      }
    }
    return(p5)
  }
  
  CD8_T_clusters <- levels(CD8_T$SubCelltype)
  CD8_T_clusters <- data.frame(original = CD8_T_clusters,
                               rename = CD8_T_clusters)
  CD8_T_clusters[2, "rename"] <- "CD8 Prolif.-\nMki67"
  CD8_T_clusters[3, "rename"] <- "CD8 Prolif.-\nCdc20"
  CD8_T_clusters[5, "rename"] <- "CD8 IFN\nResponse"
  rownames(CD8_T_clusters) <- CD8_T_clusters$original
  CD8_T_sub_diff_rename <- CD8_T_sub_diff
  CD8_T_sub_diff_rename$cluster <- CD8_T_clusters[CD8_T_sub_diff_rename$cluster,
                                                  2]
  CD8_T_colors_rename <- CD8_T_colors
  names(CD8_T_colors_rename) <- CD8_T_clusters$rename
  pdf("./绘图/19.CD8_T_不同亚群_最后一组_vs_第一组.pdf",
      width = 10, height = 5.5)
  jjVolcano(diffData = CD8_T_sub_diff_rename,
            log2FC.cutoff = 0.25, 
            base_size = 14, pSize = 0.75, #设置点的大小
            fontface = 'italic', #设置字体形式
            aesCol = c("#0099CC","#CC3333"), #设置点的颜色
            tile.col = CD8_T_colors_rename, #设置cluster的颜色
            adjustP.cutoff = 1,
            cluster.order = CD8_T_clusters$rename,
            col.type = "updown",
            topGeneN = 5 
  ) + labs(x = "")
  dev.off()
  
  Vio_plot <- function(diff_gene = Epithelial_diff, logfc_thr = 0.25,
                       cols = c("#FF4500", "#4169E1"),
                       p_thr = 0.05, group = c("TACSTD2+", "TACSTD2-"),
                       title = "TACSTD2+ VS TACSTD2-",
                       genes = NULL) {
    diff_gene <- diff_gene[order(diff_gene$avg_log2FC, decreasing=T),]
    DOWN_row <- intersect(which(diff_gene$p_val < p_thr),
                          which(diff_gene$avg_log2FC < -logfc_thr))
    UP_row <- intersect(which(diff_gene$p_val < p_thr),
                        which(diff_gene$avg_log2FC > logfc_thr))
    DOWN <- diff_gene[DOWN_row,]
    UP <- diff_gene[UP_row,]
    DOWN <- DOWN[order(DOWN$avg_log2FC,decreasing = F),]
    UP <- UP[order(UP$avg_log2FC,decreasing = T),]
    color <- rep("gray", length(diff_gene$avg_log2FC))
    cex <- rep(1, length(diff_gene$avg_log2FC))
    cex[DOWN_row] <- 1.5
    cex[UP_row] <- 1.5
    cex <- as.character(cex)
    color[DOWN_row] <- cols[2]
    color[UP_row] <- cols[1]
    color <- factor(color, levels = c(cols[1], cols[2], "gray"))
    color_value <- c(cols[1], cols[2], "gray")
    names(color_value) <- color_value
    if (!is.null(genes)) {
      genes <- genes[genes %in% rownames(diff_gene)]
      sigdiff <- diff_gene[genes,]
    } else {
      sigdiff <- rbind(DOWN[1:5,], UP[1:5,])
    }
    vio_plot <- ggplot(data = diff_gene, aes(x = avg_log2FC, y = -log10(p_val),
                                             size = cex, colour = color,
                                             alpha = "1")) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            panel.grid = element_blank(),
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 10, color = "black")) +
      labs(x = "log2 FoldChange", y = "log10 P value") +
      geom_point() +
      geom_hline(aes(yintercept=-log10(p_thr)), colour = "black", linetype="dashed") +
      geom_vline(aes(xintercept= -logfc_thr), colour="black", linetype="dashed") +
      geom_vline(aes(xintercept= logfc_thr), colour="black", linetype="dashed") +
      scale_color_manual(name = "",
                         values = color_value,
                         labels = c(paste0("\n",nrow(UP)," UP genes for\n",group[1]),
                                    paste0("\n",nrow(DOWN)," UP genes for\n",group[2]),
                                    "Not significant")) +
      scale_size_manual(values = c('1' = 1, "1.5" = 1.5)) +
      scale_alpha_manual(values = c("0.3" = 0.3, "1" = 1)) +
      guides(size = "none", alpha = "none",
             color = guide_legend(override.aes = list(size = 2))) +
      geom_text_repel(data = sigdiff, aes(x = avg_log2FC, 
                                          y = -log10(p_val), 
                                          label = gene,
                                          alpha = rep("1",nrow(sigdiff))),
                      size = 3, fontface = "italic",
                      box.padding = unit(0.8, "lines"),
                      point.padding = unit(0, "lines"), 
                      min.segment.length = 0,
                      segment.color = "black",
                      colour="#000000",
                      show.legend = FALSE,
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) 
    vio_plot
  }
  Idents(CD8_T) <- CD8_T$Sample_rename
  temp_diff <- FindMarkers(CD8_T, logfc.threshold = 0,
                           ident.1 = "OT-I/FX + CD40 mAb",
                           ident.2 = "OT-I")
  temp_diff$gene <- rownames(temp_diff)
  p <- Vio_plot(diff_gene = temp_diff, logfc_thr = 0.25,
                cols = c("#CC3333", "#0099CC"),
                p_thr = 0.05,
                group = c("OT-I/FX + CD40 mAb", "OT-I"),
                title = "CD8 T",
                genes = c("Itgb1", "Fos", "Klf2", "Gzma", "S1pr1",
                          "Tcf7", "Itga4", "Tsc22d3", "Il7r",
                          "Plac8", "Gzmb", "Ifi27l2a", "Ly6a",
                          "Tubb5", "Cxcr6", "Gzmk", "Ccl3", "Cd160",
                          "Mki67", "Tnfrsf9", "Ifng", "Pdcd1"))
  pdf(paste0("./绘图/19.CD8_T_所有亚型_最后一组_vs_第一组_火山图.pdf"),
      width = 6.2, height = 5)
  print(p)
  dev.off()
  
  temp_diff$gene <- rownames(temp_diff)
  temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
  temp_diff <- temp_diff[abs(temp_diff$avg_log2FC) > 0.25,]
  temp_diff <- temp_diff[order(temp_diff$avg_log2FC, decreasing = T),]
  openxlsx::write.xlsx(temp_diff,
                       "./绘图/19.CD8_T_所有亚型_最后一组_vs_第一组_差异表达基因.xlsx")
  
  library(ggrepel)
  for (c in CD8_SubCelltype) {
    temp <- subset(CD8_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp, logfc.threshold = 0,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff$gene <- rownames(temp_diff)
    p <- Vio_plot(diff_gene = temp_diff, logfc_thr = 0.25,
             cols = c("#CC3333", "#0099CC"),
             p_thr = 0.05,
             group = c("OT-I/FX + CD40 mAb", "OT-I"),
             title = c)
    pdf(paste0("./绘图/19.CD8_T_", gsub(pattern = "/",
                                      replacement = "_",
                                      fixed = T, x = c), "_最后一组_vs_第一组_火山图.pdf"),
        width = 6.2, height = 5)
    print(p)
    dev.off()
    
  }
  
  library(clusterProfiler)
  library(createKEGGdb)
  # createKEGGdb::create_kegg_db("mmu")
  # install.packages("./KEGG.db_1.0.tar.gz",type="source")
  library(KEGG.db)
  KEGG_DATA <- clusterProfiler::download_KEGG("mmu")
  KEGG_df <- KEGG_DATA[["KEGGPATHID2EXTID"]]
  KEGG_df_2 <- KEGG_DATA[["KEGGPATHID2NAME"]]
  rownames(KEGG_df_2) <- KEGG_df_2$from
  KEGG_df$from <- KEGG_df_2[KEGG_df$from,2]
  KEGG_df$from <- unlist(lapply(KEGG_df$from, function(x){
    unlist(strsplit(x, " - Mus musculus (house mouse)",
                    fixed = T))[1]
  }))
  library(org.Mm.eg.db)
  id <- bitr(geneID = unique(KEGG_df$to),
       fromType = "ENTREZID",
       toType = "SYMBOL",
       OrgDb = org.Mm.eg.db)
  rownames(id) <- as.character(id$ENTREZID)
  KEGG_df$to <- id[KEGG_df$to, "SYMBOL"]
  colnames(KEGG_df) <- c("term", "gene")
  KEGG_df$term <- paste0("KEGG ", KEGG_df$term)
  GeneSet1 <- read.gmt("./绘图/m5.go.bp.v2023.2.Mm.symbols.gmt")
  GeneSet1$term <- as.character(GeneSet1$term)
  GeneSet1$term <- unlist(lapply(GeneSet1$term, function(x){
    temp <- unlist(strsplit(x, "_"))
    temp <- temp[-1]
    temp <- unlist(lapply(temp, function(y){
      y <- tolower(y)
      y <- Hmisc::capitalize(y)
    }))
    paste0(temp, collapse = " ")
  }))
  GeneSet1$term <- paste0("GO ", GeneSet1$term)
  
  GeneSet2 <- read.gmt("./绘图/mh.all.v2023.2.Mm.symbols.gmt")
  GeneSet2$term <- as.character(GeneSet2$term)
  GeneSet2$term <- unlist(lapply(GeneSet2$term, function(x){
    temp <- unlist(strsplit(x, "_"))
    temp <- temp[-1]
    temp <- unlist(lapply(temp, function(y){
      y <- tolower(y)
      y <- Hmisc::capitalize(y)
    }))
    paste0(temp, collapse = " ")
  }))
  GeneSet2$term <- paste0("HALLMARK ", GeneSet2$term)

  GeneSet <- as.data.frame(rbind(GeneSet1, GeneSet2))
  GeneSet <- as.data.frame(rbind(GeneSet,
                                 KEGG_df))
  GeneSet_list <- split.data.frame(GeneSet,
                                   f = factor(GeneSet$term))
  GeneSet_list_length <- unlist(lapply(GeneSet_list, function(x){
    nrow(x)
  }))
  GeneSet_list <- GeneSet_list[GeneSet_list_length >= 5]
  
  TERM2GENE <- rbindlist(GeneSet_list)
  library(enrichplot)
  Idents(CD8_T) <- CD8_T$Sample_rename
  CD8_T_diff <- FindMarkers(CD8_T,
                            ident.1 = "OT-I/FX + CD40 mAb",
                            ident.2 = "OT-I",
                            logfc.threshold = 0)
  CD8_T_diff <- CD8_T_diff[order(CD8_T_diff$avg_log2FC,
                                 decreasing = T),]
  diff_list <- CD8_T_diff$avg_log2FC
  names(diff_list) <- rownames(CD8_T_diff)
  set.seed(123)
  CD8T_gsea <- GSEA(diff_list, TERM2GENE = TERM2GENE,
                     verbose = T, eps = 0, pvalueCutoff = 1)
  CD8T_gsea_result <- CD8T_gsea@result
  CD8T_gsea_result <- CD8T_gsea_result[CD8T_gsea_result$pvalue < 0.05,]
  CD8T_gsea_result <- CD8T_gsea_result[order(CD8T_gsea_result$NES,
                                             decreasing = T),]
  CD8T_gsea@result <- CD8T_gsea_result
  saveRDS(CD8T_gsea, "./绘图/CD8T不分亚型_第四组_vs_第一组_gsea_pathway.rds")
  openxlsx::write.xlsx(CD8T_gsea_result,
                       "./绘图/CD8T不分亚型_第四组_vs_第一组_gsea_result.xlsx")
  
  temp <- grep(pattern = "myc|e2f|oxidative|mtorc1|fatty|peroxisome|pyruvate|stat5|glycolysis|Il2",
       ignore.case = TRUE, x = CD8T_gsea@result$Description)
  temp <- CD8T_gsea@result[temp,]
  temp <- temp[-c(3,7,12),]
  library(GseaVis)
  pdf("./绘图/21.CD8_T_不分亚型_第四组_vs_第一组_gsea_plot.pdf",
      width = 7, height = 4)
  gseaNb(object = CD8T_gsea,
         geneSetID = temp$ID,
         addPval = T, subPlot = 2,
         # curveCol = scales::hue_pal()(7),
         curveCol = c("#89C75F","#F37B7D","#D24B27", "#FF8C00",
                      "#3BBCA8","#6E4B9E","#0C727C","#D8A767",
                      "#272E6A"),
         htHeight = 0.6)#,
         # subPlot = 2, addPval = TRUE,
         # pvalX = 0.95, pvalY = 0.8,
         # pvalSize = 4, pDigit = 1)
  dev.off()
  
  for (c in CD8_SubCelltype) {
    temp <- subset(CD8_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp, logfc.threshold = 0,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[order(temp_diff$avg_log2FC,
                                 decreasing = T),]
    temp_diff_list <- temp_diff$avg_log2FC
    names(temp_diff_list) <- rownames(temp_diff)
    set.seed(123)
    temp_gsea <- GSEA(temp_diff_list, TERM2GENE = TERM2GENE,
                      verbose = T, eps = 0, pvalueCutoff = 1)
    temp_gsea_result <- temp_gsea@result
    temp_gsea_result <- temp_gsea_result[temp_gsea_result$pvalue < 0.05,]
    temp_gsea_result <- temp_gsea_result[order(temp_gsea_result$NES,
                                               decreasing = T),]
    temp_gsea@result <- temp_gsea_result
    openxlsx::write.xlsx(temp_gsea_result,
                         paste0("./绘图/",
                                "CD8T_", gsub(pattern = "/",
                                              replacement = "_",
                                              x = c), "_第四组_vs_第一组_gsea_result.xlsx"))
    saveRDS(temp_gsea,
            paste0("./绘图/",
                   "CD8T_", gsub(pattern = "/",
                                 replacement = "_",
                                 x = c), "_第四组_vs_第一组_gsea_pathway.rds"))

  }
  
  
  CD8_T_Ratio <- openxlsx::read.xlsx("./绘图/CD8_T_Ratio.xlsx",
                                     check.names = F)
  colnames(CD8_T_Ratio) <- c("Celltype", "OT-I", "OT-I/XCL1", "OT-I/FX",
                             "OT-I/FX + CD40 mAb")
  CD8_T_Ratio <- reshape::melt.data.frame(CD8_T_Ratio,
                                          id.vars = "Celltype")
  unique(CD8_T_Ratio$Celltype)
  pdf("./绘图/22.CD8_T_non exhasution_exhasution.pdf", width = 3.5, height = 6)
  ggplot(data = CD8_T_Ratio[CD8_T_Ratio$Celltype == "non exhasution/exhasution",],
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "non exhasution/exhasution") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  pdf("./绘图/22.CD8_T_exhasution_Tpex.pdf", width = 3.5, height = 6)
  ggplot(data = CD8_T_Ratio[CD8_T_Ratio$Celltype == "exhasution/Tpex",],
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "exhasution/Tpex") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  pdf("./绘图/22.CD8_T_exhasution_eff.pdf", width = 3.5, height = 6)
  ggplot(data = CD8_T_Ratio[CD8_T_Ratio$Celltype == "exhasution/eff",],
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "exhasution/eff") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  pdf("./绘图/22.CD8_T_echaustion_cm_em.pdf", width = 3.5, height = 6)
  ggplot(data = CD8_T_Ratio[CD8_T_Ratio$Celltype == "echaustion/cm em",],
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "echaustion/cm em") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
}


# CD4_T
{
  CD4_T <- readRDS("./注释/CD4+/CD4_T_annotated.rds")
  table(CD4_T$SubCelltype)
  Idents(CD4_T) <- CD4_T$SubCelltype
  CD4_T <- RenameIdents(CD4_T,
                        "CD4 Tex" = "CD4 Tex",
                        "CD4 Tpex" = "CD4 Tpex",
                        "CD4 IFN Response" = "CD4 IFN response",
                        "CD4 Tem" = "CD4 Tem",
                        "CD4 Teff" = "CD4 Teff",
                        "CD4 Texpanding" = "CD4 T expanding",
                        "CD4 Treg" = "Treg",
                        "CD4 Naive/Memory" = "CD4 Naïve/Tcm"
                        )
  CD4_T$SubCelltype <- factor(as.character(Idents(CD4_T)),
                              levels = c("CD4 Tex", "CD4 Tpex", "CD4 IFN response", "CD4 Tem",
                                         "CD4 Teff", "CD4 T expanding", "Treg", "CD4 Naïve/Tcm"))
  # CD4_T$SubCelltype <- factor(CD4_T$SubCelltype,
  #                             levels = c("CD4 Tex", "CD4 Tpex", "CD4 IFN Response", "CD4 Tem",
  #                                        "CD4 Teff", "CD4 Texpanding", "CD4 Treg", "CD4 Naive/Memory"))
  CD4_T_colors <- ArchR::ArchRPalettes$paired
  CD4_T_colors <- CD4_T_colors[1:length(unique(CD4_T$SubCelltype))]
  names(CD4_T_colors) <- levels(CD4_T$SubCelltype)
  pdf("./绘图/23.CD4_T_亚群UMAP.pdf",
      height = 4.5, width = 6.2)
  DimPlot(CD4_T, group.by = "SubCelltype", cols = CD4_T_colors) +
    ggtitle("")
  dev.off()
  
  SubCelltype <- data.frame(Sub = CD4_T$SubCelltype,
                            Cluster = CD4_T$RNA_snn_res.0.5)
  SubCelltype <- SubCelltype[!duplicated(SubCelltype),]
  rownames(SubCelltype) <- SubCelltype$Sub
  SubCelltype <- SubCelltype[names(CD4_T_colors),]
  cluster_color <- CD4_T_colors
  names(cluster_color) <- SubCelltype$Cluster
  pdf("./绘图/24.CD4_T_亚群UMAP_cluster.pdf",
      height = 4.5, width = 4.8)
  DimPlot(CD4_T, group.by = "RNA_snn_res.0.5", label = T,
          cols = cluster_color) +
    ggtitle("")
  dev.off()
  
  Idents(CD4_T) <- CD4_T$SubCelltype
  pdf("./绘图/25.CD4_T_dotplot.pdf", height = 3, width = 10)
  DotPlot(CD4_T, features = c("Cd4", "Lef1", "Sell", "Ccr7", "Tcf7", "Slamf6", "Xcl1", "Foxp3",
                              "Tnfrsf4", "Ctla4", "Gzma", "Gzmb", "Prf1", "Ifng", "Tnfsf11", "Il7r",
                              "Fos", "Cd40lg", "Ifit1", "Isg15", "Ccl5", "Tox", "Pdcd1",
                              "Lag3", "Havcr2", "Cd74", "H2-Ab1")) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10, colour = "black"),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank(),
          legend.position = "right") +
    scale_color_gradient2(low = "#1E90FF", high = "#EE2C2C") +
    labs(y = "", x = "")
  dev.off()
  
  
  genes <- c("Cd8a", "Tcf7", "Slamf6", "Foxp3", "Cd74", "Ifng",
             "Il7r", "Ifit3", "Xcl1", "Tox", "Pdcd1", "Havcr2")
  feature_plot_list <- lapply(genes, function(x){
    FeaturePlot(CD4_T, features = c(x), order = T,
                min.cutoff = 0) +
      scale_color_gradientn(colours = c("gray", "#124A9B")) +
      theme_classic() +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            plot.title = element_text(size = 14, colour = "black",
                                      face = "bold", hjust = 0.5),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid = element_blank())
  })
  pdf("./绘图/26.CD4_T_细胞marker_featurePlot.pdf",
      width = 10, height = 18)
  patchwork::wrap_plots(feature_plot_list,
                        byrow = F, nrow = 6, ncol = 3)
  dev.off()
  
  CD4_T$Sample_rename <- All_seu@meta.data[colnames(CD4_T), "Sample_rename"]
  pdf("./绘图/27.CD4_T_细胞类型比例.pdf",
      width = 6, height = 5)
  percent_bar(sce = CD4_T, Ident = "SubCelltype",
              Group = "Sample_rename",
              fill_color = CD4_T_colors)
  dev.off()
  temp <- as.data.frame.array(table(CD4_T$SubCelltype, CD4_T$Sample_rename))
  for (i in 1:ncol(temp)) {
    temp[,i] <- temp[,i] / colSums(temp)[i]
  }
  temp <- data.frame(Celltype = rownames(temp),
                     temp)
  openxlsx::write.xlsx(temp,
                       "./绘图/27.CD4_T_细胞类型比例.xlsx")
  
  library(monocle)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  {
    data <- as.matrix(CD4_T@assays$RNA@counts)
    data <- as(data, 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = CD4_T@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)
    mycds <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily = negbinomial.size())
    
    mycds <- estimateSizeFactors(mycds)
    mycds <- estimateDispersions(mycds, cores = 4)
    disp_table <- dispersionTable(mycds)
    disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    mycds <- setOrderingFilter(mycds, VariableFeatures(CD4_T))
    plot_ordering_genes(mycds)
    mycds <- reduceDimension(mycds, max_components = 2,
                             reduction_method = "DDRTree")
    mycds <- orderCells(mycds)
    pdf("./绘图/28.CD4_T_monocle轨迹.pdf",
        width = 13, height = 6)
    plot_cell_trajectory(mycds, color_by = "SubCelltype",
                         cell_size = 1, show_branch_points = FALSE) +
      scale_color_manual(values = CD4_T_colors) +
      theme_classic() +
      facet_wrap(".~SubCelltype", ncol = 4) +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            strip.text.x.top = element_text(size = 10, colour = "black"),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid = element_blank(),
            legend.position = "right") +
      guides(color = guide_legend(override.aes = list(size = 2)))
    dev.off()
    
    mycds <- orderCells(mycds, root_state = "3")
    pdf("./绘图/28.CD4_T_monocle轨迹_pseudotime.pdf",
        width = 6, height = 4.5)
    plot_cell_trajectory(mycds, color_by = "Pseudotime",
                         show_branch_points = FALSE) +
      theme_classic() +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            strip.text.x.top = element_text(size = 10, colour = "black"),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid = element_blank(),
            legend.position = "right")
    dev.off()
  }
  
  
  Cell_group_heatmap <- function(seu,
           group.by = "Celltype_rename",
           group.color = Celltype_color,
           group.label = "Celltype",
           genes = unique(top5$gene),
           label_genes = unique(top5$gene)[unique(top5$gene) %in% unlist(lapply(Celltype_markers_list, FUN = as.character))],
           limit = c(-2, 2)){
    exp <- seu@assays$RNA@data[genes,]
    exp <- t(as.matrix(exp))
    exp <- aggregate.data.frame(exp,
                                by = list(seu@meta.data[,group.by]),
                                FUN = mean)
    rownames(exp) <- exp$Group.1
    exp <- exp[,-1]
    for (i in 1:ncol(exp)) {
      exp[,i] <- (exp[,i] - mean(exp[,i])) / sd(exp[,i])
    }
    exp <- as.matrix(t(exp))
    # exp <- exp[genes_order,]
    exp[exp > limit[2]] <- limit[2]
    exp[exp < limit[1]] <- limit[1]
    neg_length <- 0 - min(exp)
    pos_length <- max(exp) - 0
    sum <- neg_length + pos_length
    neg_length <- 100 * (neg_length / sum)
    neg_length <- round(neg_length)
    pos_length <- 100 * (pos_length / sum)
    pos_length <- round(pos_length)
    neg_color <- colorRampPalette(c("#174C9E", "white"))(neg_length)
    pos_color <- colorRampPalette(c("white", "#E8402B"))(pos_length)
    color_bar <- c(neg_color, pos_color)
    color_map <- circlize::colorRamp2(seq(from = min(exp), to = max(exp), length = length(color_bar)),
                                      color_bar)
    group <- apply(exp, 1, function(x){which.max(x)})
    exp_list <- split.data.frame(exp, f = group, drop = F)
    for (i in 1:length(exp_list)) {
      exp_list[[i]] <- exp_list[[i]][order(exp_list[[i]][,i],
                                           decreasing = T),,drop = F]
    }
    exp <- yulab.utils::rbindlist(exp_list)
    col_anno <- HeatmapAnnotation(df = data.frame(Group = colnames(exp)),
                                  col = list(Group = group.color),
                                  annotation_label = group.label,
                                  show_annotation_name = FALSE)
    if (!is.null(label_genes)) {
      row_label <- rowAnnotation(foo = anno_mark(at = match(label_genes, rownames(exp)),
                                                 labels = label_genes,
                                                 labels_gp = gpar(fontsize = 10),
                                                 padding = unit(0, "mm")))
      
      Heatmap(exp, col = color_map, color_space = "RGB",
              cluster_rows = F, cluster_columns = F,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              show_row_names = F,
              right_annotation = row_label, top_annotation = col_anno,
              heatmap_legend_param = list(title = "Scaled\nExpression"#,
                                          # at = c(round(min(exp),1), 0, max(exp)),
                                          # labels = c(round(min(exp),1), 0, max(exp)))
              ))
    } else {
      Heatmap(exp, col = color_map, color_space = "RGB",
              cluster_rows = F, cluster_columns = F,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              top_annotation = col_anno,
              heatmap_legend_param = list(title = "Scaled\nExpression"
                                          # at = c(round(min(exp),1), 0, max(exp)),
                                          # labels = c(round(min(exp),1), 0, max(exp))))
              ))
    }
    
  }
  
  CD4_T_DEGs <- FindAllMarkers(CD4_T, only.pos = T)
  CD4_T_DEGs %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  top20 <- as.data.frame(top20)
  top20 <- top20[!duplicated(top20$gene),]
  
  CD4_T_DEGs %>%
    group_by(cluster) %>%
    slice_head(n = 3) %>%
    ungroup() -> top3
  top3 <- as.data.frame(top3)
  top3 <- top3[!duplicated(top3$gene),]
  library(ComplexHeatmap)
  pdf("./绘图/29.CD4_T_SubCelltype_top20_差异基因_heatmap.pdf", height = 8, width = 6)
  Cell_group_heatmap(seu = CD4_T,
                     group.by = "SubCelltype",
                     group.color = CD4_T_colors,
                     group.label = "SubCelltype",
                     genes = unique(top20$gene),
                     limit = c(-2, 2),
                     label_genes = unique(top3$gene))
  dev.off()
  
  Idents(CD4_T) <- CD4_T$Sample_rename
  CD4_T_Group_DEGs <- FindAllMarkers(CD4_T, only.pos = T)
  CD4_T_Group_DEGs <- CD4_T_Group_DEGs[order(CD4_T_Group_DEGs$cluster,
                                             -CD4_T_Group_DEGs$avg_log2FC),]
  CD4_T_Group_DEGs_list <- split.data.frame(CD4_T_Group_DEGs,
                                            f = list(CD4_T_Group_DEGs$cluster))
  names(CD4_T_Group_DEGs_list) <- gsub(pattern = "/", replacement = "_",
                                       fixed = T, x = names(CD4_T_Group_DEGs_list))
  openxlsx::write.xlsx(CD4_T_Group_DEGs_list,
                       "./绘图/30.CD4_T_所有亚型_4个组别的高表达基因.xlsx")
  CD4_T_Group_DEGs <- c()
  for (i in c("OT-I/FX + CD40 mAb",
              "OT-I/FX", "OT-I/XCL1")) {
    temp_diff <- FindMarkers(CD4_T, only.pos = T,
                             ident.1 = i,
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[order(temp_diff$avg_log2FC,
                                 decreasing = T),]
    temp_diff$cluster <- i
    temp_diff$gene <- rownames(temp_diff)
    CD4_T_Group_DEGs <- c(CD4_T_Group_DEGs,
                          list(temp_diff))
  }
  names(CD4_T_Group_DEGs) <- c("OT-I_FX + CD40 mAb",
                               "OT-I_FX", "OT-I_XCL1")
  openxlsx::write.xlsx(CD4_T_Group_DEGs,
                       "./绘图/30.CD4_T_所有亚型_3个处理组_分别vs_Vector_高表达基因.xlsx")
  
  
  
  Cell_group_heatmap2 <- function(seu,
                                  group.by = "Celltype_rename",
                                  group.color = Celltype_color,
                                  group.label = "Celltype",
                                  genes = unique(top5$gene),
                                  label_genes = unique(top5$gene)[unique(top5$gene) %in% unlist(lapply(Celltype_markers_list, FUN = as.character))],
                                  limit = c(-2, 2),
                                  flip = F,
                                  show_row_names = T,
                                  show_column_names = T,
                                  cluster_columns = T,
                                  cluster_rows = T,
                                  row_names_side = "left",
                                  column_names_side = "bottom",
                                  column_names_rot = 90,
                                  row_names_rot = 0,
                                  body_width = unit(8, "cm"),
                                  body_height = unit(8, "cm"),
                                  cols = c("#174C9E", "#E8402B"),
                                  genes_order = NULL
  ){
    exp <- seu@assays$RNA@data[genes,]
    exp <- t(as.matrix(exp))
    exp <- aggregate.data.frame(exp,
                                by = list(seu@meta.data[,group.by]),
                                FUN = mean)
    rownames(exp) <- exp$Group.1
    exp <- exp[,-1]
    exp <- exp[,colSums(exp) != 0]
    for (i in 1:ncol(exp)) {
      exp[,i] <- (exp[,i] - mean(exp[,i])) / sd(exp[,i])
    }
    exp <- as.matrix(t(exp))
    if (!is.null(genes_order)) {
      exp <- exp[genes_order,]
    }
    exp[exp > limit[2]] <- limit[2]
    exp[exp < limit[1]] <- limit[1]
    neg_length <- 0 - min(exp)
    pos_length <- max(exp) - 0
    sum <- neg_length + pos_length
    neg_length <- 100 * (neg_length / sum)
    neg_length <- round(neg_length)
    pos_length <- 100 * (pos_length / sum)
    pos_length <- round(pos_length)
    neg_color <- colorRampPalette(c(cols[1], "white"))(neg_length)
    pos_color <- colorRampPalette(c("white", cols[2]))(pos_length)
    color_bar <- c(neg_color, pos_color)
    color_map <- circlize::colorRamp2(seq(from = min(exp), to = max(exp), length = length(color_bar)),
                                      color_bar)
    if (flip == T) {
      exp <- as.matrix(data.frame(t(exp),
                                  check.rows = F,
                                  check.names = F))
    }
    if (!is.null(genes_order)) {
      cluster_rows = FALSE
      cluster_columns = FALSE
    }
    Heatmap(exp, col = color_map, color_space = "RGB",
            cluster_rows = cluster_rows, cluster_columns = cluster_columns,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            # top_annotation = col_anno,
            show_row_names = show_row_names,
            show_column_names = show_column_names,
            column_dend_height = unit(0, "mm"),
            row_dend_width = unit(0, "mm"),
            row_names_side = row_names_side,
            column_names_side = column_names_side,
            show_row_dend = F, show_column_dend = F,
            column_names_rot = column_names_rot,
            row_names_rot = row_names_rot,
            width = body_width, height = body_height,
            heatmap_legend_param = list(title = "Scaled\nExpression"
                                        # at = c(round(min(exp),1), 0, max(exp)),
                                        # labels = c(round(min(exp),1), 0, max(exp))))
            ))
    
  }
  pdf("./绘图/31.CD4_T_co-stimu_co-inhibi_热图.pdf",
      height = 3.5, width = 8)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      # genes = c("Tnfrsf10b", "Cd244a", "Tnfsf10", "Tnfrsf18", "Havcr2", "Ctla4", "Pdcd1", "Tnfrsf9", "Icos",
                      #           "Cd160", "Tnfrsf4", "Lag3", "Cd27", "Cd28", "Tigit"),
                      genes = c("Cd160", "Icos", "Cd27", "Cd28", "Fas", "Tnfrsf4",
                                "Tnfrsf9", "Tnfrsf10b", "Tnfrsf18",
                                "Pdcd1", "Ctla4", "Lag3", "Havcr2", "Cd244a", "Tigit"),
                      limit = c(-2, 2),
                      flip = T,
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = T,
                      show_column_names = T,
                      cluster_columns = T,
                      cluster_rows = F,
                      column_names_rot = 45,
                      row_names_rot = 0,
                      body_width = unit(15 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"),
                      cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/31.CD4_T_TCR_signaling_热图.pdf",
      height = 3.5, width = 8)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      # genes = c("Nr4a1", "Nfatc1", "Nfkbia", "Nfatc2",
                      #           "Nfkbib", "Nfkb1", "Nfatc3", "Jun", "Fos"),
                      genes = c("Ifng", "Nr4a1", "Nr4a2", "Nfatc1", "Nfatc2",
                                "Nfkbib", "Nfkbia", "Nfatc4", "Nfat5", "Nfkb1",
                                "Nfatc3", "Jun", "Fos"),
                      genes_order = c("Nfatc2", "Nr4a1",
                                      "Nfat5", "Nfkbia", "Nr4a2", "Ifng", "Nfkbib",
                                      "Nfatc1", "Jun", "Fos", "Nfatc3", "Nfkb1"),
                      limit = c(-2, 2),
                      flip = T,
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = T,
                      show_column_names = T,
                      cluster_columns = T,
                      cluster_rows = F,
                      column_names_rot = 45,
                      row_names_rot = 0,
                      body_width = unit(13 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"),
                      cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  pdf("./绘图/31.CD4_T_Activation_Differentiation_热图.pdf",
      height = 3.5, width = 13)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      # genes = c("Cd44", "Tbx21", "Cd69", "Klrg1", "Runx3", "Cd101",
                      #           "Entpd1", "Il2rb", "Il2ra", "Klrc1", "Klrd1", "Mki67",
                      #           "Id2", "Cd38", "Eomes", "Lef1", "Bcl2", "Bach2", "Tcf7",
                      #           "Id3", "Tox", "Slamf6", "Sell", "Il7r", "Bcl6", "Runx1"),
                      genes = c("Id2", "Klrg1", "Klrc1", "Klrd1", "Nkg7", "Klf2",
                                "Foxp1", "Tbx21", "Batf", "Tox", "Il2ra", "Il2rb", "Il21r",
                                "Stat5a", "Stat5b", "Entpd1", "Cd101", "Cd38", "Eomes",
                                "Mki67", "Sell", "Cd69", "Bcl6", "Id3", "Il7r", "Lef1", "Bcl2",
                                "Bach2", "Slamf6", "Runx1", "Runx3"),
                      genes_order = c("Tbx21", "Klrc1",
                                      "Klrd1", "Id2",
                                      "Stat5a", "Il2ra", "Entpd1", "Cd101", "Cd38", "Il21r", "Il2rb",
                                      "Eomes", "Batf", "Nkg7", "Sell", "Mki67", "Cd69", "Bcl6",
                                      "Klf2", "Slamf6", "Il7r", "Bcl2", "Bach2", "Foxp1",
                                      "Lef1", "Runx1", "Stat5b", "Tox", "Klrg1", "Id3", "Runx3"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(31 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/31.CD4_T_IFN_signaling_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      # genes = c("Ifitm3", "Ifitm1", "Ifitm2", "Ifitm5", "Stat1"),
                      genes = c("Ifng", "Jak1", "Jak2", "Stat1", "Stat3", "Ifitm1",
                                "Isg15", "Ifitm3", "Ifitm2", "Ifitm5"),
                      genes_order = c("Isg15", "Stat3", "Stat1", "Ifng", "Ifitm1", "Ifitm3",
                                      "Jak2", "Ifitm2", "Jak1"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(10 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/31.CD4_T_Effector_molecules_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      genes = c("Ifng", "Gzmb", "Prf1", "Tnf", "Gzmd",
                                "Gzme", "Gzmk", "Il2", "Gzma", "Lamp1"),
                      genes_order = c("Tnf", "Gzmb", "Ifng",
                                      "Gzmk", "Prf1", "Il2", "Gzma", "Gzmd", "Lamp1", "Gzme"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = T,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(10 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  
  pdf("./绘图/31.CD4_T_Chemokine_热图.pdf",
      height = 3.5, width = 10)
  Cell_group_heatmap2(seu = CD4_T,
                      group.by = "Sample_rename",
                      group.color = CD4_T_colors,
                      group.label = "Sample",
                      # genes = c("Cx3cr1", "Ccr5", "Ccr2", "Cxcr3",
                      #           "Cxcr6", "Cxcr5", "Ccr7"),
                      genes = c("Ccl3", "Xcl1", "Ccr7", "Ccl4", "Ccl5",
                                "Ccr5", "Ccr2", "Cxcr3", "Cxcr6", "Cxcr5",
                                "Cx3cr1"),
                      genes_order = c("Xcl1", "Ccr2", "Cxcr3", "Ccr5", "Ccl3", "Ccl4", "Ccl5",
                                      "Cxcr6", "Cxcr5", "Ccr7", "Cx3cr1"),
                      limit = c(-2, 2),
                      flip = T, row_names_side = "left",
                      column_names_side = "top", show_row_names = T,
                      show_column_names = T, cluster_columns = F,
                      cluster_rows = F, column_names_rot = 45,
                      row_names_rot = 0, body_width = unit(11 * 0.7, "cm"),
                      body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
  dev.off()
  
  library(scRNAtoolVis)
  CD4_T_sub_diff <- c()
  table(CD4_T$SubCelltype, CD4_T$Sample_rename)
  CD4_SubCelltype <- unique(CD4_T$SubCelltype)
  for (c in CD4_SubCelltype) {
    temp <- subset(CD4_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
    temp_diff$cluster <- c
    temp_diff$gene <- rownames(temp_diff)
    CD4_T_sub_diff <- as.data.frame(rbind(CD4_T_sub_diff,
                                          temp_diff))
  }
  jjVolcano <- function(diffData = NULL, myMarkers = NULL, order.by = c("avg_log2FC"), 
                        log2FC.cutoff = 0.25, pvalue.cutoff = 0.05, adjustP.cutoff = 0.01, 
                        topGeneN = 5, col.type = "updown", back.col = "grey93", pSize = 0.75, 
                        aesCol = c("#0099CC", "#CC3333"), legend.position = c(0.7, 0.9), base_size = 10,
                        tile.col = jjAnno::useMyCol("paired", n = 9), 
                        cluster.order = NULL, polar = FALSE, expand = c(-1, 1), flip = FALSE, ...) 
  {
    diff.marker <- diffData %>% dplyr::filter(abs(avg_log2FC) >= 
                                                log2FC.cutoff & p_val < pvalue.cutoff)
    diff.marker <- diff.marker %>% dplyr::mutate(type = ifelse(avg_log2FC >= 
                                                                 log2FC.cutoff, "sigUp", "sigDown")) %>% dplyr::mutate(type2 = ifelse(p_val_adj < 
                                                                                                                                        adjustP.cutoff, paste("adjust Pvalue < ", adjustP.cutoff, 
                                                                                                                                                              sep = ""), paste("adjust Pvalue >= ", adjustP.cutoff, 
                                                                                                                                                                               sep = "")))
    if (!is.null(cluster.order)) {
      diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
    }
    back.data <- purrr::map_df(unique(diff.marker$cluster), function(x) {
      tmp <- diff.marker %>% dplyr::filter(cluster == x)
      new.tmp <- data.frame(cluster = x, min = min(tmp$avg_log2FC) - 
                              0.2, max = max(tmp$avg_log2FC) + 0.2)
      return(new.tmp)
    })
    top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
    top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = topGeneN, 
                                                          order_by = get(order.by))
    top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = topGeneN, 
                                                          order_by = get(order.by))
    top.marker <- rbind(top.marker.max, top.marker.min)
    if (!is.null(myMarkers)) {
      top.marker <- diff.marker %>% dplyr::filter(gene %in% 
                                                    myMarkers)
    }
    else {
      top.marker <- top.marker
    }
    p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, 
                                                    y = avg_log2FC)) + ggplot2::geom_col(data = back.data, 
                                                                                         ggplot2::aes(x = cluster, y = min), fill = back.col) + 
      ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, 
                                                       y = max), fill = back.col)
    if (col.type == "updown") {
      p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), 
                                      size = pSize) + ggplot2::scale_color_manual(values = c(sigDown = aesCol[1], 
                                                                                             sigUp = aesCol[2]))
    }
    else if (col.type == "adjustP") {
      p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2), 
                                      size = pSize) + ggplot2::scale_color_manual(values = c(aesCol[2], 
                                                                                             aesCol[1]))
    }
    p3 <- p2 + ggplot2::scale_y_continuous(n.breaks = 6) + ggplot2::theme_classic(base_size = base_size) + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                     legend.position = legend.position, legend.title = ggplot2::element_blank(), 
                     legend.background = ggplot2::element_blank()) + ggplot2::xlab("Clusters") + 
      ggplot2::ylab("Average log2FoldChange") + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
    p4 <- p3 + ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, 
                                               fill = cluster), color = "black", height = log2FC.cutoff * 
                                    2, alpha = 0.3, show.legend = FALSE) + ggplot2::scale_fill_manual(values = tile.col) + 
      ggrepel::geom_text_repel(data = top.marker, ggplot2::aes(x = cluster, 
                                                               y = avg_log2FC, label = gene), max.overlaps = 50, 
                               ...)
    if (polar == TRUE) {
      p5 <- p4 + geomtextpath::geom_textpath(ggplot2::aes(x = cluster, 
                                                          y = 0, label = cluster)) + ggplot2::scale_y_continuous(n.breaks = 6, 
                                                                                                                 expand = ggplot2::expansion(mult = expand)) + ggplot2::theme_void(base_size = base_size) + 
        ggplot2::theme(legend.position = legend.position, 
                       legend.title = ggplot2::element_blank()) + ggplot2::coord_polar(clip = "off", 
                                                                                       theta = "x")
    }
    else {
      if (flip == TRUE) {
        p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
          ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                           y = 0, label = cluster), color = "black", size = 3.3) + ggplot2::theme(axis.line.y = ggplot2::element_blank(), 
                                                                                                                  axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
          ggplot2::coord_flip()
      }
      else {
        p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
          ggplot2::geom_text(ggplot2::aes(x = cluster, 
                                          y = 0, label = cluster), color = "black", size = 3.3) + ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                                                                                                                 axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
      }
    }
    return(p5)
  }
  
  # CD4_T_clusters <- levels(CD4_T$SubCelltype)
  # CD4_T_clusters <- data.frame(original = CD4_T_clusters,
  #                              rename = CD4_T_clusters)
  # CD4_T_clusters[3, "rename"] <- "CD4 IFN\nresponse"
  # CD4_T_clusters[6, "rename"] <- "CD4 T\nexpanding"
  # rownames(CD4_T_clusters) <- CD4_T_clusters$original
  # CD4_T_sub_diff_rename <- CD4_T_sub_diff
  # CD4_T_sub_diff_rename$cluster <- CD4_T_clusters[CD4_T_sub_diff_rename$cluster,
  #                                                 2]
  # CD4_T_colors_rename <- CD4_T_colors
  # names(CD4_T_colors_rename) <- CD4_T_clusters$rename
  pdf("./绘图/32.CD4_T_不同亚群_最后一组_vs_第一组.pdf",
      width = 8, height = 5.5)
  jjVolcano(diffData = CD4_T_sub_diff,
            log2FC.cutoff = 0.25, 
            base_size = 14, pSize = 0.75, #设置点的大小
            fontface = 'italic', #设置字体形式
            aesCol = c("#0099CC","#CC3333"), #设置点的颜色
            tile.col = CD4_T_colors, #设置cluster的颜色
            adjustP.cutoff = 1,
            cluster.order = levels(CD4_T$SubCelltype),
            col.type = "updown",
            topGeneN = 5, legend.position = c(0.85, 0.9)
  ) + labs(x = "")
  dev.off()
  
  Vio_plot <- function(diff_gene = Epithelial_diff, logfc_thr = 0.25,
                       cols = c("#FF4500", "#4169E1"),
                       p_thr = 0.05, group = c("TACSTD2+", "TACSTD2-"),
                       title = "TACSTD2+ VS TACSTD2-",
                       genes = NULL) {
    diff_gene <- diff_gene[order(diff_gene$avg_log2FC, decreasing=T),]
    DOWN_row <- intersect(which(diff_gene$p_val < p_thr),
                          which(diff_gene$avg_log2FC < -logfc_thr))
    UP_row <- intersect(which(diff_gene$p_val < p_thr),
                        which(diff_gene$avg_log2FC > logfc_thr))
    DOWN <- diff_gene[DOWN_row,]
    UP <- diff_gene[UP_row,]
    DOWN <- DOWN[order(DOWN$avg_log2FC,decreasing = F),]
    UP <- UP[order(UP$avg_log2FC,decreasing = T),]
    color <- rep("gray", length(diff_gene$avg_log2FC))
    cex <- rep(1, length(diff_gene$avg_log2FC))
    cex[DOWN_row] <- 1.5
    cex[UP_row] <- 1.5
    cex <- as.character(cex)
    color[DOWN_row] <- cols[2]
    color[UP_row] <- cols[1]
    color <- factor(color, levels = c(cols[1], cols[2], "gray"))
    color_value <- c(cols[1], cols[2], "gray")
    names(color_value) <- color_value
    if (!is.null(genes)) {
      genes <- genes[genes %in% rownames(diff_gene)]
      sigdiff <- diff_gene[genes,]
    } else {
      sigdiff <- rbind(DOWN[1:5,], UP[1:5,])
    }
    vio_plot <- ggplot(data = diff_gene, aes(x = avg_log2FC, y = -log10(p_val),
                                             size = cex, colour = color,
                                             alpha = "1")) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            panel.grid = element_blank(),
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 10, color = "black")) +
      labs(x = "log2 FoldChange", y = "log10 P value") +
      geom_point() +
      geom_hline(aes(yintercept=-log10(p_thr)), colour = "black", linetype="dashed") +
      geom_vline(aes(xintercept= -logfc_thr), colour="black", linetype="dashed") +
      geom_vline(aes(xintercept= logfc_thr), colour="black", linetype="dashed") +
      scale_color_manual(name = "",
                         values = color_value,
                         labels = c(paste0("\n",nrow(UP)," UP genes for\n",group[1]),
                                    paste0("\n",nrow(DOWN)," UP genes for\n",group[2]),
                                    "Not significant")) +
      scale_size_manual(values = c('1' = 1, "1.5" = 1.5)) +
      scale_alpha_manual(values = c("0.3" = 0.3, "1" = 1)) +
      guides(size = "none", alpha = "none",
             color = guide_legend(override.aes = list(size = 2))) +
      geom_text_repel(data = sigdiff, aes(x = avg_log2FC, 
                                          y = -log10(p_val), 
                                          label = gene,
                                          alpha = rep("1",nrow(sigdiff))),
                      size = 3, fontface = "italic",
                      box.padding = unit(0.8, "lines"),
                      point.padding = unit(0, "lines"), 
                      min.segment.length = 0,
                      segment.color = "black",
                      colour="#000000",
                      show.legend = FALSE,
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) 
    vio_plot
  }
  Idents(CD4_T) <- CD4_T$Sample_rename
  temp_diff <- FindMarkers(CD4_T, logfc.threshold = 0,
                           ident.1 = "OT-I/FX + CD40 mAb",
                           ident.2 = "OT-I")
  temp_diff$gene <- rownames(temp_diff)
  p <- Vio_plot(diff_gene = temp_diff, logfc_thr = 0.25,
                cols = c("#CC3333", "#0099CC"),
                p_thr = 0.05,
                group = c("OT-I/FX + CD40 mAb", "OT-I"),
                title = "CD4 T",
                genes = c("Gzmb", "Tnfrsf9", "Klrg1", "Nkg7", "Tnfrsf4", "Ccl4", "Icos", "Snx9", "Ccr5",
                          "Irf1", "Cxcl10", "Il2ra", "Gzmk", "Ccl5", "Ccl3", "Xcl1", "Cxcr6", "Entpd1",
                          "Jun", "Ifng", "Tcf7", "Rgs10", "Itgb3", "Itgb1", "Itga4", "Fos", "Klf2", "Tsc22d3",
                          "Il7r", "S1pr1", "Tsc22d3", "Ythdc1", "Lef1", "Slamf6", "Cd69", "C1qb"))
  pdf(paste0("./绘图/33.CD4_T_所有亚型_最后一组_vs_第一组_火山图.pdf"),
      width = 6.2, height = 5)
  print(p)
  dev.off()
  
  temp_diff$gene <- rownames(temp_diff)
  temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
  temp_diff <- temp_diff[abs(temp_diff$avg_log2FC) > 0.25,]
  temp_diff <- temp_diff[order(temp_diff$avg_log2FC, decreasing = T),]
  openxlsx::write.xlsx(temp_diff,
                       "./绘图/33.CD4_T_所有亚型_最后一组_vs_第一组_差异表达基因.xlsx")
  
  library(ggrepel)
  for (c in levels(CD4_SubCelltype)) {
    if (c == "CD4 Teff") {
      next
    }
    temp <- subset(CD4_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp, logfc.threshold = 0,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff$gene <- rownames(temp_diff)
    p <- Vio_plot(diff_gene = temp_diff, logfc_thr = 0.25,
                  cols = c("#CC3333", "#0099CC"),
                  p_thr = 0.05,
                  group = c("OT-I/FX + CD40 mAb", "OT-I"),
                  title = c)
    pdf(paste0("./绘图/34.CD4_T_", gsub(pattern = "/",
                                      replacement = "_",
                                      fixed = T, x = c), "_最后一组_vs_第一组_火山图.pdf"),
        width = 6.2, height = 5)
    print(p)
    dev.off()
    
  }
  
  library(clusterProfiler)
  library(createKEGGdb)
  # createKEGGdb::create_kegg_db("mmu")
  # install.packages("./KEGG.db_1.0.tar.gz",type="source")
  library(KEGG.db)
  KEGG_DATA <- clusterProfiler::download_KEGG("mmu")
  KEGG_df <- KEGG_DATA[["KEGGPATHID2EXTID"]]
  KEGG_df_2 <- KEGG_DATA[["KEGGPATHID2NAME"]]
  rownames(KEGG_df_2) <- KEGG_df_2$from
  KEGG_df$from <- KEGG_df_2[KEGG_df$from,2]
  KEGG_df$from <- unlist(lapply(KEGG_df$from, function(x){
    unlist(strsplit(x, " - Mus musculus (house mouse)",
                    fixed = T))[1]
  }))
  library(org.Mm.eg.db)
  id <- bitr(geneID = unique(KEGG_df$to),
             fromType = "ENTREZID",
             toType = "SYMBOL",
             OrgDb = org.Mm.eg.db)
  rownames(id) <- as.character(id$ENTREZID)
  KEGG_df$to <- id[KEGG_df$to, "SYMBOL"]
  colnames(KEGG_df) <- c("term", "gene")
  KEGG_df$term <- paste0("KEGG ", KEGG_df$term)
  GeneSet1 <- read.gmt("./绘图/m5.go.bp.v2023.2.Mm.symbols.gmt")
  GeneSet1$term <- as.character(GeneSet1$term)
  GeneSet1$term <- unlist(lapply(GeneSet1$term, function(x){
    temp <- unlist(strsplit(x, "_"))
    temp <- temp[-1]
    temp <- unlist(lapply(temp, function(y){
      y <- tolower(y)
      y <- Hmisc::capitalize(y)
    }))
    paste0(temp, collapse = " ")
  }))
  GeneSet1$term <- paste0("GO ", GeneSet1$term)
  
  GeneSet2 <- read.gmt("./绘图/mh.all.v2023.2.Mm.symbols.gmt")
  GeneSet2$term <- as.character(GeneSet2$term)
  GeneSet2$term <- unlist(lapply(GeneSet2$term, function(x){
    temp <- unlist(strsplit(x, "_"))
    temp <- temp[-1]
    temp <- unlist(lapply(temp, function(y){
      y <- tolower(y)
      y <- Hmisc::capitalize(y)
    }))
    paste0(temp, collapse = " ")
  }))
  GeneSet2$term <- paste0("HALLMARK ", GeneSet2$term)
  
  GeneSet <- as.data.frame(rbind(GeneSet1, GeneSet2))
  GeneSet <- as.data.frame(rbind(GeneSet,
                                 KEGG_df))
  GeneSet_list <- split.data.frame(GeneSet,
                                   f = factor(GeneSet$term))
  GeneSet_list_length <- unlist(lapply(GeneSet_list, function(x){
    nrow(x)
  }))
  GeneSet_list <- GeneSet_list[GeneSet_list_length >= 5]
  
  TERM2GENE <- rbindlist(GeneSet_list)
  library(enrichplot)
  Idents(CD4_T) <- CD4_T$Sample_rename
  CD4_T_diff <- FindMarkers(CD4_T,
                            ident.1 = "OT-I/FX + CD40 mAb",
                            ident.2 = "OT-I",
                            logfc.threshold = 0)
  CD4_T_diff <- CD4_T_diff[order(CD4_T_diff$avg_log2FC,
                                 decreasing = T),]
  diff_list <- CD4_T_diff$avg_log2FC
  names(diff_list) <- rownames(CD4_T_diff)
  set.seed(123)
  CD4T_gsea <- GSEA(diff_list, TERM2GENE = TERM2GENE,
                    verbose = T, eps = 0, pvalueCutoff = 1)
  CD4T_gsea_result <- CD4T_gsea@result
  CD4T_gsea_result <- CD4T_gsea_result[CD4T_gsea_result$pvalue < 0.05,]
  CD4T_gsea_result <- CD4T_gsea_result[order(CD4T_gsea_result$NES,
                                             decreasing = T),]
  CD4T_gsea@result <- CD4T_gsea_result
  saveRDS(CD4T_gsea, "./绘图/CD4T不分亚型_第四组_vs_第一组_gsea_pathway.rds")
  openxlsx::write.xlsx(CD4T_gsea_result,
                       "./绘图/CD4T不分亚型_第四组_vs_第一组_gsea_result.xlsx")
  
  temp <- grep(pattern = "myc|e2f|oxidative|mtorc1|fatty|peroxisome|pyruvate|stat5|glycolysis|Il2",
               ignore.case = TRUE, x = CD4T_gsea@result$Description)
  temp <- CD4T_gsea@result[temp,]
  # temp <- temp[-c(3,7,12),]
  library(GseaVis)
  pdf("./绘图/35.CD4_T_不分亚型_第四组_vs_第一组_gsea_plot.pdf",
      width = 7, height = 4)
  gseaNb(object = CD4T_gsea,
         geneSetID = temp$ID,
         addPval = T, subPlot = 2,
         # curveCol = scales::hue_pal()(7),
         curveCol = c("#89C75F","#F37B7D","#D24B27", "#FF8C00",
                      "#3BBCA8","#6E4B9E","#0C727C","#D8A767",
                      "#272E6A"),
         htHeight = 0.6)#,
  # subPlot = 2, addPval = TRUE,
  # pvalX = 0.95, pvalY = 0.8,
  # pvalSize = 4, pDigit = 1)
  dev.off()
  
  for (c in levels(CD4_SubCelltype)) {
    if (c == "CD4 Teff") {
      next
    }
    temp <- subset(CD4_T, subset = SubCelltype == as.character(c))
    Idents(temp) <- temp$Sample_rename
    temp_diff <- FindMarkers(temp, logfc.threshold = 0,
                             ident.1 = "OT-I/FX + CD40 mAb",
                             ident.2 = "OT-I")
    temp_diff <- temp_diff[order(temp_diff$avg_log2FC,
                                 decreasing = T),]
    temp_diff_list <- temp_diff$avg_log2FC
    names(temp_diff_list) <- rownames(temp_diff)
    set.seed(123)
    temp_gsea <- GSEA(temp_diff_list, TERM2GENE = TERM2GENE,
                      verbose = T, eps = 0, pvalueCutoff = 1)
    temp_gsea_result <- temp_gsea@result
    temp_gsea_result <- temp_gsea_result[temp_gsea_result$pvalue < 0.05,]
    temp_gsea_result <- temp_gsea_result[order(temp_gsea_result$NES,
                                               decreasing = T),]
    temp_gsea@result <- temp_gsea_result
    openxlsx::write.xlsx(temp_gsea_result,
                         paste0("./绘图/",
                                "CD4T_", gsub(pattern = "/",
                                              replacement = "_",
                                              x = c), "_第四组_vs_第一组_gsea_result.xlsx"))
    saveRDS(temp_gsea,
            paste0("./绘图/",
                   "CD4T_", gsub(pattern = "/",
                                 replacement = "_",
                                 x = c), "_第四组_vs_第一组_gsea_pathway.rds"))
    
  }
  
  
  # RNA Velocity
  {
    library(velocyto.R)
    library(Seurat)
    library(dplyr)
    library(SeuratWrappers)
    library(pagoda2)
    library(harmony)
    library(SeuratDisk)

    setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Velocyto")
    features <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/1/CD40-T_1/outs/per_sample_outs/CD40-T_1/count/sample_filtered_feature_bc_matrix/features.tsv",
                         header = F, sep = "\t")
    sum(duplicated(features$V2)) # 37 个重复的SYMBOL
    row_dup <- which(duplicated(features$V2))
    features$V2[row_dup] <- paste0(features$V2[row_dup],
                                   "--",
                                   features$V1[row_dup])
    rownames(features) <- features$V1
    
    looms <- c("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Velocyto/CD40-T_2/velocyto/CD40-T_2.loom",
               "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Velocyto/Fx-T_2/velocyto/Fx-T_2.loom",
               "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Velocyto/Vector-T_2/velocyto/Vector-T_2.loom",
               "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/Velocyto/XCL1-T_2/velocyto/XCL1-T_2.loom")
    names(looms) <- c("CD40-T_2", "Fx-T_2", "Vector-T_2", "XCL1-T_2")
    looms_list <- c()
    prefix <- c("CD40-T_2:", "Fx-T_2:", "Vector-T_2:", "XCL1-T_2:")
    Sample_suffix <- c("CD40-T_2_", "Fx-T_2_", "Vector-T_2_", "XCL1-T_2_")
    for (loom in 1:4) {
      test <- read.loom.matrices(looms[loom])
      for (i in 1:3) {
        temp <- test[[i]]
        temp <- temp[rownames(temp) %in% features[,1],]
        rownames(temp) <- features[rownames(temp),2]
        colnames(temp) <- unlist(lapply(colnames(temp), function(x){
          temp_1 <- unlist(strsplit(x, prefix[loom], fixed = T))[2]
          temp_1 <- unlist(strsplit(temp_1, "x", fixed = T))[1]
          temp_1 <- paste0(Sample_suffix[loom], temp_1, "-1")
          temp_1
        }))
        temp <- temp[,colnames(temp) %in% colnames(seu)]
        test[[i]] <- temp
      }
      looms_list <- c(looms_list,
                      list(test))
    }
    names(looms_list) <- names(looms)
    
    for (loom in 1:length(looms_list)) {
      loom_seu <- as.Seurat(looms_list[[loom]])
      looms_list[[loom]] <- loom_seu
    }
    looms_seu <- merge(x = looms_list[[1]], y = looms_list[-1])
    
    ### CD4+ T
    DimPlot(CD4_T, group.by = "SubCelltype")
    CD4_looms_seu <- looms_seu[,colnames(CD4_T)]
    CD4_looms_seu$percent.mt <- CD4_T@meta.data[colnames(CD4_looms_seu), "percent.mt"]
    CD4_looms_seu$SubCelltype <- CD4_T@meta.data[colnames(CD4_looms_seu), "SubCelltype"]
    CD4_looms_seu$SubCelltype <- as.character(CD4_looms_seu$SubCelltype)
    CD4_looms_seu[["RNA"]] <- CD4_looms_seu[["spliced"]]
    CD4_looms_seu %>%
      NormalizeData(normalization.method = "LogNormalize") %>%
      FindVariableFeatures() %>%
      ScaleData(vars.to.regress = "percent.mt") %>%
      RunPCA() %>%
      RunUMAP(dims = 1:30, reduction = "pca") -> CD4_looms_seu
    DefaultAssay(CD4_looms_seu) <- "RNA"
    cells <- rownames(CD4_looms_seu@reductions[["umap"]]@cell.embeddings)
    CD4_looms_seu@reductions[["umap"]]@cell.embeddings <- CD4_T@reductions[["umap"]]@cell.embeddings[cells,]
    SaveH5Seurat(CD4_looms_seu, "CD4_Velocyto.h5Seurat")
    Convert("CD4_Velocyto.h5Seurat", dest = "h5ad")
  }
  
  Treg <- subset(CD4_T, subset = SubCelltype == "Treg")
  
  VlnPlot2 <- function(seu = Treg, group.by = "Sample_rename",
                       features = "Foxp3", cols = sample_color,
                       group_order = NULL,
                       comparisons = list(c("OT-I", "OT-I/FX"),
                                          c("OT-I", "OT-I/FX + CD40 mAb"))) {
    temp <- seu@assays$RNA@data[features,]
    temp <- as.data.frame(temp)
    colnames(temp) <- features
    temp$group <- seu@meta.data[,group.by]
    if (!is.null(group_order)) {
      temp$group <- factor(as.character(temp$group),
                           levels = group_order)
    }
    ggplot(data = temp, aes(x = temp[,"group"], y = temp[,features], fill = temp[,"group"])) +
      geom_violin() +
      geom_jitter(size = 0.1) +
      ggpubr::stat_compare_means(comparisons = comparisons) +
      scale_fill_manual(values = sample_color) +
      theme_classic() +
      ggtitle(features) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            panel.grid = element_blank(),
            legend.position = "none",
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10, color = "black")) +
      labs(x = "", y = "Expression level")
                                 
  }
  VlnPlots <- lapply(c("Foxp3", "Nrp1", "Il10", "Rorc",
                       "Nt5e", "Ifngr1", "Gzmb", "Il23r", "Satb1", "Gata3",
                       "Entpd1", "Ctla4", "Il2ra", "Pdcd1", "Icos"),
                     function(x){
                       VlnPlot2(seu = Treg, group.by = "Sample_rename",
                                features = x, cols = sample_color,
                                group_order = NULL,
                                comparisons = list(c("OT-I", "OT-I/FX"),
                                                   c("OT-I", "OT-I/FX + CD40 mAb")))
                     })
  pdf("./绘图/37.CD4_T_Treg_genes_不同组vlnplot_V2.pdf",
      height = 20, width = 10)
  cowplot::plot_grid(plotlist = VlnPlots)
  dev.off()
  
  CD4_T_Ratio <- openxlsx::read.xlsx("./绘图/CD4_T_Ratio.xlsx",
                                     check.names = F)
  colnames(CD4_T_Ratio) <- c("Celltype", "OT-I", "OT-I/XCL1", "OT-I/FX",
                             "OT-I/FX + CD40 mAb")
  CD4_T_Ratio <- reshape::melt.data.frame(CD4_T_Ratio,
                                          id.vars = "Celltype")
  unique(CD4_T_Ratio$Celltype)
  pdf("./绘图/38.CD4_T_Tpex_Tex_lineplot.pdf", width = 3.5, height = 6)
  ggplot() +
    geom_line(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "Tpex/Tex",],
              aes(x = variable, y = value),
              group = 1) +
    geom_point(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "Tpex/Tex",],
               aes(x = variable, y = value, color = variable),
               size = 4, shape = 15) +
    scale_color_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "Tpex/Tex") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  unique(CD4_T_Ratio$Celltype)
  pdf("./绘图/38.CD4_T_Tex_memory_lineplot.pdf", width = 3.5, height = 6)
  ggplot() +
    geom_line(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "Tex/memory",],
              aes(x = variable, y = value),
              group = 1) +
    geom_point(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "Tex/memory",],
               aes(x = variable, y = value, color = variable),
               size = 4, shape = 15) +
    scale_color_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "Tex/memory") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  pdf("./绘图/38.CD4_T_CD4_treg_lineplot.pdf", width = 3.5, height = 6)
  ggplot() +
    geom_line(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "CD4/treg",],
              aes(x = variable, y = value),
              group = 1) +
    geom_point(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "CD4/treg",],
               aes(x = variable, y = value, color = variable),
               size = 4, shape = 15) +
    scale_color_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "CD4/treg") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  
  CD4_T_Ratio <- openxlsx::read.xlsx("./绘图/CD4_T_Ratio.xlsx",
                                     check.names = F)
  colnames(CD4_T_Ratio) <- c("Celltype", "OT-I", "OT-I/XCL1", "OT-I/FX",
                             "OT-I/FX + CD40 mAb")
  rownames(CD4_T_Ratio) <- CD4_T_Ratio$Celltype
  CD4_T_Ratio <- CD4_T_Ratio[,-1]
  CD4_T_Ratio <- as.data.frame(t(CD4_T_Ratio))
  
  pdf("./绘图/38.CD4_T_Tpex_Tex_pieplot.pdf",
      height = 6, width = 6)
  pie(x = CD4_T_Ratio[,"Tpex/Tex"],
      labels = paste0(rownames(CD4_T_Ratio),
                      "\n", round(CD4_T_Ratio[,"Tpex/Tex"], 2)),
      col = sample_color[rownames(CD4_T_Ratio)],
      main = "Tpex/Tex")
  dev.off()
  
  pdf("./绘图/38.CD4_T_Tex_memory_pieplot.pdf",
      height = 6, width = 6)
  pie(x = CD4_T_Ratio[,"Tex/memory"],
      labels = paste0(rownames(CD4_T_Ratio),
                      "\n", round(CD4_T_Ratio[,"Tex/memory"], 2)),
      col = sample_color[rownames(CD4_T_Ratio)],
      main = "Tex/memory")
  dev.off()
  
  pdf("./绘图/38.CD4_T_CD4_treg_pieplot.pdf",
      height = 6, width = 6)
  pie(x = CD4_T_Ratio[,"CD4/treg"],
      labels = paste0(rownames(CD4_T_Ratio),
                      "\n", round(CD4_T_Ratio[,"CD4/treg"], 2)),
      col = sample_color[rownames(CD4_T_Ratio)],
      main = "CD4/treg")
  dev.off()
  
  CD8_T
  table(All_seu$Celltype_rename)
  All_seu$CD8_Treg <- "Other"
  All_seu@meta.data[rownames(CD8_T@meta.data),
                    "CD8_Treg"] <- "CD8"
  All_seu@meta.data[rownames(Treg@meta.data),
                    "CD8_Treg"] <- "Treg"
  temp <- as.data.frame.array(table(All_seu$Sample_rename,
                                    All_seu$CD8_Treg))
  for (i in 1:nrow(temp)) {
    temp[i,] <- temp[i,] / sum(temp[i,])
  }
  temp$CD8 / temp$Treg
  pdf("./绘图/38.CD4_T_CD8_Treg_pieplot.pdf",
      height = 6, width = 6)
  pie(x = temp$CD8 / temp$Treg,
      labels = paste0(rownames(CD4_T_Ratio),
                      "\n", round(temp$CD8 / temp$Treg, 2)),
      col = sample_color[rownames(CD4_T_Ratio)],
      main = "CD8/Treg")
  dev.off()
  
  CD4_T_Ratio <- data.frame(Celltype = "CD8/Treg",
                            variable = rownames(temp),
                            value = temp$CD8 / temp$Treg)
  CD4_T_Ratio$variable <- factor(CD4_T_Ratio$variable,
                                 levels = CD4_T_Ratio$variable)
  pdf("./绘图/38.CD4_T_CD8_Treg_lineplot.pdf", width = 3.5, height = 6)
  ggplot() +
    geom_line(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "CD8/Treg",],
              aes(x = variable, y = value),
              group = 1) +
    geom_point(data = CD4_T_Ratio[CD4_T_Ratio$Celltype == "CD8/Treg",],
               aes(x = variable, y = value, color = variable),
               size = 4, shape = 15) +
    scale_color_manual(values = sample_color) +
    theme_classic() +
    labs(x = "", y = "CD8/Treg") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
}


### Macro_Mono
Macro_Mono <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/注释/Macro_Mono/Macro_Mono_annotated.rds")
Macro_Mono <- subset(Macro_Mono,
                     subset = SubCelltype != "Unidentified")
DimPlot(Macro_Mono)
Macro_Mono %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt"),
            model.use = c("negbinom")) %>%
  RunPCA() -> Macro_Mono

Macro_Mono <- RunUMAP(Macro_Mono,
                      reduction = "pca",
                      dims = 1:30, n.neighbors = 50,
                      min.dist = 0.8, n.epochs = 1000)
set.seed(124)
Idents(Macro_Mono)
Macro_Mono <- RenameIdents(Macro_Mono,
                           "Mac-C1qb" = "Mac-C1qb",
                           "Mac-Mif" = "Mac-Plin2",
                           "Mac-mt" = "Mac-Mt",
                           "Mac-Plin2" = "Mac-Plin2",
                           "Mac-prolif" = "Mac-prolif.",
                           "Mono-Ace" = "Mono-Ace",
                           "Mono-Chil3" = "Mono-Chil3"
                           )
Macro_Mono$SubCelltype <- factor(as.character(Idents(Macro_Mono)),
                                 levels = c("Mac-C1qb", "Mac-Mt", "Mac-Plin2", "Mac-prolif.",
                                            "Mono-Ace", "Mono-Chil3"))
# Macro_Mono_color <- sample(ArchR::ArchRPalettes$circus,
#                            size = 6, replace = FALSE)
Macro_Mono_color <- c("#C92326", "#3083B8", "#621314", "#ED776F", "#BD5EA0", "#2E247E")
names(Macro_Mono_color) <- NULL
names(Macro_Mono_color) <- levels(Macro_Mono$SubCelltype)
cluster_meta <- Embeddings(Macro_Mono, reduction = "umap")
cluster_meta <- as.data.frame(cluster_meta)
cluster_meta$Cluster <- Macro_Mono@meta.data[rownames(cluster_meta),"SubCelltype"]
cluster_meta <- aggregate.data.frame(cluster_meta[,1:2],
                                     by = list(cluster_meta$Cluster),
                                     FUN = median)
cluster_meta$Group.1 <- 0:5
Macro_Mono_color_2 <- Macro_Mono_color
names(Macro_Mono_color_2) <- NULL
Idents(Macro_Mono) <- Macro_Mono$SubCelltype
pdf("./绘图/39.Macro_Mono_UMAP(with clusters).pdf",
    height = 4.5, width = 6)
DimPlot(Macro_Mono, cols = Macro_Mono_color_2) +
  geom_point(data = cluster_meta,
             aes(x = UMAP_1, y = UMAP_2),
             size = 7, alpha = 0.7, color = "#EEE9E9") +
  geom_text(data = cluster_meta,
            aes(x = UMAP_1, y = UMAP_2, label = Group.1))
dev.off()

pdf("./绘图/39.Macro_Mono_UMAP.pdf",
    height = 4.5, width = 6)
DimPlot(Macro_Mono, cols = Macro_Mono_color)
dev.off()

bottom_colorbar_df <- data.frame(X = levels(Idents(Macro_Mono)), Y = 1)
bottom_colorbar <- ggplot(data = bottom_colorbar_df,
                          aes(x = X, y = Y, fill = X)) +
  geom_tile() +
  labs(fill = "SubCluster") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = Macro_Mono_color)
p <- DotPlot(Macro_Mono, features = c("C1qa", "C1qb", "Apoe", "Ms4a7", "C5ar1",
                                 "mt-Co1", "mt-Co3",
                                 "Mmp12", "Cd68", "Hmox1", "Plin2",
                                 "Birc5", "Pclaf", "Stmn1",
                                 "Ace", "Lyz2", "Csf1r",
                                 "Ly6c2", "Chil3"),
             cols = c("gray", "#CD3333")) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = "") +
  coord_flip()
pdf("./绘图/40.Macro_Mono_marker_dotplot.pdf",
    width = 5, height = 5)
p %>%
aplot::insert_top(bottom_colorbar,
                  height = 0.04)
dev.off()

Macro_Mono <- AddModuleScore(Macro_Mono,
                             features = list(c("H2-Aa", "H2-Ab1", "Il1b",
                                               "Il15", "Il18", "Tnf", "Cd40",
                                               "Cd74", "Ptgs2", "Cxcl9", "Cxcl10")))
colnames(Macro_Mono@meta.data)[ncol(Macro_Mono@meta.data)] <- "Antitumoral"
Macro_Mono <- AddModuleScore(Macro_Mono,
                             features = list(c("Pf4", "Apoe", "Arg1", "Mrc1", "Il10", "Spp1", "Ccl9",
                                               "C1qa", "C1qb")))
colnames(Macro_Mono@meta.data)[ncol(Macro_Mono@meta.data)] <- "Protumoral"


AUCell_score <- AUCell::AUCell_run(exprMat = Macro_Mono@assays$RNA@data,
                                   geneSets = list("Antitumoral" = c("H2-Aa", "H2-Ab1", "Il1b",
                                                                     "Il15", "Il18", "Tnf", "Cd40",
                                                                     "Cd74", "Ptgs2", "Cxcl9", "Cxcl10"),
                                                   "Protumoral" = c("Pf4", "Apoe", "Arg1", "Mrc1", "Il10", "Spp1", "Ccl9",
                                                                    "C1qa", "C1qb")))
AUCell_score <- as.data.frame(t(AUCell_score@assays@data@listData[["AUC"]]))
AUCell_score$SubCelltype <- Macro_Mono@meta.data[rownames(AUCell_score), "SubCelltype"]
Data <- reshape::melt.data.frame(AUCell_score, id.vars = "SubCelltype")
Data <- Data[Data$SubCelltype != "Mac-Mt",]

# Data <- Macro_Mono@meta.data[, c("Antitumoral", "Protumoral")]
# Data$SubCelltype <- Macro_Mono$SubCelltype
# Data <- reshape::melt.data.frame(Data, id.vars = "SubCelltype")
# Data <- Data[Data$SubCelltype != "Mac-Mt",]
# pdf("./绘图/Macro_Mono_补充_Pro_Anti_Tumoral.pdf",
#     width = 6, height = 3.5)
# ggplot(data = Data, aes(x = SubCelltype, y = value, fill = variable)) +
#   geom_violin(aes(fill = variable), scale = "width") +
#   ggforce::geom_sina(aes(colour = variable), scale = "width",
#                      size = 1, alpha = 0.5) +
#   geom_boxplot(position = position_dodge(0.9),
#                width = 0.15, outlier.shape = NA,
#                linewidth = 1, color = "black") +
#   ggpubr::stat_compare_means(label = "p.signif") +
#   theme_bw() +
#   labs(x = "", y = "AddMoudleScore", fill = "", color = "") +
#   theme(panel.grid = element_blank(),
#         axis.text = element_text(size = 10, colour = "black"),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.text = element_text(size = 10)) +
#   scale_fill_manual(values = c("#5F9EA0", "#FF6347")) +
#   scale_color_manual(values = c("#5F9EA0", "#FF6347"))
# dev.off()

pdf("./绘图/Macro_Mono_补充_Pro_Anti_Tumoral_AUCell.pdf",
    width = 6, height = 3.5)
ggplot(data = Data, aes(x = SubCelltype, y = value, fill = variable)) +
  geom_violin(aes(fill = variable), scale = "width") +
  ggforce::geom_sina(aes(colour = variable), scale = "width",
                     size = 1, alpha = 0.5) +
  geom_boxplot(position = position_dodge(0.9),
               width = 0.15, outlier.shape = NA,
               linewidth = 1, color = "black") +
  ggpubr::stat_compare_means(label = "p.signif") +
  theme_bw() +
  labs(x = "", y = "AUCell score", fill = "", color = "") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#1FA673", "#D65513")) +
  scale_color_manual(values = c("#1FA673", "#D65513"))
dev.off()

AUCell_score$Anti_vs_Pro <- (AUCell_score$Antitumoral + 0.001) / (AUCell_score$Protumoral + 0.001)
AUCell_score$Anti_vs_Pro_bi <- "Other"
AUCell_score$Anti_vs_Pro_bi[AUCell_score$Anti_vs_Pro > 1] <- "Antitumoral cell"
AUCell_score$Anti_vs_Pro_bi[AUCell_score$Anti_vs_Pro < 1] <- "Protumoral cell"
AUCell_score$Sample_rename <- Macro_Mono@meta.data[rownames(AUCell_score), "Sample_rename"]
pdf("./绘图/Macro_Mono_补充_Pro_Anti_Tumoral_percentage.pdf",
    width = 4.5, height = 4.5)
ggplot(data = AUCell_score[AUCell_score$Anti_vs_Pro_bi != "Other",],
       aes(x = Sample_rename, fill = Anti_vs_Pro_bi)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = c("#1FA673", "#D65513")) +
  theme_bw() +
  labs(x = "", y = "Relative abundance", fill = "", color = "") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(labels = scales::percent)
dev.off()

AUCell_score$Sample_rename <- Macro_Mono@meta.data[rownames(AUCell_score), c("Sample_rename")]
ggplot(data = AUCell_score, aes(x = Sample_rename, y = Antitumoral)) +
  geom_violin(aes(fill = Sample_rename), scale = "width") +
  ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                     size = 1, alpha = 0.5) +
  geom_boxplot(position = position_dodge(0.9),
               width = 0.15, outlier.shape = NA,
               linewidth = 1, color = "black") +
  facet_wrap(.~SubCelltype)

ggplot(data = AUCell_score, aes(x = Sample_rename, y = Protumoral)) +
  geom_violin(aes(fill = Sample_rename), scale = "width") +
  ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                     size = 1, alpha = 0.5) +
  geom_boxplot(position = position_dodge(0.9),
               width = 0.15, outlier.shape = NA,
               linewidth = 1, color = "black") +
  facet_wrap(.~SubCelltype)

{
  unique(Macro_Mono$SubCelltype)
  Idents(Macro_Mono) <- as.character(Macro_Mono$SubCelltype)
  df <- data.frame(Original = c("Mac-Plin2", "Mac-C1qb", "Mono-Chil3",
                                "Mac-Mt", "Mono-Ace", "Mac-prolif."),
                   Rename = c("Protumoral", "Protumoral", "Antitumoral",
                              "Other", "Antitumoral", "Protumoral"))
  rownames(df) <- df$Original
  df2 <- Macro_Mono@meta.data[, c("Sample_rename",
                                  "SubCelltype")]
  df2$Rename <- df[df2$SubCelltype, "Rename"]
  df2 <- df2[df2$Rename != "Other",]
  ggplot(data = df2, aes(x = Sample_rename, fill = Rename)) +
    geom_bar(stat = "count", position = "fill") +
    scale_fill_manual(values = c("#5F9EA0", "#FF6347"))
  
  Data <- Macro_Mono@meta.data[, c("SubCelltype", "Sample_rename", "Antitumoral", "Protumoral")]
  # Data <- reshape::melt.data.frame(Data, id.vars = "SubCelltype")
  Data <- Data[Data$SubCelltype != "Mac-Mt",]
  ggplot(data = Data, aes(x = Sample_rename, y = Protumoral)) +
    geom_violin(aes(fill = Sample_rename), scale = "width") +
    ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                       size = 1, alpha = 0.5) +
    geom_boxplot(position = position_dodge(0.9),
                 width = 0.15, outlier.shape = NA,
                 linewidth = 1, color = "black") +
    facet_wrap(.~SubCelltype, scales = "free")
  ggplot(data = Data, aes(x = Sample_rename, y = Antitumoral)) +
    geom_violin(aes(fill = Sample_rename), scale = "width") +
    ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                       size = 1, alpha = 0.5) +
    geom_boxplot(position = position_dodge(0.9),
                 width = 0.15, outlier.shape = NA,
                 linewidth = 1, color = "black") +
    facet_wrap(.~SubCelltype, scales = "free")
  
  ggplot(data = Data, aes(x = Sample_rename, y = Antitumoral)) +
    geom_violin(aes(fill = Sample_rename), scale = "width") +
    ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                       size = 1, alpha = 0.5) +
    geom_boxplot(position = position_dodge(0.9),
                 width = 0.15, outlier.shape = NA,
                 linewidth = 1, color = "black")
  ggplot(data = Data, aes(x = Sample_rename, y = Protumoral)) +
    geom_violin(aes(fill = Sample_rename), scale = "width") +
    ggforce::geom_sina(aes(colour = Sample_rename), scale = "width",
                       size = 1, alpha = 0.5) +
    geom_boxplot(position = position_dodge(0.9),
                 width = 0.15, outlier.shape = NA,
                 linewidth = 1, color = "black")
}




library(homologene)
s.genes <- cc.genes$s.genes
s.genes <- homologene(s.genes, inTax = 9606, outTax = 10090)
s.genes <- s.genes$`10090`
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- homologene(g2m.genes, inTax = 9606, outTax = 10090)
g2m.genes <- g2m.genes$`10090`
Macro_Mono <- CellCycleScoring(Macro_Mono, s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = FALSE)
Macro_Mono$Phase <- factor(Macro_Mono$Phase,
                           levels = c("G1", "S", "G2M"))
table(Macro_Mono$SubCelltype, Macro_Mono$Phase)
left_colorbar_df <- data.frame(X = 1, Y = levels(Idents(Macro_Mono)))
left_colorbar <- ggplot(data = left_colorbar_df,
                          aes(x = X, y = Y, fill = Y)) +
  geom_tile() +
  labs(fill = "SubCluster") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = Macro_Mono_color)
p <- percent_bar(sce = Macro_Mono, Ident = "Phase",
            Group = "SubCelltype") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(fill = "Phase") +
  coord_flip() +
  scale_fill_manual(values = c("#6959CD", "#6B8E23", "#BDB76B"))
pdf("./绘图/41.Macro_Mono_Cellcyle.pdf",
    height = 4, width = 7)
p %>%
  aplot::insert_left(left_colorbar,
                     width = 0.04)
dev.off()

Idents(Macro_Mono) <- Macro_Mono$Sample
Macro_Mono <- RenameIdents(Macro_Mono,
                           "CD40-T_2" = "OT-I/FX + CD40 mAb",
                           "Fx-T_2" = "OT-I/FX",
                           "XCL1-T_2" = "OT-I/XCL1",
                           "Vector-T_2" = "OT-I"
                           )
Macro_Mono$Sample_rename <- factor(as.character(Idents(Macro_Mono)),
                                   levels = c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb"))
sample_color <- c("#1E699D", "#8F94C0", "#D87177", "#B42225")
names(sample_color) <- c("OT-I", "OT-I/XCL1", "OT-I/FX", "OT-I/FX + CD40 mAb")

pdf("./绘图/42.Macro_Mono_细胞类型比例.pdf",
    width = 5, height = 5)
percent_bar(sce = Macro_Mono, Ident = "SubCelltype",
            Group = "Sample_rename",
            fill_color = Macro_Mono_color)
dev.off()
temp <- as.data.frame.array(table(Macro_Mono$SubCelltype,
                                  Macro_Mono$Sample_rename))
for (i in 1:ncol(temp)) {
  temp[,i] <- temp[,i] / colSums(temp)[i]
}
temp <- data.frame(Celltype = rownames(temp),
                   temp, check.rows = F, check.names = F)
openxlsx::write.xlsx(temp,
                     "./绘图/42.Macro_Mono_细胞类型比例.xlsx")

# Idents(Macro_Mono) <- factor(Macro_Mono$SubCelltype,
#                              levels = c("Mac-C1qb", "Mac-Mif", "Mac-mt", "Mac-Plin2",
#                                         "Mac-prolif", "Mono-Ace", "Mono-Chil3"))
Idents(Macro_Mono) <- Macro_Mono$SubCelltype
Macro_Mono_DEGs <- FindAllMarkers(Macro_Mono, only.pos = T)
Macro_Mono_DEGs %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
top20 <- as.data.frame(top20)
Macro_Mono <- ScaleData(Macro_Mono, features = rownames(Macro_Mono))
pdf("./绘图/43.Macro_Mono_细胞类型_top20_deg_heatmap.pdf",
    width = 7, height = 15)
DoHeatmap(Macro_Mono, features = top20$gene,
          group.colors = Macro_Mono_color,
          label = FALSE)
dev.off()


genes <- c("C1qb", "mt-Co1", "Plin2",
           "Stmn1", "Ace", "Chil3")
feature_plot_list <- lapply(genes, function(x){
  FeaturePlot(Macro_Mono, features = c(x), order = T,
              min.cutoff = 0) +
    scale_color_gradientn(colours = c("#00868B", "#F5DEB3", "#CDCD00", "#FF0000")) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          plot.title = element_text(size = 14, colour = "black",
                                    face = "bold", hjust = 0.5),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank())
})
pdf("./绘图/45.Macro_Mono_细胞marker_featurePlot.pdf",
    width = 12, height = 10)
patchwork::wrap_plots(feature_plot_list,
                      byrow = T, nrow = 3, ncol = 3)
dev.off()

Macro_Mono_sub <- subset(Macro_Mono, subset = SubCelltype %in% c("Mac-C1qb", "Mono-Chil3"))
pdf("./绘图/46.Macro_Mono_sub_dimplot_splited_by_group.pdf",
    height = 3.2, width = 10)
DimPlot(Macro_Mono_sub, split.by = "Sample_rename",
        cols = Macro_Mono_color)
dev.off()

temp <- as.data.frame.array(table(Macro_Mono$SubCelltype,
                                  Macro_Mono$Sample_rename))
for (i in 1:ncol(temp)) {
  temp[,i] <- temp[,i] / sum(temp[,i])
}
temp <- data.frame(Celltype = rownames(temp),
                   temp)
colnames(temp) <- c("Celltype", "OT-I", "OT-I/XCL1",
                    "OT-I/FX", "OT-I/FX + CD40 mAb")
temp <- reshape::melt.data.frame(temp,
                                 id.vars = "Celltype")
unique(temp$Celltype)
temp <- temp[temp$Celltype %in% c("Mac-C1qb", "Mono-Chil3"),]
pdf("./绘图/47.Macro_Mono_不同样本中细胞类型比例.pdf",
    width = 4.8, height = 4.5)
ggplot(data = temp,
       aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_color) +
  theme_classic() +
  facet_wrap(".~Celltype", nrow = 1) +
  labs(x = "", y = "Percentage of cell type") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        strip.text = element_text(size = 13, colour = "black"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

pdf("./绘图/48.Macro_Mono_Interferon-related genes_vlnplot.pdf",
    height = 9, width = 10)
VlnPlot(Macro_Mono, cols = Macro_Mono_color, pt.size = 0,
        features = c("Cxcl16", "Slamf7", "Irf4", "Irf7",
                     "Cxcl2", "Ifitm1", "Ifitm3", "Ifitm6",
                     "Isg15", "Arg1"))
dev.off()


pdf("./绘图/48.Macro_Mono_Mac-Plin2_genes_vlnplot.pdf",
    height = 9, width = 10)
VlnPlot(Macro_Mono, cols = Macro_Mono_color, pt.size = 0,
        features = c("Cxcl16", "Slamf7", "Irf4", "Irf7",
                     "Cxcl2", "Ifitm1", "Ifitm3", "Ifitm6",
                     "Isg15", "Arg1"))
dev.off()

pdf("./绘图/48.Macro_Mono_Mac-C1qb_genes_vlnplot.pdf",
    height = 12, width = 10)
VlnPlot(Macro_Mono, cols = Macro_Mono_color, pt.size = 0,
        features = c("Apoe", "C1qa", "C1qb", "C1qc",
                     "Lpl", "Ctsb", "Mrc1", "Pf4", "Cd63",
                     "Wfdc17", "Ccl7", "Trem2", "Dab2", "Ebi3"))
dev.off()

pdf("./绘图/48.Macro_Mono_Mac-Chil3_genes_vlnplot.pdf",
    height = 12, width = 10)
VlnPlot(Macro_Mono, cols = Macro_Mono_color, pt.size = 0,
        features = c("Cxcl2", "Cxcl10", "Ifitm6", "Tnf", "Il1b",
                     "Plac8", "Irf1", "Irf7", "Cd14", "Hp", "Mgst1",
                     "Wfdc17", "Anxa2", "Fn1", "Ccl4", "Cybb"))
dev.off()

library(clusterProfiler)
library(org.Mm.eg.db)
Macro_Mono_DEGs_list <- split.data.frame(Macro_Mono_DEGs, f = list(Macro_Mono_DEGs$cluster))
for (i in 1:length(Macro_Mono_DEGs_list)) {
  temp <- Macro_Mono_DEGs_list[[i]]
  temp <- temp[order(temp$avg_log2FC, decreasing = T),]
  Macro_Mono_DEGs_list[[i]] <- temp[1:50,]
}
for (i in 1:length(Macro_Mono_DEGs_list)) {
  temp <- Macro_Mono_DEGs_list[[i]]
  temp <- bitr(geneID = temp$gene, fromType = "SYMBOL",
               toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  temp <- temp$ENTREZID
  Macro_Mono_DEGs_list[[i]] <- temp
}

Cluster_GO <- compareCluster(Macro_Mono_DEGs_list,
                             fun = "enrichGO", ont = "BP",
                             readable = TRUE, OrgDb = org.Mm.eg.db,
                             pvalueCutoff = 1, qvalueCutoff = 1)
Cluster_GO@compareClusterResult <- Cluster_GO@compareClusterResult[Cluster_GO@compareClusterResult$pvalue < 0.05,]
Cluster_GO@compareClusterResult$Description <- Hmisc::capitalize(Cluster_GO@compareClusterResult$Description)
temp <- Cluster_GO@compareClusterResult
temp_list <- split.data.frame(temp, f = list(temp$Cluster))
openxlsx::write.xlsx(temp_list, "./绘图/44.Macro_Mono_subtype_GO.xlsx")
temp <- Cluster_GO@compareClusterResult
# dup_term <- temp$Description[duplicated(temp$Description)]
# temp <- temp[-which(temp$Description %in% dup_term),]
temp$log10P <- -log10(temp$pvalue)
# temp %>%
#   group_by(Cluster) %>%
#   slice_head(n = 5) %>%
#   ungroup() -> top5
top5 <- temp[temp$ID %in% c("GO:0002495", "GO:0002443", "GO:0030595", "GO:0050900", "GO:0034341",
                            "GO:0050851", "GO:0002757", "GO:0050852", "GO:0002253", "GO:0001909",
                            "GO:0006007", "GO:0006735", "GO:0006090", "GO:0030595", "GO:0060326",
                            "GO:0000280", "GO:0007051", "GO:0007059", "GO:0044772", "GO:0045787",
                            "GO:0050900", "GO:0007159", "GO:1901652", "GO:0002253", "GO:0006909",
                            "GO:0019221", "GO:0031349", "GO:0002221", "GO:0002758", "GO:0050727"),]
sum(duplicated(top5$Description))
pathway <- c("GO:0002495", "GO:0002443", "GO:0030595", "GO:0050900", "GO:0034341",
             "GO:0050851", "GO:0002757", "GO:0050852", "GO:0002253", "GO:0001909",
             "GO:0006007", "GO:0006735", "GO:0006090", "GO:0030595", "GO:0060326",
             "GO:0000280", "GO:0007051", "GO:0007059", "GO:0044772", "GO:0045787",
             "GO:0050900", "GO:0007159", "GO:1901652", "GO:0002253", "GO:0006909",
             "GO:0019221", "GO:0031349", "GO:0002221", "GO:0002758", "GO:0050727")
pathway <- pathway[!duplicated(pathway)]
pathway_df <- top5[top5$ID %in% pathway,]
pathway_df <- pathway_df[,2:3]
pathway_df <- pathway_df[!duplicated(pathway_df),]
top5$Description <- factor(top5$Description,
                           levels = rev(pathway_df$Description))
GeneRatio2 <- c()
for (i in 1:nrow(top5)) {
  temp2 <- unlist(strsplit(top5$GeneRatio[i], "/", fixed = T))
  GeneRatio2 <- c(GeneRatio2,
                  as.numeric(temp2[1]) / as.numeric(temp2[2]))
}
top5$GeneRatio2 <- GeneRatio2

p <- ggplot() +
  geom_point(data = top5, aes(x = Cluster, y = Description,
                              fill = log10P, size = GeneRatio2),
             color = "black", shape = 21) +
  scale_fill_gradient(low = "#27408B", high = "#CD2626") +
  theme_bw() +
  # theme(axis.text.y = element_text(size = 12, colour = rev(Macro_Mono_color[top5$Cluster])),
  theme(axis.text.y = element_text(size = 12, colour = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)) +
  labs(x = "", fill = "-log10(P-value)",
       size = "GeneRatio", y = "")

top_colorbar_df <- data.frame(X = levels(Idents(Macro_Mono)),
                              Y = 1)
top_colorbar <- ggplot(data = top_colorbar_df,
                        aes(x = X, y = Y, fill = X)) +
  geom_tile() +
  labs(fill = "SubCluster") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = Macro_Mono_color)

pdf("./绘图/44.Macro_Mono_Subtype_Specific_GO.pdf",
    height = 7, width = 10)
p %>%
  aplot::insert_top(top_colorbar,
                    height = 0.02)
dev.off()


# bottom_colorbar_df <- data.frame(X = levels(Idents(Macro_Mono)), Y = 1)
# bottom_colorbar <- ggplot(data = bottom_colorbar_df,
#                           aes(x = X, y = Y, fill = X)) +
#   geom_tile() +
#   labs(fill = "SubCluster") +
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12)) +
#   scale_fill_manual(values = Macro_Mono_color)
p <- DotPlot(Macro_Mono, features = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Nos2", "Il1b",
                                      "Il12a", "Il12b", "Il15", "Il18", "Tnf", "Cd40",
                                      "Cd74", "Cd80", "Cd86", "Ptgs2", "Cxcl9", "Cxcl10",
                                      "Stat1", "Ciita", "Sod2", "Acod1", "Cd274", "Plac8",
                                      "Tgfbi", "Chil3", "Lyz2", "Ccl3", "Ccl4", "Ifitm1",
                                      "Ifitm6", "Irf1", "Siglec1", "Cd209a")) +
  theme_bw() +
  # scale_color_gradientn(colours = c("#008B8B", "#3CB371", "#FFA500", "#FFB90F", "#FF1493", "#FF0000", "#CD2626")) +
  scale_color_gradientn(colours = c("#008B8B", "#00CD66", "#FFA500", "#CD0000")) +
  scale_y_discrete(position = "right") + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
  labs(x = "", y = "") +
  coord_flip()
pdf("./绘图/49.Macro_Mono_genes1_dotplot.pdf",
    width = 4.5, height = 7)
p #%>%
  # aplot::insert_top(bottom_colorbar,
  #                   height = 0.04)
dev.off()

p <- DotPlot(Macro_Mono, features = c("Fn1", "Pf4", "Apoe", "Arg1", "Mrc1", "Il1r2", "Il10",
                                      "Spp1", "C5ar1", "Ccl9", "C1qa", "C1qb", "C1qc", "Trem2",
                                      "Cd24a", "Vegfb", "Vegfa", "Csf3r", "Ccl2", "Il4ra", "Il34",
                                      "Csf1", "Cd276", "Pdcd1", "Stat6", "Cd163", "Wfdc17",
                                      "Tgfb1", "Egr2", "Socs3", "Retnla", "Hp"),
             cols = c("gray", "#CD3333")) +
  theme_bw() +
  # scale_color_gradientn(colours = c("#008B8B", "#3CB371", "#FFA500", "#FFB90F", "#FF1493", "#FF0000", "#CD2626")) +
  scale_color_gradientn(colours = c("#008B8B", "#00CD66", "#FFA500", "#CD0000")) +
  scale_y_discrete(position = "right") + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
  labs(x = "", y = "") +
  coord_flip()
pdf("./绘图/49.Macro_Mono_genes2_dotplot.pdf",
    width = 4.5, height = 6.8)
p # %>%
  # aplot::insert_top(bottom_colorbar,
  #                   height = 0.04)
dev.off()


diff_genes_num <- c()
for (i in as.character(unique(Macro_Mono$Sample_rename))[1:3]) {
  Idents(Macro_Mono) <- Macro_Mono$Sample_rename
  temp_diff <- FindMarkers(Macro_Mono, ident.1 = i,
                           ident.2 = "OT-I")
  temp_diff <- temp_diff[temp_diff$p_val < 0.05,]
  up <- sum(temp_diff$avg_log2FC > 0)
  down <- sum(temp_diff$avg_log2FC < 0)
  diff_genes_num <- as.data.frame(rbind(diff_genes_num,
                                        data.frame(Sample = i,
                                                   Gene_num = c(up, -down))))
}
diff_genes_num$Sample <- factor(diff_genes_num$Sample,
                                levels = c("OT-I/FX + CD40 mAb",
                                           "OT-I/FX", "OT-I/XCL1"))

pdf("./绘图/50.Macro_Mono_差异表达基因数目.pdf", 
    width = 7, height = 2.5)
ggplot() +
  geom_bar(data = diff_genes_num, aes(y = Sample, x = Gene_num, fill = Sample),
           stat = "identity", width = 0.7) +
  geom_text(data = diff_genes_num[diff_genes_num$Gene_num > 0,], aes(y = Sample, x = Gene_num,
                                                 fill = Sample, label = abs(Gene_num)),
            hjust = -0.2) +
  geom_text(data = diff_genes_num[diff_genes_num$Gene_num < 0,], aes(y = Sample, x = Gene_num,
                                                   fill = Sample, label = abs(Gene_num)),
            hjust = 1.2) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  annotate("text", x = -350, y = 3, label = "Down regulated", 
           size = 4, hjust = 0) +
  annotate("text", x = 300, y = 3, label = "Up regulated", 
           size = 4, hjust = 0) +
  scale_x_continuous(limits = c(-400, 550)) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(colour = "black", size = 10),
        plot.title = element_text(colour = "black", size = 12,
                                  face = "bold", hjust = 0.5),
        legend.position = "none") +
  labs(x = "Gene number", y = "") +
  scale_fill_manual(values = sample_color)
dev.off()

temp_diff <- FindMarkers(Macro_Mono, ident.1 = "OT-I/FX + CD40 mAb",
                         ident.2 = "OT-I", logfc.threshold = 0)
temp_diff$gene <- rownames(temp_diff)
library(ggrepel)
p <- Vio_plot(diff_gene = temp_diff, logfc_thr = 0.25,
              cols = c("#B42225", "#1E699D"),
              p_thr = 0.05,
              group = c("OT-I/FX + CD40 mAb", "OT-I"),
              title = c("Macro/Mono"),
              genes = c("C1qb", "C1qc", "C1qa", "Ctsb", "Apoe",
                        "Ccl7", "Wfdc17", "Lpl", "Tmed10", "Cd300ld3",
                        "C5ar1",
                        "Ly6i", "Il1b", "Ly6c2", "Cxcl2", "Acod1", "Cd14", "H2-Q6", "Plac8",
                        "Tnf", "Irf1", "Cd274", "Ifitm6", "Cxcl10", "Upp1", "Egr1"))
pdf("./绘图/51.Macro_Mono_最后一组_第一组_火山图.pdf",
    width = 6.2, height = 5)
p
dev.off()

temp_diff$Group <- ifelse(temp_diff$avg_log2FC > 0,
                          "OT-I_FX + CD40 mAb", "OT-I")
openxlsx::write.xlsx(split.data.frame(temp_diff,
                                      f = list(temp_diff$Group)),
                     "./绘图/51.Macro_Mono_最后一组_第一组_差异基因.xlsx")

library(ComplexHeatmap)
pdf("./绘图/52.Macro_Mono_geneset1_热图.pdf",
    height = 3.5, width = 10)
Cell_group_heatmap2(seu = Macro_Mono,
                    group.by = "Sample_rename",
                    group.color = sample_color,
                    group.label = "Sample",
                    genes = c("Ccl7", "Mrc1", "Timp2", "Trem2",
                              "Ms4a7", "Ccl2", "Dab2", "Serpinb6a",
                              "Lgals1", "Pf4", "Cd63", "Ctsd",
                              "Igf1", "Apoe", "Lpl", "H2-Q7", "Polr2l",
                              "Irf1", "Malat1", "Ccl5", "Btg1", "Hk2", "Tmsb10",
                              "AW112010"),
                    genes_order = c("Ccl7", "Mrc1", "Timp2", "Trem2",
                                    "Ms4a7", "Ccl2", "Dab2", "Serpinb6a",
                                    "Lgals1", "Pf4", "Cd63", "Ctsd",
                                    "Igf1", "Apoe", "Lpl", "H2-Q7", "Polr2l",
                                    "Irf1", "Malat1", "Ccl5", "Btg1", "Hk2", "Tmsb10",
                                    "AW112010"),
                    limit = c(-2, 2),
                    flip = T, row_names_side = "left",
                    column_names_side = "top", show_row_names = T,
                    show_column_names = T, cluster_columns = T,
                    cluster_rows = F, column_names_rot = 45,
                    row_names_rot = 0, body_width = unit(24 * 0.7, "cm"),
                    body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
dev.off()


pdf("./绘图/52.Macro_Mono_geneset2_热图.pdf",
    height = 3.5, width = 10)
Cell_group_heatmap2(seu = Macro_Mono,
                    group.by = "Sample_rename",
                    group.color = sample_color,
                    group.label = "Sample",
                    genes = c("H2-K1", "H2-Eb1", "H2-Aa", "Cd74",
                              "Ly6a", "Gzma", 
                              "Rpl6", "Rps27", "Rps29", "Rps28", "Gm10076"),
                    genes_order = c("H2-K1", "H2-Eb1", "H2-Aa", "Cd74",
                                    "Ly6a", "Gzma", 
                                    "Rpl6", "Rps27", "Rps29", "Rps28", "Gm10076"),
                    limit = c(-2, 2),
                    flip = T, row_names_side = "left",
                    column_names_side = "top", show_row_names = T,
                    show_column_names = T, cluster_columns = T,
                    cluster_rows = F, column_names_rot = 45,
                    row_names_rot = 0, body_width = unit(11 * 0.7, "cm"),
                    body_height = unit(4 * 0.7, "cm"), cols = c("#174C9E", "#E8402B"))
dev.off()

GSE150970_meta <- read.csv("GSE150970_All_CD45_Cell_Anno_Coordinates.txt",
                           sep = "\t", header = T, row.names = 1)
setwd("./GSE150970_RAW/")
files <- list.files(pattern = ".gz")
samples <- unlist(lapply(files, function(x){
  temp <- unlist(strsplit(x,"_"))[2]
  gsub(pattern = ".", replacement = "-", fixed = T, x = temp)
}))
file_sample <- data.frame(files = files,
                          samples = samples)
# file_sample <- file_sample[file_sample$samples %in% samples_sc,]
lapply(unique(file_sample$samples),function(x){
  dir.create(x, recursive = T) #为每个样本创建子文件夹
  original <- file_sample[which(file_sample$samples %in% x),1]
  after <- c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz" )
  file.rename(from = original, to = paste(x, after, sep = "/")) # 移动文件并重命名
})
GSE150970 <- c()
for (i in unique(file_sample$samples)) {
  print(i)
  temp <- Read10X(i)
  colnames(temp) <- paste0(i, "_", colnames(temp))
  temp <- CreateSeuratObject(temp)
  GSE150970 <- c(GSE150970, list(temp))
}
GSE150970 <- merge(x = GSE150970[[1]],
                   y = GSE150970[-1])
GSE150970$Cellname <- as.character(colnames(GSE150970))
GSE150970$Cellname <- unlist(lapply(GSE150970$Cellname,
                                    function(x){
                                      unlist(strsplit(x, "-1", fixed = T))[1]
                                    }))
GSE150970 <- subset(GSE150970, subset = Cellname %in% rownames(GSE150970_meta))
GSE150970$Celltype <- GSE150970_meta[GSE150970$Cellname, "cluster"]
sort(unique(GSE150970$Celltype))
GSE150970 <- NormalizeData(GSE150970)
GSE150970_M <- subset(GSE150970, subset = Celltype %in% c("M_1", "M_2", "M_3", "M_4", "M_5_Mono", "M_6", "M_7_Pro", "M_8", "M_9"))
Idents(GSE150970_M) <- GSE150970_M$Celltype
GSE150970_M_diff <- FindAllMarkers(GSE150970_M, only.pos = T, logfc.threshold = 0.25)
M3_diff <- GSE150970_M_diff[GSE150970_M_diff$cluster == "M_3",]

Idents(Macro_Mono) <- Macro_Mono$SubCelltype
subcluster_diff <- FindAllMarkers(Macro_Mono, only.pos = T)
subcluster_diff_list <- split.data.frame(subcluster_diff, f = subcluster_diff$cluster)
for (i in 1:length(subcluster_diff_list)) {
  subcluster_diff_list[[i]] <- subcluster_diff_list[[i]][order(subcluster_diff_list[[i]]$avg_log2FC,
                                                               decreasing = T),]
}
for (i in 1:length(subcluster_diff_list)) {
  subcluster_diff_list[[i]] <- subcluster_diff_list[[i]][1:100,]
}
for (i in 1:length(subcluster_diff_list)) {
  print(sum(subcluster_diff_list[[i]]$gene %in% M3_diff$gene))
}

setwd("../")
openxlsx::write.xlsx(M3_diff, "./绘图/53.M3_diff_list.xlsx")


##### TCR
CD8_T <- readRDS("./注释/CD8+/CD8_T_annotated.rds")
table(CD8_T$SubCelltype)
CD8_T$SubCelltype <- factor(CD8_T$SubCelltype,
                            levels = c("CD8 Tex", "CD8 Proliferating-Mki67",
                                       "CD8 Proliferating-Cdc20", "CD8 Tpex",
                                       "CD8 IFN Response", "CD8 Tem", "CD8 Teff",
                                       "CD8 Texpanding", "CD8 Naive/Tcm"))
CD8_T_colors <- ArchR::ArchRPalettes$circus
CD8_T_colors <- CD8_T_colors[1:length(unique(CD8_T$SubCelltype))]
names(CD8_T_colors) <- levels(CD8_T$SubCelltype)

CD4_T <- readRDS("./注释/CD4+/CD4_T_annotated.rds")
table(CD4_T$SubCelltype)
Idents(CD4_T) <- CD4_T$SubCelltype
CD4_T <- RenameIdents(CD4_T,
                      "CD4 Tex" = "CD4 Tex",
                      "CD4 Tpex" = "CD4 Tpex",
                      "CD4 IFN Response" = "CD4 IFN response",
                      "CD4 Tem" = "CD4 Tem",
                      "CD4 Teff" = "CD4 Teff",
                      "CD4 Texpanding" = "CD4 T expanding",
                      "CD4 Treg" = "Treg",
                      "CD4 Naive/Memory" = "CD4 Naïve/Tcm"
)
CD4_T$SubCelltype <- factor(as.character(Idents(CD4_T)),
                            levels = c("CD4 Tex", "CD4 Tpex", "CD4 IFN response", "CD4 Tem",
                                       "CD4 Teff", "CD4 T expanding", "Treg", "CD4 Naïve/Tcm"))
CD4_T_colors <- ArchR::ArchRPalettes$paired
CD4_T_colors <- CD4_T_colors[1:length(unique(CD4_T$SubCelltype))]
names(CD4_T_colors) <- levels(CD4_T$SubCelltype)

All_seu <- readRDS("第一次注释_4.Data2_Annotated_with_CD45_classification_renamed.rds")
unique(All_seu$Sample_rename)
table(All_seu$Sample, All_seu$Sample_rename)

# CD8 TCR
{
  CD40_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/CD40-T_2/outs/per_sample_outs/CD40-T_2/vdj_t/filtered_contig_annotations.csv")
  CD40_T_2_TCR$barcode <- paste0("CD40-T_2_", CD40_T_2_TCR$barcode)
  CD40_T_2_TCR$Sample <- "OT-I/FX + CD40 mAb"
  Vector_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/Vector-T_2/outs/per_sample_outs/Vector-T_2/vdj_t/filtered_contig_annotations.csv")
  Vector_T_2_TCR$barcode <- paste0("Vector-T_2_", Vector_T_2_TCR$barcode)
  Vector_T_2_TCR$Sample <- "OT-I"
  XCL1_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/XCL-T_2/outs/per_sample_outs/XCL-T_2/vdj_t/filtered_contig_annotations.csv")
  XCL1_T_2_TCR$barcode <- paste0("XCL1-T_2_", XCL1_T_2_TCR$barcode)
  XCL1_T_2_TCR$Sample <- "OT-I/XCL1"
  Fx_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/Fx-T_2/outs/per_sample_outs/Fx-T_2/vdj_t/filtered_contig_annotations.csv")
  Fx_T_2_TCR$barcode <- paste0("Fx-T_2_", Fx_T_2_TCR$barcode)
  Fx_T_2_TCR$Sample <- "OT-I/FX"
  All_TCR <- as.data.frame(rbind(CD40_T_2_TCR,
                                 Vector_T_2_TCR,
                                 XCL1_T_2_TCR,
                                 Fx_T_2_TCR))
  All_TCR_list <- split.data.frame(All_TCR,
                                   f = All_TCR$barcode)
  All_TCR_type <- c()
  for (i in 1:length(All_TCR_list)) {
    temp <- All_TCR_list[[i]]
    temp <- unique(temp$chain)
    if (length(temp) == 2) {
      All_TCR_type <- c(All_TCR_type,
                        "TRB and TRA")
      next
    }
    if (temp == "TRB") {
      All_TCR_type <- c(All_TCR_type,
                        "TRB Only")
      next
    }
    if (temp == "TRA") {
      All_TCR_type <- c(All_TCR_type,
                        "TRA Only")
      next
    }
  }
  table(All_TCR_type)
  All_TCR_type <- data.frame(Cells = names(All_TCR_list),
                             TCR_type = All_TCR_type)
  rownames(All_TCR_type) <- All_TCR_type$Cells
  All_TCR_remove_index <- c()
  All_TCR_list_filter_after <- c()
  for (i in 1:length(All_TCR_list)) {
    temp <- All_TCR_list[[i]]
    if (!(c("TRB") %in% temp$chain)) {
      All_TCR_remove_index <- c(All_TCR_remove_index,
                                i)
      next
    } else {
      temp <- filter(.data = temp,
                     umis == max(umis),
                     # reads == max(reads) & umis == max(umis),
                     .by = "chain")
      temp <- filter(.data = temp,
                     reads == max(reads),
                     # reads == max(reads) & umis == max(umis),
                     .by = "chain")
      temp <- as.data.frame(temp)
      All_TCR_list_filter_after <- c(All_TCR_list_filter_after,
                                     list(temp))
    }
  }
  All_TCR_filter_after <- data.table::rbindlist(All_TCR_list_filter_after)
  All_TCR_filter_after_list <- split.data.frame(All_TCR_filter_after,
                                                f = list(All_TCR_filter_after$barcode))
  for (i in 1:length(All_TCR_filter_after_list)) {
    temp <- All_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    All_TCR_filter_after_list[[i]] <- temp2
  }
  All_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(All_TCR_filter_after_list))
  All_TCR_filter_after_list_df$barcode <- unlist(lapply(All_TCR_filter_after_list_df$barcode,
                                                        function(x){
                                                          unlist(strsplit(x, ";"))[1]
                                                        }))
  rownames(All_TCR_filter_after_list_df) <- All_TCR_filter_after_list_df$barcode
  
  CD8_T$TCR_type <- All_TCR_type[colnames(CD8_T), 2]
  CD8_T$TCR_type[is.na(CD8_T$TCR_type)] <- "No recovery"
  CD8_T$TCR_type <- factor(CD8_T$TCR_type,
                           levels = c("TRB Only", "TRA Only",
                                      "TRB and TRA", "No recovery"))
  table(CD8_T$TCR_type)
  
  pdf("./绘图/54.CD8_T_TCR_type.pdf",
      height = 4.5, width = 6.2)
  DimPlot(CD8_T, group.by = "TCR_type", pt.size = 0.7) +
    scale_color_manual(values = c("#228B22", "#20B2AA",
                                  "#BA55D3", "#B5B5B5"))
  dev.off()
  
  CD8_T$Sample_rename <- All_seu@meta.data[colnames(CD8_T), "Sample_rename"]
  pdf("./绘图/55.CD8_T_TCR_type_percent_bar.pdf",
      height = 5.5, width = 6)
  percent_bar(CD8_T, Ident = "TCR_type", Group = "Sample_rename",
              fill_color = c("#228B22", "#20B2AA",
                             "#BA55D3", "#B5B5B5"))
  dev.off()
  
  CD8_T <- AddModuleScore(CD8_T, features = list(c("Gzma", "Gzmb", "Gzmc",
                                                   "Gzmd", "Gzme", "Gzmf",
                                                   "Gzmg", "Gzmk", "Gzmm")))
  colnames(CD8_T@meta.data)[59] <- "Cytotoxicity Score"
  
  pdf("./绘图/56.CD8_T_毒性打分_没有额外基因.pdf",
      width = 4, height = 6.5)
  ggplot(data = CD8_T@meta.data[,c("Sample_rename",
                                   "Cytotoxicity Score")],
         aes(x = Sample_rename, y = `Cytotoxicity Score`,
             fill = Sample_rename)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1,
                 color = "black", fill = "white") +
    theme_classic() +
    scale_fill_manual(values = sample_color) +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    labs(x = "") +
    ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                  c("OT-I", "OT-I/FX"),
                                                  c("OT-I", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/XCL1", "OT-I/FX"),
                                                  c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/FX", "OT-I/FX + CD40 mAb")
    ))
  dev.off()
  
  CD8_T <- AddModuleScore(CD8_T, features = list(c("Gzma", "Gzmb", "Gzmc",
                                                   "Gzmd", "Gzme", "Gzmf",
                                                   "Gzmg", "Gzmk", "Gzmm",
                                                   "Ifng", "Tnf", "Prf1")),
                          search = TRUE)
  colnames(CD8_T@meta.data)[60] <- "Cytotoxicity Score_1"
  
  pdf("./绘图/57.CD8_T_毒性打分_额外基因.pdf",
      width = 4, height = 6.5)
  ggplot(data = CD8_T@meta.data[,c("Sample_rename",
                                   "Cytotoxicity Score_1")],
         aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
             fill = Sample_rename)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1,
                 color = "black", fill = "white") +
    theme_classic() +
    scale_fill_manual(values = sample_color) +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    labs(x = "", y = "Cytotoxicity Score") +
    ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                  c("OT-I", "OT-I/FX"),
                                                  c("OT-I", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/XCL1", "OT-I/FX"),
                                                  c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/FX", "OT-I/FX + CD40 mAb")
    ))
  dev.off()
  
  CD8_T_TCR_filter <- subset(CD8_T, subset = TCR_type %in% c("TRB and TRA"))
  CD8_T_TCR_filter$Chain <- All_TCR_filter_after_list_df[colnames(CD8_T_TCR_filter), "chain"]
  table(CD8_T_TCR_filter$Chain)
  CD8_T_TCR_filter$raw_clonotype_id <- All_TCR_filter_after_list_df[colnames(CD8_T_TCR_filter), "raw_clonotype_id"]
  CD8_T_TCR_filter$raw_clonotype_id <- unlist(lapply(CD8_T_TCR_filter$raw_clonotype_id,
                                                     function(x){
                                                       unique(unlist(strsplit(x, ";")))
                                                     }))
  clonotype_summary <- table(CD8_T_TCR_filter$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  clonotype_summary$Clonotype_group <- ""
  clonotype_summary$Clonotype_group[clonotype_summary$Freq == 1] <- "Single (0< X ≤1)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 1 & clonotype_summary$Freq <= 5] <- "Small (1< X ≤5)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 5 & clonotype_summary$Freq <= 10] <- "Medium (5< X ≤10)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 10 & clonotype_summary$Freq <= 30] <- "Large (10< X ≤30)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 30] <- "Hyperexpanded (30 < X)"
  table(clonotype_summary$Clonotype_group)
  rownames(clonotype_summary) <- clonotype_summary$Var1
  clonotype_summary$Clonotype
  CD8_T_TCR_filter$Clonotype <- clonotype_summary[CD8_T_TCR_filter$raw_clonotype_id, "Clonotype"]
  CD8_T_TCR_filter$Clonotype_group <- clonotype_summary[CD8_T_TCR_filter$raw_clonotype_id, "Clonotype_group"]
  unique(CD8_T_TCR_filter$Clonotype_group)
  CD8_T_TCR_filter$Clonotype_group <- factor(as.character(CD8_T_TCR_filter$Clonotype_group),
                                             levels = rev(c("Single (0< X ≤1)", "Small (1< X ≤5)",
                                                            "Medium (5< X ≤10)", "Large (10< X ≤30)",
                                                            "Hyperexpanded (30 < X)")))
  pdf("./绘图/58.CD8_T_Clonotype_group_percent_bar.pdf",
      height = 5.5, width = 6)
  percent_bar(CD8_T_TCR_filter, Ident = "Clonotype_group", Group = "Sample_rename",
              fill_color = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                 "#BF5968", "#580C16"))) +
    labs(y = "Clonotype frequency")
  dev.off()
  
  pdf("./绘图/59.CD8_T_Clonotype_group_UMAP.pdf",
      height = 5, width = 7.3)
  DimPlot(CD8_T_TCR_filter, group.by = "Clonotype_group", pt.size = 0.7) +
    scale_color_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                      "#BF5968", "#580C16"))) +
    ggtitle("CD8 T")
  dev.off()
  pdf("./绘图/60.CD8_T_Clonotype_group_UMAP_splited.pdf",
      height = 4.7, width = 15)
  DimPlot(CD8_T_TCR_filter, group.by = "Clonotype_group",
          split.by = "Sample_rename", pt.size = 0.7) +
    scale_color_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                      "#BF5968", "#580C16"))) +
    ggtitle("CD8 T")
  dev.off()
  
  CD8_T_TCR_filter_Naive <- subset(CD8_T_TCR_filter,
                                   subset = SubCelltype == "CD8 Naive/Tcm")
  meta <- CD8_T_TCR_filter_Naive@meta.data
  meta <- data.frame(Cells = rownames(meta),
                     Clonotype = meta$Clonotype,
                     Clonotype_group = meta$Clonotype_group,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta$Frequency <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- All_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$CDR3 <- All_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta <- meta[order(meta$Frequency,
                     decreasing = T),]
  openxlsx::write.xlsx(meta, "./绘图/61.CD8_T_naive_TCR.xlsx")
  meta <- meta[,-1]
  meta <- meta[!duplicated(meta),]
  openxlsx::write.xlsx(meta, "./绘图/62.CD8_T_naive_TCR_clonotype.xlsx")
  
  CD8_T_TCR_filter_expanded <- subset(CD8_T_TCR_filter,
                                      subset = Clonotype_group == "Hyperexpanded (30 < X)")
  pdf("./绘图/63.CD8_T_expanded_亚群UMAP.pdf",
      height = 4.5, width = 6.7)
  DimPlot(CD8_T_TCR_filter_expanded,
          group.by = "SubCelltype",
          cols = CD8_T_colors, pt.size = 0.7) +
    ggtitle("")
  dev.off()
  
  pdf("./绘图/64.CD8_T_Clonotype_group_percent_bar_splited.pdf",
      height = 5.3, width = 16)
  p <- ggplot(data = CD8_T_TCR_filter@meta.data,
              aes(x = SubCelltype, fill = Clonotype_group)) +
    geom_bar(stat = "count", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(".~Sample_rename", nrow = 1) +
    theme_bw() +
    labs(x = "", y = "Clonotype frequency", fill = "") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                     "#BF5968", "#580C16")))
  p
  dev.off()
  
  temp <- table(CD8_T_TCR_filter$Clonotype, CD8_T_TCR_filter$Sample_rename)
  temp <- as.data.frame.array(temp)
  temp_bi <- temp > 0
  shared_clonotype <- which(rowSums(temp_bi) == 4)
  shared_clonotype <- names(shared_clonotype)
  CD8_T_TCR_filter$Shared_clonotype <- ifelse(CD8_T_TCR_filter$Clonotype %in% shared_clonotype,
                                              "Shared clones", "Different clones")
  temp <- as.data.frame.array(table(CD8_T_TCR_filter$Sample_rename,
                                    CD8_T_TCR_filter$Shared_clonotype))
  temp <- as.data.frame(t(temp))
  df_all <- c()
  for (i in 1:4) {
    df <- data.frame(category = c("Different clones","Shared clones"),
                     count = temp[,i])
    df$fraction <- df$count / sum(df$count)
    df$ymax <- cumsum(df$fraction)
    df$ymin <- c(0, head(df$ymax, n = -1))
    df$Sample <- colnames(temp)[i]
    df_all <- as.data.frame(rbind(df_all,
                                  df))
  }
  
  df_all$category <- factor(df_all$category,
                            levels = c("Shared clones",
                                       "Different clones"))
  df_all$Sample <- factor(df_all$Sample,
                          levels = unique(df_all$Sample))
  pdf("./绘图/65.CD8_T_Shared_Clonotype_percent.pdf",
      width = 15, height = 4)
  ggplot(df_all, aes(ymax = ymax, ymin = ymin,
                     xmax = 4, xmin = 3)) +
    geom_rect(aes(fill = category),
              color = "black") +
    theme_bw() +
    xlim(2, 4) +
    facet_wrap(".~Sample", nrow = 1) +
    scale_fill_manual(values = c("#CD0000", "gray")) +
    coord_polar(theta = "y") +
    labs(fill = "") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold"))
  dev.off()
  
  clonotype_summary_2 <- clonotype_summary[clonotype_summary$Clonotype %in% shared_clonotype,]
  clonotype_summary_2 <- clonotype_summary_2[order(clonotype_summary_2$Freq,
                                                   decreasing = T),]
  first_clono <- clonotype_summary_2$Clonotype[1]
  second_clono <- clonotype_summary_2$Clonotype[2]
  third_clono <- clonotype_summary_2$Clonotype[3]
  
  CD8_T_TCR_filter$Top3 <- "Other"
  CD8_T_TCR_filter$Top3[CD8_T_TCR_filter$Clonotype %in% first_clono] <- "1st"
  CD8_T_TCR_filter$Top3[CD8_T_TCR_filter$Clonotype %in% second_clono] <- "2nd"
  CD8_T_TCR_filter$Top3[CD8_T_TCR_filter$Clonotype %in% third_clono] <- "3rd"
  pdf("./绘图/66.CD8_T_top3_frequent_shared_clonotypes.pdf",
      height = 4.7, width = 15)
  DimPlot(CD8_T_TCR_filter, group.by = "Top3",
          split.by = "Sample_rename", pt.size = 0.7) +
    scale_color_manual(values = c("#580C16", "#BF5968", "#F7BB7E",
                                  "gray")) +
    ggtitle("CD8 T") +
    labs(color = "Most frequent\nshared clonotypes")
  dev.off()
  
  meta <- CD8_T_TCR_filter@meta.data
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     "Clone Group" = meta$Clonotype_group,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- All_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- All_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- All_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- All_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- All_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- All_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- All_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/67.CD8_T_TCR_num.xlsx")
  
  meta$CDR3A_p13E <- 0
  meta$CDR3A_p13E[grep(pattern = "DYSNNRLT", x = meta$CDR3A)] <- 1
  meta$CDR3B_p13E <- 0
  meta$CDR3B_p13E[grep(pattern = "LELGG", x = meta$CDR3B)] <- 1
  meta$p13E <- meta$CDR3A_p13E + meta$CDR3B_p13E
  meta$p13E <- ifelse(meta$p13E > 0, "p13E-specific", "Unknown specificity")
  rownames(meta) <- meta$Cells
  CD8_T_TCR_filter$p13E <- meta[colnames(CD8_T_TCR_filter), "p13E"]
  
  meta$CDR3B_OT_I <- 0
  meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = meta$CDR3B)] <- 1
  meta$OT_I <- ifelse(meta$CDR3B_OT_I > 0, "OT-I specific", "Unknown specificity")
  CD8_T_TCR_filter$OT_I <- meta[colnames(CD8_T_TCR_filter), "OT_I"]
  {
    pdf("./绘图/68.CD8_T_毒性打分_没有额外基因_p13E-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$p13E == "p13E-specific",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/68.CD8_T_毒性打分_没有额外基因_p13E-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$p13E == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/69.CD8_T_毒性打分_额外基因_p13E-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$p13E == "p13E-specific",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/69.CD8_T_毒性打分_额外基因_p13E-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$p13E == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
  }
  
  {
    pdf("./绘图/70.CD8_T_毒性打分_没有额外基因_OT_I-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$OT_I == "OT-I specific",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/71.CD8_T_毒性打分_没有额外基因_OT_I-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$OT_I == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/72.CD8_T_毒性打分_额外基因_OT_I-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$OT_I == "OT-I specific",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/73.CD8_T_毒性打分_额外基因_OT_I-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD8_T_TCR_filter@meta.data[CD8_T_TCR_filter$OT_I == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
  }
}

# CD4 TCR
{
  CD40_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/CD40-T_2/outs/per_sample_outs/CD40-T_2/vdj_t/filtered_contig_annotations.csv")
  CD40_T_2_TCR$barcode <- paste0("CD40-T_2_", CD40_T_2_TCR$barcode)
  CD40_T_2_TCR$Sample <- "OT-I/FX + CD40 mAb"
  Vector_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/Vector-T_2/outs/per_sample_outs/Vector-T_2/vdj_t/filtered_contig_annotations.csv")
  Vector_T_2_TCR$barcode <- paste0("Vector-T_2_", Vector_T_2_TCR$barcode)
  Vector_T_2_TCR$Sample <- "OT-I"
  XCL1_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/XCL-T_2/outs/per_sample_outs/XCL-T_2/vdj_t/filtered_contig_annotations.csv")
  XCL1_T_2_TCR$barcode <- paste0("XCL1-T_2_", XCL1_T_2_TCR$barcode)
  XCL1_T_2_TCR$Sample <- "OT-I/XCL1"
  Fx_T_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/Fx-T_2/outs/per_sample_outs/Fx-T_2/vdj_t/filtered_contig_annotations.csv")
  Fx_T_2_TCR$barcode <- paste0("Fx-T_2_", Fx_T_2_TCR$barcode)
  Fx_T_2_TCR$Sample <- "OT-I/FX"
  All_TCR <- as.data.frame(rbind(CD40_T_2_TCR,
                                 Vector_T_2_TCR,
                                 XCL1_T_2_TCR,
                                 Fx_T_2_TCR))
  All_TCR_list <- split.data.frame(All_TCR,
                                   f = All_TCR$barcode)
  All_TCR_type <- c()
  for (i in 1:length(All_TCR_list)) {
    temp <- All_TCR_list[[i]]
    temp <- unique(temp$chain)
    if (length(temp) == 2) {
      All_TCR_type <- c(All_TCR_type,
                        "TRB and TRA")
      next
    }
    if (temp == "TRB") {
      All_TCR_type <- c(All_TCR_type,
                        "TRB Only")
      next
    }
    if (temp == "TRA") {
      All_TCR_type <- c(All_TCR_type,
                        "TRA Only")
      next
    }
  }
  table(All_TCR_type)
  All_TCR_type <- data.frame(Cells = names(All_TCR_list),
                             TCR_type = All_TCR_type)
  rownames(All_TCR_type) <- All_TCR_type$Cells
  All_TCR_remove_index <- c()
  All_TCR_list_filter_after <- c()
  for (i in 1:length(All_TCR_list)) {
    temp <- All_TCR_list[[i]]
    if (!(c("TRB") %in% temp$chain)) {
      All_TCR_remove_index <- c(All_TCR_remove_index,
                                i)
      next
    } else {
      temp <- filter(.data = temp,
                     umis == max(umis),
                     # reads == max(reads) & umis == max(umis),
                     .by = "chain")
      temp <- filter(.data = temp,
                     reads == max(reads),
                     # reads == max(reads) & umis == max(umis),
                     .by = "chain")
      temp <- as.data.frame(temp)
      All_TCR_list_filter_after <- c(All_TCR_list_filter_after,
                                     list(temp))
    }
  }
  All_TCR_filter_after <- data.table::rbindlist(All_TCR_list_filter_after)
  All_TCR_filter_after_list <- split.data.frame(All_TCR_filter_after,
                                                f = list(All_TCR_filter_after$barcode))
  for (i in 1:length(All_TCR_filter_after_list)) {
    temp <- All_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    All_TCR_filter_after_list[[i]] <- temp2
  }
  All_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(All_TCR_filter_after_list))
  All_TCR_filter_after_list_df$barcode <- unlist(lapply(All_TCR_filter_after_list_df$barcode,
                                                        function(x){
                                                          unlist(strsplit(x, ";"))[1]
                                                        }))
  rownames(All_TCR_filter_after_list_df) <- All_TCR_filter_after_list_df$barcode
  
  CD4_T$TCR_type <- All_TCR_type[colnames(CD4_T), 2]
  CD4_T$TCR_type[is.na(CD4_T$TCR_type)] <- "No recovery"
  CD4_T$TCR_type <- factor(CD4_T$TCR_type,
                           levels = c("TRB Only", "TRA Only",
                                      "TRB and TRA", "No recovery"))
  table(CD4_T$TCR_type)
  
  pdf("./绘图/74.CD4_T_TCR_type.pdf",
      height = 4.5, width = 6.2)
  DimPlot(CD4_T, group.by = "TCR_type", pt.size = 0.7) +
    scale_color_manual(values = c("#228B22", "#20B2AA",
                                  "#BA55D3", "#B5B5B5"))
  dev.off()
  
  CD4_T$Sample_rename <- All_seu@meta.data[colnames(CD4_T), "Sample_rename"]
  pdf("./绘图/75.CD4_T_TCR_type_percent_bar.pdf",
      height = 5.5, width = 6)
  percent_bar(CD4_T, Ident = "TCR_type", Group = "Sample_rename",
              fill_color = c("#228B22", "#20B2AA",
                             "#BA55D3", "#B5B5B5"))
  dev.off()
  
  # CD4_T <- AddModuleScore(CD4_T, features = list(c("Gzma", "Gzmb", "Gzmc",
  #                                                  "Gzmd", "Gzme", "Gzmf",
  #                                                  "Gzmg", "Gzmk", "Gzmm")))
  # colnames(CD4_T@meta.data)[59] <- "Cytotoxicity Score"
  
  pdf("./绘图/76.CD4_T_毒性打分_没有额外基因.pdf",
      width = 4, height = 6.5)
  ggplot(data = CD4_T@meta.data[,c("Sample_rename",
                                   "Cytotoxicity Score")],
         aes(x = Sample_rename, y = `Cytotoxicity Score`,
             fill = Sample_rename)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1,
                 color = "black", fill = "white") +
    theme_classic() +
    scale_fill_manual(values = sample_color) +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    labs(x = "") +
    ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                  c("OT-I", "OT-I/FX"),
                                                  c("OT-I", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/XCL1", "OT-I/FX"),
                                                  c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/FX", "OT-I/FX + CD40 mAb")
    ))
  dev.off()
  
  # CD4_T <- AddModuleScore(CD4_T, features = list(c("Gzma", "Gzmb", "Gzmc",
  #                                                  "Gzmd", "Gzme", "Gzmf",
  #                                                  "Gzmg", "Gzmk", "Gzmm",
  #                                                  "Ifng", "Tnf", "Prf1")),
  #                         search = TRUE)
  # colnames(CD4_T@meta.data)[60] <- "Cytotoxicity Score_1"
  
  pdf("./绘图/77.CD4_T_毒性打分_额外基因.pdf",
      width = 4, height = 6.5)
  ggplot(data = CD4_T@meta.data[,c("Sample_rename",
                                   "Cytotoxicity Score_1")],
         aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
             fill = Sample_rename)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1,
                 color = "black", fill = "white") +
    theme_classic() +
    scale_fill_manual(values = sample_color) +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    labs(x = "", y = "Cytotoxicity Score") +
    ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                  c("OT-I", "OT-I/FX"),
                                                  c("OT-I", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/XCL1", "OT-I/FX"),
                                                  c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                  c("OT-I/FX", "OT-I/FX + CD40 mAb")
    ))
  dev.off()
  
  CD4_T_TCR_filter <- subset(CD4_T, subset = TCR_type %in% c("TRB and TRA"))
  CD4_T_TCR_filter$Chain <- All_TCR_filter_after_list_df[colnames(CD4_T_TCR_filter), "chain"]
  table(CD4_T_TCR_filter$Chain)
  CD4_T_TCR_filter$raw_clonotype_id <- All_TCR_filter_after_list_df[colnames(CD4_T_TCR_filter), "raw_clonotype_id"]
  CD4_T_TCR_filter$raw_clonotype_id <- unlist(lapply(CD4_T_TCR_filter$raw_clonotype_id,
                                                     function(x){
                                                       unique(unlist(strsplit(x, ";")))
                                                     }))
  clonotype_summary <- table(CD4_T_TCR_filter$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  clonotype_summary$Clonotype_group <- ""
  clonotype_summary$Clonotype_group[clonotype_summary$Freq == 1] <- "Single (0< X ≤1)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 1 & clonotype_summary$Freq <= 5] <- "Small (1< X ≤5)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 5 & clonotype_summary$Freq <= 10] <- "Medium (5< X ≤10)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 10 & clonotype_summary$Freq <= 30] <- "Large (10< X ≤30)"
  clonotype_summary$Clonotype_group[clonotype_summary$Freq > 30] <- "Hyperexpanded (30 < X)"
  table(clonotype_summary$Clonotype_group)
  rownames(clonotype_summary) <- clonotype_summary$Var1
  clonotype_summary$Clonotype
  CD4_T_TCR_filter$Clonotype <- clonotype_summary[CD4_T_TCR_filter$raw_clonotype_id, "Clonotype"]
  CD4_T_TCR_filter$Clonotype_group <- clonotype_summary[CD4_T_TCR_filter$raw_clonotype_id, "Clonotype_group"]
  unique(CD4_T_TCR_filter$Clonotype_group)
  CD4_T_TCR_filter$Clonotype_group <- factor(as.character(CD4_T_TCR_filter$Clonotype_group),
                                             levels = rev(c("Single (0< X ≤1)", "Small (1< X ≤5)",
                                                            "Medium (5< X ≤10)", "Large (10< X ≤30)")))
  pdf("./绘图/78.CD4_T_Clonotype_group_percent_bar.pdf",
      height = 5.5, width = 6)
  percent_bar(CD4_T_TCR_filter, Ident = "Clonotype_group", Group = "Sample_rename",
              fill_color = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                 "#BF5968"))) +
    labs(y = "Clonotype frequency")
  dev.off()
  
  pdf("./绘图/79.CD4_T_Clonotype_group_UMAP.pdf",
      height = 5, width = 7.3)
  DimPlot(CD4_T_TCR_filter, group.by = "Clonotype_group", pt.size = 0.7) +
    scale_color_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                      "#BF5968"))) +
    ggtitle("CD4 T")
  dev.off()
  pdf("./绘图/80.CD4_T_Clonotype_group_UMAP_splited.pdf",
      height = 4.7, width = 15)
  DimPlot(CD4_T_TCR_filter, group.by = "Clonotype_group",
          split.by = "Sample_rename", pt.size = 0.7) +
    scale_color_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                      "#BF5968"))) +
    ggtitle("CD4 T")
  dev.off()
  
  pdf("./绘图/81.CD4_T_Clonotype_group_percent_bar_splited.pdf",
      height = 5.3, width = 16)
  p <- ggplot(data = CD4_T_TCR_filter@meta.data,
              aes(x = SubCelltype, fill = Clonotype_group)) +
    geom_bar(stat = "count", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(".~Sample_rename", nrow = 1) +
    theme_bw() +
    labs(x = "", y = "Clonotype frequency", fill = "") +
    theme(axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = rev(c("#154A7C", "#759DC5", "#F7BB7E",
                                     "#BF5968")))
  p
  dev.off()
  
  temp <- table(CD4_T_TCR_filter$Clonotype, CD4_T_TCR_filter$Sample_rename)
  temp <- as.data.frame.array(temp)
  temp_bi <- temp > 0
  shared_clonotype <- which(rowSums(temp_bi) == 4)
  shared_clonotype <- names(shared_clonotype)
  CD4_T_TCR_filter$Shared_clonotype <- ifelse(CD4_T_TCR_filter$Clonotype %in% shared_clonotype,
                                              "Shared clones", "Different clones")
  temp <- as.data.frame.array(table(CD4_T_TCR_filter$Sample_rename,
                                    CD4_T_TCR_filter$Shared_clonotype))
  temp <- as.data.frame(t(temp))
  df_all <- c()
  for (i in 1:4) {
    df <- data.frame(category = c("Different clones","Shared clones"),
                     count = temp[,i])
    df$fraction <- df$count / sum(df$count)
    df$ymax <- cumsum(df$fraction)
    df$ymin <- c(0, head(df$ymax, n = -1))
    df$Sample <- colnames(temp)[i]
    df_all <- as.data.frame(rbind(df_all,
                                  df))
  }
  
  df_all$category <- factor(df_all$category,
                            levels = c("Shared clones",
                                       "Different clones"))
  df_all$Sample <- factor(df_all$Sample,
                          levels = unique(df_all$Sample))
  pdf("./绘图/82.CD4_T_Shared_Clonotype_percent.pdf",
      width = 15, height = 4)
  ggplot(df_all, aes(ymax = ymax, ymin = ymin,
                     xmax = 4, xmin = 3)) +
    geom_rect(aes(fill = category),
              color = "black") +
    theme_bw() +
    xlim(2, 4) +
    facet_wrap(".~Sample", nrow = 1) +
    scale_fill_manual(values = c("#CD0000", "gray")) +
    coord_polar(theta = "y") +
    labs(fill = "") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.title = element_text(size = 13, color = "black"),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold"))
  dev.off()
  
  clonotype_summary_2 <- clonotype_summary[clonotype_summary$Clonotype %in% shared_clonotype,]
  clonotype_summary_2 <- clonotype_summary_2[order(clonotype_summary_2$Freq,
                                                   decreasing = T),]
  first_clono <- clonotype_summary_2$Clonotype[1]
  second_clono <- clonotype_summary_2$Clonotype[2]
  third_clono <- clonotype_summary_2$Clonotype[3]
  
  CD4_T_TCR_filter$Top3 <- "Other"
  CD4_T_TCR_filter$Top3[CD4_T_TCR_filter$Clonotype %in% first_clono] <- "1st"
  CD4_T_TCR_filter$Top3[CD4_T_TCR_filter$Clonotype %in% second_clono] <- "2nd"
  CD4_T_TCR_filter$Top3[CD4_T_TCR_filter$Clonotype %in% third_clono] <- "3rd"
  pdf("./绘图/83.CD4_T_top3_frequent_shared_clonotypes.pdf",
      height = 4.7, width = 15)
  DimPlot(CD4_T_TCR_filter, group.by = "Top3",
          split.by = "Sample_rename", pt.size = 0.7) +
    scale_color_manual(values = c("#580C16", "#BF5968", "#F7BB7E",
                                  "gray")) +
    ggtitle("CD4 T") +
    labs(color = "Most frequent\nshared clonotypes")
  dev.off()
  
  meta <- CD4_T_TCR_filter@meta.data
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     "Clone Group" = meta$Clonotype_group,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- All_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- All_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- All_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- All_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- All_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- All_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- All_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/84.CD4_T_TCR_num.xlsx")
  
  meta$CDR3A_p13E <- 0
  meta$CDR3A_p13E[grep(pattern = "DYSNNRLT", x = meta$CDR3A)] <- 1
  meta$CDR3B_p13E <- 0
  meta$CDR3B_p13E[grep(pattern = "LELGG", x = meta$CDR3B)] <- 1
  meta$p13E <- meta$CDR3A_p13E + meta$CDR3B_p13E
  meta$p13E <- ifelse(meta$p13E > 0, "p13E-specific", "Unknown specificity")
  rownames(meta) <- meta$Cells
  CD4_T_TCR_filter$p13E <- meta[colnames(CD4_T_TCR_filter), "p13E"]
  
  meta$CDR3B_OT_I <- 0
  meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = meta$CDR3B)] <- 1
  meta$OT_I <- ifelse(meta$CDR3B_OT_I > 0, "OT-I specific", "Unknown specificity")
  CD4_T_TCR_filter$OT_I <- meta[colnames(CD4_T_TCR_filter), "OT_I"]
  {
    pdf("./绘图/85.CD4_T_毒性打分_没有额外基因_p13E-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$p13E == "p13E-specific",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/86.CD4_T_毒性打分_没有额外基因_p13E-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$p13E == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/87.CD4_T_毒性打分_额外基因_p13E-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$p13E == "p13E-specific",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/88.CD4_T_毒性打分_额外基因_p13E-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$p13E == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
  }
  
  {
    pdf("./绘图/89.CD4_T_毒性打分_没有额外基因_OT_I-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$OT_I == "OT-I specific",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/90.CD4_T_毒性打分_没有额外基因_OT_I-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$OT_I == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score")],
           aes(x = Sample_rename, y = `Cytotoxicity Score`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/91.CD4_T_毒性打分_额外基因_OT_I-specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$OT_I == "OT-I specific",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
    
    pdf("./绘图/92.CD4_T_毒性打分_额外基因_OT_I-not_specific.pdf",
        width = 4, height = 6.5)
    ggplot(data = CD4_T_TCR_filter@meta.data[CD4_T_TCR_filter$OT_I == "Unknown specificity",
                                             c("Sample_rename", "Cytotoxicity Score_1")],
           aes(x = Sample_rename, y = `Cytotoxicity Score_1`,
               fill = Sample_rename)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1,
                   color = "black", fill = "white") +
      theme_classic() +
      scale_fill_manual(values = sample_color) +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Cytotoxicity Score") +
      ggpubr::stat_compare_means(comparisons = list(c("OT-I", "OT-I/XCL1"),
                                                    c("OT-I", "OT-I/FX"),
                                                    c("OT-I", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/XCL1", "OT-I/FX"),
                                                    c("OT-I/XCL1", "OT-I/FX + CD40 mAb"),
                                                    c("OT-I/FX", "OT-I/FX + CD40 mAb")
      ))
    dev.off()
  }
}

#############
NK <- readRDS("./注释/NK/NK_annotated.rds")
NK_colors <- ArchR::ArchRPalettes$kelly
names(NK_colors) <- NULL
NK <- subset(NK, subset = SubCelltype != "NK-antigen presentation")

pdf("./绘图/93.NK_UMAP.pdf",
    height = 4.5, width = 6)
DimPlot(NK, cols = NK_colors)
dev.off()

cluster_meta <- Embeddings(NK, reduction = "umap")
cluster_meta <- as.data.frame(cluster_meta)
cluster_meta$Cluster <- NK@meta.data[rownames(cluster_meta),"SubCelltype"]
cluster_meta <- aggregate.data.frame(cluster_meta[,1:2],
                                     by = list(cluster_meta$Cluster),
                                     FUN = median)
cluster_meta$Group.1 <- 0:6
Idents(Macro_Mono) <- Macro_Mono$SubCelltype
pdf("./绘图/94.NK_UMAP(with clusters).pdf",
    height = 4.5, width = 6)
DimPlot(NK, cols = NK_colors) +
  geom_point(data = cluster_meta,
             aes(x = UMAP_1, y = UMAP_2),
             size = 7, alpha = 0.7, color = "#EEE9E9") +
  geom_text(data = cluster_meta,
            aes(x = UMAP_1, y = UMAP_2, label = Group.1))
dev.off()

Idents(NK) <- NK$Sample
NK <- RenameIdents(NK,
                   "CD40-T_2" = "OT-I/FX + CD40 mAb",
                   "Fx-T_2" = "OT-I/FX",
                   "Vector-T_2" = "OT-I",
                   "XCL1-T_2" = "OT-I/XCL1")
NK$Sample_rename <- factor(as.character(Idents(NK)),
                           levels = c("OT-I", "OT-I/XCL1",
                                      "OT-I/FX", "OT-I/FX + CD40 mAb"))
pdf("./绘图/95.NK_细胞类型比例_perSample.pdf",
    width = 5, height = 5)
percent_bar(sce = NK, Ident = "SubCelltype",
            Group = "Sample_rename",
            fill_color = NK_colors)
dev.off()

pdf("./绘图/96.NK_细胞类型比例_perCluster.pdf",
    width = 7, height = 5)
percent_bar(sce = NK, Ident = "Sample_rename",
            Group = "SubCelltype",
            fill_color = sample_color)
dev.off()

Idents(NK) <- factor(NK$SubCelltype,
                     levels = sort(unique(NK$SubCelltype)))
features <- list("NK-1" = c("Plac8", "Ctla2a", "Fosl2", "Vegfa", "Sult2b1", "Tigit", "Ier5l"),
                 "NK-2" = c("Ltb", "Ly6e", "Emb", "Ly6c2"),
                 "NK-3" = c("Hsp90ab1", "Rps17", "Rps20", "Rpl12", "Rps2"),
                 "NK-4" = c("Ccl4", "Ccl3", "Nr4a1", "Icam1", "Ifng", "Nfkbid", "Nfkbiz", "Nfkbia"),
                 "NK-5" = c("Ccl5", "S100a6", "Lgals1", "Klrg1", "Gzma"),
                 "NKT" = c("Il7r", "Ifi27l2a", "Ly6a"))
p <- DotPlot(NK, features = features) +
  theme_classic() +
  # scale_color_gradientn(colours = c("#008B8B", "#3CB371", "#FFA500", "#FFB90F", "#FF1493", "#FF0000", "#CD2626")) +
  scale_color_gradientn(colours = c("#008B8B", "#00CD66", "#FFA500", "#CD0000")) +
  # scale_y_discrete(position = "right") + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, colour = "black")) +
  labs(x = "", y = "")
pdf("./绘图/97.NK_markers_dotplot.pdf",
    width = 11, height = 3.5)
p #%>%
# aplot::insert_top(bottom_colorbar,
#                   height = 0.04)
dev.off()

library(dplyr)
NK_DEGs <- FindAllMarkers(NK, only.pos = T)
NK_DEGs %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
top20 <- as.data.frame(top20)
NK <- ScaleData(NK, features = rownames(NK))
pdf("./绘图/98.NK_subtype_top20_deg_heatmap.pdf",
    width = 8.5, height = 14)
DoHeatmap(NK, features = top20$gene,
          group.colors = NK_colors,
          label = FALSE,
          disp.min = -2,
          disp.max = 2) +
  scale_fill_gradientn(colors = c("#00008B","white","#B22222"))
dev.off()

NK_DEGs %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
top10 <- as.data.frame(top10)
pdf("./绘图/99.NK_subtype_top10_deg_heatmap.pdf",
    width = 8.5, height = 9)
DoHeatmap(NK, features = top10$gene,
          group.colors = NK_colors,
          label = FALSE,
          disp.min = -2,
          disp.max = 2) +
  scale_fill_gradientn(colors = c("#00008B","white","#B22222"))
dev.off()

Biocarta <- read.gmt("m2.cp.biocarta.v2023.2.Mm.symbols.gmt")
Biocarta$term <- as.character(Biocarta$term)
Biocarta$term <- unlist(lapply(Biocarta$term, function(x){
  temp <- unlist(strsplit(x, split = "_", fixed = T))
  temp <- temp[-1]
  temp <- Hmisc::capitalize(tolower(temp))
  paste0(temp, collapse = " ")
}))

NK_DEGs_list <- split.data.frame(NK_DEGs, f = NK_DEGs$cluster)
NK_DEGs_GO <- c()
for (i in 1:length(NK_DEGs_list)) {
  temp <- NK_DEGs_list[[i]]
  temp <- enricher(gene = temp$gene, pvalueCutoff = 1,
                   qvalueCutoff = 1, TERM2GENE = Biocarta)
  temp <- temp@result
  temp <- temp[temp$pvalue < 0.05,]
  NK_DEGs_GO <- c(NK_DEGs_GO,
                  list(temp))
}
names(NK_DEGs_GO) <- names(NK_DEGs_list)
openxlsx::write.xlsx(NK_DEGs_GO,
                     "./绘图/100.NK_subtype_GO.xlsx")

for (i in 1:length(NK_DEGs_GO)) {
  NK_DEGs_GO[[i]]$Cluster <- names(NK_DEGs_GO)[i]
}
temp <- data.table::rbindlist(NK_DEGs_GO)
temp$log10P <- -log10(temp$pvalue)
temp %>%
  group_by(Cluster) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
top5 <- as.data.frame(top5)
sum(duplicated(top5$Description))
pathway <- top5$Description
pathway <- pathway[!duplicated(pathway)]
pathway_df <- top5[top5$Description %in% pathway,]
pathway_df <- pathway_df[,c("ID", "Description")]
pathway_df <- pathway_df[!duplicated(pathway_df),]
top5$Description <- factor(top5$Description,
                           levels = rev(pathway_df$Description))
GeneRatio2 <- c()
for (i in 1:nrow(top5)) {
  temp2 <- unlist(strsplit(top5$GeneRatio[i], "/", fixed = T))
  GeneRatio2 <- c(GeneRatio2,
                  as.numeric(temp2[1]) / as.numeric(temp2[2]))
}
top5$GeneRatio2 <- GeneRatio2

p <- ggplot() +
  geom_point(data = top5, aes(x = Cluster, y = Description,
                              fill = log10P, size = GeneRatio2),
             color = "black", shape = 21) +
  scale_fill_gradient(low = "#EE9A00", high = "#CD2626") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black",
                                   angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "", fill = "-log10(P-value)",
       size = "GeneRatio", y = "")

pdf("./绘图/101.NK_subtype_GO.pdf",
    height = 6, width = 6)
p
dev.off()




#############
XCL1_tdLN_2_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/2/Fx_tdLN_2/outs/per_sample_outs/Fx_tdLN_2/vdj_t/filtered_contig_annotations.csv")
XCL1_tdLN_2_TCR$barcode <- paste0("XCL1-tdLN_2_", XCL1_tdLN_2_TCR$barcode)
XCL1_tdLN_2_TCR$Sample <- "OT_I/XCL1 tdLN"
XCL1_tdLN_2_TCR$DataSet <- "Data 2"
XCL1_tdLN_2_TCR_list <- split.data.frame(XCL1_tdLN_2_TCR,
                                       f = XCL1_tdLN_2_TCR$barcode)
XCL1_tdLN_2_remove_index <- c()
XCL1_tdLN_2_TCR_list_filter_after <- c()
for (i in 1:length(XCL1_tdLN_2_TCR_list)) {
  temp <- XCL1_tdLN_2_TCR_list[[i]]
  if (sum(c("TRB","TRA") %in% temp$chain) < 2) {
    XCL1_tdLN_2_remove_index <- c(XCL1_tdLN_2_remove_index,
                                  i)
    next
  } else {
    temp <- dplyr::filter(.data = temp,
                   umis == max(umis),
                   # reads == max(reads) & umis == max(umis),
                   .by = "chain")
    temp <- dplyr::filter(.data = temp,
                   reads == max(reads),
                   # reads == max(reads) & umis == max(umis),
                   .by = "chain")
    temp <- as.data.frame(temp)
    XCL1_tdLN_2_TCR_list_filter_after <- c(XCL1_tdLN_2_TCR_list_filter_after,
                                         list(temp))
  }
}
XCL1_tdLN_2_TCR_filter_after <- data.table::rbindlist(XCL1_tdLN_2_TCR_list_filter_after)

FX_tdLN_1_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/1/FX-TCR_1/outs/per_sample_outs/FX-TCR_1/vdj_t/filtered_contig_annotations.csv")
FX_tdLN_1_TCR$barcode <- paste0("FX-tdLN_1_", FX_tdLN_1_TCR$barcode)
FX_tdLN_1_TCR$Sample <- "OT-I/FX tdLN"
FX_tdLN_1_TCR$DataSet <- "Data 1"
FX_tdLN_1_TCR_list <- split.data.frame(FX_tdLN_1_TCR,
                                       f = FX_tdLN_1_TCR$barcode)
FX_tdLN_1_remove_index <- c()
FX_tdLN_1_TCR_list_filter_after <- c()
for (i in 1:length(FX_tdLN_1_TCR_list)) {
  temp <- FX_tdLN_1_TCR_list[[i]]
  if (sum(c("TRB","TRA") %in% temp$chain) < 2) {
    FX_tdLN_1_remove_index <- c(FX_tdLN_1_remove_index,
                                i)
    next
  } else {
    temp <- dplyr::filter(.data = temp,
                          umis == max(umis),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- dplyr::filter(.data = temp,
                          reads == max(reads),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- as.data.frame(temp)
    FX_tdLN_1_TCR_list_filter_after <- c(FX_tdLN_1_TCR_list_filter_after,
                                           list(temp))
  }
}
FX_tdLN_1_TCR_filter_after <- data.table::rbindlist(FX_tdLN_1_TCR_list_filter_after)

Vector_tdLN_1_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/1/Vector-TCR_1/outs/per_sample_outs/Vector-TCR_1/vdj_t/filtered_contig_annotations.csv")
Vector_tdLN_1_TCR$barcode <- paste0("Vector-tdLN_1_", Vector_tdLN_1_TCR$barcode)
Vector_tdLN_1_TCR$Sample <- "OT-I tdLN"
Vector_tdLN_1_TCR$DataSet <- "Data 1"
Vector_tdLN_1_TCR_list <- split.data.frame(Vector_tdLN_1_TCR,
                                       f = Vector_tdLN_1_TCR$barcode)
Vector_tdLN_1_remove_index <- c()
Vector_tdLN_1_TCR_list_filter_after <- c()
for (i in 1:length(Vector_tdLN_1_TCR_list)) {
  temp <- Vector_tdLN_1_TCR_list[[i]]
  if (sum(c("TRB","TRA") %in% temp$chain) < 2) {
    Vector_tdLN_1_remove_index <- c(Vector_tdLN_1_remove_index,
                                    i)
    next
  } else {
    temp <- dplyr::filter(.data = temp,
                          umis == max(umis),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- dplyr::filter(.data = temp,
                          reads == max(reads),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- as.data.frame(temp)
    Vector_tdLN_1_TCR_list_filter_after <- c(Vector_tdLN_1_TCR_list_filter_after,
                                         list(temp))
  }
}
Vector_tdLN_1_TCR_filter_after <- data.table::rbindlist(Vector_tdLN_1_TCR_list_filter_after)

CD40_tdLN_1_TCR <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR/1.Cellranger/1/CD40-TCR_1/outs/per_sample_outs/CD40-TCR_1/vdj_t/filtered_contig_annotations.csv")
CD40_tdLN_1_TCR$barcode <- paste0("CD40-tdLN_1_", CD40_tdLN_1_TCR$barcode)
CD40_tdLN_1_TCR$Sample <- "OT-I/FX + CD40 mAb tdLN"
CD40_tdLN_1_TCR$DataSet <- "Data 1"
CD40_tdLN_1_TCR_list <- split.data.frame(CD40_tdLN_1_TCR,
                                           f = CD40_tdLN_1_TCR$barcode)
CD40_tdLN_1_remove_index <- c()
CD40_tdLN_1_TCR_list_filter_after <- c()
for (i in 1:length(CD40_tdLN_1_TCR_list)) {
  temp <- CD40_tdLN_1_TCR_list[[i]]
  if (sum(c("TRB","TRA") %in% temp$chain) < 2) {
    CD40_tdLN_1_remove_index <- c(CD40_tdLN_1_remove_index,
                                  i)
    next
  } else {
    temp <- dplyr::filter(.data = temp,
                          umis == max(umis),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- dplyr::filter(.data = temp,
                          reads == max(reads),
                          # reads == max(reads) & umis == max(umis),
                          .by = "chain")
    temp <- as.data.frame(temp)
    CD40_tdLN_1_TCR_list_filter_after <- c(CD40_tdLN_1_TCR_list_filter_after,
                                           list(temp))
  }
}
CD40_tdLN_1_TCR_filter_after <- data.table::rbindlist(CD40_tdLN_1_TCR_list_filter_after)


{
  CD40_tdLN_1_TCR_filter_after_list <- split.data.frame(CD40_tdLN_1_TCR_filter_after,
                                                        f = list(CD40_tdLN_1_TCR_filter_after$barcode))
  for (i in 1:length(CD40_tdLN_1_TCR_filter_after_list)) {
    temp <- CD40_tdLN_1_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    CD40_tdLN_1_TCR_filter_after_list[[i]] <- temp2
  }
  CD40_tdLN_1_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(CD40_tdLN_1_TCR_filter_after_list))
  CD40_tdLN_1_TCR_filter_after_list_df$barcode <- unlist(lapply(CD40_tdLN_1_TCR_filter_after_list_df$barcode,
                                                                function(x){
                                                                  unlist(strsplit(x, ";"))[1]
                                                                }))
  rownames(CD40_tdLN_1_TCR_filter_after_list_df) <- CD40_tdLN_1_TCR_filter_after_list_df$barcode
}

{
  Vector_tdLN_1_TCR_filter_after_list <- split.data.frame(Vector_tdLN_1_TCR_filter_after,
                                                        f = list(Vector_tdLN_1_TCR_filter_after$barcode))
  for (i in 1:length(Vector_tdLN_1_TCR_filter_after_list)) {
    temp <- Vector_tdLN_1_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    Vector_tdLN_1_TCR_filter_after_list[[i]] <- temp2
  }
  Vector_tdLN_1_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(Vector_tdLN_1_TCR_filter_after_list))
  Vector_tdLN_1_TCR_filter_after_list_df$barcode <- unlist(lapply(Vector_tdLN_1_TCR_filter_after_list_df$barcode,
                                                                function(x){
                                                                  unlist(strsplit(x, ";"))[1]
                                                                }))
  rownames(Vector_tdLN_1_TCR_filter_after_list_df) <- Vector_tdLN_1_TCR_filter_after_list_df$barcode
}

{
  FX_tdLN_1_TCR_filter_after_list <- split.data.frame(FX_tdLN_1_TCR_filter_after,
                                                          f = list(FX_tdLN_1_TCR_filter_after$barcode))
  for (i in 1:length(FX_tdLN_1_TCR_filter_after_list)) {
    temp <- FX_tdLN_1_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    FX_tdLN_1_TCR_filter_after_list[[i]] <- temp2
  }
  FX_tdLN_1_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(FX_tdLN_1_TCR_filter_after_list))
  FX_tdLN_1_TCR_filter_after_list_df$barcode <- unlist(lapply(FX_tdLN_1_TCR_filter_after_list_df$barcode,
                                                                  function(x){
                                                                    unlist(strsplit(x, ";"))[1]
                                                                  }))
  rownames(FX_tdLN_1_TCR_filter_after_list_df) <- FX_tdLN_1_TCR_filter_after_list_df$barcode
}

{
  XCL1_tdLN_2_TCR_filter_after_list <- split.data.frame(XCL1_tdLN_2_TCR_filter_after,
                                                      f = list(XCL1_tdLN_2_TCR_filter_after$barcode))
  for (i in 1:length(XCL1_tdLN_2_TCR_filter_after_list)) {
    temp <- XCL1_tdLN_2_TCR_filter_after_list[[i]]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$chain, decreasing = T),]
    temp2 <- c()
    for (k in 1:ncol(temp)) {
      temp2 <- c(temp2,
                 paste0(temp[,k], collapse = ";"))
    }
    names(temp2) <- colnames(temp)
    temp2 <- as.data.frame(temp2)
    temp2 <- as.data.frame(t(temp2))
    XCL1_tdLN_2_TCR_filter_after_list[[i]] <- temp2
  }
  XCL1_tdLN_2_TCR_filter_after_list_df <- as.data.frame(data.table::rbindlist(XCL1_tdLN_2_TCR_filter_after_list))
  XCL1_tdLN_2_TCR_filter_after_list_df$barcode <- unlist(lapply(XCL1_tdLN_2_TCR_filter_after_list_df$barcode,
                                                              function(x){
                                                                unlist(strsplit(x, ";"))[1]
                                                              }))
  rownames(XCL1_tdLN_2_TCR_filter_after_list_df) <- XCL1_tdLN_2_TCR_filter_after_list_df$barcode
}


{
  clonotype_summary <- table(XCL1_tdLN_2_TCR_filter_after_list_df$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  rownames(clonotype_summary) <- clonotype_summary$Var1
  XCL1_tdLN_2_TCR_filter_after_list_df$Clonotype <- clonotype_summary[XCL1_tdLN_2_TCR_filter_after_list_df$raw_clonotype_id, "Clonotype"]
  meta <- XCL1_tdLN_2_TCR_filter_after_list_df
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- XCL1_tdLN_2_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/102.XCL1_tdLN_DataSet2_TCR_num.xlsx")
  
}

{
  clonotype_summary <- table(FX_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  rownames(clonotype_summary) <- clonotype_summary$Var1
  FX_tdLN_1_TCR_filter_after_list_df$Clonotype <- clonotype_summary[FX_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id, "Clonotype"]
  meta <- FX_tdLN_1_TCR_filter_after_list_df
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- FX_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/102.FX_tdLN_DataSet1_TCR_num.xlsx")
  
}

{
  clonotype_summary <- table(Vector_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  rownames(clonotype_summary) <- clonotype_summary$Var1
  Vector_tdLN_1_TCR_filter_after_list_df$Clonotype <- clonotype_summary[Vector_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id, "Clonotype"]
  meta <- Vector_tdLN_1_TCR_filter_after_list_df
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- Vector_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/102.Vector_tdLN_DataSet1_TCR_num.xlsx")
}

{
  clonotype_summary <- table(CD40_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id)
  clonotype_summary <- as.data.frame(clonotype_summary)
  clonotype_summary <- clonotype_summary[order(clonotype_summary$Freq,
                                               decreasing = T),]
  clonotype_summary$Clonotype <- paste0("Clonotype_", 1:nrow(clonotype_summary))
  rownames(clonotype_summary) <- clonotype_summary$Var1
  CD40_tdLN_1_TCR_filter_after_list_df$Clonotype <- clonotype_summary[CD40_tdLN_1_TCR_filter_after_list_df$raw_clonotype_id, "Clonotype"]
  meta <- CD40_tdLN_1_TCR_filter_after_list_df
  meta <- data.frame(Cells = rownames(meta),
                     "Clone Name" = meta$Clonotype,
                     raw_clonotype_id = meta$raw_clonotype_id)
  meta[,"Clone Size"] <- clonotype_summary[meta$raw_clonotype_id, "Freq"]
  meta$Chain <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "chain"]
  meta$TRAV <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRAV <- unlist(lapply(meta$TRAV, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$CDR3A <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3A <- unlist(lapply(meta$CDR3A, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRAJ <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRAJ <- unlist(lapply(meta$TRAJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[2]
  }))
  meta$TRBV <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "v_gene"]
  meta$TRBV <- unlist(lapply(meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$CDR3B <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "cdr3"]
  meta$CDR3B <- unlist(lapply(meta$CDR3B, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta$TRBJ <- CD40_tdLN_1_TCR_filter_after_list_df[meta$Cells, "j_gene"]
  meta$TRBJ <- unlist(lapply(meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  meta <- meta[order(meta$`Clone Size`,
                     decreasing = T),]
  meta$index <- paste0(meta$TRAV,"|",meta$CDR3A,"|",meta$TRAJ,"|",
                       meta$TRBV,"|",meta$CDR3B,"|",meta$TRBJ)
  meta$index_num <- table(meta$index)[meta$index]
  meta$index_num <- as.numeric(meta$index_num)
  meta_2 <- meta[!duplicated(meta[,"index"]),]
  meta_2 <- meta_2[order(meta_2$index_num,
                         decreasing = TRUE),]
  meta_2 <- data.frame("Clone Name" = meta_2$`Clone.Name`,
                       "Clone Size" = meta_2$index_num,
                       "TRAV" = meta_2$TRAV,
                       "CDR3A" = meta_2$CDR3A,
                       "TRAJ" = meta_2$TRAJ,
                       "TRBV" = meta_2$TRBV,
                       "CDR3B" = meta_2$CDR3B,
                       "TRBJ" = meta_2$TRBJ,
                       check.rows = F, check.names = F)
  openxlsx::write.xlsx(meta_2, "./绘图/102.CD40_tdLN_DataSet1_TCR_num.xlsx")
}

## OT-I percent
CD40_tdLN_1_meta <- CD40_tdLN_1_TCR_filter_after_list_df
CD40_tdLN_1_meta$CDR3B <- CD40_tdLN_1_meta$cdr3
CD40_tdLN_1_meta$CDR3B <- unlist(lapply(CD40_tdLN_1_meta$CDR3B, function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
}))
CD40_tdLN_1_meta$CDR3B_OT_I <- 0
CD40_tdLN_1_meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = CD40_tdLN_1_meta$CDR3B)] <- 1
CD40_tdLN_1_meta$OT_I <- ifelse(CD40_tdLN_1_meta$CDR3B_OT_I > 0, "OT-I TCRs", "Host (non-OT-I) TCRs")
table(CD40_tdLN_1_meta$OT_I)

Vector_tdLN_1_meta <- Vector_tdLN_1_TCR_filter_after_list_df
Vector_tdLN_1_meta$CDR3B <- Vector_tdLN_1_meta$cdr3
Vector_tdLN_1_meta$CDR3B <- unlist(lapply(Vector_tdLN_1_meta$CDR3B, function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
}))
Vector_tdLN_1_meta$CDR3B_OT_I <- 0
Vector_tdLN_1_meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = Vector_tdLN_1_meta$CDR3B)] <- 1
Vector_tdLN_1_meta$OT_I <- ifelse(Vector_tdLN_1_meta$CDR3B_OT_I > 0, "OT-I TCRs", "Host (non-OT-I) TCRs")
table(Vector_tdLN_1_meta$OT_I)

FX_tdLN_1_meta <- FX_tdLN_1_TCR_filter_after_list_df
FX_tdLN_1_meta$CDR3B <- FX_tdLN_1_meta$cdr3
FX_tdLN_1_meta$CDR3B <- unlist(lapply(FX_tdLN_1_meta$CDR3B, function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
}))
FX_tdLN_1_meta$CDR3B_OT_I <- 0
FX_tdLN_1_meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = FX_tdLN_1_meta$CDR3B)] <- 1
FX_tdLN_1_meta$OT_I <- ifelse(FX_tdLN_1_meta$CDR3B_OT_I > 0, "OT-I TCRs", "Host (non-OT-I) TCRs")
table(FX_tdLN_1_meta$OT_I)

XCL1_tdLN_2_meta <- XCL1_tdLN_2_TCR_filter_after_list_df
XCL1_tdLN_2_meta$CDR3B <- XCL1_tdLN_2_meta$cdr3
XCL1_tdLN_2_meta$CDR3B <- unlist(lapply(XCL1_tdLN_2_meta$CDR3B, function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
}))
XCL1_tdLN_2_meta$CDR3B_OT_I <- 0
XCL1_tdLN_2_meta$CDR3B_OT_I[grep(pattern = "CASSRANYEQYF", x = XCL1_tdLN_2_meta$CDR3B)] <- 1
XCL1_tdLN_2_meta$OT_I <- ifelse(XCL1_tdLN_2_meta$CDR3B_OT_I > 0, "OT-I TCRs", "Host (non-OT-I) TCRs")
table(XCL1_tdLN_2_meta$OT_I)

tdLN <- as.data.frame(rbind(CD40_tdLN_1_meta,
                            Vector_tdLN_1_meta,
                            FX_tdLN_1_meta,
                            XCL1_tdLN_2_meta))

temp <- table(tdLN$Sample, tdLN$OT_I)
temp <- as.data.frame.array(temp)
temp <- temp / rowSums(temp)
temp <- data.frame(Sample = unlist(lapply(rownames(temp), function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
})), temp, check.rows = F, check.names = F)
temp <- reshape::melt.data.frame(temp, id.vars = "Sample")
temp$Sample <- factor(temp$Sample, levels = c("OT-I tdLN", "OT_I/XCL1 tdLN",
                                              "OT-I/FX tdLN", "OT-I/FX + CD40 mAb tdLN"))
temp$variable <- factor(temp$variable,
                        levels = c("OT-I TCRs", "Host (non-OT-I) TCRs"))
p <- ggplot(data = temp, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(x = "", y = "Percentage", fill = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#000080", "gray"))
pdf("./绘图/103.tdLN_sample_OT-I_percent_bar.pdf",
    height = 4.5, width = 5.5)
p
dev.off()

Samples_tdLN <- split.data.frame(tdLN, f = list(tdLN$Sample))
CD40_tdLN_1_meta <- Samples_tdLN$`OT-I/FX + CD40 mAb tdLN`
# CD40_tdLN_1_meta <- CD40_tdLN_1_meta[CD40_tdLN_1_meta$OT_I == "Host (non-OT-I) TCRs",]
temp <- names(table(CD40_tdLN_1_meta$Clonotype))[table(CD40_tdLN_1_meta$Clonotype) > 1]
CD40_tdLN_1_meta <- CD40_tdLN_1_meta[CD40_tdLN_1_meta$Clonotype %in% temp,]
sort(table(CD40_tdLN_1_meta$Clonotype), decreasing = T)


# CD40_tdLN_1_meta <- CD40_tdLN_1_meta[1:floor(nrow(CD40_tdLN_1_meta)*0.15),]
Vector_tdLN_1_meta <- Samples_tdLN$`OT-I tdLN`
# Vector_tdLN_1_meta <- Vector_tdLN_1_meta[Vector_tdLN_1_meta$OT_I == "Host (non-OT-I) TCRs",]
# temp <- names(table(Vector_tdLN_1_meta$Clonotype))[table(Vector_tdLN_1_meta$Clonotype) > 1]
Vector_tdLN_1_meta <- Vector_tdLN_1_meta[Vector_tdLN_1_meta$Clonotype %in% temp,]

FX_tdLN_1_meta <- Samples_tdLN$`OT-I/FX tdLN`
FX_tdLN_1_meta <- FX_tdLN_1_meta[FX_tdLN_1_meta$OT_I == "Host (non-OT-I) TCRs",]
temp <- names(table(FX_tdLN_1_meta$Clonotype))[table(FX_tdLN_1_meta$Clonotype) > 1]
FX_tdLN_1_meta <- FX_tdLN_1_meta[FX_tdLN_1_meta$Clonotype %in% temp,]
sort(table(FX_tdLN_1_meta$Clonotype), decreasing = T)

XCL1_tdLN_2_meta <- Samples_tdLN$`OT_I/XCL1 tdLN`
XCL1_tdLN_2_meta <- XCL1_tdLN_2_meta[XCL1_tdLN_2_meta$OT_I == "Host (non-OT-I) TCRs",]
temp <- names(table(XCL1_tdLN_2_meta$Clonotype))[table(XCL1_tdLN_2_meta$Clonotype) > 1]
XCL1_tdLN_2_meta <- XCL1_tdLN_2_meta[XCL1_tdLN_2_meta$Clonotype %in% temp,]
sort(table(XCL1_tdLN_2_meta$Clonotype), decreasing = T)



# cumulative_sum <- cumsum(sort(table(CD40_tdLN_1_meta$Clonotype), decreasing = T))
# cumulative_proportion <- cumulative_sum / sum(sort(table(CD40_tdLN_1_meta$Clonotype), decreasing = T))
# CD40_tdLN_1_meta <- CD40_tdLN_1_meta[CD40_tdLN_1_meta$Clonotype %in% names(cumulative_proportion)[cumulative_proportion <= 0.15],]
CD40_tdLN_1_meta$Top10 <- CD40_tdLN_1_meta$Clonotype
index <- !(CD40_tdLN_1_meta$Top10 %in% c("Clonotype_2", "Clonotype_3", "Clonotype_4",
                                             "Clonotype_5", "Clonotype_6", "Clonotype_7",
                                             "Clonotype_8", "Clonotype_10", "Clonotype_11",
                                             "Clonotype_12"))
CD40_tdLN_1_meta$Top10[index] <- "Other"
CD40_tdLN_1_meta$Top10 <- factor(as.character(CD40_tdLN_1_meta$Top10),
                                 levels = c("Clonotype_2", "Clonotype_3", "Clonotype_4",
                                            "Clonotype_5", "Clonotype_6", "Clonotype_7",
                                            "Clonotype_8", "Clonotype_10", "Clonotype_11",
                                            "Clonotype_12", "Other"))

Vector_tdLN_1_meta$Top10 <- Vector_tdLN_1_meta$Clonotype
index <- !(Vector_tdLN_1_meta$Top10 %in% c("Clonotype_1", "Clonotype_2", "Clonotype_3",
                                         "Clonotype_4", "Clonotype_5", "Clonotype_6",
                                         "Clonotype_7", "Clonotype_8", "Clonotype_9",
                                         "Clonotype_10"))
Vector_tdLN_1_meta$Top10[index] <- "Other"

FX_tdLN_1_meta$Top10 <- FX_tdLN_1_meta$Clonotype
index <- !(FX_tdLN_1_meta$Top10 %in% c("Clonotype_1", "Clonotype_2", "Clonotype_3",
                                           "Clonotype_4", "Clonotype_5", "Clonotype_6",
                                           "Clonotype_7", "Clonotype_8", "Clonotype_9",
                                           "Clonotype_10"))
FX_tdLN_1_meta$Top10[index] <- "Other"

XCL1_tdLN_2_meta$Top10 <- XCL1_tdLN_2_meta$Clonotype
index <- !(XCL1_tdLN_2_meta$Top10 %in% c("Clonotype_1", "Clonotype_2", "Clonotype_3",
                                       "Clonotype_4", "Clonotype_5", "Clonotype_6",
                                       "Clonotype_7", "Clonotype_8", "Clonotype_9",
                                       "Clonotype_10"))
XCL1_tdLN_2_meta$Top10[index] <- "Other"

tdLN <- as.data.frame(rbind(CD40_tdLN_1_meta,
                            Vector_tdLN_1_meta,
                            FX_tdLN_1_meta,
                            XCL1_tdLN_2_meta))
tdLN$Sample <- unlist(lapply(tdLN$Sample, function(x){
  unlist(strsplit(x, ";", fixed = T))[1]
}))
tdLN$Sample <- factor(tdLN$Sample, levels = c("OT-I tdLN", "OT_I/XCL1 tdLN",
                                              "OT-I/FX tdLN", "OT-I/FX + CD40 mAb tdLN"))
tdLN$OT_I <- factor(tdLN$OT_I,
                        levels = c("OT-I TCRs", "Host (non-OT-I) TCRs"))
tdLN$Top10 <- factor(tdLN$Top10,
                     levels = c("Clonotype_1", "Clonotype_2", "Clonotype_3",
                                "Clonotype_4", "Clonotype_5", "Clonotype_6",
                                "Clonotype_7", "Clonotype_8", "Clonotype_9",
                                "Clonotype_10", "Other"))
p <- ggplot(data = tdLN, aes(x = Sample, fill = Top10)) +
  geom_bar(stat = "count", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(x = "", y = "Percentage", fill = "") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c(scales::hue_pal()(10), "gray"))
pdf("./绘图/104.tdLN_sample_top_10_percent_bar.pdf",
    height = 4.5, width = 4.5)
p
dev.off()


CD40_tdLN_1_meta
Vector_tdLN_1_meta
FX_tdLN_1_meta
XCL1_tdLN_2_meta

{
  CD40_tdLN_1_meta$TRBV <- CD40_tdLN_1_meta$v_gene
  CD40_tdLN_1_meta$TRBV <- unlist(lapply(CD40_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  CD40_tdLN_1_meta$TRBJ <- CD40_tdLN_1_meta$j_gene
  CD40_tdLN_1_meta$TRBJ <- unlist(lapply(CD40_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  TRBV <- unlist(lapply(CD40_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBV <- sort(unique(TRBV))
  TRBJ <- unlist(lapply(CD40_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBJ <- sort(unique(TRBJ))
  temp <- matrix(0, nrow = length(TRBV), ncol = length(TRBJ))
  rownames(temp) <- TRBV
  colnames(temp) <- TRBJ
  for (i in 1:nrow(CD40_tdLN_1_meta)) {
    V <- CD40_tdLN_1_meta$TRBV[i]
    V <- unlist(strsplit(V, "+", fixed = T))
    J <- CD40_tdLN_1_meta$TRBJ[i]
    J <- unlist(strsplit(J, "+", fixed = T))
    for (v in V) {
      for (j in J) {
        temp[v, j] <- temp[v, j] + 1
      }
    }
  }
  library(pheatmap)
  pdf("./绘图/105.CD40_tdLN_V_J_heatmap.pdf", height = 8, width = 6)
  pheatmap(temp, cluster_rows = F, cluster_cols = F,
           border_color = NA, angle_col = 90,
           cellwidth = 16, cellheight = 16,
           color = colorRampPalette(c("white", "#6699cc", "#cc3333"))(100))
  dev.off()
}
{
  Vector_tdLN_1_meta$TRBV <- Vector_tdLN_1_meta$v_gene
  Vector_tdLN_1_meta$TRBV <- unlist(lapply(Vector_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  Vector_tdLN_1_meta$TRBJ <- Vector_tdLN_1_meta$j_gene
  Vector_tdLN_1_meta$TRBJ <- unlist(lapply(Vector_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  TRBV <- unlist(lapply(Vector_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBV <- sort(unique(TRBV))
  TRBJ <- unlist(lapply(Vector_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBJ <- sort(unique(TRBJ))
  temp <- matrix(0, nrow = length(TRBV), ncol = length(TRBJ))
  rownames(temp) <- TRBV
  colnames(temp) <- TRBJ
  for (i in 1:nrow(Vector_tdLN_1_meta)) {
    V <- Vector_tdLN_1_meta$TRBV[i]
    V <- unlist(strsplit(V, "+", fixed = T))
    J <- Vector_tdLN_1_meta$TRBJ[i]
    J <- unlist(strsplit(J, "+", fixed = T))
    for (v in V) {
      for (j in J) {
        temp[v, j] <- temp[v, j] + 1
      }
    }
  }
  library(pheatmap)
  pdf("./绘图/105.Vector_tdLN_V_J_heatmap.pdf", height = 8, width = 6)
  pheatmap(temp, cluster_rows = F, cluster_cols = F,
           border_color = NA, angle_col = 90,
           cellwidth = 16, cellheight = 16,
           color = colorRampPalette(c("white", "#6699cc", "#cc3333"))(100))
  dev.off()
}
{
  FX_tdLN_1_meta$TRBV <- FX_tdLN_1_meta$v_gene
  FX_tdLN_1_meta$TRBV <- unlist(lapply(FX_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  FX_tdLN_1_meta$TRBJ <- FX_tdLN_1_meta$j_gene
  FX_tdLN_1_meta$TRBJ <- unlist(lapply(FX_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  TRBV <- unlist(lapply(FX_tdLN_1_meta$TRBV, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBV <- sort(unique(TRBV))
  TRBJ <- unlist(lapply(FX_tdLN_1_meta$TRBJ, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBJ <- sort(unique(TRBJ))
  temp <- matrix(0, nrow = length(TRBV), ncol = length(TRBJ))
  rownames(temp) <- TRBV
  colnames(temp) <- TRBJ
  for (i in 1:nrow(FX_tdLN_1_meta)) {
    V <- FX_tdLN_1_meta$TRBV[i]
    V <- unlist(strsplit(V, "+", fixed = T))
    J <- FX_tdLN_1_meta$TRBJ[i]
    J <- unlist(strsplit(J, "+", fixed = T))
    for (v in V) {
      for (j in J) {
        temp[v, j] <- temp[v, j] + 1
      }
    }
  }
  library(pheatmap)
  pdf("./绘图/105.FX_tdLN_V_J_heatmap.pdf", height = 8, width = 6)
  pheatmap(temp, cluster_rows = F, cluster_cols = F,
           border_color = NA, angle_col = 90,
           cellwidth = 16, cellheight = 16,
           color = colorRampPalette(c("white", "#6699cc", "#cc3333"))(100))
  dev.off()
}
{
  XCL1_tdLN_2_meta$TRBV <- XCL1_tdLN_2_meta$v_gene
  XCL1_tdLN_2_meta$TRBV <- unlist(lapply(XCL1_tdLN_2_meta$TRBV, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  XCL1_tdLN_2_meta$TRBJ <- XCL1_tdLN_2_meta$j_gene
  XCL1_tdLN_2_meta$TRBJ <- unlist(lapply(XCL1_tdLN_2_meta$TRBJ, function(x){
    unlist(strsplit(x, ";", fixed = T))[1]
  }))
  TRBV <- unlist(lapply(XCL1_tdLN_2_meta$TRBV, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBV <- sort(unique(TRBV))
  TRBJ <- unlist(lapply(XCL1_tdLN_2_meta$TRBJ, function(x){
    unlist(strsplit(x, "+", fixed = T))
  }))
  TRBJ <- sort(unique(TRBJ))
  temp <- matrix(0, nrow = length(TRBV), ncol = length(TRBJ))
  rownames(temp) <- TRBV
  colnames(temp) <- TRBJ
  for (i in 1:nrow(XCL1_tdLN_2_meta)) {
    V <- XCL1_tdLN_2_meta$TRBV[i]
    V <- unlist(strsplit(V, "+", fixed = T))
    J <- XCL1_tdLN_2_meta$TRBJ[i]
    J <- unlist(strsplit(J, "+", fixed = T))
    for (v in V) {
      for (j in J) {
        temp[v, j] <- temp[v, j] + 1
      }
    }
  }
  library(pheatmap)
  pdf("./绘图/105.XCL1_tdLN_V_J_heatmap.pdf", height = 8, width = 6)
  pheatmap(temp, cluster_rows = F, cluster_cols = F,
           border_color = NA, angle_col = 90,
           cellwidth = 16, cellheight = 16,
           color = colorRampPalette(c("white", "#6699cc", "#cc3333"))(100))
  dev.off()
}


