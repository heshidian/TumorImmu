setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/5.RNA+TCR-BCR")
library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)

CD8 <- readRDS("./注释/CD8+/CD8_T_annotated.rds")
DimPlot(CD8)
CD8$Celltype_2 <- "CD8_T"
unique(CD8$Sample)

CD4 <- readRDS("./注释/CD4+/CD4_T_annotated.rds")
DimPlot(CD4)
CD4$Celltype_2 <- "CD4_T"
unique(CD4$Sample)

DC <- readRDS("./注释/DC/DC_annotated.rds")
DimPlot(DC)
unique(DC$SubCelltype)
DC <- subset(DC, subset = SubCelltype != "Unidentified")
DC$Celltype_2 <- "DC"
unique(DC$Sample)

CD8_CD4_DC <- merge(x = CD8, y = list(CD4, DC))
table(CD8_CD4_DC$Celltype_2, CD8_CD4_DC$Sample)

CD8_CD4_DC <- NormalizeData(CD8_CD4_DC)
CD8_CD4_DC_Sample <- SplitObject(CD8_CD4_DC, split.by = "Sample")
CD8_CD4_DC_Sample <- lapply(CD8_CD4_DC_Sample,
                            function(x){NormalizeData(x)})
samples <- names(CD8_CD4_DC_Sample)

human_mouse <- read.csv("./cellphoneDB/ensembl_human_mouse_symbol.txt",
                        sep = "\t", header = T)
human_mouse <- human_mouse[!duplicated(human_mouse$Gene.name),]
human_mouse <- human_mouse[human_mouse$Gene.name != "",]
rownames(human_mouse) <- human_mouse$Gene.name
human_mouse <- human_mouse[-which(human_mouse$Human.gene.name == ""),]

for (i in c("Fx-T_2", "XCL1-T_2", "Vector-T_2")) {
  temp <- CD8_CD4_DC_Sample[[i]]
  temp_data <- temp@assays$RNA@data
  temp_data <- temp_data[rownames(temp_data) %in% human_mouse$Gene.name,]
  temp_data <- temp_data[rowSums(temp_data) > 0,]
  rownames(temp_data) <- human_mouse[rownames(temp_data), "Human.gene.name"]
  # Save normalised counts - NOT scaled!
  if (!dir.exists(paste0("./cellphoneDB/", i, "/Exp"))) {
    dir.create(paste0("./cellphoneDB/", i, "/Exp"))
  }
  writeMM(temp_data, file = paste0("./cellphoneDB/", i, "/Exp/matrix.mtx"))
  # save gene and cell names
  write(x = rownames(temp_data), file = paste0("./cellphoneDB/", i, "/Exp/features.tsv"))
  write(x = colnames(temp_data), file = paste0("./cellphoneDB/", i, "/Exp/barcodes.tsv"))
  meta <- data.frame(Cell = colnames(temp), cell_type = temp$Celltype_2)
  write.csv(meta, paste0("./cellphoneDB/", i, "/meta.csv"),
            row.names = F)
}

direction <- c("DC|CD8_T", "DC|CD4_T", "CD8_T|DC", "CD4_T|DC")
cellphoneDB_result <- c()
for (i in c("Fx-T_2", "XCL1-T_2", "Vector-T_2")) {
  files <- list.files(paste0("./cellphoneDB/", i),
                      full.names = TRUE)
  p_value <- files[grep(pattern = "statistical_analysis_pvalues", x = files)]
  mean <- files[grep(pattern = "statistical_analysis_means", x = files)]
  p_value <- read.csv(p_value, sep = "\t", check.names = F)
  mean <- read.csv(mean, sep = "\t", check.names = F)
  rownames(p_value) <- p_value$interacting_pair
  rownames(mean) <- mean$interacting_pair
  for (j in direction) {
    temp_p <- p_value[,c(colnames(p_value)[1:13], j)]
    temp_mean <- mean[,c(colnames(mean)[1:13], j)]
    temp_mean <- temp_mean[rownames(temp_p),]
    temp_result <- data.frame(p_value[,1:13],
                              mean = temp_mean[,14],
                              p = temp_p[,14],
                              interaction = gsub(pattern = "|", replacement = "->", x = j, fixed = T),
                              sample = i)
    cellphoneDB_result <- as.data.frame(rbind(cellphoneDB_result,
                                              temp_result))
  }
}
# cellphoneDB_result <- cellphoneDB_result[cellphoneDB_result$p < 0.05,]
for (i in 1:nrow(cellphoneDB_result)) {
  if (cellphoneDB_result$gene_a[i] == "") {
    if (cellphoneDB_result$gene_b[i] == "") {
      ligand <- unlist(strsplit(cellphoneDB_result$partner_a[i],
                                "complex:", fixed = T))[2]
      receptor <- unlist(strsplit(cellphoneDB_result$partner_b[i],
                                  "complex:", fixed = T))[2]
    } else {
      receptor <- cellphoneDB_result$gene_b[i]
      ligand <- unlist(strsplit(cellphoneDB_result$interacting_pair[i],
                                paste0("_", receptor), fixed = TRUE))[1]
    }
  } else {
    ligand <- cellphoneDB_result$gene_a[i]
    index <- nchar(ligand)
    receptor <- substr(cellphoneDB_result$interacting_pair[i],
                       index+2, nchar(cellphoneDB_result$interacting_pair[i]))
  }
  cellphoneDB_result$interacting_pair[i] <- paste0(ligand, "-", receptor)
}

sample_rename <- data.frame(original = c("Vector-T_2", "XCL1-T_2", "Fx-T_2"),
                            rename = c("OT-I", "OT-I/XCL1", "OT-I/FX"),
                            row.names = c("Vector-T_2", "XCL1-T_2", "Fx-T_2"))
cellphoneDB_result$sample <- sample_rename[cellphoneDB_result$sample,
                                           "rename"]
cellphoneDB_result$sample <- factor(cellphoneDB_result$sample,
                                    levels = c("OT-I", "OT-I/XCL1", "OT-I/FX"))
# cellphoneDB_result$mean[cellphoneDB_result$mean > 1.5] <- 1.5
cellphoneDB_result$sig <- ifelse(cellphoneDB_result$p < 0.05,
                                 "P < 0.05", "P >= 0.05")


cellphoneDB_result_1 <- cellphoneDB_result[cellphoneDB_result$interaction %in% c("DC->CD8_T", "DC->CD4_T"),]
cellphoneDB_result_1$interacting_pair[grep(pattern = "ENTP", x = cellphoneDB_result_1$interacting_pair)]
cellphoneDB_result_1 <- cellphoneDB_result_1[cellphoneDB_result_1$interacting_pair %in% c("IL1B-IL1_receptor_inhibitor",
                                                                                          "IL18-IL18_receptor",
                                                                                          "IL15-IL15_receptor",
                                                                                          "IL12-IL12_receptor",
                                                                                          "IL27-IL27_receptor",
                                                                                          "CXCL9-CXCR3", "CXCL16-CXCR6",
                                                                                          "CXCL10-DPP4", "CXCL10-CXCR3",
                                                                                          "CCL22-DPP4", "CD86-CD28",
                                                                                          "CD80-CD28", "CD274-PDCD1",
                                                                                          "CD80-CTLA4", "TNFSF9-TNFRSF9",
                                                                                          "TNF-TNFRSF1B", "TNF-TNFRSF1A",
                                                                                          "PDCD1LG2-PDCD1", "ENTPD1-ADORA2A"),]
quantile(cellphoneDB_result_1$mean)
cellphoneDB_result_1$mean[cellphoneDB_result_1$mean > 1.5] <- 1.5
cellphoneDB_result_1$interacting_pair <- factor(cellphoneDB_result_1$interacting_pair,
                                                levels = rev(c("IL1B-IL1_receptor_inhibitor",
                                                               "IL18-IL18_receptor",
                                                               "IL15-IL15_receptor",
                                                               "IL12-IL12_receptor",
                                                               "IL27-IL27_receptor",
                                                               "CXCL9-CXCR3", "CXCL16-CXCR6",
                                                               "CXCL10-DPP4", "CXCL10-CXCR3",
                                                               "CCL22-DPP4", "CD86-CD28",
                                                               "CD80-CD28", "CD274-PDCD1",
                                                               "CD80-CTLA4", "TNFSF9-TNFRSF9",
                                                               "TNF-TNFRSF1B", "TNF-TNFRSF1A",
                                                               "PDCD1LG2-PDCD1", "ENTPD1-ADORA2A")))
pdf("./cellphoneDB/1.DC_T.pdf", height = 5, width = 6)
ggplot(data = cellphoneDB_result_1,
       aes(x = sample, y = interacting_pair)) +
  geom_point(data = cellphoneDB_result_1,
             aes(x = sample, y = interacting_pair,
                 size = mean, fill = mean),
             shape = 21, colour = "black") +
  geom_point(data = cellphoneDB_result_1,
             aes(x = sample, y = interacting_pair,
                 size = mean, color = sig), fill = NA,
             shape = 21) +
  facet_wrap(.~interaction) +
  scale_size(breaks = c(0.01, 0.375, 0.75, 1.125, 1.5),  # 指定值
             labels = c(0.01, 0.375, 0.75, 1.125, 1.5),
             range = c(0, 8)) +
  scale_fill_gradientn(colours = c("blue", "white", "red", "red"),
                       values = c(0, 0.375, 0.75, 1.125, 1.5)) +
  scale_color_manual(values = c("black", "gray")) +
  labs(size = "cellphoneDB\nmean value",
       fill = "cellphoneDB\nmean value",
       color = "P-value", x = "", y = "Ligand-Receptor") +
  # theme_bw() +
  theme(# panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        strip.text.x.top = element_text(size = 12, colour = "black"))
dev.off()


cellphoneDB_result_2 <- cellphoneDB_result[cellphoneDB_result$interaction %in% c("CD8_T->DC", "CD4_T->DC"),]
cellphoneDB_result$interacting_pair[grep(pattern = "CD47", x = cellphoneDB_result$interacting_pair)]
cellphoneDB_result_2 <- cellphoneDB_result_2[cellphoneDB_result_2$interacting_pair %in% c("XCL1-XCR1",
                                                                                          "CCL5-CCR5",
                                                                                          "TNFSF13B-TNFRSF13B",
                                                                                          "TNFSF11-TNFRSF11A",
                                                                                          "TNF-TNFRSF1A",
                                                                                          "SEMA4D-CD72", "SELPLG-SELL",
                                                                                          "LTB-LTBR", "IFNG-Type_II_IFNR",
                                                                                          "ICAM2-CD209", "FLT3LG-FLT3",
                                                                                          "CD40LG-CD40", "TGFB1-TGFbeta_receptor1",
                                                                                          "PPIA-BSG", "CD48-CD244", "CD47-SIRB1_complex"),]

quantile(cellphoneDB_result_2$mean)
cellphoneDB_result_2$mean[cellphoneDB_result_2$mean > 1.5] <- 1.5
cellphoneDB_result_2$interacting_pair <- factor(cellphoneDB_result_2$interacting_pair,
                                                levels = rev(c("XCL1-XCR1",
                                                               "CCL5-CCR5",
                                                               "TNFSF13B-TNFRSF13B",
                                                               "TNFSF11-TNFRSF11A",
                                                               "TNF-TNFRSF1A",
                                                               "SEMA4D-CD72", "SELPLG-SELL",
                                                               "LTB-LTBR", "IFNG-Type_II_IFNR",
                                                               "ICAM2-CD209", "FLT3LG-FLT3",
                                                               "CD40LG-CD40", "TGFB1-TGFbeta_receptor1",
                                                               "PPIA-BSG", "CD48-CD244", "CD47-SIRB1_complex")))
pdf("./cellphoneDB/2.T_DC.pdf", height = 4.5, width = 6)
ggplot(data = cellphoneDB_result_2,
       aes(x = sample, y = interacting_pair)) +
  geom_point(data = cellphoneDB_result_2,
             aes(x = sample, y = interacting_pair,
                 size = mean, fill = mean),
             shape = 21, colour = "black") +
  geom_point(data = cellphoneDB_result_2,
             aes(x = sample, y = interacting_pair,
                 size = mean, color = sig), fill = NA,
             shape = 21) +
  facet_wrap(.~interaction) +
  scale_size(breaks = c(0.01, 0.375, 0.75, 1.125, 1.5),  # 指定值
             labels = c(0.01, 0.375, 0.75, 1.125, 1.5),
             range = c(0, 8)) +
  scale_fill_gradientn(colours = c("blue", "white", "red", "red"),
                       values = c(0, 0.375, 0.75, 1.125, 1.5)) +
  scale_color_manual(values = c("black", "gray")) +
  labs(size = "cellphoneDB\nmean value",
       fill = "cellphoneDB\nmean value",
       color = "P-value", x = "", y = "Ligand-Receptor") +
  # theme_bw() +
  theme(# panel.grid = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 10, colour = "black"),
    legend.title = element_text(size = 12, colour = "black"),
    strip.text.x.top = element_text(size = 12, colour = "black"))
dev.off()


