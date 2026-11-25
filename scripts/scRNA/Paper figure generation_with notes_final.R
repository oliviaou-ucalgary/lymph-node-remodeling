########### Basic processing #########
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggsci)
library(ggplot2)
library(MAST)

#Load dataset and create Seurat objects
CTRL.data <- Read10X_h5(filename = "./naive/filtered_feature_bc_matrix.h5")
CTRL <- CreateSeuratObject(counts = CTRL.data, project = "CTRL FRC", min.cells = 3, min.features = 200)
CTRL_FRC <- subset(CTRL, subset = Ptprc ==0 & Cdh5 ==0 & Prox1 ==0 & Lyve1 ==0 & Pdgfra>1)
CTRL_PvC <- subset(CTRL, subset = Prox1 ==0 & Acta2>30 & Itga7>-1)
name_CTRL_FRC <- colnames(CTRL_FRC)
name_CTRL_PvC <- colnames(CTRL_PvC)
shared_CTRL <- intersect(x = name_CTRL_FRC, y = name_CTRL_PvC)
CTRL_PvC <- subset (CTRL_PvC, cells = shared_CTRL, invert = T)
CTRL_FRCs <- merge(CTRL_FRC,CTRL_PvC)

OXd4.data <- Read10X_h5(filename = "./OXd4/filtered_feature_bc_matrix.h5")
OXd4 <- CreateSeuratObject(counts = OXd4.data, project = "OXd4 FRC", min.cells = 3, min.features = 200)
OXd4_FRC <- subset(OXd4, subset = Ptprc ==0 & Cdh5 ==0 & Prox1 ==0 & Lyve1 ==0 & Pdgfra>1)
OXd4_PvC <- subset(OXd4, subset = Prox1 ==0 & Acta2>30 & Itga7>-1)
name_OXd4_FRC <- colnames(OXd4_FRC)
name_OXd4_PvC <- colnames(OXd4_PvC)
shared_OXd4 <- intersect(x = name_OXd4_FRC, y = name_OXd4_PvC)
OXd4_PvC <- subset (OXd4_PvC, cells = shared_OXd4, invert = T)
OXd4_FRCs <- merge(OXd4_FRC,OXd4_PvC)


#Look for matrix for quality control
all_data_for_filt <- merge(CTRL_FRCs, OXd4_FRCs, add.cell.ids = c("CTRL", "OXd4"))
all_data_for_filt <- PercentageFeatureSet(all_data_for_filt, "^mt-", col.name = "percent_mito")
all_data_for_filt <- PercentageFeatureSet(all_data_for_filt, "^Rp[sl]", col.name = "percent_ribo")
all_data_for_filt <- PercentageFeatureSet(all_data_for_filt, "^Hb[^(p)]", col.name = "percent_hb")
all_data_for_filt <- PercentageFeatureSet(all_data_for_filt, "Pecam1|Pf4", col.name = "percent_plat")
feats <- c("percent_mito", "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(all_data_for_filt, group.by = "orig.ident", features = feats, pt.size = 0.1)
par(mar = c(4, 8, 2, 1))
C <- all_data_for_filt@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)


#Data integration#
all_sample.list <- c(CTRL_FRCs, OXd4_FRCs)

all_sample.list <- lapply(X = all_sample.list, FUN = function (x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x[["percent.ribo"]] <- PercentageFeatureSet(x, pattern = "^Rp[sl]")
  x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 & percent.ribo > 7)
  x <- x[!grepl("Mala1", rownames(x)), ]
  #x <- x[!grepl('^RP[SL]', rownames(x)), ]
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = all_sample.list)
FRCs_anchors <- FindIntegrationAnchors(object.list = all_sample.list, anchor.features = features)
FRCs_combined <- IntegrateData(anchorset = FRCs_anchors)

DefaultAssay(FRCs_combined) <- "RNA"

Cell_cycle_Regev_BLOOD <- read.delim("./cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
s.genes <- Cell_cycle_Regev_BLOOD[1:43,1]
g2m.genes <- Cell_cycle_Regev_BLOOD[44:97,1]
FRCs_combined <- CellCycleScoring(object = FRCs_combined, s.features = s.genes, 
                                                g2m.features = g2m.genes, set.ident = TRUE)


DefaultAssay(FRCs_combined) <- "integrated"
FRCs_combined_regressed <- ScaleData(FRCs_combined, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "percent.ribo"))
FRCs_combined_regressed <- RunPCA(FRCs_combined_regressed)
FRCs_combined_regressed_30 <- FindNeighbors(FRCs_combined_regressed, reduction = "pca", dims = 1:30)
FRCs_combined_regressed_30 <- FindClusters(FRCs_combined_regressed_30, resolution = 0.8)
FRCs_combined_regressed_30 <- RunUMAP(FRCs_combined_regressed_30, reduction = "pca", dims = 1:30)

# Assign a new name and duplicate the final object
FRC_paper <- FRCs_combined_regressed_30
DefaultAssay(FRC_paper) <- "RNA"

#Find out the signatures between all clusters 
FRC_diff.markers <- FindAllMarkers(FRC_paper, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
FRC_diff.markers_top_8 <- FRC_diff.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
DotPlot(FRC_paper, features = unique (FRC_diff.markers_top_8$gene), col.min = 0, col.max = 5, dot.scale = 4, 
        split.by = "orig.ident", cols = c("lightblue2", "red3")) + RotatedAxis() + FontSize(8)

Idents(FRC_paper) <- "seurat_clusters"
FRC_paper <- subset(FRC_paper, subset = seurat_clusters == "11", invert = T) #Remove the excessive and tiny proliferating population

#Rename Idents based on top 8 signature markers
FRC_paper <- RenameIdents(FRC_paper, `0` = "Ccl19lo TRC", `1` = "Ccl19lo TRC", `2` = "MRC", `3` = "Ccl19hi TRC", 
                          `4` = "MedRC", `5` = "Inmt+ SC", `6` = "Cxcl13inter RC", `7` = "dividing_cells", `8` = "FDC", 
                          `9` = "Cd34+ SC", `10` = "PvC", `12` = "Cxcl9+ TRC", `13` = "Cd34+ SC", `14` = "Cd34+ SC")
FRC_palette <- c("#71D0F5DD", "#46732EDD", "#075149AA", "#1A9993", "#D2AF81CC", "#FED439CC", "#8A9197CC", 
                  "#FD7446AA", "#197EC0DD", "#F05C3BAA", "#C80813CC")
DimPlot(FRC_paper, reduction = "umap", split.by = "orig.ident", pt.size = 1,
        cols = FRC_palette, label = F)
ggsave("umap_FRC_paper_no_label.png", width = 15, height = 8, dpi = 600)
ggsave("umap_FRC_paper_no_label.pdf", width = 15, height = 8, dpi = 600)

featureplot_list <- c("Ccl21a", "Ccl19", "Grem1", "Cxcl9", "Inmt", "Ch25h", "Rarres1",
                      "Cr2", "Cxcl13", "Tnfsf11", "Cxcl12", "Nr4a1", "Fbln1",
                      "Cd34", "Ccl2", "Ccl7", "Acta2", "Madcam1")
for (feature_gene_i in featureplot_list) {
    filename_i <- paste("FeaturePlot_", feature_gene_i, ".png", sep="")
    filename_ii <- paste("FeaturePlot_", feature_gene_i, ".pdf", sep="")
    FeaturePlot(FRC_paper, features = feature_gene_i, min.cutoff = "q10", ncol = 2, split.by = "orig.ident", 
                pt.size = 0.05, cols = c("lightgrey","darkgreen")) & theme(title = element_text(size = 5), 
                                                           axis.title = element_text(size = 7))
    ggsave(filename_i, width = 6, height = 3, dpi = 450)
    ggsave(filename_ii, width = 6, height = 3, dpi = 450)
}

FeaturePlot(FRC_paper, features = "Lyve1", min.cutoff = "q10", ncol = 2, split.by = "orig.ident", 
            pt.size = 0.05, cols = c("lightgrey","darkgreen")) & theme(title = element_text(size = 5), 
                                                                       axis.title = element_text(size = 7))
#Cell Cycle
FRC_paper_Cell_Cycle <- CellCycleScoring(object = FRC_paper, s.features = s.genes, g2m.features = g2m.genes, 
                                                set.ident = TRUE)
DimPlot(FRC_paper_Cell_Cycle, reduction = "umap", split.by = "orig.ident", cols = c("purple", "grey", "orange2"),
        pt.size = 0.5)
ggsave("umap_FRC_cell_cycle.png", width = 7.5, height = 4, dpi = 450)
ggsave("umap_FRC_cell_cycle.pdf", width = 7.5, height = 4, dpi = 450)

#DGE_all_RNA
FRC_diff.markers <- FindAllMarkers(FRC_paper, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
FRC_diff.markers_top_8 <- FRC_diff.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
DotPlot(FRC_paper, features = unique (FRC_diff.markers_top_8$gene), col.min = 0, col.max = 5, dot.scale = 4, 
        split.by = "orig.ident", cols = c("lightblue2", "red3")) + RotatedAxis() + FontSize(8)
ggsave("top8featurePlot_split.png", width = 20, height = 10, dpi = 600)
ggsave("top8featurePlot_split.pdf", width = 20, height = 10, dpi = 600)

DotPlot(FRC_paper, features = unique (FRC_diff.markers_top_8$gene), col.min = 0, col.max = 5, dot.scale = 4, 
        cols = c("lightgrey","blue1")) + RotatedAxis() + FontSize(8)
ggsave("top8featurePlot.png", width = 20, height = 5, dpi = 600)
ggsave("top8featurePlot.pdf", width = 20, height = 5, dpi = 600)

#DGE_all_integrated
DefaultAssay(FRC_paper) <- "integrated"
FRC_DGE_integrated <- FindAllMarkers(FRC_paper, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
FRC_DGE_integrated_top_8 <- FRC_DGE_integrated %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
DotPlot(FRC_paper, features = unique (FRC_DGE_integrated_top_8$gene), col.min = 0, col.max = 5, dot.scale = 4, 
        split.by = "orig.ident", cols = c("lightblue2", "red3")) + RotatedAxis() + FontSize(8)
ggsave("top8featurePlot_integrated.png", width = 20, height = 10, dpi = 600)
ggsave("top8featurePlot_integrated.pdf", width = 20, height = 10, dpi = 600)

DotPlot(FRC_paper, features = unique(FRC_DGE_integrated_top_8$gene), col.min = 0, col.max = 5, dot.scale = 4, 
        cols = c("lightgrey","blue1")) + RotatedAxis() + FontSize(8)
ggsave("top8featurePlot_integrated.png", width = 20, height = 5, dpi = 600)
ggsave("top8featurePlot_integrated.pdf", width = 20, height = 5, dpi = 600)
DefaultAssay(FRC_paper) <- "RNA"

#VlnPlots

TRC_list <- c("Pdpn", "Ccl21a", "Ccl19")
Cxcl13_list <- c("Pdpn", "Cxcl13")
for (feature_gene_i in TRC_list) {
  filename_i <- paste("VlnPlot_", feature_gene_i, ".png", sep="")
  filename_ii <- paste("VlnPlot_", feature_gene_i, ".pdf", sep="")
  VlnPlot(FRC_paper,idents = c("Ccl19hi TRC", "Ccl19lo TRC"), features = feature_gene_i, split.by = "orig.ident", 
              pt.size = 0.1, cols = c("lightblue2", "red3")) + theme(legend.position = 'right', 
                                                                     axis.text.x = element_text(angle = 0, hjust=0.5))
  ggsave(filename_i, width = 5, height = 4, dpi = 450)
  ggsave(filename_ii, width = 5, height = 4, dpi = 450)
}

for (feature_gene_i in Cxcl13_list) {
  filename_i <- paste("VlnPlot_Cxcl13_RC_", feature_gene_i, ".png", sep="")
  filename_ii <- paste("VlnPlot_Cxcl13_RC_", feature_gene_i, ".pdf", sep="")
  VlnPlot(FRC_paper,idents = c("Cxcl13inter RC"), features = feature_gene_i, split.by = "orig.ident", 
          pt.size = 0.1, cols = c("lightblue2", "red3")) + theme(legend.position = 'right', 
                                                                 axis.text.x = element_text(angle = 0, hjust=0.5))
  ggsave(filename_i, width = 4, height = 4, dpi = 450)
  ggsave(filename_ii, width = 4, height = 4, dpi = 450)
}

#Calculate % per cluster across age
Treatment <- split(FRC_paper@meta.data, FRC_paper@meta.data$orig.ident)
cell_per_treatment <- lapply(Treatment, function(x) {
  #x <- l_condition[[1]]
  x_cells <- unlist(lapply(split(x, x$seurat_clusters), nrow))
  x_cells
})



###### Convert for scVelo ######
CTRL_FRCs <- subset(FRC_paper, subset = orig.ident == "CTRL FRC")
SaveH5Seurat(CTRL_FRCs, filename = "CTRL_FRC_regressed.h5Seurat")
Convert("CTRL_FRC_regressed.h5Seurat", dest = "h5ad")
FRC_palette <- c("#71D0F5DD", "#46732EDD", "#075149AA", "#1A9993", "#D2AF81CC", "#FED439CC", "#8A9197CC", 
                 "#FD7446AA", "#197EC0DD", "#F05C3BAA", "#C80813CC")
ctrl_umap <- DimPlot(CTRL_FRCs, reduction = "umap", split.by = "orig.ident", pt.size = 1,
                     cols = FRC_palette, label = F)
ctrl_build <- ggplot_build(ctrl_umap)
ctrl_data <- ctrl_build$data[[1]]
write.csv(ctrl_data, file = "ctrl_cluster_color.csv")
write.csv(Embeddings(CTRL_FRCs, reduction = "umap"), file = "ctrl_cell_embeddings.csv")

OXd4_FRCs <- subset(FRC_paper, subset = orig.ident == "OXd4 FRC")
SaveH5Seurat(OXd4_FRCs, filename = "OXd4_FRC_regressed.h5Seurat")
Convert("OXd4_FRC_regressed.h5Seurat", dest = "h5ad")
OXd4_umap <- DimPlot(OXd4_FRCs, reduction = "umap", cols = pal_simpsons("springfield", alpha = 0.8)(15))
OXd4_build <- ggplot_build(OXd4_umap)
OXd4_data <- OXd4_build$data[[1]]
write.csv(OXd4_data, file = "OXd4_cluster_color.csv")
write.csv(Embeddings(OXd4_FRCs, reduction = "umap"), file = "OXd4_cell_embeddings.csv")


########## Enrichment analysis ###########
library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationHub")
library("ReactomePA")
library("enrichplot")


pval_FRC_diff.markers <- subset(FRC_diff.markers, rowSums(FRC_diff.markers[5] < 0.05) > 0)
df <- pval_FRC_diff.markers[,7:6]
dfsample <- split(df$gene, df$cluster)
dfsample$`Ccl19hi TRC` = bitr(dfsample$`Ccl19hi TRC`, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = "org.Mm.eg.db")
dfsample$`Ccl19lo TRC` = bitr(dfsample$`Ccl19lo TRC`, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = "org.Mm.eg.db")
dfsample$`Cxcl9+ TRC` = bitr(dfsample$`Cxcl9+ TRC`, fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = "org.Mm.eg.db")
dfsample$`Cxcl13inter RC` = bitr(dfsample$`Cxcl13inter RC`, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = "org.Mm.eg.db")
dfsample$`dividing_cells` = bitr(dfsample$`dividing_cells`, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = "org.Mm.eg.db")
dfsample$`MRC` = bitr(dfsample$`MRC`, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Mm.eg.db")
dfsample$`FDC` = bitr(dfsample$`FDC`, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Mm.eg.db")
dfsample$`MedRC` = bitr(dfsample$`MedRC`, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = "org.Mm.eg.db")
dfsample$`Inmt+ SC` = bitr(dfsample$`Inmt+ SC`, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db")
dfsample$`Cd34+ SC` = bitr(dfsample$`Cd34+ SC`, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db")
dfsample$`PvC` = bitr(dfsample$`PvC`, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Mm.eg.db")
genelist <- list("Ccl19hi TRC" = dfsample$`Ccl19hi TRC`$ENTREZID,
                 "Ccl19lo TRC" = dfsample$`Ccl19lo TRC`$ENTREZID,
                 "Cxcl9+ TRC" = dfsample$`Cxcl9+ TRC`$ENTREZID,
                 "Cxcl13inter RC" = dfsample$`Cxcl13inter RC`$ENTREZID,
                 "dividing_cells" = dfsample$`dividing_cells`$ENTREZID,
                 "MRC" = dfsample$`MRC`$ENTREZID,
                 "FDC" = dfsample$`FDC`$ENTREZID,
                 "MedRC" = dfsample$`MedRC`$ENTREZID,
                 "Inmt+ SC" = dfsample$`Inmt+ SC`$ENTREZID,
                 "Cd34+ SC" = dfsample$`Cd34+ SC`$ENTREZID,
                 "PvC" = dfsample$`PvC`$ENTREZID)

pval_FRC_DGE_integrated <- subset(FRC_DGE_integrated, rowSums(FRC_DGE_integrated[5] < 0.05) > 0)
df_int <- pval_FRC_DGE_integrated[,7:6]
dfsample_int <- split(df_int$gene, df_int$cluster)
dfsample_int$`Ccl19hi TRC` = bitr(dfsample_int$`Ccl19hi TRC`, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = "org.Mm.eg.db")
dfsample_int$`Ccl19lo TRC` = bitr(dfsample_int$`Ccl19lo TRC`, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = "org.Mm.eg.db")
dfsample_int$`Cxcl13inter RC` = bitr(dfsample_int$`Cxcl13inter RC`, fromType = "SYMBOL", toType = "ENTREZID",
                                     OrgDb = "org.Mm.eg.db")
int_genelist <- list("Ccl19hi TRC" = dfsample_int$`Ccl19hi TRC`$SYMBOL,
                     "Ccl19lo TRC" = dfsample_int$`Ccl19lo TRC`$SYMBOL,
                     "Cxcl13inter RC" = dfsample_int$`Cxcl13inter RC`$SYMBOL)
GO_BP_3clusters_int <- compareCluster(geneClusters = int_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                                      ont = "BP", pAdjustMethod="BH", keyType = "SYMBOL")
GO_CC_3clusters_int <- compareCluster(geneClusters = int_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                                      ont = "CC", pAdjustMethod="BH", keyType = "SYMBOL")
GO_MF_3clusters_int <- compareCluster(geneClusters = int_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                                      ont = "MF", pAdjustMethod="BH", keyType = "SYMBOL")

cnet_plot_selected <- c("extracellular matrix organization", "collagen fibril organization", "cell chemotaxis",
                        "B cell activation", "cellular response to fibroblast growth factor stimulus")

cnetplot(GO_BP_3clusters_int, 
         showCategory = cnet_plot_selected, 
         includeAll = T, 
         colorEdge= FALSE,
         cex.params = list(gene_label = 1.5, 
                           category_label = 3)) + scale_fill_manual(values = c("#71D0F5DD", "#075149AA", "#FED439CC"))
ggsave("BP_Go_cnetplot.png", width = 12, height = 10, dpi = 450)


GOCluster_MF <- compareCluster(geneClusters = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                               ont = "MF", pAdjustMethod="BH")
selected_pathway <- GOCluster_MF@compareClusterResult[grep('(collagen|extracellular matrix)', 
                                                           GOCluster_MF@compareClusterResult$Description), ]
selected_pathway <- sample(selected_pathway$Description)
dotplot(GOCluster_MF, showCategory = selected_pathway, includeAll = T, 
        title = "GO_Molecular Functions") + scale_size_area(max_size = 12)

GOCluster_CC <- compareCluster(geneClusters = cnet_plot_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                               ont = "CC", pAdjustMethod="BH")
selected_pathway <- GOCluster_CC@compareClusterResult[grep('(collagen|extracellular matrix)', 
                                                           GOCluster_CC@compareClusterResult$Description), ]
selected_pathway <- sample(selected_pathway$Description)
dotplot(GOCluster_CC, showCategory = selected_pathway, includeAll = T, 
        title = "GO_Cellular Components") + scale_size_area(max_size = 12)

GOCluster_BP <- compareCluster(geneClusters = cnet_plot_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                               ont = "BP", pAdjustMethod="BH")
selected_pathway <- GOCluster_BP@compareClusterResult[grep('(collagen|extracellular matrix)', 
                                                           GOCluster_BP@compareClusterResult$Description), ]
selected_pathway <- sample(selected_pathway$Description)
dotplot(GOCluster_BP, showCategory = selected_pathway, includeAll = T, 
        title = "GO_Biological Process") + scale_size_area(max_size = 12)



#Cxcl13inter RCs#
df_FRC_rps_removed_DEGs <- subset(pval_FRC_rps_removed_DEGs, rowSums(pval_FRC_rps_removed_DEGs[5] < 0.05) > 0)
df <- df_FRC_rps_removed_DEGs[,7:6]
dfsample <- split(df$gene, df$cluster)
dfsample$`Cxcl13inter RC` = bitr(dfsample$`Cxcl13inter RC`, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = "org.Mm.eg.db")
Cxcl13int_BP <- enrichGO(dfsample$`Cxcl13inter RC`$SYMBOL, OrgDb = "org.Mm.eg.db", ont = "BP", pAdjustMethod="BH",
                         keyType = "SYMBOL")
barplot(Cxcl13int_BP)

Cxcl13int_RCs["gene"] <- rownames(Cxcl13int_RCs)
Cxcl13int_RCs <- subset(Cxcl13int_RCs, rowSums(Cxcl13int_RCs[5] < 0.05 & Cxcl13int_RCs[2] > 0) > 0)
pval_Cxcl13int_RCs <- Cxcl13int_RCs[ ,6]
pval_Cxcl13int_RCs = bitr(pval_Cxcl13int_RCs, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Mm.eg.db")
Cxcl13int_BP_2 <- enrichGO(pval_Cxcl13int_RCs$SYMBOL, OrgDb = "org.Mm.eg.db", ont = "BP", pAdjustMethod="BH",
                         keyType = "SYMBOL")
barplot(Cxcl13int_BP_2)

df_Cxcl13_vs_Ccl19lo_RCs <- split(df$gene, df$cluster)
df_Cxcl13_vs_Ccl19lo_RCs$`Cxcl13inter RC` = bitr(df_Cxcl13_vs_Ccl19lo_RCs$`Cxcl13inter RC`, 
                                                 fromType = "SYMBOL", 
                                                 toType = "ENTREZID",
                                                 OrgDb = "org.Mm.eg.db")
df_Cxcl13_vs_Ccl19lo_RCs$`Ccl19lo TRC` = bitr(df_Cxcl13_vs_Ccl19lo_RCs$`Ccl19lo TRC`, 
                                                 fromType = "SYMBOL", 
                                                 toType = "ENTREZID",
                                                 OrgDb = "org.Mm.eg.db")
Cxcl13_vs_Ccl19lo_RCs_genelist <- list("Cxcl13inter RC" = df_Cxcl13_vs_Ccl19lo_RCs$`Cxcl13inter RC`$SYMBOL,
                                       "Ccl19lo TRC" = df_Cxcl13_vs_Ccl19lo_RCs$`Ccl19lo TRC`$SYMBOL)

BP_Cxcl13_vs_Ccl19lo_RCs <- compareCluster(Cxcl13_vs_Ccl19lo_RCs_genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                                           ont = "BP", pAdjustMethod="BH", keyType = "SYMBOL")
dotplot(BP_Cxcl13_vs_Ccl19lo_RCs, showCategory = 5, includeAll = T, 
        title = "GO_Biological Process") + scale_size_area(max_size = 12)

###### Monocle analysis ######
## Naive Monocle Trajectory Analysis

name_of_clusters <- c(`0` = "Ccl19lo TRC", `1` = "Ccl19lo TRC", `2` = "MRC", `3` = "Ccl19hi TRC", 
                      `4` = "MedRC", `5` = "Inmt+ SC", `6` = "Cxcl13inter RC", `7` = "dividing_cells", `8` = "FDC", 
                      `9` = "Cd34+ SC", `10` = "PvC", `12` = "Cxcl9+ TRC", `13` = "Cd34+ SC", `14` = "Cd34+ SC")
mono_palette <- c("#71D0F5DD","#71D0F5DD", "#46732EDD", "#075149AA", "#1A9993", "#D2AF81CC", "#FED439CC", "#8A9197CC", 
                 "#FD7446AA", "#197EC0DD", "#F05C3BAA", "#C80813CC", "#197EC0DD", "#197EC0DD")


monocle_CTRL <- subset(FRC_paper, subset = orig.ident == "CTRL FRC")
monocle_CTRL_deg <- FindAllMarkers(monocle_CTRL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
monocle_CTRL <- as.CellDataSet(monocle_CTRL)
monocle_CTRL <- estimateSizeFactors(monocle_CTRL)
monocle_CTRL <- estimateDispersions(monocle_CTRL)

## select top 1000 genes as filter
deg_CTRL_1000 <- row.names(monocle_CTRL_deg)[order(monocle_CTRL_deg$p_val_adj)][1:1000]
monocle_CTRL_1000 <- setOrderingFilter(monocle_CTRL, deg_CTRL_1000)
monocle_CTRL_1000 <- reduceDimension(monocle_CTRL_1000, method = 'DDRTree')
monocle_CTRL_1000 <- orderCells(monocle_CTRL_1000)

## ordering cells by assigning root nodes
GM_state <- function(monocle_CTRL_1000){
  if (length(unique(monocle_CTRL_1000$State)) > 1){
    T0_counts <- table(monocle_CTRL_1000$State, monocle_CTRL_1000$seurat_clusters)[,"13"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
monocle_CTRL_1000 <- orderCells(monocle_CTRL_1000, root_state = GM_state(monocle_CTRL_1000))

plot_cell_trajectory(monocle_CTRL_1000)
ggsave("trajectory_CTRL_state.png", dpi = 600, width = 4.5, height = 5)
ggsave("trajectory_CTRL_state.pdf", dpi = 600, width = 4.5, height = 5)

plot_cell_trajectory(monocle_CTRL_1000, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters) + scale_color_simpsons("springfield", alpha = 0.8) + facet_grid(~seurat_clusters, labeller = as_labeller(name_of_clusters))
ggsave("trajectory_CTRL_splited_by_cluster.png", dpi = 600, width = 20, height = 3.5)
ggsave("trajectory_CTRL_splited_by_cluster.pdf", dpi = 600, width = 20, height = 3.5)

plot_cell_trajectory(monocle_CTRL_1000, color_by = "Pseudotime")
ggsave("trajectory_CTRL_psudotime.png", dpi = 600, width = 4.5, height = 5)
ggsave("trajectory_CTRL_psudotime.pdf", dpi = 600, width = 4.5, height = 5)

monocle_OXd4 <- subset(FRC_paper, subset = orig.ident == "OXd4 FRC")
monocle_OXd4_deg <- FindAllMarkers(monocle_OXd4, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
monocle_OXd4 <- as.CellDataSet(monocle_OXd4)
monocle_OXd4 <- estimateSizeFactors(monocle_OXd4)
monocle_OXd4 <- estimateDispersions(monocle_OXd4)

## select top 1000 genes as filter
deg_OXd4_1000 <- row.names(monocle_OXd4_deg)[order(monocle_OXd4_deg$p_val_adj)][1:1000]
monocle_OXd4_1000 <- setOrderingFilter(monocle_OXd4, deg_OXd4_1000)
monocle_OXd4_1000 <- reduceDimension(monocle_OXd4_1000, method = 'DDRTree')
monocle_OXd4_1000 <- orderCells(monocle_OXd4_1000)

## ordering cells by assigning root nodes
GM_state <- function(monocle_OXd4_1000){
  if (length(unique(monocle_OXd4_1000$State)) > 1){
    T0_counts <- table(monocle_OXd4_1000$State, monocle_OXd4_1000$seurat_clusters)[,"13"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
monocle_OXd4_1000 <- orderCells(monocle_OXd4_1000, root_state = GM_state(monocle_OXd4_1000))

plot_cell_trajectory(monocle_OXd4_1000)
ggsave("trajectory_OXd4_state.png", dpi = 600, width = 4.5, height = 5)
ggsave("trajectory_OXd4_state.pdf", dpi = 600, width = 4.5, height = 5)

plot_cell_trajectory(monocle_OXd4_1000, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters) + facet_grid(~seurat_clusters, labeller = as_labeller(name_of_clusters)) + scale_color_manual(values = mono_palette)
ggsave("trajectory_OXd4_splited_by_cluster.png", dpi = 600, width = 20, height = 3.5)
ggsave("trajectory_OXd4_splited_by_cluster.pdf", dpi = 600, width = 20, height = 3.5)

plot_cell_trajectory(monocle_OXd4_1000, color_by = "Pseudotime")
ggsave("trajectory_OXd4_psudotime.png", dpi = 600, width = 4.5, height = 5)
ggsave("trajectory_OXd4_psudotime.pdf", dpi = 600, width = 4.5, height = 5)

plot_cell_trajectory(monocle_CTRL_1000, markers = "Cd34", 
                     use_color_gradient = T) + scale_color_gradient2 (low="grey", high="purple3")
ggsave("trajectory_CTRL_Cd34.pdf", dpi = 600, width = 4.5, height = 5)

plot_cell_trajectory(monocle_OXd4_1000, markers = "Cd34", 
                     use_color_gradient = T) + scale_color_gradient2 (low="grey", high="purple3")
ggsave("trajectory_OX_Cd34.pdf", dpi = 600, width = 4.5, height = 5)


###### DEG Heatmaps ######
FRC_rps_removed <- FRC_paper[grep("^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa", rownames(FRC_paper), value = T, invert = T), ]
FRC_rps_removed[["percent.ribo"]] <- PercentageFeatureSet(FRC_rps_removed, pattern = "^Rp[sl]")

for (cluster_i in cluster_list) {
  filename_i <- paste("BP_", cluster_i, ".png", sep = "")
  filename_ii <- paste("BP_", cluster_i, ".pdf", sep = "")
  csv_i <- paste("BP_", cluster_i, ".csv", sep = "")
  
  DE_i <- subset(FRC_rps_removed, idents = cluster_i)
  Idents(DE_i) <- "orig.ident"
  OXd4_vs_ctrl_i <- FindAllMarkers(object = DE_i, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, 
                                   test.use = "MAST")
  pval_i <- subset(OXd4_vs_ctrl_i, rowSums(OXd4_vs_ctrl_i[5] < 0.05) > 0)
  df_i <- pval_i[,7:6]
  dfsample_i <- split(df_i$gene, df_i$cluster)
  dfsample_i$`CTRL FRC` = bitr(dfsample_i$`CTRL FRC`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  dfsample_i$`OXd4 FRC` = bitr(dfsample_i$`OXd4 FRC`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  genelist_i <- list("CTRL FRC" = dfsample_i$`CTRL FRC`$ENTREZID,
                     "OXd4 FRC" = dfsample_i$`OXd4 FRC`$ENTREZID)
  
  GO_BP_i <- compareCluster(geneClusters = genelist_i, fun = "enrichGO", OrgDb = "org.Mm.eg.db",
                                 ont = "BP", pAdjustMethod="BH")
  write.csv(GO_BP_i@compareClusterResult, file = csv_i)
  title_i <- paste(cluster_i, sep = "")
  dotplot(GO_BP_i, showCategory = 5, includeAll = T, title = title_i) + scale_size_area(max_size = 12)
  ggsave(filename_i, dpi= 450, width = 6, height = 8)
  ggsave(filename_ii, dpi= 450, width = 6, height = 8)
}


chemokines <- c("Cxcl12", "Cxcl13", "Ccl19", "Ccl21a", "Cxcl1", "Cxcl16", "Tnfsf13b", "Il6", "Il7")
ECM_components <- c("Ecm1", "Ecm2", "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col6a1")
ECM_remodeling <- c("Mmp2", "Mmp9", "Mmp14", "Mmp23", "Timp1", "Timp2", "Timp3")

Cxcl13int_AvgExp <- AggregateExpression(subset(FRC_rps_removed, ident = "Cxcl13inter RC"), 
                                        normalization.method = "LogNormalize")
Dividing_AvgExp <- AggregateExpression(subset(FRC_rps_removed, ident = "dividing_cells"), 
                                       normalization.method = "LogNormalize")

cluster_list <- c("Ccl19lo TRC", "Ccl19hi TRC", "Cxcl9+ TRC", "MRC", "FDC",
                  "MedRC", "Inmt+ SC", "Cd34+ SC", "PvC")

chemokine_all_clusters <- list()
for (cluster_i in cluster_list) {
  AvgExp_i <- AggregateExpression(subset(FRC_rps_removed, ident = cluster_i), group.by = "orig.ident")
  chemokine_all_clusters[[cluster_i]] <- as.matrix(AvgExp_i$RNA[rownames(AvgExp_i$RNA) %in% chemokines, ])
}

chemokine_all_clusters[["Cxcl13inter RC"]] <- as.matrix(Cxcl13int_AvgExp$RNA[rownames(Cxcl13int_AvgExp$RNA) %in% chemokines, ])
chemokine_all_clusters[["dividing_cells"]] <- as.matrix(Dividing_AvgExp$RNA[rownames(Dividing_AvgExp$RNA) %in% chemokines, ])

chemokine_matrix <- data.frame(chemokine_all_clusters)
chemokine_matrix <- chemokine_matrix[order(chemokine_matrix[,1], decreasing = T),]
pheatmap(chemokine_matrix, color = colorRampPalette(c("blue", "white", "red"), space="rgb")(128),
         cluster_cols = F, cluster_rows = F, cellwidth = 15, cellheigth = 15, fontsize_col = 12, fontsize = 15,
         fontsize_row = 12, fontsize_number = 8, main = "Chemokines, cytokines & survival factors", angle_col = "45",
         scale = "row", legend_breaks = -3:3, labels_col = c("", "Ccl19lo TRC"," ", "Ccl19hi TRC","", "Cxcl9+ TRC", "",
                                                             "MRC","", "FDC","", "MedRC","", "Inmt+ SC", " ",
                                                             "Cd34+ SC", "", "PvC", "Cxcl13inter RC", "Dividing cells"), 
         gaps_col = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 19))

ECM_components_all <- list()
for (cluster_i in cluster_list) {
  AvgExp_i <- AggregateExpression(subset(FRC_rps_removed, ident = cluster_i), group.by = "orig.ident")
  ECM_components_all[[cluster_i]] <- as.matrix(AvgExp_i$RNA[rownames(AvgExp_i$RNA) %in% ECM_components, ])
}

ECM_components_all[["Cxcl13inter RC"]] <- as.matrix(Cxcl13int_AvgExp$RNA[rownames(Cxcl13int_AvgExp$RNA) %in% ECM_components, ])
ECM_components_all[["dividing_cells"]] <- as.matrix(Dividing_AvgExp$RNA[rownames(Dividing_AvgExp$RNA) %in% ECM_components, ])

ECM_components_all_matrix <- data.frame(ECM_components_all)
ECM_components_all_matrix <- ECM_components_all_matrix[order(ECM_components_all_matrix[,1], decreasing = T),]
pheatmap(ECM_components_all_matrix, color = colorRampPalette(c("blue", "white", "red"), space="rgb")(128),
         cluster_cols = F, cluster_rows = F, cellwidth = 15, cellheigth = 15, fontsize_col = 12, fontsize = 15,
         fontsize_row = 12, fontsize_number = 8, main = "Conduit Components", angle_col = "45",
         scale = "row", legend_breaks = -3:3, labels_col = c("", "Ccl19lo TRC"," ", "Ccl19hi TRC","", "Cxcl9+ TRC", "",
                                                             "MRC","", "FDC","", "MedRC","", "Inmt+ SC", " ",
                                                             "Cd34+ SC", "", "PvC", "Cxcl13inter RC", "Dividing cells"), 
         gaps_col = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 19))

ECM_remodel_all <- list()
for (cluster_i in cluster_list) {
  AvgExp_i <- AggregateExpression(subset(FRC_rps_removed, ident = cluster_i), group.by = "orig.ident")
  ECM_remodel_all[[cluster_i]] <- as.matrix(AvgExp_i$RNA[rownames(AvgExp_i$RNA) %in% ECM_remodeling, ])
}

ECM_remodel_all[["Cxcl13inter RC"]] <- as.matrix(Cxcl13int_AvgExp$RNA[rownames(Cxcl13int_AvgExp$RNA) %in% ECM_remodeling, ])
ECM_remodel_all[["dividing_cells"]] <- as.matrix(Dividing_AvgExp$RNA[rownames(Dividing_AvgExp$RNA) %in% ECM_remodeling, ])


ECM_remodel_all_matrix <- data.frame(ECM_remodel_all)
pheatmap(ECM_remodel_all_matrix, color = colorRampPalette(c("blue", "white", "red"), space="rgb")(128),
         cluster_cols = F, cluster_rows = F, cellwidth = 15, cellheigth = 15, fontsize_col = 12, fontsize = 15,
         fontsize_row = 12, fontsize_number = 8, main = "ECM Remodeling", angle_col = "45",
         scale = "row", legend_breaks = -3:3, labels_col = c("", "Ccl19lo TRC"," ", "Ccl19hi TRC","", "Cxcl9+ TRC", "",
                                                             "MRC","", "FDC","", "MedRC","", "Inmt+ SC", " ",
                                                             "Cd34+ SC", "", "PvC", "Cxcl13inter RC", "Dividing cells"), 
         gaps_col = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 19))
