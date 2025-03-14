library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(stringr)
library(harmony)

# Import data
txg <- Read10X_h5("../../benchmark_large_data/gex/txg_3p_filtered_feature_bc_matrix.h5")
pip <- Read10X("../../benchmark_large_data/gex/pip_sens4/")

# Define a common feature set
featureset <- intersect(rownames(txg), rownames(pip))
so <- CreateSeuratObject(cbind(txg[featureset,],pip[featureset,]))
so$technology <- c(
  rep("3' GEM-X", dim(txg)[2]),
  rep("PIP-seq V", dim(pip)[2])
)

data.frame(
  featureset,
  txg = rowSums(txg[featureset,])/sum(txg[featureset,])*1000000,
  pip = rowSums(pip[featureset,])/sum(pip[featureset,])*1000000
) -> count_df

count_df %>%
  ggplot(aes(x = txg + 1, y = pip + 1)) +
  geom_point(alpha = 0.1, color = "dodgerblue4") + 
  scale_x_log10() + scale_y_log10() +
  theme_bw()+
  labs(x = " GEM-X 3p expression (CPM, log10 scaled)", y = "PIP-seq V expression (CPM, log10 scaled)") +
  ggtitle("Pseudobulk gene expression")

cor(log1p(count_df$txg), log1p(count_df$pip))

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(., c("technology")) %>%
  RunUMAP(dims = 1:30, reduction = "harmony")

so <- so %>% FindNeighbors(reduction = "harmony")
so <- so %>% FindClusters(resolution = 0.7)
DimPlot(so, label = FALSE, group.by = c( "technology"), shuffle = TRUE) &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_color_manual(values = c("dodgerblue2", "firebrick"))

DimPlot(so, group.by = "seurat_clusters", label = TRUE)
FindMarkers(so, "18", "1", only.pos = TRUE)
FindMarkers(so, "3' GEM-X",
            "PIP-seq V",
            group.by = "technology")

cells.order.random <- sample(Cells(so))

FeaturePlot(so, features = c("CD3D", "CD8A", "MS4A1", "FCGR3A", "NKG7", "PPBP"),  max.cutoff = "q90", ncol = 3) &
  labs(x = "UMAP1", y = "UMAP2") & 
  scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))

FeaturePlot(so, features = c("MATR3", "TAB2","DHFR", "RBFOX2"),  max.cutoff = "q90", split.by = "technology") &
  labs(x = "UMAP1", y = "UMAP2") & 
  scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))

