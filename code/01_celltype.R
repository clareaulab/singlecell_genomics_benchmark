library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat("../../benchmark_large_data/pbmc_multimodal.h5seurat")

txg <- Read10X_h5("../../benchmark_large_data/gex/txg_3p_filtered_feature_bc_matrix.h5")
pip <- Read10X("../../benchmark_large_data/gex/pip_sens4/")

project <- function(raw){
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  return(df)
}

df_pip <- project(pip)
df_txg <- project(txg)

mdf <- merge(
  df_pip %>% group_by(predicted.celltype.l2) %>% 
    summarize(count = n()) %>% mutate(prop_pip = count/sum(count)*100),
  df_txg %>% group_by(predicted.celltype.l2) %>% 
    summarize(count = n()) %>% mutate(prop_txg = count/sum(count)*100), by = "predicted.celltype.l2", all = TRUE
) 
saveRDS(mdf ,file = "../output/both_classifications.rds")
mdf <- readRDS("../output/both_classifications.rds")

library(ggrepel)
mdf[complete.cases(mdf), ] %>%
  ggplot(aes(x = prop_txg, y = prop_pip, label = predicted.celltype.l2)) +
  geom_text_repel() + geom_point() + theme_bw() +
  ggtitle("Azimuth L2 Cell types") + 
  labs(x = "% cells with annotation, 3' GEM-X", y= "% cells with annotation, PIPseq V")
