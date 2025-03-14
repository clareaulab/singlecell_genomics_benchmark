library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(stringr)
library(rhdf5)

# https://github.com/slowkow/saturation/blob/main/saturation.R
# This doesn't do any of the filtering
downsample_analysis <- function(pcr_dups_per_molecule_cells, feature_per_molecule_cells, gex_barcode_per_molecule_cells,
                                prob = 0.1, seed_value = 2020, what = "experiment") {
  set.seed(seed_value)
  my_read_counts <- rbinom(n = length(pcr_dups_per_molecule_cells), size = pcr_dups_per_molecule_cells, prob = prob)
  keep <- my_read_counts > 0
  sdt <- data.table(
    pcr_counts = my_read_counts[keep],
    barcodes = gex_barcode_per_molecule_cells[keep],
    gene = feature_per_molecule_cells[keep]
  )
  per_barcode_nUmis <- sdt[, .(N = sum(pcr_counts >= 1), readsPerCell = sum(pcr_counts)), by = list(barcodes)] 
  per_barcode_nGenes <- sdt[, .(N_both = sum(pcr_counts >= 1)), by = list(barcodes,gene)][, .(N = sum(N_both >= 1)), by = list(barcodes)] 
  
  data.frame(
    prob = prob,
    seed = seed_value,
    what = what,
    mean_reads_per_cell = mean(per_barcode_nUmis$readsPerCell),
    median_reads_per_cell = median(per_barcode_nUmis$readsPerCell),
    mean_umis_per_cell = mean(per_barcode_nUmis$N),
    median_umis_per_cell = median(per_barcode_nUmis$N),
    mean_genes_per_cell = mean(per_barcode_nGenes$N),
    median_genes_per_cell = median(per_barcode_nGenes$N)
  )
}

prob_vec <- c(0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 1 , 2.5, 5, 6.69, 10)/10

# Helper function to process stuff
process_sample_10x <- function(molecule_info_file, cell_barcodes, what = "what", probs = prob_vec){

  #Import molecular counts
  pcr_dups_per_molecule <- h5read(molecule_info_file, "count")
  barcode_idx_per_molecule <- h5read(molecule_info_file, "barcode_idx")
  feature_per_molecule <- h5read(molecule_info_file, "feature_idx")
  features <- h5read(molecule_info_file, "features")[["name"]]
  gex_barcode_per_molecule <- h5read(molecule_info_file, "barcodes")[barcode_idx_per_molecule + 1]
  total_reads <- sum(pcr_dups_per_molecule)
  
  # Filter vectors per cell
  in_cells <- gex_barcode_per_molecule %in% cell_barcodes
  feature_per_molecule_cells <- feature_per_molecule[in_cells]
  gex_barcode_per_molecule_cells <- gex_barcode_per_molecule[in_cells]
  pcr_dups_per_molecule_cells <- pcr_dups_per_molecule[in_cells]
  
  # loop over probabilities
  ds_df <- lapply(probs, function(prob){
    print(prob)
    downsample_analysis(pcr_dups_per_molecule_cells, feature_per_molecule_cells, gex_barcode_per_molecule_cells,
                        prob = prob, seed_value = 2020, what = what) %>%
      mutate(raw_reads_per_cell = total_reads*prob/length(cell_barcodes))
    
  }) %>% rbindlist() %>% data.frame()
  ds_df
}

process_sample_pipseq <- function(molecule_info_file, cell_barcodes, what = "what", probs = prob_vec){
  
  #Import molecular counts
  pcr_dups_per_molecule <- h5read(molecule_info_file, "counts")
  barcodes_per_molecule <- h5read(molecule_info_file, "barcodes")
  feature_per_molecule <- h5read(molecule_info_file, "genes")
  features <- h5read(molecule_info_file, "gene_list")
  cell_barcode_per_molecule <-  h5read(molecule_info_file, "barcode_list")[barcodes_per_molecule + 1]
  gex_name_per_molecule <- h5read(molecule_info_file, "gene_list")[feature_per_molecule + 1]
  total_reads <- sum(pcr_dups_per_molecule)
  
  # Filter vectors per cell
  in_cells <- cell_barcode_per_molecule %in% cell_barcodes
  feature_per_molecule_cells <- gex_name_per_molecule[in_cells]
  gex_barcode_per_molecule_cells <- cell_barcode_per_molecule[in_cells]
  pcr_dups_per_molecule_cells <- pcr_dups_per_molecule[in_cells]
  
  # loop over probabilities
  ds_df <- lapply(probs, function(prob){
    print(prob)
    downsample_analysis(pcr_dups_per_molecule_cells, feature_per_molecule_cells, gex_barcode_per_molecule_cells,
                        prob = prob, seed_value = 2020, what = what) %>%
      mutate(raw_reads_per_cell = total_reads*prob/length(cell_barcodes))
  }) %>% rbindlist() %>% data.frame()
  ds_df
}

# Import data 
cell_barcodes2 <- gsub("-1", "", fread("../barcodes/3p_barcodes.tsv.gz", header = FALSE)[[1]])
molecule_info_file2 <- "../../benchmark_large_data/molecule_info/3p_molecule_info.h5"

cell_barcodes3 <- fread("../barcodes/pipseq_sens4_barcodes.tsv.gz", header = FALSE)[[1]]
molecule_info_file3 <- "../../benchmark_large_data/molecule_info/pipseq_molecule_info.h5"

# Process each one
ps_3p <- process_sample_10x(molecule_info_file2, cell_barcodes2, what = "10x_GEMX_3p")
ps_pip <- process_sample_pipseq(molecule_info_file3, cell_barcodes3, what = "PIPseq_v5")

full_ds_df <- rbind(
  ps_3p,ps_pip
)

pA <- ggplot(full_ds_df %>% filter(raw_reads_per_cell < 30000), aes(x = raw_reads_per_cell, y = median_umis_per_cell, color = what)) +
  geom_point(size = 3) + geom_line() +
  pretty_plot(fontsize = 14) + labs(x = 'Mean reads / cell', y = "Median UMIs per cell", color = "technology") + L_border() + 
  scale_color_manual(values = jdb_palette("corona")[c(1:3)]) 
pA
cowplot::ggsave2(pA, file = "../plots/pipseq_gemx_benchmarking.pdf", width = 8 ,height = 5)

# compute downsampling probabilities to lowest comparator library
full_ds_df %>% group_by(what) %>% top_n(1, prob) %>%
  ungroup() %>%
  mutate(prop = min(raw_reads_per_cell)/raw_reads_per_cell) %>% data.frame()

rbind(
  process_sample_10x(molecule_info_file2, cell_barcodes2, what = "10x_GEMX_3p", 1),
  process_sample_pipseq(molecule_info_file3, cell_barcodes3, what = "PIPseq_v5", 0.6697392)
) %>% data.frame()

