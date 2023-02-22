
#### running through ENIGMA testing



library(ggplot2)
library(ENIGMA)
library(stringr)


## load in bulk data

Data <- read.table('~/Bioinformatics/bsn9/bulk_bsn9_counts_protein_coding_0702622.tsv')

Data <- Data[apply(Data, 1, sd) > 0,] ### remove invariant genes

cell_types <- unique(str_split_fixed(colnames(Data), 'r', 2)[,1]) ## get a list of cell types


## load in single cell prop2count data

seurat_pseudobulk_tissue <- read.table('~/Bioinformatics/single_cell_data/Tissue_pseudobulk_prop2count_082022.tsv')


sc_meta = data.frame(row.names = colnames(seurat_pseudobulk_tissue),
                     sample = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,2],
                     cell_type = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,1])


seurat_pseudobulk_tissue_CPM <-  sweep(seurat_pseudobulk_tissue, MARGIN = 2, STATS = colSums(seurat_pseudobulk_tissue), FUN = '/') * 1000000

seurat_pseudobulk_tissue_ave_CPM <- sapply(unique(sc_meta$cell_type), function(x){
  tmp <- seurat_pseudobulk_tissue_CPM[, str_split_fixed(colnames(seurat_pseudobulk_tissue_CPM), '__', 2)[,1]==x]
  return(rowMeans(tmp))
})



#### subset to common genes

enigma.common.genes <- intersect(rownames(Data), rownames(seurat_pseudobulk_tissue))

Data <- Data[enigma.common.genes,]


#### make the ENIGMA object
egm_aggre = ENIGMA::create_ENIGMA(bulk = as.matrix(Data), ref = as.matrix(seurat_pseudobulk_tissue_ave_CPM), ref_type = "aggre")

#egm_aggre = batch_correct(egm_aggre) ## batch correction gives poor AUROC compared to no batch correction

egm_aggre = get_cell_proportion(egm_aggre, method = "RLR")

#egmL2 = ENIGMA_L2_max_norm(egm_aggre, preprocess = 'sqrt') ## sqrt preprocessing does not perform as well as log2(+1) on gene detection
egmL2_log = ENIGMA_L2_max_norm(egm_aggre, preprocess = 'log')

#saveRDS(egmL2_log, 'bulk_bsn9_enigma_L2_norm_log_deconv_082022.rds')


## extract the neuron specific columns
egmL2_log_neuron <- egmL2_log@result_CSE@assays@data$counts[, egmL2_log@result_CSE$cell_type=='Neuron']
colnames(egmL2_log_neuron) <- str_split_fixed(colnames(egmL2_log_neuron), ':', 2)[,1]


write.table(egmL2_log_neuron, 'bulk_bsn9_enigma_L2_norm_log_NeuronalCounts_082022.tsv', sep = '\t')
