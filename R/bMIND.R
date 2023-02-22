## load libraries
library(MIND)
library(data.table)
library(pbapply)
library(dplyr)
library(tidyverse)


## load in bulk data

Data <- read.table('~/project/Hammarlund/singlecell/bulk_bsn9_counts_protein_coding_0702622.tsv')

## load in proportions estimations

proportions <- read.table('~/project/Hammarlund/singlecell/deconvolution_algorithms/bMIND/NNLS_average_across_100_bootstraps.log1p.30Cells.072922.tsv')

## load pseudobulk data

seurat_pseudobulk_tissue <- read.table('~/project/Hammarlund/singlecell/deconvolution_algorithms/bMIND/seurat_pseudobulk_tissue_072922.tsv')


sc_meta = data.frame(row.names = colnames(seurat_pseudobulk_tissue),
                     sample = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,2],
                     cell_type = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,1])

bMIND.common.genes <- intersect(rownames(Data), rownames(seurat_pseudobulk_tissue))

counts_prior <- get_prior(sc = seurat_pseudobulk_tissue, meta_sc = sc_meta, filter_pd = F)


bMIND_genes <- list()
counter = 0
for(gene in bMIND.common.genes){
  
  counter <- counter + 1
  if(counter == 432){return(Null)}else{
    bMIND_genes[[gene]] <-
      tryCatch(
        bMIND(bulk = log2(1+Data[gene,]),
              frac = t(proportions[colnames(counts_prior$profile),]),
              profile = counts_prior$profile),
        error=function(e) NULL)}
  
  print(counter)
  
}

saveRDS(bMIND_genes, 'bMIND_bsn9_deconvolution_list_using_prior_only_080522.RDS')
