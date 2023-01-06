#### libraries

library(SCADIE)
library(dplyr)
library(stringr)
library(edgeR)
library(data.table)
library(NMF)



##### put in the modified functions (modified just to be more verbose)

update_H <- function(W_input,Y){
  
  output<-NMF::fcnnls(x=W_input,y =Y )
  output_H <-output$x%*%diag(1/apply(output$x,MARGIN = 2,sum))
  return( output_H )
}

update_H_dwls <- function(W_input,Y){
  
  output_H <- matrix(NA,nrow=ncol(W_input),ncol=ncol(Y))
  for (i in 1:ncol(output_H)){
    output_H[,i] <- solveDampenedWLS(S =as.matrix(W_input),B=as.matrix(Y[,i],ncol=1) )
  }
  
  return( output_H )
}


update_W <-function(H_input,Y){
  
  if ( (length(which( duplicated(Y, MARGIN = 2) ))>0)|(length(which( duplicated(Y, MARGIN = 1) ))>0)|length(which(apply(Y,2,sd)==0))>0| length(which(apply(Y,1,sd)==0)) ){
    stop("Duplicated or constant cols/rows in Y matrix, please check!")
  }
  
  output<-NMF::fcnnls(x=t(H_input),y =t(Y))
  
  return( t(output$x) )
}

update_W_bs <-function(H_input,Y){
  
  try(output<-NMF::fcnnls(x=t(H_input),y =t(Y)), silent=TRUE)
  while (exists("output") == FALSE){
    n = nrow(t(Y)); p=ncol(t(Y)); yz = t(Y) + matrix(rnorm(n*p,0,0.1), nrow=n)
    try(output<-NMF::fcnnls(x=t(H_input),y =yz),silent=TRUE)
  }
  return( t(output$x) )
}

SIR_itr_general_mod <- function (ini_H_adjusted, ini_W_adjusted, bulk_expr_sub, bulk_expr_full, 
                                 n_ct, itr = 200, H_update_method = "NNLS", H_update_gene = "all", 
                                 signature_gene_row_index, duplicated_rows = F) {
  H_tmp <- ini_H_adjusted
  W_tmp <- ini_W_adjusted
  for (z in 2:itr) {
    print(z)
    if (duplicated_rows == T) {
      W_tmp <- update_W_bs(H_input = H_tmp, Y = bulk_expr_full)
    }
    else {
      W_tmp <- update_W(H_input = H_tmp, Y = bulk_expr_full)
    }
    if (H_update_method == "NNLS") {
      if (H_update_gene == "all") {
        H_tmp <- update_H(W_input = W_tmp, Y = bulk_expr_full)
      }
      else if (H_update_gene == "signature") {
        H_tmp <- update_H(W_input = W_tmp[signature_gene_row_index, 
        ], Y = bulk_expr_full[signature_gene_row_index, 
        ])
      }
      stopifnot(H_update_gene == "all" | H_update_gene == 
                  "signature")
    }
    else if (H_update_method == "DWLS") {
      if (H_update_gene == "all") {
        H_tmp <- update_H_dwls(W_input = W_tmp, Y = bulk_expr_full)
      }
      else if (H_update_gene == "signature") {
        H_tmp <- update_H_dwls(W_input = W_tmp[signature_gene_row_index, 
        ], Y = bulk_expr_full[signature_gene_row_index, 
        ])
      }
      stopifnot(H_update_gene == "all" | H_update_gene == 
                  "signature")
    }
    stopifnot(H_update_method == "NNLS" | H_update_method == 
                "DWLS")
  }
  return(list(W_end = W_tmp, H_end = H_tmp, W_ini = ini_W_adjusted, 
              H_ini = ini_H_adjusted))
}



### load the bulk data


Data <- read.table('bulk_bsn9_counts_protein_coding_0702622.tsv')
unique(str_split_fixed(colnames(Data), 'r', 2)[,1])

## load the proportions estimates

proportions <- read.table('NNLS_average_across_100_bootstraps.log1p.30Cells.072922.tsv')


cell_types <- unique(str_split_fixed(colnames(Data), 'r', 2)[, 1])

Data_sub <- Data[-which(apply(Data,1,sd)==0),]

Data_sub <- log1p(Data_sub)

SCADIE_NNLS <- sapply(cell_types, function(x){
  
  
  SCADIE_bulk <- list()
  
  ## add the bulk matrix for cell to the list 
  cell_index <- str_split_fixed(colnames(Data_sub), 'r', 2)[, 1] == x
  cat(x)
  cat('\n')
  SCADIE_bulk$bulk_full <- as.matrix(Data_sub[,cell_index])
  
  SCADIE_bulk$bulk_sub <- SCADIE_bulk$bulk_full[apply(SCADIE_bulk$bulk_full , 1, sd) > 0,]
  
  cat('removed invariant genes\n')
  
  ## Next, we include the initial $H$s for each group, again, make sure they are matrices:
  SCADIE_bulk$initial_H <- as.matrix(proportions[,cell_index])
  
  cat('loaded proportions estimate\n')
  ## calculate initial Ws
  SCADIE_bulk$initial_W <- t((fcnnls( x = t(SCADIE_bulk$initial_H),y= t(SCADIE_bulk$bulk_sub)  ) )$x) 
  cat('calcuated initial W\n')
  
  NNLS_tmp <- SIR_itr_general_mod(ini_H_adjusted = SCADIE_bulk$initial_H, 
                                  ini_W_adjusted = SCADIE_bulk$initial_W, 
                                  bulk_expr_full = SCADIE_bulk$bulk_sub,
                                  bulk_expr_sub =  SCADIE_bulk$bulk_sub,
                                  n_ct = ncol(SCADIE_bulk$initial_W), 
                                  itr = 100, 
                                  H_update_method = 'NNLS', 
                                  H_update_gene = 'all',
                                  signature_gene_row_index = NA,
                                  duplicated_rows = T)
  
  return(NNLS_tmp$W_end)
  
})

saveRDS(SCADIE_NNLS, 'bulk_bsn9_SCADIE_NNLS_log1p_080822.rds')



### convert to a single matrix for neuronal data

seurat_pseudobulk_tissue <- read.table('~/Bioinformatics/single_cell_data/seurat_pseudobulk_tissue_072922.tsv')


sc_meta = data.frame(row.names = colnames(seurat_pseudobulk_tissue),
                     sample = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,2],
                     cell_type = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,1])


SCADIE_NNLS_neuron <- sapply(SCADIE_NNLS, function(x){
  
  neuron <- x[,'Neuron']
  return(neuron)
})



SCADIE_NNLS.df <- data.frame(row.names = rownames(Data),
                             matrix(data = 0,
                                    nrow = nrow(Data),
                                    ncol = length(cell_types)))
colnames(SCADIE_NNLS.df) <- cell_types

for(cell in cell_types){
  tmp <- SCADIE_NNLS_neuron[[cell]]
  SCADIE_NNLS.df[names(tmp), cell] <- SCADIE_NNLS_neuron[[cell]]
}

SCADIE_NNLS.df <- ecpm1(SCADIE_NNLS.df)
SCADIE_NNLS.df[SCADIE_NNLS.df < 0] <- 0


write.table(SCADIE_NNLS.df, 'bulk_bsn9_SCADIE_NNLS_neuronal_counts_082022.tsv', sep = '\t')
