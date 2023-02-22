### subtraction using proportion transformed bulk data, and actual proportions from single cell data

# Initializations ----
#library(Biobase)
#library(BisqueRNA)
library(Seurat)
library(stringr)
library(Matrix)
#library(abind)
library(tidyverse)
library(dplyr)
library(bayestestR)
library(nnls)
library(pbapply)
#library(glmnet)
print('loaded libraries')


#### ~define functions~
calc_bulk_auc_one_sample_rownames <- function(bulk_vector, sample_name, training_set_bulk, training_genes, bulk_matrix){

  # log transform the bulk data, then calculate the AUC for the ROC curve.

  ## bulk_vector = genes x samples matrix, needs to be named with gene names
  ## training_set_bulk = genes x bulk_cell_types matrix
  ## training_genes = list of genes


  training_genes <- training_genes[order(training_genes)]
  cell_type <- str_split_fixed(sample_name, 'r', 2)[,1]
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- rownames(bulk_matrix)
  dcnt <- dcnt[training_genes]


  train_bulk <- training_set_bulk[training_genes,]


  train_bulk_single <- train_bulk[,cell_type]
  names(train_bulk_single) <- rownames(train_bulk)
  diags_d <- tibble(threshold = c(seq(0,0.05, 0.0001),seq(0.1,1,0.05)),
                    TPR = map_dbl(threshold, ~get_tpr(dcnt, train_bulk_single, .x)),
                    FPR = map_dbl(threshold, ~get_fpr(dcnt, train_bulk_single, .x)))

  bulk_auc <- auc(diags_d$FPR, diags_d$TPR , method = 'trap')
  return(bulk_auc)
}


### calculate TPR value for a particular threshold
get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}



calc_bulk_auc_one_sample <- function(bulk_vector, sample_name, training_set_bulk, training_genes){

  # log transform the bulk data, then calculate the AUC for the ROC curve.

  ## bulk_vector = genes x samples matrix, needs to be named with gene names
  ## training_set_bulk = genes x bulk_cell_types matrix
  ## training_genes = list of genes

  cell_type <- str_split_fixed(sample_name, 'r', 2)[,1]
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- names(bulk_vector)
  dcnt <- dcnt[training_genes]


  train_bulk <- training_set_bulk[training_genes,]


  train_bulk_single <- train_bulk[,cell_type]
  names(train_bulk_single) <- rownames(train_bulk)
  diags_d <- tibble(threshold = c(0,seq(0,15,0.1)),
                    TPR = map_dbl(threshold, ~get_tpr(dcnt, train_bulk_single, .x)),
                    FPR = map_dbl(threshold, ~get_fpr(dcnt, train_bulk_single, .x)))

  bulk_auc <- auc(rev(diags_d$FPR), rev(diags_d$TPR) , method = 'trap')
  return(bulk_auc)
}

calc_bulk_auc_one_sample_PR <- function(bulk_vector, sample_name, training_set_bulk, training_genes){

  # log transform the bulk data, then calculate the AUC for the ROC curve.

  ## bulk_vector = genes x samples matrix, needs to be named with gene names
  ## training_set_bulk = genes x bulk_cell_types matrix
  ## training_genes = list of genes

  cell_type <- str_split_fixed(sample_name, 'r', 2)[,1]
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- names(bulk_vector)
  dcnt <- dcnt[training_genes]


  train_bulk <- training_set_bulk[training_genes,]


  train_bulk_single <- train_bulk[,cell_type]
  names(train_bulk_single) <- rownames(train_bulk)
  diags_d <- tibble(threshold = c(0,seq(0,15,0.1)),
                    TPR = map_dbl(threshold, ~get_tpr(dcnt, train_bulk_single, .x)),
                    FDR = map_dbl(threshold, ~get_fdr(dcnt, train_bulk_single, .x)))

  bulk_auc <- auc((1-diags_d$FDR), (diags_d$TPR) , method = 'trap')
  return(bulk_auc)
}

sample_AUC <- function(sample_name, sample_matrix, ground_truth){
  cell <- str_split_fixed(sample_name, 'r', 2)[,1]
  ground_truth <- ground_truth[rownames(ground_truth) %in% rownames(sample_matrix), ]
  sample_matrix <- sample_matrix[rownames(ground_truth), ]
  bulk_vector <- sample_matrix[, sample_name]
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- rownames(ground_truth)


  ground_truth_single <- ground_truth[,cell]
  names(ground_truth_single) <- rownames(ground_truth)
  diags_d <- tibble(threshold = c(0,seq(0,15,0.1)),
                    TPR = map_dbl(threshold, ~get_tpr(dcnt, ground_truth_single, .x)),
                    FPR = map_dbl(threshold, ~get_fpr(dcnt, ground_truth_single, .x)))

  bulk_auc <- bayestestR::auc(diags_d$FPR, diags_d$TPR , method = 'trap')
  return(bulk_auc)
}

fHg <- function(x){
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x) !=0)
      {
        p <- x / sum(x)
        res <- -sum(p*log2(p), na.rm=TRUE)
        res <- 1 - (res/log2(length(p))) #Modification: To bring to normalized scale
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}
fPem <- function(x)
{
	if(!all(is.na(x)))
 	{
 		x <- as.matrix(x)
 		x[x<0] <- NA
 		x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
		x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
 		x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score

 		x[x<1] <- 1
 		x <- log10(x)

 		x<- abs(x)
 		#res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max) #choose only the maximal score for each gene
 		res <- x/max(x, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
 	} else {
 		res <- NA
 		print("No data avalable.")
 	}
  	return(res)
}
fTau <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(max(x)!=0)
 			{
 				x <- (1-(x/max(x)))
 				res <- sum(x, na.rm=TRUE)
 				res <- res/(length(x)-1)
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
fSpm <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x) !=0)
 			{
 				spm <- x^2/(x%*%x)
 				res <- max(spm) #Modification:To bring to normalized scale. Choose max
        #res <- spm
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
Spm <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x) !=0)
 			{
 				spm <- x^2/(x%*%x)
 				#res <- max(spm) #Modification:To bring to normalized scale. Choose max
        res <- spm
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}

subtract_single_log <- function(learning_rate, bulk, pseudobulk, proportions_vector, anti_identity_vector, specificity_score){

  ## learning_rate = scalar from 0 to 1
  ## bulk = genes x samples matrix
  ## pseudbulk = genes x sc_cell_types matrix
  ## proportions_table = sc_cell_types x samples matrix
  #pseudobulk <- sc_pseudobulk_use
  sc.log <- log(pseudobulk[,cells_in_use])
  sc.log[sc.log<0]=0
  sc.log <- sc.log * specificity_score

  proportions_vector <- proportions_vector * anti_identity_vector

  subtraction_vector <- (learning_rate * ( as.matrix(sc.log) %*%  as.matrix(proportions_vector)))

  bulk_deconv <- log(bulk) - subtraction_vector
  bulk_deconv <- exp(bulk_deconv)

  return(bulk_deconv)
}

get_best_LR_single_log <- function(bulk_vector, sample_name, pseudobulk, training_set_bulk, training_genes, proportions_table, LR_list, specificity_score){

  # take a list of learning rates, and find the best one as defined by the one with the highest AUC value
  ## bulk = vector of gene values for a single sample
  ## pseudbulk = genes x sc_cell_types matrix
  ## training_set_bulk = genes x bulk_cell_types matrix
  ## training_genes = list of genes
  ## proportions table = sc_cell_types estimates for a single sample
  ## LR_list = list of values to use as the learning rates: ex: 1/(2**seq(0,10,1)) returns a list of 1, 1/2, 1/4, 1/8, etc...

  if((0 %in% LR_list) == F){LR_list <- c(0,LR_list)}
  bulk_deconv_list <- lapply(LR_list, function(LR){
    #print(LR)
    bulk_subtract_df <- subtract_single_log(LR, bulk_vector, pseudobulk, proportions_table, anti_identity_vector, specificity_score)
    bulk_subtract <- bulk_subtract_df[,1]
    names(bulk_subtract) <- rownames(bulk_subtract_df)
    bulk_subtract <- bulk_subtract[order(names(bulk_subtract))]
    return(bulk_subtract)

  })

  auc_list <- lapply(bulk_deconv_list, function(bulk_subtracted){ ## for all elements in bulk_deconv_list, calculate their AUC
    #bulk_subtracted <- bulk_subtracted/sum(bulk_subtracted)
    #bulk_subtracted <- bulk_subtracted * 1000000
    auc_add <- calc_bulk_auc_one_sample(bulk_subtracted, sample_name, training_set_bulk, training_genes)
    #print(auc_add)
    return(auc_add)
  })
  best_LR <- which(unlist(auc_list) == max(unlist(auc_list)), arr.ind = TRUE)
  print(best_LR)
  if(1 %in% best_LR){
      best_LR <- 1
  }
  else { best_LR <- tail(best_LR, n = 1) }
  return(list(bulk_deconv_list[best_LR][[1]], LR_list[best_LR])) ## return a matrix with the highest AUC
  #return(bulk_deconv_list)
}

sc_object <- readRDS('100720_L4_all_cells_Seurat.rds')

sc_object <- sc_object[,sc_object$Tissue != 'Unknown' & sc_object$Tissue != 'Unannotated']

sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DB01')] <- 'DB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DA9')] <- 'DA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VC_4_5')] <- 'VC'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VB01', 'VB02')] <- 'VB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RMD_DV', 'RMD_LR')] <- 'RMD'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RME_DV', 'RME_LR')] <- 'RME'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VA12')] <- 'VA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('IL2_DV', 'IL2_LR')] <- 'IL2'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('AWC_ON', 'AWC_OFF')] <- 'AWC'
sc_object <- sc_object[,!(sc_object$Cell.type %in% c('RIV_stressed', 'SMD_stressed'))]
print('merged single cell clusters')
non_neuronal_list <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')

sc_object$neuron_level <- sc_object$Tissue

sc_object$neuron_level[sc_object$neuron_level=='Neuron'] <- sc_object$Cell.type[sc_object$neuron_level=='Neuron']


prop_by_type <- pbsapply(unique(sc_object$neuron_level), function(cell_type){
  temp <- sc_object@assays$RNA@counts[,sc_object$neuron_level==cell_type]
  rowSums(temp > 0)/ncol(temp)})
prop_by_type <- prop_by_type[order(rownames(prop_by_type)), order(colnames(prop_by_type))]

CeNGEN_TPM <- pbsapply(unique(sc_object$neuron_level), function(cell_type){
  temp <- sc_object@assays$RNA@counts[,sc_object$neuron_level==cell_type]
  temp <- sweep(temp, 2, colSums(temp), '/') * 1000000
  rowMeans(temp)})
CeNGEN_TPM <- CeNGEN_TPM[order(rownames(CeNGEN_TPM)), order(colnames(CeNGEN_TPM))]


all_values_matched <- data.frame(
  TPM = reshape::melt(CeNGEN_TPM[rownames(prop_by_type), colnames(prop_by_type)])[,'value'],
  log_TPM = log(reshape::melt(CeNGEN_TPM[rownames(prop_by_type), colnames(prop_by_type)])[,'value']),
  proportion = reshape::melt(prop_by_type)[,'value'])

all_values_matched[all_values_matched==-Inf] <- NA
all_values_matched[all_values_matched$proportion == 1, ] <- NA
all_values_matched <- na.omit(all_values_matched)


fit3 <- nls(log_TPM ~ ((a*log(proportion/(1-proportion))) + b), data = all_values_matched, start = list(a = 1.5, b = 13))

sum_fit3 <- summary(fit3)

logit_scalar <- sum_fit3$parameters[,1][['a']]
logit_shifter <- sum_fit3$parameters[,1][['b']]


#bulk_data <- read.table('~/Bioinformatics/bsn5/bsn5_counts_110121.tsv', sep = '\t')
bulk_load <- read.table('bsn9_featureCoounts_ab2828_112821.txt', header = T)

colnames(bulk_load)[7:ncol(bulk_load)] <- str_split_fixed(colnames(bulk_load)[7:ncol(bulk_load)], '\\.',6)[,5]


## cut metadata
rownames(bulk_load) <- bulk_load$Geneid

bulk_meta <- bulk_load[, 1:6]

bulk_load <- bulk_load[, 7:ncol(bulk_load)]
bulk_load <- bulk_load[order(rownames(bulk_load)), order(colnames(bulk_load))]

### collapse technical replicates


colnames(bulk_load) <-str_split_fixed(colnames(bulk_load),"t",2)[,1]
bulk_load <- data.frame(vapply(unique(colnames(bulk_load)), function(x)
  rowSums(bulk_load[,colnames(bulk_load)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(bulk_load)) ))
dim(bulk_load)

## cut reference samples and outliers
bulk_load <- bulk_load[,!(str_split_fixed(colnames(bulk_load), 'r', 2)[,1] %in% c('Ref'))]


bulk_load$RICr133 <- NULL
bulk_load$PVMr122 <- NULL
bulk_load$ADFr99 <- NULL
bulk_load$M4r117 <- NULL
bulk_load$AVKr113 <- NULL

ws281 <- data.frame(wbData::wb_load_gene_ids('281'))
rownames(ws281) <- ws281$gene_id
genes_keep <- ws281 %>% dplyr::filter(biotype != 'rRNA_gene' & biotype != 'miRNA_gene' &
                                 biotype != 'piRNA_gene' & biotype != 'piRNA_gene' &
                                 biotype != 'transposable_element_gene' & biotype != 'antisense_lncRNA_gene' &
                                 biotype != 'scRNA_gene' & biotype != 'gene') %>% rownames()

bulk_data <- bulk_load[intersect(genes_keep, rownames(bulk_load)),]


sc_counts_from_prop <- logit_scalar * log(prop_by_type[,colnames(prop_by_type)!='gene_name']/(1-prop_by_type[,colnames(prop_by_type)!='gene_name'])) + logit_shifter

sc_counts_from_prop[sc_counts_from_prop==Inf] <- max(sc_counts_from_prop[sc_counts_from_prop!=Inf])
sc_counts_from_prop <- exp(sc_counts_from_prop)


sc_pseudobulk <- data.frame(sc_counts_from_prop)
sc_pseudobulk <- data.frame(sc_pseudobulk)
sc_pseudobulk <- data.frame(CeNGEN_TPM)

### subset down to just the genes that are in the bulk dataset AND the single cell dataset
common.genes <- base::intersect(rownames(bulk_data), rownames(sc_pseudobulk))

bulk_data <- bulk_data[common.genes,]
sc_pseudobulk <- sc_pseudobulk[common.genes,]



#bulk_prop_est <- 1/( 1+ (exp((logit_shifter - log(bulk_data))/(logit_scalar)) ))

cell_types <- unique(str_split_fixed(colnames(bulk_data), 'r', 2)[,1])

cell_types <- unique(sc_object$neuron_level[sc_object$Tissue == 'Neuron'])
cell_types <- cell_types[order(cell_types)]

### gt genes

sc_gene_thresholds <- readRDS('211028_genes_categorized_by_pattern.rds')

ubiquitous_genes <- sc_gene_thresholds[['ubiquitous']]
non_neuronal_genes <- read.table('nn_genes.tsv')[,1]
#non_neuronal_genes <- sc_gene_thresholds[['nonneuronal']]


gt1 <- data.frame(matrix(data = 1, nrow = length(ubiquitous_genes), ncol = length(cell_types)))
rownames(gt1) <- ubiquitous_genes
colnames(gt1) <- cell_types

gt0 <- data.frame(matrix(data = 0, nrow = length(non_neuronal_genes), ncol = length(cell_types)))
rownames(gt0) <- non_neuronal_genes
colnames(gt0) <- cell_types

gt <- rbind(gt1, gt0)


write.table(gt, 'ubituiqtous_and_nonNeuronal_gt_genes_matrix_042222.tsv')
### calculate pre subtraction AUC

training_genes <- intersect(rownames(gt), common.genes)

training_genes <- sample(training_genes, size = length(training_genes)*0.7, replace=F)


bulk_deconv_use <- bulk_data
#bulk_deconv_cpm <- sweep(bulk_deconv, MARGIN=2, STATS=colSums(bulk_deconv), FUN='/') * 1000000
bulk_deconv_use <- lapply(bulk_deconv_use, as.numeric)
bulk_deconv_use <- do.call(cbind, bulk_deconv_use)
rownames(bulk_deconv_use) <- rownames(bulk_data)
colnames(bulk_deconv_use) <- colnames(bulk_data)

bulk_deconv_use <- bulk_deconv_use[order(rownames(bulk_deconv_use)), order(colnames(bulk_deconv_use))]
bulk_deconv <- bulk_deconv_use
bulk_deconv <- bulk_deconv[order(rownames(bulk_deconv)), order(colnames(bulk_deconv))]

sc_pseudobulk_use <- sc_pseudobulk
sc_pseudobulk_use <- lapply(sc_pseudobulk_use, as.numeric)
sc_pseudobulk_use <- do.call(cbind, sc_pseudobulk_use)
rownames(sc_pseudobulk_use) <- rownames(sc_pseudobulk)
colnames(sc_pseudobulk_use) <- colnames(sc_pseudobulk)
sc_pseudobulk_use <- data.frame(sc_pseudobulk_use)
sc_pseudobulk_use$VD <- sc_pseudobulk_use$VD_DD
sc_pseudobulk_use$DD <- sc_pseudobulk_use$VD_DD
sc_pseudobulk_use$VD_DD <- NULL

sc_pseudobulk_use <- sc_pseudobulk_use[order(rownames(sc_pseudobulk_use)), order(colnames(sc_pseudobulk_use))]


## plan

cell_types <- unique(str_split_fixed(colnames(bulk_deconv), 'r', 2)[,1], )
cell_types <- cell_types[order(cell_types)]

training_set_bulk <- gt[training_genes,]
contaminant_tissues <- c("Intestine", "Hypodermis", "Reproductive", "Muscle_mesoderm", "Excretory", "Glia", "Pharynx")

sc_pseudobulk_use_tmp <- log1p(sc_pseudobulk)
#sc_pseudobulk_use_tmp[sc_pseudobulk_use_tmp < 0] <- 0

#specificity_score <- apply(sc_pseudobulk_use_tmp, 1, fHg)
#specificity_score <- apply(sc_pseudobulk_use_tmp, 1, fTau)
specificity_score <- apply(sc_pseudobulk_use_tmp, 1, fSpm)

summary(specificity_score[intersect(ubiquitous_genes, common.genes)])

summary(specificity_score[intersect(non_neuronal_genes, common.genes)])



total_rounds_per_sample <- numeric(0)
total_rounds_per_sample <- c(total_rounds_per_sample, 1:ncol(bulk_deconv))
names(total_rounds_per_sample) <- colnames(bulk_deconv)
total_improvement_per_sample <- numeric(0)
total_improvement_per_sample <- c(total_improvement_per_sample, 1:ncol(bulk_deconv))
names(total_improvement_per_sample) <- colnames(bulk_deconv)


step_sizes.df <- data.frame(matrix(nrow = 50, ncol = ncol(bulk_deconv)))
colnames(step_sizes.df) <- colnames(bulk_deconv)
step_sizes.df[is.na(step_sizes.df)] <- 0

auc_before <- pbsapply(colnames(bulk_data), function(sample_name){
  calc_bulk_auc_one_sample_rownames(bulk_data[,sample_name], sample_name, training_set_bulk, training_genes, bulk_data)
})
auc_after <- pbsapply(colnames(bulk_data), function(sample_name){
  calc_bulk_auc_one_sample_rownames(bulk_deconv[,sample_name], sample_name, training_set_bulk, training_genes, bulk_deconv)
})
auc_afterhg <- pbsapply(colnames(bulk_data), function(sample_name){
  calc_bulk_auc_one_sample_rownames(bulk_deconv[,sample_name], sample_name, training_set_bulk, training_genes, bulk_deconv)
})


bulk_deconv <- bulk_data
bulk_deconv <- bulk_deconv[order(rownames(bulk_deconv)), order(colnames(bulk_deconv))]


for(sample_1 in colnames(bulk_deconv)){
    print(sample_1)
    sample_step_sizes <- list()
    for(i in seq(0,100,1)){

          if(i==0){bulk_deconv_target <- bulk_deconv[,sample_1] ## some steps required a dataframe
          names(bulk_deconv_target) <- rownames(bulk_deconv)
          starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target, sample_1, training_set_bulk, training_genes)
          pre_auc <- starting_auc
          print(c('starting AUC', starting_auc))
          cell <- strsplit(sample_1, 'r')[[1]][1]}

          else{
            print(c('iteration',i))
            cells_in_use <- c(cell, contaminant_tissues)
            sc_pseudobulk_use_tmp <- sc_pseudobulk_use[,cells_in_use]
            sc_pseudobulk_use_tmp <- sweep(sc_pseudobulk_use_tmp, 2, colSums(sc_pseudobulk_use_tmp), '/')
            sc_pseudobulk_use_tmp <- sc_pseudobulk_use_tmp * sum(bulk_deconv_target)
            estimates <- nnls( A = log1p(as.matrix(sc_pseudobulk_use_tmp * specificity_score )),
                               b = log1p(bulk_deconv_target * specificity_score ) )$x
            names(estimates) <- cells_in_use
            estimates <- estimates/sum(estimates)
            print(estimates)

            print('estimate done')


            #sc_pseudobulk_use_tmplog <- log(sc_pseudobulk_use[,cells_in_use])
            #sc_pseudobulk_use_tmplog[sc_pseudobulk_use_tmplog < 0] <- 0

            #specificity_score <- data.frame(fPem(sc_pseudobulk_use_tmplog))
            specificity_score_use <- specificity_score[rownames(sc_pseudobulk_use_tmp)]


            anti_identity_vector <- (names(estimates) %in% contaminant_tissues) * 1

            bulk_deconv_target <- get_best_LR_single_log(bulk_deconv_target, sample_1, sc_pseudobulk_use_tmp, training_set_bulk,training_genes, data.frame(estimates), 40/(2**seq(0,20,1)), specificity_score_use)

            sample_step_sizes[[i]] <- bulk_deconv_target[[2]]
            bulk_deconv_target <- bulk_deconv_target[[1]]
            print(bulk_deconv_target[1:5])
            print(c('subtracted', calc_bulk_auc_one_sample(bulk_deconv_target, sample_1, training_set_bulk, training_genes)))




            post_auc <- calc_bulk_auc_one_sample(bulk_deconv_target, sample_1, training_set_bulk, training_genes)

            print(c(i, post_auc))


            if(post_auc == pre_auc){
                print(c('percent improvement -->',(post_auc - starting_auc)*100))
                bulk_deconv[,sample_1] <- bulk_deconv_target
                #total_rounds_per_sample[sample_1] <- i
                #total_improvement_per_sample[sample_1] <- ((post_auc - starting_auc)*100)
                #step_sizes[,]
                for(g in seq(1, length(sample_step_sizes), 1)){step_sizes.df[g,sample_1] <- sample_step_sizes[[g]] }
                #step_sizes.df[[sample_1]] <- sample_step_sizes

                sc_pseudobulk_use_tmp <- sc_pseudobulk_use[,cells_in_use]
                sc_pseudobulk_use_tmp <- sweep(sc_pseudobulk_use_tmp, 2, colSums(sc_pseudobulk_use_tmp), '/')
                sc_pseudobulk_use_tmp <- sc_pseudobulk_use_tmp * sum(bulk_deconv_target)
                estimates <- nnls( A=log1p(as.matrix(sc_pseudobulk_use_tmp * specificity_score )),
                                   b=log1p(bulk_deconv_target * specificity_score ) )$x
                names(estimates) <- cells_in_use
                estimates <- estimates/sum(estimates)
                print(estimates)

                break}
            pre_auc <- post_auc

            cat('\n\n\n\n\n\n\n\n\n') }
    }
}


write.table(bulk_deconv, 'bulk_ROC_subtracted_prop_count_SPM_log_042122.tsv', sep='\t')
