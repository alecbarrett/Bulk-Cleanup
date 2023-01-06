#### code for viz package

# Goal: load genexsamples or genesxcell_types matrixes, and compare them to a ground truth, making a PR curve, and rank order bar plots for a few target cell types

# Some manual processing required, but not much

# Main function (Run_GT_viz) takes you from GenesxCell_types matrices to a pptx file

### libraries
library("officer")
library("tidyverse")
library("officer")
library("rvg")
library("purrr")
library("edgeR")


### load in the functions
source('./GT_binary_visualization_source.R')

 



### example usage


### ground truth

all_genes_gt_bulk <- read.csv('bulk_all_ground_truth_042021.csv', row.names = 1)



### load data 


unaltered_bulk_data <- read.table('bsn9_bulk_data_052522.tsv') 

subtracted_bulk_data <- read.table('bulk_ROC_subtracted_prop_count_SPM_log_042122.tsv')

proportions <- read.table('SingleCell_proportions_Bulk_annotations.tsv')

prop2TPM <- read.table('sc_TPM_from_prop_052522.tsv')

dynamicprop2TPM <- read.table('sc_TPM_from_dynamic_prop_052522.tsv')

dynamic_proportions <- dynamic_proportions_normalization(proportions)


## normalize the bulk data

unaltered_bulk_data_TMM <- TMM_normalize_data(unaltered_bulk_data)
aggr_unaltered_bulk_data_TMM <- average_within_cellType(unaltered_bulk_data_TMM)


subtracted_bulk_data_TMM <- TMM_normalize_data(subtracted_bulk_data)
aggr_subtracted_bulk_data_TMM <- average_within_cellType(subtracted_bulk_data_TMM)



## integrate the data
static_integrated <- integrate_by_cellType(aggr_subtracted_bulk_data_TMM, prop2TPM)
dynamic_integrated <- integrate_by_cellType(aggr_subtracted_bulk_data_TMM, dynamicprop2TPM)




dataset_list <- list(aggr_unaltered_bulk_data_TMM, aggr_subtracted_bulk_data_TMM, 
                     static_integrated, dynamic_integrated, proportions, dynamic_proportions) ## 
dataset_names <- c('unaltered bulk', 'subtracted bulk', 
                   'static_integrated', 'dynamic_integrated', 
                   'proportions', 'dynamic_proportions') ### give them names! 


pptx_object <- Run_GT_viz(list_of_datasets = dataset_list, 
                           dataset_names = dataset_names,
                           ground_truth = all_genes_gt_bulk
)

print(pptx_object, 'testing_full_script.pptx')
