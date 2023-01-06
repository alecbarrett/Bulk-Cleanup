#### functions for ground truth vizualization package


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

######### more functions
integrate_by_cellType <- function(dataset_1, dataset_2, pseudocount = 1){
  
  common.genes <- intersect(rownames(dataset_1), rownames(dataset_2))
  common.cells <- intersect(colnames(dataset_1), colnames(dataset_2))
  
  dataset_1 <- dataset_1[common.genes, common.cells]
  dataset_2 <- dataset_2[common.genes, common.cells]
  
  integrated <- exp( ( log(dataset_1 + pseudocount) + log(dataset_2 + pseudocount) )/2) - pseudocount
  
}


dynamic_proportions_normalization <- function(proportions, low_hard = 0.02, high_hard = 0.01){
  
  proportions[apply(proportions, 1, min) > high_hard, ] <- 1 
  proportions[apply(proportions, 1, max) < low_hard, ] <- 0 
  mid_range <- (apply(proportions, 1, max) > low_hard) & (apply(proportions, 1, min) < high_hard)
  proportions[mid_range, ] <- 
    proportions[mid_range, ]/apply(proportions[mid_range,], 1, max)
  
  return(proportions)
  
}

TMM_normalize_data <- function(data, sep='r'){
  group <- str_split_fixed(colnames(data), sep, 2)[,1]
  TMM <- edgeR::DGEList(data, group = group)
  TMM <- edgeR::calcNormFactors(TMM)
  TMM <- edgeR::cpm(TMM)
  return(TMM)
}

average_within_cellType <- function(data, sep='r'){
  
  tmp <- data
  colnames(tmp) <-str_split_fixed(colnames(tmp), sep, 2)[,1]
  tmp <- data.frame(vapply(unique(colnames(tmp)), function(x) 
    rowMeans(tmp[,colnames(tmp)== x,drop=FALSE], na.rm=TRUE),
    numeric(nrow(tmp)) ))
  
  return(tmp)
  
}

Conform_datasets <- function(list_of_datasets, dataset_names){
  
  total_datasets <- length(list_of_datasets)
  all_genes <- unlist(lapply(list_of_datasets, rownames))
  common.genes <- names(table(all_genes)[table(all_genes)==total_datasets])
  
  all_cells <- unlist(lapply(list_of_datasets, colnames))
  common.cellTypes <- names(table(all_cells)[table(all_cells)==total_datasets])
  
  new_list <- lapply(list_of_datasets, function(data){
    
    data <- data[common.genes, common.cellTypes]
    
  })
  names(new_list) <- dataset_names
  
  return(new_list)
  
}


return_gt_metrics <- function(cell_level_data, ground_truth, source_name){
  
  common.genes <- intersect(rownames(cell_level_data), rownames(ground_truth))
  
  neurons <- intersect(colnames(cell_level_data), colnames(ground_truth))
  
  data_plot <- cell_level_data[common.genes, neurons]
  
  gt <- ground_truth[common.genes, neurons]
  
  metrics_tibble <- tibble(threshold = c(0,2**seq(-10,12,0.1)),
                           TPR = map_dbl(threshold, ~get_tpr(data_plot, gt, .x)),
                           FPR = map_dbl(threshold, ~get_fpr(data_plot, gt, .x)),
                           FDR = map_dbl(threshold, ~get_fdr(data_plot, gt, .x)),
                           source_name = source_name)
  return(metrics_tibble)
  
}

Plot_PR_curve <- function(list_of_gt_metrics){
  
  auc_list <- lapply(names(list_of_gt_metrics), function(source_name){
    tib <- list_of_gt_metrics[[source_name]]
    auc <- bayestestR::auc(1-tib$FDR, tib$TPR)
    auc <- round(auc, digits = 3)
    return(paste(source_name, auc, sep = ': '))
  })
  
  auc_list <- paste(auc_list, collapse = '\n')
  
  merged_list <- do.call(rbind, list_of_gt_metrics)
  
  p <- ggplot(data = merged_list, aes(x = 1-FDR, y=TPR, color= source_name)) +
    geom_path(size = 3) +
    annotate("text", label = auc_list, x = (1-max(merged_list$FDR))+.25, y = .25, fontface = 'bold', size = 6) +
    theme_classic(base_size = 20) 
  
  return(p) 
}

return_gt_by_cellType_plot <- function(cell_level_data, ground_truth, source_name, cell_name){
  
  if(max(cell_level_data > 10)){ cell_level_data <- log1p(cell_level_data)}
  
  data.frame(value = cell_level_data[rownames(ground_truth),cell_name], 
             gt = ground_truth[rownames(ground_truth), cell_name],
             gene = rownames(ground_truth)) |>
    arrange(value) |> 
    mutate(gene=factor(gene, levels=gene)) |>
    ggplot(aes(y = gene, x = value, fill = gt, color = gt, label = gt)) +
    
    geom_bar(stat = 'identity') +
    geom_text() +
    
    annotate(geom = 'text', label = paste(sum(cell_level_data[rownames(ground_truth),cell_name] == 0 & ground_truth[,cell_name] == 1), 'of',
                                          sum(ground_truth[,cell_name]), 'genes not detected'), 
             x = max(cell_level_data[rownames(ground_truth),cell_name])/2, y = 10, hjust = 1, vjust = 1) +
    
    labs(title = paste(cell_name, source_name)) +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold', size = 25))
  
}

Plot_gt_by_cellType <- function(list_of_datasets, ground_truth, list_of_cells){
  plot_list <- lapply(list_of_cells, function(cell_name){
    lapply(names(list_of_datasets), function(source_name){
      
      cell_level_data <- list_of_datasets[[source_name]]
      
      return_gt_by_cellType_plot(cell_level_data=cell_level_data, ground_truth, source_name, cell_name)
      
      
    })
  })
  
  return(do.call(c, plot_list))
  
}


make_pptx_object <- function(list_of_plots){
  
  myppt <- read_pptx()
  mylab <- layout_summary(myppt)[[1]] # Slide Layout Name
  mytmp <- layout_summary(myppt)[[2]][1]
  
  for(plt in list_of_plots){
    myppt <- myppt |>
      add_slide(x = _, layout="Title Slide", master=mytmp) |>
      ph_with(value=plt, location = ph_location_fullsize())
    
  }
  
  return(myppt)
}




### the actual function to run
Run_GT_viz <- function(list_of_datasets,
                       dataset_names,
                       ground_truth,
                       list_of_cells = c('ASG', 'OLQ', 'PVD'),
                       pptx_name = paste(Sys.Date(), '.pptx')){
  
  list_of_datasets <- Conform_datasets(list_of_datasets, dataset_names,  ground_truth)
  
  list_of_gt_metrics <- lapply(dataset_names, function(source_name){
    cell_level_data <- list_of_datasets[[source_name]]
    return_gt_metrics(cell_level_data, ground_truth, source_name)
  })
  names(list_of_gt_metrics) <- dataset_names
  
  PR_plot <- Plot_PR_curve(list_of_gt_metrics)
  
  gt_by_cellType <- Plot_gt_by_cellType(list_of_datasets, ground_truth, list_of_cells)
  
  
  all_plots <- c(list(PR_plot), gt_by_cellType)
  
  make_pptx_object(all_plots)
}

