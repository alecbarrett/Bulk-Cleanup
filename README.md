# Bulk-Cleanup

An in progress tool for cleaning FACS generated bulk RNA-seq counts profiles. The current iteration is specific to my dataset, I am currently working on an R library implementation.

## Inputs
|Data|Rows|Columns|
|----|-----|-------|
|Bulk|genes|Samples|
|single-cell|genes|cell types|
|ground truth|genes|cell types|

The user must specify a cell type name that corresponds to the sample's target (i.e. glutamatergic neurons), and also provide a list of putative contaminating cell types.

The ground truth dataset needs to be a binary matrix, showing where each gene is expressed. In my trials I have focused on cleaning neuronal bulk RNA-seq samples, and for the ground truth matrix I have used ubiquitously expressed genes (1s for all cell types) and genes known to be expressed only outside the nervous system (0s for all cell types).

All proportion estimates, gene level weights, and subtractions, are performed after a log transformation.

This code iteratively cleans each bulk sample by:
1. Generating a gene level weighting score [0-1] based on its expression pattern.
    * Genes expressed in one cell type and no others get a high weight, genes expressed in a small handful of cells get mid weights, and genes expressed in many or all cells get weights close to 0.
2. Running through an algorithm to clean the bulk dataset:
    1. Estimate the cell type composition of the bulk samples using a non-negative least squares regression, using the gene level weights.
    2. Subtract the contaminant profiles using the single cell counts data, with an additional learning rate constraint (sampling across possible learning rates) (bulk_ = bulk - (single_cell * composition_estimates * learning_rate))
        * If 50 learning rates are selected, then there will be 50 bulk_ profiles created, one for each learning rate
    3. compare the original & subtracted bulk profiles to the ground truth expression matrix, calculate AUROC
    5. Return the profile with the highest AUROC
        * If multiple learning rates gave the highest AUROC, pick the lowest learning rate (if given the option, subtract less)
        * If none of the subtracted profiles have a higher AUROC than the original bulk sample, halt the loop!
        
        
This algorithm uses two metrics as a guide to decide how much to subtract, the composition estimate, and the correspondence to ground truth. If the estimated composition is 100% the target cell type, then no subtraction will occur, and if no subtractions can improve the correspondence to grouth truth, then no subtraction will occur. This provides two influences to ensure that the algorithm does not go overboard in its subtraction, and remove real signal.

