# ARGO
differential expression testing along spatio-temporal trajectories using gaussian processes


We present ARGO, a statistical tool for performing differential testing along continuous spatial and (pseudo)-temporal trajectories with single cell gene expression datasets. Details about how the contrast tests are performed, test statistics computed and significance levels derived will be added to this repository in the near future. For now, we can treat the algorithm as a black-box with the following attributes:

*inputs* : 
 - gene expression matrix from a single-cell or spatial transcriptomics dataset
 - an ordering of cells along some continuous underlying trajectory, such as their inferred pseudotime values.
 - design matrix for the experiment -- what samples belong to which conditions and contrast vector that must be tested 

*output* : 
 - p-value for a global hypothesis test for each gene that expression does not vary significantly between conditions anywhere along the trajectory
 - a level-$\alpha$ confidence band along the trajectory, that allows us to infer intervals where significant deviation between conditions is identified. 


The model internally uses gaussian processes to fit expression data for each gene, and uses the intrinsic co-regionalization model to derive test statistics for the hypothesis tests. 

This repository presents a reproducible example of benchmarking performance of ARGO with existing differential gene expression analysis tools using simulated single-cell counts data using scDesign3[1]. We work with the sample pancreatic endocrinogenesis single cell trajectory data presented in the package vignette[2] to generate ground-truth perturbation data between two conditions, and evaluate performance of the different algorithms in this data. The repository contains a Makefile to run this minimal example, where we generate the perturbed counts datasets, run differential testing across seven different methods on this data and compute precision-recall and specificity / sensitivity metrics and finally make plots of AUPRC and F1 scores for different perturbation parameters. 

All scripts and results for this analysis are contained in the testing/ folder, with output figures stored in testing/figs/ 

The algorithm can also be used to perform differential abundance analysis along the trajectories, and these applications will be added to the repository in the near future. 
  
The necessary computational environment can be created by running make setup and the described analysis can be run by calling make all. 

The repository structure is captured by this tree:

.
├── data
│   ├── raw
│   ├── sheets
│   └── temp
├── figs
├── output
│   ├── enrichr
│   │   ├── diffx
│   │   └── pam_shape_cluster
│   └── simulations
├── references
├── renv
├── src
│   ├── argo
│   ├── benchmark
│   ├── plotting
│   ├── preprocessing
│   ├── R
│   └── tools
└── testing
    ├── data
    │   ├── raw
    │   └── temp
    ├── figs
    ├── output
    │   └── PANCREAS
    └── src
        ├── R
        ├── scripts
        └── tests




## References
[1] https://www.nature.com/articles/s41587-023-01772-1
[2] https://songdongyuan1994.github.io/scDesign3/docs/articles/scDesign3-introduction-vignette.html 
