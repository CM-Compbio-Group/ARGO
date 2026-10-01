Differential testing along spatio-temporal trajectories for sequencing data

Visweswaran Ravikumar and Cong Ma

The Gilbert S. Omenn Department of Computational Medicine and Bioinformatics, University of Michigan, Ann Arbor

Modern high-throughput sequencing technologies enable biological processes to be profiled at high resolution, producing measurements of gene expression and other biological state variables along spatial or (pseudo-)temporal trajectories. Although many statistical tests have been developed for comparing bulk and single-cell gene expression and cell type abundance data, they assume discrete experimental units such as treatment conditions or cell-type clusters. Modeling these variables as continuous signals along spatio-temporal trajectories, and performing differential testing on them, remains a difficult, unresolved problem. Here we present Argo, our statistical framework for flexibly modeling biological measurements along continuous trajectories and performing differential testing across multiple conditions. In contrast with standard differential testing tools, our algorithm is capable of both identifying significantly varying genes or cell-types as well as precisely identifying intervals along the trajectory where the deviance occurs. In modeling spatial datasets, this enables identifying particular tissue domains where expression of marker genes or a canonical cell type is varying between conditions.  


## Brief Description of Project :   
 
We present ARGO, a statistical tool for performing differential testing along continuous spatial and (pseudo)-temporal trajectories with single cell gene expression datasets. Details about how the contrast tests are performed, test statistics computed and significance levels derived will be added to this repository in the near future. For now, we can treat the algorithm as a black-box with the following attributes:

**inputs** : 
 - gene expression matrix from a single-cell or spatial transcriptomics dataset
 - an ordering of cells along some continuous underlying trajectory, such as their inferred pseudotime values.
 - design matrix for the experiment -- what samples belong to which conditions and contrast vector that must be tested 

**outputs** : 
 - p-value for a global hypothesis test for each gene that expression does not vary significantly between conditions anywhere along the trajectory
 - a level-$\alpha$ confidence band along the trajectory, that allows us to infer intervals where significant deviation between conditions is identified. 


The model internally uses gaussian processes to fit expression data for each gene, and uses the intrinsic co-regionalization model to derive test statistics for the hypothesis tests. The algorithm can also be used to perform differential abundance analysis along the trajectories, and these applications will be added to the repository in the near future.  


## Example run on simulated dataset

This repository presents a reproducible example of benchmarking performance of ARGO with existing differential gene expression analysis tools using simulated single-cell counts data using scDesign3. We work with the sample pancreatic endocrinogenesis single cell trajectory data presented in the package vignette to generate ground-truth perturbation data between two conditions, and evaluate performance of the different algorithms in this data. The repository contains a Makefile to run this minimal example, where we generate the perturbed counts datasets, run differential testing across seven different methods on this data and compute precision-recall and specificity / sensitivity metrics and finally make plots of AUPRC and F1 scores for different perturbation parameters. 

All scripts and results for this analysis are contained in the testing/ folder, with output figures stored in *testing/figs/* 

 
## Setting up Compute Environment
 
The necessary computational environment can be created by running `make setup` and the described analysis can be run by calling `make all`. 


## Project Structure
The repository structure is captured by this tree:

```text
.  
├── data     
│   ├── raw     
│   ├── sheets     
│   └── temp     
├── figs     
├── output    
│   ├── enrichr     
│   │   ├── diffx     
│   │   └── pam_shape_cluster     
│   └── simulations     
├── references    
├── renv    
├── src     
│   ├── argo     
│   ├── benchmark     
│   ├── plotting     
│   ├── preprocessing     
│   ├── R    
│   └── tools    
└── testing     
    ├── data    
    │   ├── raw    
    │   └── temp    
    ├── figs    
    ├── output    
    │   └── PANCREAS    
    └── src    
        ├── R    
        ├── scripts     
        └── tests

```

## References
[1] https://www.nature.com/articles/s41587-023-01772-1  
[2] https://songdongyuan1994.github.io/scDesign3/docs/articles/scDesign3-introduction-vignette.html

