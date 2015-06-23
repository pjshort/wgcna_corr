# wgcna_corr
Using Weighted Gene Coexpression Network Analysis (WGNCA) to analyze microarray data for correlations between gene modules (groups of genes with similar expression profiles) and phenotypes.

# Soft Thresholding Parameter
This script should be used to set the hyperparameters to build the gene co-expression network. Essentially, a several soft thresholding parameters $k$ are tested to determine which parameters cause the resulting network to best approximate a scale-free network. It has been shown that biological networks (and many other networks such as social networks) exhibit a scale-free (also known as power law) topology.

