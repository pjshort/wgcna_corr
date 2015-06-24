# wgcna_corr
Using Weighted Gene Coexpression Network Analysis (WGNCA) to analyze microarray data for correlations between gene modules (groups of genes with similar expression profiles) and phenotypes.

# Soft Thresholding Parameter
```bash
Rscript thresholding_parameters.R --cel_path=/path/to/cel/files --out=/directory/to/save/plots/and/RData
```
This script should be used to set the hyperparameters to build the gene co-expression network. Essentially, a several soft thresholding parameters $k$ are tested to determine which parameters cause the resulting network to best approximate a scale-free network. It has been shown that biological networks (and many other networks such as social networks) exhibit a scale-free (also known as power law) topology.

_Important output to keep:_ Thresholding parameter $k$ and ExpressionSet object saved as an RData file. These will be used in downstream analyses.

# Defining Gene Modules

```bash
Rscript generate_modules.R --eset=/path/to/expressionset/Rdata/file --soft_thresh_k=integer --out_dir=/directory/to/save/plots/and/module/assignment
```

This script will generate gene -> gene module assignment by building a network, clustering into a dendrogram, and using DynamicTreeCut to slice the dendrogram branches into module assignments. A tab-delimited text file will be produced that maps genes to fine-grained modules and to larger modules (by combining modules with eigengene correlation >80%).

PDF plots showing correlation between modules as well as module colors in the context of the dendrogram of genes will be produced in the specified out_dir.

_Important output to keep:_ Gene module assignments and ExpressionSet object (from thresholding_parameters.R). These will be used to determine module eigengenes to correlate with phenotypes.

# Correlating Module Eigengenes with Phenotypes

```bash
Rscript eigengene_phenotype_correlation.R --eset=/path/to/expressionset/Rdata/file --gene_modules=/path/to/gene/modules/text/file --out_dir=/directory/to/save/correlation/matrices --discrete_phenotypes=csv_file --continuous_phenotyes=csv_file
```

This is the final data generation step of the workflow. After this, multiple sets (i.e. discovery and validation) can be analyzed together to look for modular overlaps. 

_Important output to keep:_ Correlation and p-value matrices. Together with gene module assignment text file, two different sets of _de novo_ module to phenotype correlations can be cross-referenced.

# Putting it all together

Running a discovery set from top to bottom could look like this:
```bash
Rscript thresholding_parameters.R --cel_path="../data/TransplantCELs/discovery/" --out="../results/discovery/"

# determined optimal k = 10

Rscript generate_modules.R --eset="../results/discovery/expression_set.RData" --soft_thresh_k=10 --out_dir="../results/discovery/"

Rscript eigengene_phenotype_correlation.R --eset="../results/discovery/expression_set.RData" --gene_modules="../results/discovery/gene_modules.txt" --out_dir="../results/discovery/" --discrete_phenotypes="../data/transplant_discrete_phenotypes.csv" --continuous_phenotyes="../data/transplant_continuous_phenotypes.csv"
```

We could then run the identical pipeline for a validation set. At the end, we will have gene module to phenotype correlations for both discovery and validation sets. When there is real biological signal, we should see overlap in gene membership for modules in the discovery and validation set which are both correlated with the phenotype of interest.

# Discovery and Validation
An example R Markdown (tool for reproducible code embedded in html document) is in the analysis directory and outlines one method for comparing discovery and validation results that have been produced following the "putting it all together" pipeline in the section above.

Essentially, any gene modules in the discovery and validation set that correlate with the same clinical outcome will be grouped together. A heatmap showing the overlap in gene membership for modules in the discovery and validation set will provide a visual interpretation of the corroboration of discovery signals by the validation set.

