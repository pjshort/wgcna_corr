# generate gene module assignments for a cohort of expression sets and soft
# thresholding parameter $k$

# the module assignment used here will be as fine-grained as possible
# many of the modules will be correlated with one another and plots of
# module dendrograms and inter-module correlation will be created

# gene modules will be saved as tab-delimited text files with columns
# gene_name    module_fine    module_coarse
# where module_coarse will be module_fine merged at 80% correlation level

# inputs:
#   CEL files
#   Soft Thresholding Parameter

# outputs:
#   Module Assignment text file

### dependencies
library(optparse)
library(oligo)
library(Biobase)
library(genefilter)
library(WGCNA)
library(stringr)
library(hugene10sttranscriptcluster.db)
library(genefilter)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("../R/eset_tools.R")
### command line options
option_list <- list(
  make_option("--cel_path", default="../data/TransplantCELs/discovery/",
              help="Top level directory to start recursive search for *.CEL files."),
  make_option("--eset", default=NULL, help = "Optionally pass an ExpressionSet object saved as .RData - 
              this will save time filtering for gene annotations. Object MUST have gene names as rows."),
  make_option("--soft_thresh_k", default=1, help = "Set soft thresholding parameter used to determine
              connections in the gene coexpression network. Appropriate parameter can be picked with the 
              help of thresholding_parameters.R"),
  make_option("--out_dir", default="../results/gene_modules/",
              help="Location to save the ExpressionSet object. Defaults to ./eset.Rdata"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

#test conditions
args$soft_thresh_k = 10
args$eset = "../results/discovery//expression_set.RData"
args$out_dir = "../results/discovery/"


# TODO - edit script to take ExpressionSet object instead of CEL files if user specifies
if (is.null(args$eset)){
  write("No ExpressionSet RData found - calculating expression set from CEL files.", stderr())
  cel_files = list.files(path = args$cel_path, pattern = ".*CEL", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  eset = rma_eset(CEL_list = cel_files)
  
  if (!file.exists("../data/HuGene-1_0-st-v1.na35.hg19.transcript.csv")){
write("Did not find HuGene probe annotation file - download  from 
http://www.affymetrix.com/support/technical/annotationfilesmain.affx and 
put in the data subdirectory", stderr())}
  affy_annotation = "../data/HuGene-1_0-st-v1.na35.hg19.transcript.csv"
  
  # filter for genes/probes with entrez ID - this will have the same set of genes as eset23
  eset <- gene_filter_eset(eset, affy_annotation)
} else {
  load(args$eset)  # should load an ExpressionSet object
}

expr_data <- t(exprs(eset)) # flip genes to columns and samples to rows

# change sample names to three digit Tritan IDs (in the 600s)
rownames(expr_data) <- str_extract(rownames(expr_data), "[6][0-9][0-9]") # matches 601 - 699
colnames(eset) <- rownames(expr_data)

#n_genes = ncol(expr_data)
#n_samples = nrow(expr_data)

if (args$soft_thresh_k == 1){
  warning("Soft threshold not set (or intentionally set to k = 1) - the data will likely not represent
           a scale free topology under this setting. See the ReadMe and run thresholding_parameters.R
           to determine a reasonable k.")
}

soft_thresh_power = args$soft_thresh_k

write("Calculating dissimilarity from Topological Overlap Matrix.", stderr())
dissTOM = 1 - TOMsimilarityFromExpr(expr_data, power = soft_thresh_power)

# create dendrogram with hclust
write("Creating dendrogram from gene-gene TOM.", stderr())
geneTree = hclust(as.dist(dissTOM), method = "average")

# module identification using dynamic tree cut
write("Running dynamic tree cut with deepSplit parameter 2 and minClusterSize of 30 - these are
       the default parameters.")
modules = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
module_colors = labels2colors(modules)

# merge modules at 80% correlation threshold
MEDissThres = 0.20
merge = mergeCloseModules(expr_data, module_colors, cutHeight = MEDissThres, verbose = 3)
merged_colors = merge$colors

gene_modules = cbind("gene_name" = rownames(eset), module_colors, merged_colors)

# create directory if it does not exist
write(sprintf("Writing output into %s", args$out_dir), stderr())
dir.create(args$out_dir)

# save gene modules to out dir
write.table(gene_modules, file = paste0(args$out_dir, "/gene_modules.txt"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# calculate module eigengenes to draw similarity dendrograms
fine_MEs = orderMEs(moduleEigengenes(expr_data, module_colors)$eigengenes)
colnames(fine_MEs) <- substring(colnames(fine_MEs), 3)
merged_MEs = orderMEs(moduleEigengenes(expr_data, merged_colors)$eigengenes)
colnames(merged_MEs) <- substring(colnames(merged_MEs), 3)


# save plots to out dir - dendrogram (fine-grained and merged) and module dendrogram
pdf(paste0(args$out_dir, "/gene_module_correlation.pdf"))
MEDiss = 1-cor(fine_MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
plot(METree, main = "Clustering Module Eigengenes",
     xlab = "", sub = "")

# cut the tree at 0.20 (80% correlation)
MEDissThres = 0.20
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

pdf(paste0(args$out_dir, "/gene_expr_dendrogram.pdf"))
plotDendroAndColors(geneTree, cbind(module_colors, merged_colors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

