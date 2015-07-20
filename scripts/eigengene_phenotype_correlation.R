# correlate gene module eigengene expression (first principle component) with phenotypes of interest

# generates a heatmap as well as RData matrix of module to phenotype correlation and p-values 

# inputs:
#   ExpressionSet (RData)
#   Gene Module Assignments

# outputs:
#   correlation matrix
#   correlation and heatmap

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
source("../R/correlation_tools.R")

### command line options
option_list <- list(
  make_option("--cel_path", default="../data/TransplantCELs/discovery/",
              help="Top level directory to start recursive search for *.CEL files."),
  make_option("--eset", default=NULL, help = "Optionally pass an ExpressionSet object saved as .RData - 
              this will save time filtering for gene annotations. Object MUST have gene names as rows."),
  make_option("--gene_modules", default="../data/gene_modules.txt", help = "Text file assigning gene names to module colors."),
  make_option("--discrete_phenotypes", default="../data/transplant_discrete_phenotypes.csv", help = "Discrete phenotypes for correlation (uses spearman rank correlation)."),
  make_option("--continuous_phenotypes", default="../data/transplant_continuous_phenotypes.csv", help = "Continuous phenotypes for correlation (uses pearson correlation)."),
  make_option("--use_merged", default = FALSE, help = "Use merged modules rather than fine-grained modules? Defaults to FALSE."),
  make_option("--out_dir", default="../results/gene_modules/",
              help="Location to save the ExpressionSet object. Defaults to ./eset.Rdata"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

args$eset = "../results/discovery_full/expression_set.RData"
args$gene_modules = "../results/discovery_full/gene_modules.txt"
args$out_dir = "../results/discovery_full/"

# load gene modules from tab delimited text file
gene_modules = read.table(args$gene_modules, header=TRUE, sep = "\t")

# load ExpressionSet
load(args$eset) # expression set should be named eset
expr_data = t(exprs(eset))
n_samples = ncol(eset)

# get gene names from rows of expression set
genes = rownames(eset)

if (args$use_merged == TRUE){
  module_colors = gene_modules[match(genes, gene_modules$gene_name), "merged_colors"]
} else {
  module_colors = gene_modules[match(genes, gene_modules$gene_name), "module_colors"]
}

# correlate and produce heat figs
continuous_pheno = read.csv(args$continuous_phenotypes)
discrete_pheno = read.csv(args$discrete_phenotypes)

# get only relevant patient IDs - should be columns of eset
patient_ids = colnames(eset)
continuous_pheno = continuous_pheno[match(patient_ids, continuous_pheno$patient_number),]
discrete_pheno = discrete_pheno[match(patient_ids, discrete_pheno$patient_number),]

# remove any phenotypes with >75% missing data
continuous_pheno <- continuous_pheno[ , colSums(is.na(continuous_pheno)) < n_samples*0.75]
discrete_pheno <- discrete_pheno[ ,colSums(is.na(discrete_pheno)) < n_samples*0.75]

# false discovery rate for multiple hypothesis testing
fdr = 0.25

if (length(table(module_colors > 30))){
  write("Too many gene modules to write names out heatmap - toggle --use_merged=TRUE to merged closely
         correlated modules.")
  colors_only = TRUE
} else {
  colors_only = FALSE
}

d = pheno_heat_fig(expr_data, module_colors, discrete_pheno, 
                               fname = paste0(args$out_dir, "/discrete_corr_heatmap.pdf"), 
                               signif_only = FALSE, fdr = fdr, colors_only = colors_only, method = "spearman")
discrete_corr = d[["corr"]]
discrete_p = d[["p"]]
rownames(discrete_corr) = substring(rownames(discrete_corr), 3)
rownames(discrete_p) = substring(rownames(discrete_p), 3)

c = pheno_heat_fig(expr_data, module_colors, continuous_pheno, 
                                fname = paste0(args$out_dir, "/continuous_corr_heatmap.pdf"), 
                                signif_only = TRUE, fdr = fdr, colors_only = colors_only, method = "pearson")
continuous_corr = c[["corr"]]
continuous_p = c[["p"]]
rownames(continuous_corr) = substring(rownames(continuous_corr), 3)
rownames(continuous_p) = substring(rownames(continuous_p), 3)

write("Saving correlation and p-value matrices (modules x phenotypes) to RData file.", stderr())
save(continuous_corr, continuous_p, discrete_corr, discrete_p, file = paste0(args$out_dir, "/module_to_phenotype_correlations.RData"))


