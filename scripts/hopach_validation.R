# HOPACH gene clustering for independent validation

# this script will operate very simply, taking in an expression set and returning
# a .RData file with the hopach object with assignments from genes (rows of expression set)
# to gene modules

# inputs:
#   ExpressionSet (RData)
#   Gene Module Assignments

# outputs:
#   correlation matrix
#   correlation and heatmap

### dependencies
library(hopach)
library(Biobase)
library(optparse)
options(stringsAsFactors = FALSE)


### command line options
option_list <- list(
  make_option("--eset", default=NULL, help = "Pass an ExpressionSet object saved as .RData - 
              this will save time filtering for gene annotations. Object MUST have gene names as rows."),
  make_option("--subset", default=0, help = "If true, should be integer defining the number of genes
               to take for analysis (these will be genes with highest variance). Otherwise, defaults to 0 (take all genes and do not subset)."),
  make_option("--greedy", action="store_true", default=FALSE, help = "If true, should be integer defining the number of genes
               to take for analysis (these will be genes with highest variance)."),
  make_option("--out_dir", default="../results/gene_modules/",
              help="Location to save the ExpressionSet object and hopach module assignments.
              Defaults to ./eset.Rdata"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

load(args$eset)  # loads a genes x patients Eset object

if (args$subset == 0){
  exp.subset = exprs(eset)
} else {
  vars = apply(exprs(eset), 1, var)
  subset = vars > quantile(vars, (nrow(exprs(eset))-args$subset)/nrow(exprs(eset)))
  exp.subset = exprs(eset)[subset,]
}

hopach_gene_names = rownames(exp.subset)
gene.dist = distancematrix(exp.subset, "cosangle")

if (args$greedy){
  gene.hobj = hopach(exp.subset, dmat=gene.dist, clusters="greedy")
} else {
  gene.hobj = hopach(exp.subset, dmat=gene.dist)
}

save(exp.subset, gene.hobj, file = paste0(args$out_dir, "hopach_modules.RData"))
