# generate graphs for estimating soft thresholding parameters

# inputs:
#   CEL files

### dependencies
library(optparse)
library(stringr)
library(oligo)
library(Biobase)
library(genefilter)
library(WGCNA)
library(flashClust)
library(hugene10sttranscriptcluster.db)
library(genefilter)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("../R/eset_tools.R")
### command line options
option_list <- list(
  make_option("--cel_path", default="../data/TransplantCELs/batch1_2009",
              help="Top level directory to start recursive search for *.CEL files."),
  make_option("--eset", default = NULL, help = "ExpressionSet object saved as an RData file (must be saved as variable named eset)."),
  make_option("--out_dir", default="../results",
              help="Location to save the ExpressionSet object. Defaults to ./eset.Rdata"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

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

colnames(eset) <- str_extract(colnames(eset), "[6][0-9][0-9]")
expr_data = t(exprs(eset))

### generate network using different soft thresholding parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)

write(sprintf("Writing output into %s", args$out_dir), stderr())
dir.create(args$out_dir)

pdf(paste0(args$out_dir, "/soft_thresholding_correlation.pdf"))
# from tutorial
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

write("The plot produced can be used to pick a soft thresholding parameter 
      for subsequent analyses. Red line is drawn at 0.90, but k should be
      chosen at the point where correlation begins to level out and where
      mean connectivity is at a reasonable value.", stderr())

write("Saving expression set in out dir.", stderr())


save(eset, file = paste0(args$out_dir, "/expression_set.RData"))
