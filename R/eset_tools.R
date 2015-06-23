# construct normalized expression sets from list of CEL files 
# uses RMA from the oligo package

library(oligo)
library(hugene10sttranscriptcluster.db)
library(genefilter)

rma_eset <- function(CEL_list){
  
  # take a list of .CEL files and make a normalized ExpressionSet object
  cels <- read.celfiles(CEL_list)
  
  # build an expression set
  eset <- oligo::rma(cels, target="core")
  
  return(eset)
}

# download this file from http://www.affymetrix.com/support/technical/annotationfilesmain.affx
#affy_annotation = "../data/HuGene-1_0-st-v1.na35.hg19.transcript.csv"

gene_filter_eset <- function(eset, affy_annotation){
  
  # take an expression set and affy annotation file and filter:
  # only genes with Entrez IDs
  # one probe per gene (largest IQR)
  
  write("Annotating expression set with gene names and Entrez IDs.", stderr())
  
  # annotate the expressionSet object
  fData <- data.frame(ID = featureNames(eset))
  rownames(fData) <- fData$ID
  
  # official Gene Name
  ann <- mget(as.character(fData$ID), hugene10sttranscriptclusterSYMBOL, ifnotfound=NA)
  ann <- lapply(ann, paste, collapse=";")
  fData$Name <- unlist(ann)
  
  # Entrez ID
  ann <- mget(as.character(rownames(fData)),hugene10sttranscriptclusterENTREZID, ifnotfound=NA)
  ann <- lapply(ann, paste, collapse=";")
  fData$EntrezID <- unlist(ann)
  
  # Chromosome
  ann <- mget(as.character(rownames(fData)), hugene10sttranscriptclusterCHR, ifnotfound=NA)
  ann <- lapply(ann, paste, collapse=";")
  fData$Chromosome <- unlist(ann)
  
  # Long Gene Name
  ann <- mget(as.character(rownames(fData)), hugene10sttranscriptclusterGENENAME, ifnotfound=NA)
  ann <- lapply(ann, paste, collapse=";")
  fData$LongName <- unlist(ann)
  
  # Affymetrix probe status annotation
  affy.annot <- read.csv(affy_annotation, skip=21, header=T)
  
  write("Annotating probes with affy probe information.", stderr())
  rownames(affy.annot) <- affy.annot$transcript_cluster_id
  fData$ProbeStatus <- factor(as.character(affy.annot[rownames(fData),"category"]))
  metadata <- data.frame(labelDescription = c("Manufacturers ID", "Official Symbol", 
                                              "EntrezID", "Chromosome", "Gene Name", 
                                              "Affy Probe Status"), row.names=c("ID", 
                                                                                "Name", "EntrezID", "Chromosome", "LongName", 
                                                                                "ProbeStatus"))
  
  colnames(fData) <- c("ID", "Name", "EntrezID", "Chromosome", "LongName", "ProbeStatus")
  
  features <- new("AnnotatedDataFrame", data = fData, varMetadata = metadata)
  
  featureData(eset) <- features
  
  
  write("Reducing eset to one probe per gene, only genes with annotated Entrez IDs.", stderr())
  
  # one gene per probe, only probes that map to an entrezID
  entrezIds <- mget(featureNames(eset), envir = hugene10sttranscriptclusterENTREZID, ifnotfound=NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  numNoEntrezId <- length(featureNames(eset)) - length(haveEntrezId) 
  
  eset <- eset[haveEntrezId, ]
  
  # make sure each probe only maps to 1 entrezID
  esIqr <- apply(exprs(eset), 1, IQR)
  uniqGenes <- findLargest(featureNames(eset), esIqr, "hugene10sttranscriptcluster")
  eset <- eset[uniqGenes, ]
  numSelected <- length(featureNames(eset))
  
  # Now make gene symbol the featureName - makes more sense to work with
  eset <- eset[fData(eset)$ProbeStatus == "main", ]
  
  # remove any duplicate gene names that have made it to here
  eset <- eset[!duplicated(fData(eset)$Name),]
  
  featureNames(eset) <- fData(eset)$Name
  
  
  return(eset)
  
}

match_eset_probes <- function(discovery_eset, validation_eset) {
  
  # slices validation expression set to include only probes used in discovery set
  
  to_keep = match(fData(discovery_eset)$ID, featureNames(validation_eset))
  validation_eset <- validation_eset[to_keep,]
  featureNames(validation_eset) <- fData(discovery_eset)$Name
  
  metadata <- data.frame(labelDescription = c("Manufacturers ID", "Official Symbol", 
                                              "EntrezID", "Chromosome", "Gene Name", 
                                              "Affy Probe Status"), row.names=c("ID", "Name", "EntrezID", "Chromosome", "LongName","ProbeStatus"))
  
  features <- new("AnnotatedDataFrame", data = fData(discovery_eset), varMetadata = metadata)
  featureData(validation_eset) <- features
  
  if (dim(discovery_eset)[1] != dim(validation_eset)[1]){
    stop("Genes/probes in discovery set not equal to validation set. Check eset_tools::match_eset_probes for troubleshooting.")
  }
  
  return(validation_eset)
  
}