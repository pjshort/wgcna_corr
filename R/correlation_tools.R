# tools for assessing and plotting correlations between gene modules and phenotypes

pheno_heat_fig <- function(expr_data, module_colors, phenotypes, 
                           fname = "pheno_heat_fig.pdf", method = "pearson", 
                           signif_only = TRUE, fdr = 0.05, colors_only = FALSE){
  
  MEs0 = moduleEigengenes(expr_data, module_colors)$eigengenes
  MEs = orderMEs(MEs0)
  n_samples = dim(expr_data)[1]  # rows in expr data
  
  # MEs will be n_samples x n_modules. values in particular column will be the eigengene expression (i.e. linear combination of module member genes that explains the most variance) for each sample 
  
  # using the pearson correlation to correlate eigengene expression with continuous variables
  corr = WGCNA::cor(MEs, data.matrix(phenotypes), use = "p", method = method)
  
  # continuous trait p values - adjust for number of tests
  p_raw = corPvalueStudent(corr, n_samples)
  
  # benjamini-hochberg error rate correction
  p = apply(p_raw, MARGIN = 2, FUN = function(p) p.adjust(p, length(corr), method = "fdr"))
  
  phenos_to_display = seq(1, dim(phenotypes)[2])
  
  if (signif_only == TRUE){
    lowest_p = apply(p, MARGIN = 2, FUN = min)
    phenos_to_display = lowest_p < fdr
    corr = corr[ , phenos_to_display]
    p = p[ , phenos_to_display]
  }
  
  textMatrix = signif(corr, 2)
  
  # mark correlations with a * that fall below FDR 0.05
  textMatrix[which(p < fdr)] = paste0(textMatrix[which(p < fdr)], "*")
  
  dim(textMatrix) = dim(corr)
  
  xlabs = names(phenotypes)[phenos_to_display]
  xlabs[colSums(p < fdr) > 0] = paste0(xlabs[colSums(p < fdr) > 0], "*")
  
  ylabs = substring(names(MEs), first = 3)
  ylabs[rowSums(p < fdr) > 0] = paste0(ylabs[rowSums(p < fdr) > 0], "*")
  
  # Display the correlation values within a heatmap plot
  pdf(fname)
  par(mar = c(8, 7.5, 3, 3));
  if (colors_only){
    labeledHeatmap(Matrix = corr, xLabels = xlabs, yLabels = names(MEs),
                   yColorLabels = TRUE, colors = blueWhiteRed(50), 
                   textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, 
                   zlim = c(-1,1), main = paste("Gene-Module to Phenotype Relationships"), 
                   cex.lab.x = 0.75)
  } else { labeledHeatmap(Matrix = corr, xLabels = xlabs, yLabels = names(MEs), 
                          ySymbols = ylabs, yColorLabels = TRUE, 
                          colors = blueWhiteRed(50), textMatrix = textMatrix, 
                          setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), 
                          main = paste("Gene-Module to Phenotype Relationships"), 
                          cex.lab.x = 0.75)
  }
  dev.off()
  
  # note that labeledHeatmap does not render within knitr - need to save and load as a chunk
  return(list("corr" = corr, "p" = p_raw))
  
}

significant_corrs <- function(MEs, phenotypes, method = "pearson", fdr = 0.05){
  
  # takes module eigengene expression levels and phenotype data matrix
  # returns dataframe of module to phenotype correlations with correlation strength and p value
  
  corr = WGCNA::cor(MEs, phenotypes, use = "p", method = method)
  
  # get p value from correlation and degrees of freedom
  p = corPvalueStudent(corr, n_samples)
  
  # benjamini-hochberg error rate correction
  p_adjust = apply(p, MARGIN = 2, FUN = function(p) p.adjust(p, length(corr), method = "fdr"))
  
  # get row and column with adjusted p < false discovery rate
  l_mat = which(p_adjust < fdr, arr.ind = TRUE)
  l_flat = which(p_adjust < fdr)
  modules = names(MEs)[l_mat[,1]]
  traits = colnames(phenotypes)[l_mat[,2]]
  cs = corr[l_flat]
  ps = p[l_flat]
  
  associations <- data.frame(cbind(modules, traits, cs, ps))
  return(associations)
  
}