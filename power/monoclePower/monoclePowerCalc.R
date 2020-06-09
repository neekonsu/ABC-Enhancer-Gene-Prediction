
library(monocle)
library(data.table)
library(GenomicRanges)

runSingleMonoclePowerCalc <- function(effectSize, cds, trt.col, empirical.pval.ref, TRIALS, use.disp.fit, adjust.disp.by.effect.size) {
  # Compute the power to detect a decrease in expression from a real monocle dataset
  # Modeled based on the DESeq power functions in src/CRISPRScreen/PowerCalculations.R
  # trt.col: must be a boolean column name from featureData representing the cells for which a decrease will be simulated
  
  print(paste0("Running Power Calc at Effect Size = ", effectSize, " for ", nrow(cds), " genes"))
  
  #Run simulation
  sim.results.list <- replicate(TRIALS, runSingleMonocleTest(cds = cds, trt.col = trt.col, effectSize = effectSize, use.disp.fit = use.disp.fit, adjust.disp.by.effect.size = adjust.disp.by.effect.size), simplify = FALSE)
  if (!is.null(empirical.pval.ref)) sim.results.list <- lapply(sim.results.list, computeSignificanceUsingReferencePvalues, empirical.pval.ref)
  sim.results <- as.data.table(rbindlist(sim.results.list))
  
  #Compute power
  if (!is.null(empirical.pval.ref)) {
    power.calc <- sim.results[,list(
                                    this.pvals = paste(pval, collapse = ","),
                                    this.empirical.pvals = paste(empirical.pvalue, collapse = ",")), by = gene_short_name]
    
    # Remove power calc, b/c need to do a global fdr
    # this.power.empirical = sum(empirical.qvalue < qval.cutoff)/TRIALS,
    # this.power.non.empirical = sum(qval < qval.cutoff)/TRIALS,
  } else {
    power.calc <- sim.results[,list(this.pvals = paste(pval, collapse = ",")), by = gene_short_name]
  }
  
  print(power.calc)
  colnames(power.calc)[colnames(power.calc) == "this.pvals"] <- paste0("NonEmpiricalPvalsAtEffectSize=", effectSize)

  return(power.calc)
}

runSingleMonocleTest <- function(cds, trt.col, effectSize, use.disp.fit, adjust.disp.by.effect.size) {
  #Do one permutation: Simulate dataset and run differential expression test

  a = Sys.time()
  sim.cds <- simulateMonocleDataSetFromRealData(cds = cds, trt.col = trt.col, effectSize = effectSize, use.disp.fit = use.disp.fit, adjust.disp.by.effect.size = adjust.disp.by.effect.size)
  diff.test <- differentialGeneTest(sim.cds, fullModelFormulaStr = paste0("~", trt.col), verbose = TRUE)
  b = Sys.time()    
  
  print(paste0("This permutation took: ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 2), " seconds"))
 
  return(diff.test)
}


simulateMonocleDataSetFromRealData <- function(cds, trt.col, effectSize, use.disp.fit, adjust.disp.by.effect.size) {
  #Given a monocle data set, simulate a new monocle data set on the same number of genes and cells using a neg.binom distribution
  #Use the dispersions and Size factors from the real dataset
  #The trt.col column indicates which cells should be simulated to have a decrease in expression
  #
  #Loosely based on src/CRISPRScreen/PowerCalculations.R ~ simulateDESeqDataSet
  
  nFeatures <- nrow(cds)
  nSamples <- ncol(cds)
  
  baseMean <- cds@dispFitInfo$blind$disp_table$mu
  sizeFactors <- cds@phenoData@data$Size_Factor
  
  # Get indices of treated cells
  treated.cells <- which(cds@phenoData@data[, trt.col]) 
  
  # Make mu matrix for simulation
  mu <- matrix( rep(baseMean, nSamples), ncol=nSamples)
  mu <- sweep(mu, 2, sizeFactors, "*")
  mu[, treated.cells] <- mu[, treated.cells] * (1 - effectSize)
  
  # Dispersion
  # Three options:
  # 1. Use observed dispersions
  # 2. Use dispersions from fit
  # 3. Use dispersions from fit AND adjust dispersion based on effect size in simulation
  if (use.disp.fit) {
    if (adjust.disp.by.effect.size) {
      dispersion <- cds@dispFitInfo$blind$disp_func(mu)
    } else {
      dispersion <- cds@dispFitInfo$blind$disp_func(baseMean) 
    }
  } else {
    dispersion <- cds@dispFitInfo$blind$disp_table$disp
  }
  
  # Simulate counts
  count.data <- matrix(rnbinom(nSamples*nFeatures, mu=mu, size=1/dispersion), ncol=nSamples)
  cds.out <- newCellDataSet(as(count.data, "sparseMatrix"),
                            phenoData = cds@phenoData,
                            featureData = cds@featureData,
                            expressionFamily = negbinomial.size())
  cds.out@dispFitInfo <- cds@dispFitInfo
  cds.out@phenoData@data$Size_Factor <- cds@phenoData@data$Size_Factor
  
  return(cds.out)
  
}

subsetMonocleCDS <- function(cds, gene.list, pheno.col.list=NULL) {
  if (!is.null(pheno.col.list)) {
    cds@phenoData@data <- cds@phenoData@data[, pheno.col.list]
  }
  
  #Make sure we only get genes that have dispersion info. This is necessary for the gasperini highmoi pilot
  gene.list <- intersect(gene.list, cds@dispFitInfo$blind$disp_table$gene_id)
  
  cds <- cds[gene.list, ]
  cds@dispFitInfo$blind$disp_table <- subset(cds@dispFitInfo$blind$disp_table, gene_id %in% gene.list)
  cds@dispFitInfo$blind$disp_table <- cds@dispFitInfo$blind$disp_table[match(gene.ids, cds@dispFitInfo$blind$disp_table$gene_id),] #match ordering
  
  #stopifnot(all(sort(cds@dispFitInfo$blind$disp_table$gene_id) == sort(cds@featureData@data$id)))
  stopifnot(all(cds@dispFitInfo$blind$disp_table$gene_id == cds@featureData@data$id))
  stopifnot(nrow(cds@dispFitInfo$blind$disp_table) == length(gene.list))

  return(cds)
}

computeSignificanceUsingReferencePvalues <- function(df, pval.ref) {
  df$empirical.pvalue <- unlist(lapply(df$pval, getEmpiricalPval, pval.ref))
  df$empirical.qvalue <- p.adjust(df$empirical.pvalue, method = "BH")
  return(df)
}

getEmpiricalPval <- function(pval, pval.ref) {
  e.pval <- (sum(pval > pval.ref) + 1) / (length(pval.ref) + 1)
  return(e.pval)
}
