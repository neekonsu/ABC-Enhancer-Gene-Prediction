###########
# Run Monocle given a CDS in RDS form. 
#
# use .r-3.4.1 

library(data.table)
library(optparse)
library(monocle)
source("/seq/lincRNA/jnasser/EP_prediction/code/LanderLab-EP-Prediction/src/CRISPRScreen/monoclePowerCalc.R")

option.list <- list(
  make_option("--cds", type="character", help="Monocle Dataset"),
  make_option("--geneList", type="character", help="List of genes to test"),
  make_option("--outDir", type="character", help="output directory"),
  make_option("--fullModelStr", default="", type="character"),
  make_option("--reducedModelStr", default="", type="character"),
  make_option("--region.id", default="", type="character", "string specifying the targeted enhancer. This is for reporting purposes only!"),
  make_option("--logicalColsToTest", type="character", default="", help="Comma delimited list of columns in cds. If provided will be ORed and added to fullModelStr"),
  make_option("--grepCol", type="character", default="", help="Make a boolean column by grepping through this column. Will be added to fullModelStr"),
  make_option("--grep", type="character", default="", help="Pattern to search for in grepCol"),
  make_option("--referencePvalues", default = "", help="File containing a single column of reference pvalues to use for significance testing"),
  make_option("--convertGeneNames", default = FALSE, action="store_true", help="Use featureData of cds to convert gene names to ids"),
  make_option("--cores", default=1, help="Number of cores for monocle differential expression test"),
  make_option("--coresForPowerCalc", default=1, help="Number of cores for parallelizing over different power effect sizes"),
  make_option("--nPermsForPowerCalc", default = 0, help="Set to 0 to skip power calc"),
  make_option("--powerEffectSizes", default = ".1,.25,.5", help="Effect sizes for which to run power calc"),
  make_option("--useDispFit", default = FALSE, action="store_true", help="In power calculation, use dispersion from fit. Otherwise use observed dispersion"),
  make_option("--adjustDispByEffectSize", default = FALSE, action="store_true", help="In power calculation, adjust dispersion estimate by modeled effect size. Only applicable if --useDispFit is set")
)
opt <- parse_args(OptionParser(option_list=option.list))

print(opt)

cds <- readRDS(opt$cds)

#get the gene names to run
if (opt$convertGeneNames) {
  genes <- fread(opt$geneList, header = FALSE)$V1
  gene.ids <- subset(cds@featureData@data, gene_short_name %in% genes)$id
} else {
  gene.ids <- fread(opt$geneList, header = FALSE)$V1
}

#Make the treatment column
if (opt$logicalColsToTest != ""){
  or.vec <- strsplit(opt$logicalColsToTest, ",")[[1]]
  
  if (length(or.vec) == 1) {
    cds@phenoData@data[, "temp.col"] <- cds@phenoData@data[, or.vec]
  } else {
    cds@phenoData@data[, "temp.col"] <- apply(cds@phenoData@data[, or.vec], 1, any)
  }
  
  full.formula <- ifelse(opt$fullModelStr == "", "temp.col", paste0("temp.col", "+", opt$fullModelStr))
  num.trt <- sum(cds@phenoData@data[, "temp.col"])
  num.ctrl <- sum(!cds@phenoData@data[, "temp.col"])
} else if (opt$grepCol != "") {
  cds@phenoData@data[, "temp.col"] <- grepl(opt$grep, cds@phenoData@data[, opt$grepCol])
  
  full.formula <- ifelse(opt$fullModelStr == "", "temp.col", paste0("temp.col", "+", opt$fullModelStr))
  num.trt <- sum(cds@phenoData@data[, "temp.col"])
  num.ctrl <- sum(!cds@phenoData@data[, "temp.col"])
} else {
  full.formula <- opt$fullModelStr
  num.trt <- NaN
  num.ctrl <- NaN
}
print(full.formula)


#Subset CDS
#This is to keep memory profile low
if (TRUE) {
  print("Subsetting CDS")
  if (opt$logicalColsToTest != "") {
    pheno.cols <- setdiff(c("temp.col", strsplit(opt$reducedModelStr, "[+]")[[1]]), "1")
  } else {
    #May not be very stable...
    pheno.cols <- union(strsplit(opt$reducedModelStr, "[+]")[[1]], strsplit(opt$fullModelStr, "[+]")[[1]])
    pheno.cols <- c("temp.col", pheno.cols)
  }
  this.cds <- subsetMonocleCDS(cds, gene.ids, c("Size_Factor", pheno.cols))
  rm(cds)
  gc()
} else {
  this.cds <- cds
}

#Need to overwrite the monocle function
#Do this using assignInNamespace...
compareModels <- function (full_models, reduced_models) {
  stopifnot(length(full_models) == length(reduced_models))
  
  # print("In compare models")
  # print(Sys.time())
  
  test_res <- mapply(function(x, y) {
    if (is.null(x) == FALSE && is.null(y) == FALSE) {
      lrt <- VGAM::lrtest(x, y)
      pval = lrt@Body["Pr(>Chisq)"][2, ]
      family = x@family@vfamily
      if (length(family) > 1) 
        family = family[1]
      data.frame(status = "OK", family = family, pval = pval, 
                 intercept = x@coefficients[[1]], model.coeff = x@coefficients[[2]])
    }
    else {
      data.frame(status = "FAIL", family = NA, pval = 1)
    }
  }, full_models, reduced_models, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method = "BH")
  test_res
}
assignInNamespace("compareModels", compareModels, "monocle")

#Run differential expression test
print("Running Test")
a = Sys.time()
diff.test <- monocle::differentialGeneTest(this.cds,
                                           fullModelFormulaStr = paste0("~", full.formula),
                                           reducedModelFormulaStr = paste0("~", opt$reducedModelStr),
                                           cores = opt$cores,
                                           verbose = TRUE)
b = Sys.time()
print(paste0("This enhancer took: ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 2), " seconds"))

#Compute significance based on reference pvalues
if (opt$referencePvalues != "") {
  print("Computing empirical pvalues based on reference!")
  pval.ref <- fread(opt$referencePvalues, header=F)$V1
  diff.test <- computeSignificanceUsingReferencePvalues(diff.test, pval.ref)
} else {
  pval.ref <- NULL
}

diff.test$full.model <- paste0("~", full.formula)
diff.test$reduced.model <- paste0("~", opt$reducedModelStr)
diff.test$region.id <- opt$region.id
diff.test$n.trt <- num.trt
diff.test$n.ctrl <- num.ctrl
write.table(diff.test, file.path(opt$outDir, "monocleDiffTestResults.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

#Run Power
if (opt$nPermsForPowerCalc > 0) {
  print("Starting power calc")
  effect.sizes <- as.numeric(strsplit(opt$powerEffectSizes, split = ",")[[1]])
  
  #Get the covariate we want to compute power for
  if (opt$logicalColsToTest != "" | opt$grepCol != "") {
    power.col <- "temp.col"
  } else {
    #May not be very stable...
    power.col <- setdiff(strsplit(full.formula, "[+]")[[1]], strsplit(opt$reducedModelStr, "[+]")[[1]])
    print(power.col)
  }
  
  if (opt$coresForPowerCalc == 1) {
    res.list <- lapply(effect.sizes, runSingleMonoclePowerCalc, 
                       cds = this.cds, 
                       trt.col = power.col, 
                       empirical.pval.ref = pval.ref, 
                       TRIALS = opt$nPermsForPowerCalc,
                       use.disp.fit = opt$useDispFit,
                       adjust.disp.by.effect.size = opt$adjustDispByEffectSize)
  } else {
    cl <- makeCluster(opt$coresForPowerCalc)
    clusterExport(cl=cl, varlist=c("this.cds", "power.col", "pval.ref", "opt"))
    res.list <- parLapply(cl, effect.sizes, runSingleMonoclePowerCalc, 
                       cds = this.cds, 
                       trt.col = power.col, 
                       empirical.pval.ref = pval.ref, 
                       TRIALS = opt$nPermsForPowerCalc,
                       use.disp.fit = opt$useDispFit,
                       adjust.disp.by.effect.size = opt$adjustDispByEffectSize)
    stopCluster(cl)
  }
  all.res <- Reduce(function(...) merge(..., all=T), res.list)
  
  diff.test.with.power <- merge(diff.test, all.res, by = "gene_short_name")
  write.table(diff.test.with.power, file.path(opt$outDir, "monocleDiffTestResultsWithPower.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
}
