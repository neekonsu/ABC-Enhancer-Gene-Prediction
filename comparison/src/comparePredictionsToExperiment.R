#####################
# Joseph Nasser
# June 4 2020
#
# Evaluate predictive model against experimental data
#
# Requirements: 
#       - R 3.4.0
#         - data.table (1.11.4)
#         - GenomicRanges (1.28.6)
#         - ROCR (1.0-7)
#         - ggplot2 (3.0.0)
#         - caTools (1.17.1)
#
# Usage:
#       Rscript comparePredictionsToExperiment.R --predictions pred.table.txt --experimentalData expt.txt --plotConfig plot.config.txt --predConfig pred.config.txt
#
# This code is intended to evaluate a various enhancer-gene prediction models against a set of experimental data
# The code will overlap the predictions and experimental data and generate appropriate evaluation metrics (eg pr curves)
#
# The plotConfig file should contain one line per PR curve plot. Scatter plots will be generated for all prediction columns defined in predConfig
#
# Configuring Prediction Columns:
# In the ideal use case each experimental element should correspond to exactly one element-gene prediction. 
# However, there are instances where a single experimental element overlaps multiple predictions (eg a large deletion) or when an experimentally tested element-gene pair is not present in the predictions files.
# The predConfig file describes how to handle these cases. 
#
# Other:
# - Can't distinguish between a missing prediction and a non-prediction
#
# To do:
# 1. Join on GeneSymbol or GeneTSS?
# 2. Validate preds 
#    - if min(pred$score) < pred.config$score$fill.val (should fail)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(caTools))

option.list <- list(
  make_option("--predictions", type="character", help="Predictions file (accepts comma delimited list)"),
  make_option("--experimentalData", type="character", help="File listing perturbational data (accepts comma delimited list)"),
  make_option("--experimentalPositiveColumn", type="character", default="Regulated", help="Column of experimentalData to consider an experimental positive"),
  make_option("--plotConfig", type="character", help="File describing which plots to make"),
  make_option("--predConfig", type="character", help="File describing how to aggregate/fill prediction columns"),
  make_option("--ignoreExptMissingPredictions", default=FALSE, action="store_true", help="Ignore EG pairs which do not have predictions. Do not fill based on predConfig"),
  make_option("--outDir", default = ".", type="character", help="Output directory"),
  make_option("--code", default = ".", type="character", help="comparison code helper functions")
)
opt <- parse_args(OptionParser(option_list=option.list))

#For Testing at Broad
# basedir <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/example/"
# opt$code <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/src/comparison.R"
# opt$predictions <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/example/input/pred.table.txt"
# opt$experimentalData <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/example/input/K562.ExperimentalData.slim.txt"
# opt$plotConfig <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/src/plot.config.txt"
# opt$predConfig <- "/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_better_comparison/ABC-Enhancer-Gene-Prediction/comparison/src/pred.config.txt"
# opt$outDir <- paste0(basedir, "out/")
# opt$ignoreExptMissingPredictions <- FALSE

source(opt$code)
dir.create(file.path(opt$outDir))
write.table(t(as.data.frame(opt)), file.path(opt$outDir, "params.txt"), sep = "\t", quote = F, col.names = F, row.names = T)

#Read input data
print("Reading input files")
pred.table <- fread(opt$predictions)
pred.list <- loadPredictions(pred.table)
expt <- loadFileString(opt$experimentalData)
expt <- subset(expt, IncludeInModel)
predConfig <- fread(opt$predConfig)
plotConfig <- fread(opt$plotConfig)

#QC Input Files
qcExpt(expt, opt)

#QC Predictions

# Merge experimental data and predictions
print("Merging experiment and predictions")
merged <- combineAllExptPred(expt = expt, 
                            pred.list = pred.list,
                            config = predConfig,
                            outdir = opt$outDir,
                            fill.missing = !opt$ignoreExptMissingPredictions)

inverse.predictors <- getInversePredictors(pred.list, predConfig)

plotCellType <- function(cellType) {

  if (cellType != 'combined') {
    this.merged <- subset(merged, CellType == cellType)
  } else {
    this.merged <- merged
  }

  # Make plots
  tryCatch({
      print(paste0("Making plots for celltype ", cellType))
      this.merged <- prepForPlotting(this.merged)
      #this.merged <- subset(this.merged, !(class %in% c("tss", "promoter")))
      dir.create(file.path(opt$outDir, cellType))
      makePlots(this.merged, plotConfig, inverse.predictors, opt$experimentalPositiveColumn, file.path(opt$outDir, cellType))
      }, error = function(e) {
        print(e)
      })
}

lapply(c('combined', unique(merged$CellType)), plotCellType)
