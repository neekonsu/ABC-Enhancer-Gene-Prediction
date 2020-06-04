library(GenomicRanges)

qcExpt <- function(expt, opt) {
  #Check for duplicate experiments
  dupe <- any(duplicated(expt[, c("CellType","GeneSymbol","chrElement","startElement","endElement")] ))
  if (dupe) {
    print("Error: The experimental data file contains duplicate experiments!")
    stop()
  }
  
  #check to make sure regulated column contains TRUE/FALSE
  # reg.vals <- sort(unique(expt[, get(opt$experimentalPositiveColumn)]))
  # if (!(identical(reg.vals, c(FALSE, TRUE)) | identical(reg.vals, c(0, 1)))) {
  #   print("Error: The experimental data column must contain exactly two distinct values: TRUE and FALSE")
  #   stop()
  # }

  #check to make sure regulated column contains TRUE/FALSE
  reg.vals <- sort(unique(expt[, get(opt$experimentalPositiveColumn)]))
  if (!(all(reg.vals %in% c(FALSE, TRUE)) | all(reg.vals %in% c(0, 1)))) {
    print("Error: The experimental data column must contain TRUE/FALSE")
    stop()
  }
  if (length(reg.vals) == 1) {
    print("Note: all values are either positives or negatives. Plotting code will fail, but merged prediction/experimental table will be output.")
  }
}

combineAllExptPred <- function(expt, pred.list, config, outdir, fill.missing) {
  merged.list <- lapply(names(pred.list), function(s) combineSingleExptPred(expt = expt, pred = pred.list[[s]], pred.name = s, config = config, outdir = outdir, fill.missing = fill.missing))
  merged <- Reduce(function(x, y) merge(x, y, all=TRUE), merged.list)
  write.table(merged, file.path(outdir, "expt.pred.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  return(merged)
}

combineSingleExptPred <- function(expt, pred, pred.name, config, outdir, fill.missing=TRUE) {
  #Subset config to columns that actuall appear. Otherwise code will fail
  config <- subset(config, pred.col %in% colnames(pred))
  
  pred.gr <- with(pred, GRanges(paste0(CellType,":",chrElement,":",GeneSymbol), IRanges(startElement, endElement)))
  expt.gr <- with(expt, GRanges(paste0(CellType,":",chrElement,":",GeneSymbol), IRanges(startElement, endElement)))
  ovl <- GenomicRanges::findOverlaps(expt.gr, pred.gr)
  
  #Merge predictions with experimental data
  merged <- cbind(expt[queryHits(ovl)], pred[subjectHits(ovl), config$pred.col, with = F])
  
  #Sometimes a perturbed element will overlap multiple model elements (eg in the case of a large deletion)
  #In these cases need to summarize, Eg sum ABC.Score across model elements overlapping the deletion
  #This requires a config file describing how each prediction column should be aggregated
  agg.cols <- c("chrElement","startElement","endElement","GeneSymbol","CellType","Significant","Regulated","EffectSize") #"class",
  merged <- collapseEnhancersOverlappingMultiplePredictions(merged, config, agg.cols)
  
  #Experimental data missing predictions
  #A tested enhancer element may not have a prediction
  #For ABC this is typically the case if the tested element does not overlap a DHS peak.
  #In this case we need to fill the predictions table
  #However, this code doesn't work if the predictor did not make a prediction for this gene.
  expt.missing.predictions <- expt[setdiff(seq(nrow(expt)), queryHits(ovl)),]
  write.table(expt.missing.predictions, file.path(outdir, "expt.missing.predictions.txt"), sep="\t", quote=F, col.names=T, row.names=F)
  print("The following experimental data is not present in predictions file: ")
  print(expt.missing.predictions[, ..agg.cols])
  if (fill.missing) {
    expt.missing.predictions <- fillMissingPredictions(expt.missing.predictions, config, agg.cols)
    cols.we.want <- c(agg.cols, config$pred.col) #'class'
    merged <- rbind(merged, expt.missing.predictions[, ..cols.we.want])
    print("Experimental data missing predictions filled. Will be considered in PR curve!")
    print(expt.missing.predictions[, ..cols.we.want])
  } else {
    print("Experimental data missing predictions ignored. Will not be considered in PR curve!")
  }

  #Rename new columns based on prediction dataset name
  colnames(merged)[colnames(merged) %in% config$pred.col] <- paste0(pred.name, ".", colnames(merged)[colnames(merged) %in% config$pred.col])
  
  return(merged)
}

fillMissingPredictions <- function(df, config, agg.cols) {
  #Fill in missing predictions as described in the config file
  for (ii in seq(nrow(config))) {
    df[, config$pred.col[[ii]]] <- config$fill.val[ii]
  }
  
  unk.cols <- setdiff(c('class', agg.cols), unique(c(colnames(df), config$pred.cols)))
  df[, unk.cols] <- "Merge:UNKNOWN"
  return(df)
}

collapseEnhancersOverlappingMultiplePredictions <- function(df, config, agg.cols) {
  
  #Summarize columns as defined in config
  list.for.agg <- as.list(df[, ..agg.cols])
  all.list <- mapply(function(pred.col, agg.func) aggregate(df[, ..pred.col], by = list.for.agg, FUN = agg.func), 
                     config$pred.col, config$agg.func, SIMPLIFY=F)
  
  #Special handling for aggregating the class column
  class.agg <- function(x) {
    if ("promoter" %in% x) {
      return("promoter")
    } else if ("tss" %in% x) {
      return("tss")
    } else if ("genic" %in% x) {
      return("genic")
    } else if ("distal" %in% x) {
      return("distal")
    } else if ("intergenic" %in% x) {
      return("intergenic")
    } else {
      return("UNKNOWN")
    }
  }

  if ('class' %in% colnames(df)) {
    class.temp <- aggregate(df$class, by = list.for.agg, FUN = class.agg)
    colnames(class.temp)[colnames(class.temp) == 'x'] <- 'class'
    all.list$class <- class.temp
  }
  
  #Merge all the aggregates together to make collapsed dataframe  
  full.result <- Reduce(function(df1, df2) merge(df1, df2, by = agg.cols), all.list)
  
  return(full.result)
}

prepForPlotting <- function(df) {
  df$scatterplot.color <- with(df, ifelse(Regulated, "Activating", ifelse(Significant, "Repressive", "Not Significant")))
  return(df)
}

makePlots <- function(merged, config, inverse.predictors, pos.col, outdir, min.sensitivity = .7) {
  #config <- fread(config)  
  pred.cols <- unique(unlist(lapply(config$pred.cols, function(s) {strsplit(s,",")[[1]]})))
  pred.cols <- intersect(pred.cols, colnames(merged))
  
  #Make scatter plots
  lapply(pred.cols, function(s) {makeScatterPlot(merged, s, "EffectSize", outdir)})
  
  merged <- as.data.table(merged) #seems to be necessary to run in interactive session on compute node
  
  #Hack for predictors where lower values are more confident (eg genomic distance, pvalue)
  #Multiply these by -1
  inverse.predictors <- intersect(inverse.predictors, colnames(merged))
  merged[, inverse.predictors] <- -1*merged[, ..inverse.predictors]

  #Compute performance objects
  pr <- sapply(pred.cols, 
               function(s) {performance(prediction(unlist(merged[, ..s]), unlist(merged[, ..pos.col])), measure="prec", x.measure="rec")})
  merged[, inverse.predictors] <- -1*merged[, ..inverse.predictors]

  pr.df <- pr2df(pr)
  write.table(pr.df, file.path(outdir, "pr.curve.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

  #write PR summary table (AUC, cutoff, etc)
  perf.summary <- makePRSummaryTable(pr, min.sensitivity, outdir)

  #Assign prediction class and write output
  merged <- addPredictionClassLabels(merged, perf.summary)
  write.table(merged, file.path(outdir, "expt.pred.annotated.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  #Make PR curve plots
  pct.pos <- sum(unlist(merged[, ..pos.col]))/nrow(merged)
  for (ii in seq(nrow(config))) {
    makePRCurvePlot(pr.df, config$plot.name[[ii]], config$pred.cols[[ii]], outdir, pct.pos = pct.pos)
  }
}

makeScatterPlot <- function(df, x.col, y.col, outdir) {
  g <- ggplot(df,
              aes(x = get(x.col),
                  y = get(y.col),
                  color = scatterplot.color)) + 
    geom_point() + 
    scale_color_manual(values=c("Activating" = "red", 
                                "Repressive" = "blue", 
                                "Not Significant" = "gray")) + 
    labs(x = x.col, y = y.col, color = "")
  
  ggsave(file.path(outdir, paste0(x.col, ".", y.col, ".scatter.pdf")), g, device = "pdf")
}

makePRCurvePlot <- function(pr.df, plot.name, col.list, outdir, pct.pos) {
  col.list <- strsplit(as.character(col.list), ",")[[1]]
  pr.df <- subset(pr.df, pred.col %in% col.list)
  
  #separate boolean predictors from continuous predictors
  pr.cutoff <- by(pr.df, pr.df$pred.col, function(df) unique(df$alpha))
  boolean.predictors <- names(pr.cutoff)[unlist(lapply(pr.cutoff, function(s) identical(s, c(Inf, 1, 0))), use.names=F)]
  
  cont.pred <- subset(pr.df, !(pred.col %in% boolean.predictors))
  bool.pred <- subset(pr.df, pred.col %in% boolean.predictors)
  bool.pred <- subset(bool.pred, alpha == 1)
  
  g <- ggplot(cont.pred,
           aes(x = recall,
               y = precision,
               color = pred.col)) + 
      geom_line() + 
      labs(title = plot.name, color = "") + 
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
      geom_hline(yintercept = pct.pos, linetype = 2, color = 'black')
    
  if (nrow(bool.pred) > 0) {
    g <- g + geom_point(data = bool.pred,
                        size = 3)
  }
      
  ggsave(file.path(outdir, paste0(plot.name, ".pr.pdf")), g, device = "pdf")
}

pr2df <- function(pr) {
  #Convert a list of ROCR performance objects into a dataframe
  
  doOne <- function(this.pr) {
    df <- as.data.frame(list(
      alpha = this.pr@alpha.values[[1]],
      precision = this.pr@y.values[[1]],
      recall = this.pr@x.values[[1]]
    ))
    return(df)
  }
  
  pr.list <- lapply(pr, doOne)
  
  for (ii in seq(length(pr.list))) {
    pr.list[[ii]]$pred.col <- names(pr.list)[ii]
  }
  
  return(rbindlist(pr.list))
}

makePRSummaryTable <- function(pr, min.sensitivity = .7, outdir) {
  auc <- lapply(pr, function(s) computeAUC(s@x.values[[1]], s@y.values[[1]]))
  cutoff <- lapply(pr, function(s) computeCutoffGivenDesiredSensitivity(s, min.sensitivity))

  perf.summary <- rbindlist(list(cutoff = as.list(as.numeric(cutoff)), 
                                 AUC = as.list(as.numeric(auc))))

  perf.summary <- t(perf.summary)
  perf.summary <- as.data.table(cbind(names(pr), min.sensitivity, perf.summary))
  colnames(perf.summary) <- c("predictor", "min.sensitivity", "cutoff", "AUPRC")

  write.table(perf.summary, file.path(outdir, "pr.summary.txt"), sep='\t', quote = F, row.names = F, col.names = T)

  return(perf.summary)
}

computeAUC <- function(x.vals, y.vals) {
  good.idx <- which(!is.na(x.vals) & !is.na(y.vals))
  return(trapz(x.vals[good.idx], y.vals[good.idx]))
}

computeCutoffGivenDesiredSensitivity <- function(pr, min.sensitivity) {
  sens <- pr@x.values[[1]]
  prec <- pr@y.values[[1]]
  cutoff.sensitivity <- min(sens[sens >= min.sensitivity])
  idx <- which.max(sens == cutoff.sensitivity)
  idx2 <- idx[which.max(prec[idx])]

  cutoff <- pr@alpha.values[[1]][idx2]
  return(cutoff)
}

addPredictionClassLabels <- function(merged, perf.summary) {
  #assign prediction class label

  for (ii in seq(nrow(perf.summary))) {
    merged <- addOneLabel(merged, as.numeric(perf.summary$cutoff[[ii]]), perf.summary$predictor[[ii]])
  }

  return(merged)
}


addOneLabel <- function(df, cutoff, score.col, pos.col = "Regulated") {
  label.name <- paste0(score.col, ".pred.class")
  df[, label.name] <- "NA"

  #browser()

  df[which(!is.na(df[, ..score.col]) & df[, ..score.col] > cutoff & df[, ..pos.col]), label.name] <- "TP"
  df[which(!is.na(df[, ..score.col]) & df[, ..score.col] <= cutoff & !df[, ..pos.col]), label.name] <- "TN"
  df[which(!is.na(df[, ..score.col]) & df[, ..score.col] > cutoff & !df[, ..pos.col]), label.name] <- "FP"
  df[which(!is.na(df[, ..score.col]) & df[, ..score.col] <= cutoff & df[, ..pos.col]), label.name] <- "FN"

  return(df)
}

getInversePredictors <- function(pred.list, pred.config) {
  inv.pred <- subset(predConfig, lowerIsMoreConfident)$pred.col
  inv.cols <- with(expand.grid(pred.table$name, inv.pred), paste0(Var1, ".", Var2))
  return(inv.cols)
}

fread_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(f, fill = TRUE, ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

fread_gz_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(paste0("gunzip -c ", f), ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

smart_fread <- function(f, ...) {
  if (summary(file(f))$class == "gzfile") {
    out <- fread_gz_ignore_empty(f, ...)
  } else {
    out <- fread_ignore_empty(f, ...)
  }
  
  #con <- file(f)
  #on.exit(close(con), add = TRUE)
  closeAllConnections()
  return(out)
}

loadDelimited <- function(file.list) {
  data.list <- lapply(file.list, smart_fread)
  return(rbindlist(data.list))
}

loadFileString <- function(file.str, delim = ",") {
  file.list <- strsplit(file.str, split = delim)[[1]]
  return(loadDelimited(file.list))
}

loadPredictions <- function(pred.table) {
  #df <- fread(pred.table)
  pred.list <- lapply(pred.table$path, function(s) loadFileString(s))
  names(pred.list) <- pred.table$name
  return(pred.list)
}
