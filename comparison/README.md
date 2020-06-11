# Evaluate Enhancer-Gene prediction models against experimental data

This codebase is designed to evaluate the performance of an enhancer-gene linking model against experimental data. This codebase supports the evaluation of multiple predictors against a single experimental data file. Comparing performance of a single predictor against multiple experiments is not natively supported.

Fundamentally, the idea is to overlap each experimentally tested element-gene-celltype tuple with a predicted element-gene-celltype tuple. Diagnostic plots such as PR curves are then produced on this overlapped dataset. Care must be taken when an experimentally tested element overlaps multiple predicted elements, or when an experimentally tested element does not overlap any predicted elements. See the configuration section below for how to handle these cases. 

Other notes

 * There must be both experimental positives and negatives in the experimental data file in order to produce PR curves. If only positives are available, the plotting module will fail but useful intermediate files will still be produced
 * Each prediction file may have multiple score columns. The names of these columns may be the same across prediction files - they will be differentiated based on dataset name (as defined by the prediction table). The combination of dataset name and prediction column should be unique. TO DO: Check for this.
 * The code currently overlaps on GeneSymbol (not on gene TSS coordinates)

## Requirements
 * Required Inputs (see below for formats)
 	* One experimental data file 
 	* At least one predictions file
 	* Configuration file describing how predictor should be aggregated
 	
### Dependencies
```
R (3.4.0)

R packages:
data.table (1.11.4)
GenomicRanges (1.28.6)
ROCR (1.0-7)
ggplot2 (3.0.0)
caTools (1.17.1)
```

## File Formats

* Experimental Data
  * <https://docs.google.com/spreadsheets/d/1Tl_fdPeAeiVkZnettWxeMJW3zL5zfFRcC2TV-vmQRLQ/edit#gid=257851280>. 
  * Not all of these columns are required. See example/input/K562.ExperimentalData.slim.txt for list of required columns

* Predictions
  * <https://docs.google.com/spreadsheets/d/1BQBFC4PzPA8v3tA_OkSp2lU1YpO74uwmkEa2TJTt-ic/edit#gid=0>
  * See example/input/K562.ABC.Predictions.AvgHiC.chrX.ENCODE.format.txt.gz for example. 

## Configuring predictions

Each predictor must have a corresponding entry in the predConfig file. The behavior of the comparison code for this predictor depends on:
 
 * pred.col: must match the name of the column in the predictions file
 * agg.func: In the case that an experimentally tested element overlaps multiple predicted elements, how should the predicted elements be aggregated to the level of the tested element. 
 * fill.val: In the case that an experimentally tested element does not overlap any predicted elements, what value should be filled in.
 * lowerIsMoreConfident: Set to TRUE if lower values of the predictor signify more confidence in the prediction. This is appropraite for predictors such as linear distance or pvalue. It is generally preferred for this column to be FALSE.

## Configuring Plotting
Each row of this file denotes a separate PR curve plot. Each predictor should be of the form: {DATASETNAME}.{SCORECOLUMN}

## Sample command
```
Rscript src/comparePredictionsToExperiment.R \
--predictions example/input/pred.table.txt \
--experimentalData example/input/K562.ExperimentalData.slim.txt \
--experimentalPositiveColumn "Regulated" \
--plotConfig src/plot.config.txt \
--predConfig src/pred.config.txt \
--code src/comparison.R \
--outDir example/out

```

