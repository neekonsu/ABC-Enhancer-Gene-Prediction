The files in this directory allows for the Comparison between ABC Predictions and Experimental Data: 

## Requirements 
       - R 3.4.0
       - data.table (1.11.4)
       - GenomicRanges (1.28.6)
       - ROCR (1.0-7)
       - ggplot2 (3.0.0)
       - caTools (1.17.1

Sample Command:
```
Rscript comparePredictionsToExperiment.R --experimentalData example_data/EnhancerPredictionsAllPutative_chr22.txt --predictions data/ExperimentalData.Gasperini.FulcoNasser.191021.txt --plotConfig CRISPR/plot.config.txt --predConfig CRISPR/pred.config.txt
```

## Description of Comparison Code 
This code is intended to evaluate a predictive model against a set of experimental data
1. Generates PR curves and scatter plots

## Specifications
plotConfig file should contain one line per PR curve plot. Scatter plots will be generated for all prediction columns defined in predConfig

## Configuring Prediction Columns:
In the ideal use case each experimental element should correspond to exactly one element-gene prediction. 
However, there are instances where a single experimental element overlaps multiple predictions (eg a large deletion) or when an experimentally tested element-gene pair is not present in the predictions files.
The predConfig file describes how to handle these cases. 

## Other:
- Assumes prediction columns are monotonic increasing! (A hack is employed for distance)
- Can't distinguish between a missing prediction and a non-prediction

