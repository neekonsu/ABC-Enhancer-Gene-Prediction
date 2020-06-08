# Evaluate Enhancer-Gene prediction models against experimental data

The following configuration files should be edited according to the specific analysis:

```
src/pred.config.txt: Defines how each predictor should handle aggregation, missingness, etc

```

```
src/plot.config.txt: Defines which PR curves to make

```

Example Command:
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
