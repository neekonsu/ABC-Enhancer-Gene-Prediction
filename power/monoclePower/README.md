Sample command:

```
Rscript /seq/lincRNA/jnasser/EP_prediction/code/LanderLab-EP-Prediction/src/CRISPRScreen//runMonocleOnCDS.R \
--outDir /seq/lincRNA/RAP/External/Gasperini2019/200304_NewRegions/PowerCalc///TargetRegions/Regions//NewRegion1/ \
--cds /seq/lincRNA/RAP/External/Gasperini2019/GEO/gasperini_processed/accessed_on_190423/GSE120861_at_scale_screen.cds.rds \
--geneList /seq/lincRNA/RAP/External/Gasperini2019/200304_NewRegions/PowerCalc///TargetRegions/Regions//NewRegion1/genes.txt \
--fullModelStr guide_count+percent.mito+prep_batch \
--reducedModelStr guide_count+percent.mito+prep_batch \
--grepCol barcode \
--grep "ACTCCTTTACAGGTTTCATG|GTACAGTTATGTAAAGATGA" \
--region.id NewRegion1 --useDispFit --adjustDispByEffectSize \
--nPermsForPowerCalc 100 --powerEffectSizes .25 \
>& /seq/lincRNA/RAP/External/Gasperini2019/200304_NewRegions/PowerCalc///TargetRegions/Regions//NewRegion1/log.qout
```
