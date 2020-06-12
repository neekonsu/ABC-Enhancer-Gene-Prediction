

conda activate /seq/lincRNA/RAP/ABC/200330_ENCODE/env/abcenv

PROJECT=/seq/lincRNA/jnasser/EP_prediction/code/ABC_code_snakemake2/; 
CODEDIR=$PROJECT/ABC-Enhancer-Gene-Prediction/; cd $CODEDIR

snakemake -s "$CODEDIR/Snakemake/Snakefile" \
        --configfile $CODEDIR/Snakemake/abc.yaml \
        --config sample_json=$CODEDIR/Snakemake/K562_chr22.json