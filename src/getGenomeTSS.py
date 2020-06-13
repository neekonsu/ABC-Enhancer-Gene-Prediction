#! bin/python3

import pandas as pd
import numpy as np
import argparse
import os
from tools import write_params, run_command
from neighborhoods import count_features_for_bed, count_single_feature_for_bed 


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    epilog = ("")
    parser = argparse.ArgumentParser(description='Selects Top 2 Most Active TSS',
            epilog=epilog,
            formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--tss_file', required=required_args, help="tss isoform file")
    parser.add_argument('--dhs', required=required_args, help="Accessibility bam file")
    parser.add_argument('--h3k27ac', required=required_args, help="H3K27ac-seq bam file")
    parser.add_argument('--chrom_sizes', required=required_args, help="File listing chromosome size annotaions")
    parser.add_argument('--gene_outf', required=required_args, help="GeneList output")
    parser.add_argument('--genetss_outf', required=required_args, help="GeneListTSS output")
    parser.add_argument('--outDir', required=required_args)
    args = parser.parse_args()
    return args 

def read_tss_file(tss_file):
    """
    Reads in TSS File
    """
    tss_df = pd.read_csv(args.tss_file, sep="\t", names=['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand', 'start_Gene', 'end_Gene', 'TargetGene'])
    return tss_df

def create_dataframes(data, len_):
    df = pd.DataFrame()
    df['chr'] = data.iloc[:len_, 0].values
    df['start'] = data.iloc[len_:(len_*2), 0].values
    df['end'] = data.iloc[(len_*2): (len_*3), 0].values
    df['TargetGene'] = data.iloc[(len_*3): (len_*4), 0].values
    df['score'] = data.iloc[(len_*4): (len_*5), 0].values
    df['strand'] = data.iloc[(len_*5):, 0].values
    return df

def filter_promoters_by_distance(promoters):
    """
    Takes in Promoter Isoforms and returns Top 2 Promoter Isoforms based on RPM 
    IF promoter start sites are 500bp from each other, otherwise, return top promoter
    """
    # ensure that promoters are at least 500bp apart
    top_promoter = promoters.iloc[0, :]
    # if no promoter isoform exists within 500bp, just pick top promoter 
    # for now, use the distance between the promoter TSS 
    start = top_promoter[1]
    try:
        temp = promoters.iloc[1:, :]
        temp['dist'] = temp[2] - start
        index = temp.loc[temp['dist'] >= 500].index.astype('int')
        second_promoter = temp.loc[index[0], :]
        concatenated_top_promoters = pd.concat([top_promoter, second_promoter])
        top_promoter = filter_promoters_by_activity(concatenated_top_promoters)
        return top_promoter
    except:
        return top_promoter

def filter_promoters_by_activity(promoters):
    activities = promoters['PromoterActivityQuantile']
    activity_foldchange = activities[1]/activities[0]
    if activity_foldchange > 0.8: 
        return promoters 
    else:
        return promoters.iloc[0, :]

def filter_expressed_df(expressed_tsscounts):
    gene_tss_df = None
    for gene in expressed_tsscounts['TargetGene']:
        tss1kb_file_subset = expressed_tsscounts.loc[expressed_tsscounts['TargetGene']==gene]
        sorted_tss1kb_file_subset = tss1kb_file_subset.sort_values(by=['PromoterActivityQuantile'], ascending=False)
        # ensure that distances between promoters are at least 500bp from each other
        top_two = filter_promoters_by_distance(sorted_tss1kb_file_subset)
        if gene_tss_df is None:
            gene_tss_df = top_two
        else:
            gene_tss_df = pd.concat([gene_tss_df, top_two])
    return gene_tss_df

def filter_nonexpressed_df(nonexpressed):
    final_df = None
    for gene in nonexpressed['TargetGene']:
        matched = nonexpressed['PromoterActivityQuantile'].loc[nonexpressed['TargetGene'].str.contains(str(gene).split("-")[0])]
        max_index = np.sort(matched).index.astype('int')
        max_match = nonexpressed.iloc[max_index, :]
        if final_df is None:
            final_df = max_match
        else:
            final_df = pd.concat([final_df, max_match])
        return final_df

def process_genome_tss(args):
    """
    Takes in ENSEMBL Gene List and outputs 1/2 Promoter Isoforms for each gene 
    Promoter_ID = {Gene_Name}_{PromoterChr:Start-End}
    """
    os.makedirs(os.path.join(args.outDir), exist_ok=True)
    write_params(args, os.path.join(args.outDir, "params_generateTSS.txt"))
    
    filebase = str(os.path.basename(args.tss_file)).split(".")[0]
    feature_name =  ["H3K27ac", "DHS"]
    feature_files = [args.h3k27ac, args.dhs]
    features = {key:value for key, value in zip(feature_name, feature_files)}
    tss_df = read_tss_file(args.tss_file)
    tss1kb_file = args.tss_file
    genome_sizes = args.chrom_sizes
    outdir = args.outDir

    tsscounts = count_features_for_bed(tss_df, tss1kb_file, genome_sizes, features, outdir, "Genes.TSS1kb", force=True, use_fast_count=True)
    for feature, feature_file in zip(feature_name, feature_files):
        output_file = os.path.join(args.outDir,"{}.{}.CountReads.bedgraph".format(filebase, feature))
        print("Taking in isoform TSS file and generating Counts")
        # Take in isoform file and count reads 
        tss_df_1 = count_single_feature_for_bed(tss_df, args.tss_file, args.chrom_sizes, feature_file, feature, args.outDir, "Genes.TSS1kb", skip_rpkm_quantile=False, force=False, use_fast_count=True)
    chrom_sizes = args.chrom_sizes
    tss_file = args.tss_file
    sort_command = "bedtools sort -faidx {chrom_sizes} -i {tss_file} > {tss_file}.sorted; mv {tss_file}.sorted {tss_file}".format(**locals())
    run_command(sort_command)
    print("Finished Sorting Gene TSS File")

    # Take top 2 promoters based on counts 
    tsscounts['PromoterActivityQuantile'] = ((0.0001+tsscounts['H3K27ac.RPKM.quantile'])*(0.0001+tsscounts['DHS.RPKM.quantile'])).rank(method='average', na_option="top", ascending=True, pct=True)
    print("Looping though all genes present to select out Top Two Promoters based on RPM")
    tsscounts.to_csv(os.path.join(args.outDir, "PromoterActivityQuantile.tsv"), sep="\t", index=False)
    # filter for expressed genes 
    # This loop only needs to run on expressed genes
    expressed_tsscounts = tsscounts.loc[tsscounts['PromoterActivityQuantile']!=0.0]
    filtered_expressed_tsscounts = filter_expressed_df(expressed_tsscounts)

    nonexpressed_dup = tsscounts.loc[tsscounts['PromoterActivityQuantile']==0.0]
    # filter for single promoter entries 
    nonexpressed_unique = filter_nonexpressed_df(nonexpressed_dup)
    gene_tss_df = pd.concat([filtered_expressed_tsscounts, nonexpressed_unique])

    print("Saving Files")
    gene_tss_df[['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand']].to_csv("ENSEMBL_GeneTSS.txt.tmp", sep="\t", index=False, header=False)
    gene_tss_df[['chr', 'start_Gene', 'end_Gene', 'TargetGene', 'score', 'strand']].to_csv("ENSEMBL_Genes.txt.tmp", sep="\t", index=False, header=False)

def concatenate_entries(args):
    data = pd.read_csv("ENSEMBL_Genes.txt.tmp", sep="\t", header=None)
    data_tss = pd.read_csv("ENSEMBL_GeneTSS.txt.tmp", sep="\t", header=None)
    print("Concatenating Entries for final output file")
    len_ = len(data.loc[data[0].str.contains('chr')])
    df = create_dataframes(data, len_)
    df.to_csv(os.path.join(args.outDir, args.gene_outf), sep="\t", header=False, index=False)
    df2 = create_dataframes(data_tss, len_)
    df2.to_csv(os.path.join(args.outDir, args.genetss_outf), sep="\t", header=False, index=False)
    
def clean_up():
    print("Removing tmp files created during run")
    os.remove("ENSEMBL_GeneTSS.txt.tmp")
    os.remove("ENSEMBL_Genes.txt.tmp")
    print("Done!")

if __name__=="__main__":
    args = parseargs()
    process_genome_tss(args)
    concatenate_entries(args)
    clean_up()
