#! bin/python3

import pandas as pd
import numpy as np
import argparse
import os
from tools import write_params, run_command
from neighborhoods import count_single_feature_for_bed 

# TODO: Grab Gene Regions from ENSEMBL File : CdsStart CdsEnd to represent Gene Bodies for respective input into run.neighborhoods.py

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    epilog = ("")
    parser = argparse.ArgumentParser(description='Selects Top 2 Most Active TSS',
            epilog=epilog,
            formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--tss_file', required=required_args, help="tss isoform file")
    parser.add_argument('--h3k27ac', required=required_args, help="H3K27ac-seq bam file")
    parser.add_argument('--chrom_sizes', required=required_args, help="File listing chromosome size annotaions")
    parser.add_argument('--outDir', required=required_args)
    args = parser.parse_args()
    return args 

#def select_promoter_regions(output_file):
def read_tss_file(tss_file):
    """
    Reads in TSS File
    """
    tss_df = pd.read_csv(args.tss_file, sep="\t", names=['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand', 'start_Gene', 'end_Gene', 'TargetGene'])
    return tss_df

def filter_promoters_by_distance(promoters):
    """
    Takes in Promoter Isoforms and returns Top 2 Promoter Isoforms based on RPM 
    IF promoter start sites are 200bp from each other, otherwise, return top promoter
    """
    # ensure that promoters are at least 200bp apart
    top_promoter = promoters.iloc[0, :]
    # if no promoter isoform exists within 200bp, just pick top promoter 
    # for now, use the distance between the promoter TSS 
    start = top_promoter[1]
    try:
        temp = promoters.iloc[1:, :]
        temp['dist'] = temp[2] - start
        index = temp.loc[temp['dist'] >= 200].index.astype('int')
        second_promoter = temp.loc[index[0], :]
        top_promoter = pd.concat([top_promoter, second_promoter])
        return top_promoter
    except:
        return top_promoter

def process_genome_tss(args):
    """
    Takes in ENSEMBL Gene List and outputs 1/2 Promoter Isoforms for each gene 
    Promoter_ID = {Gene_Name}_{PromoterChr:Start-End}
    """
    os.makedirs(os.path.join(args.outDir), exist_ok=True)
    write_params(args, os.path.join(args.outDir, "params_generateTSS.txt"))
    
    filebase = str(os.path.basename(args.tss_file)).split(".")[0]
    feature_name =  "H3K27ac."+os.path.basename(args.h3k27ac)
    output_file = os.path.join(args.outDir,"{}.{}.CountReads.bedgraph".format(filebase, feature_name))
    tss_df = read_tss_file(args.tss_file)

    # Take in isoform file and count reads 
    tss_df = count_single_feature_for_bed(tss_df, args.tss_file, args.chrom_sizes, args.h3k27ac, "H3K27ac", args.outDir, "Genes.TSS1kb", skip_rpkm_quantile=False, force=False, use_fast_count=True)
    chrom_sizes = args.chrom_sizes
    tss_file = args.tss_file
    sort_command = "bedtools sort -faidx {chrom_sizes} -i {tss_file} > {tss_file}.sorted; mv {tss_file}.sorted {tss_file}".format(**locals())
    run_command(sort_command)
    print("Finished Running code")
    # sorted tss1kb_file
    tss1kb_file = read_tss_file(args.tss_file)
    tss1kb_file['H3K27ac.RPM'] = tss_df[feature_name+'.RPKM']

    # This loop is taking a while; figure out speed up! 
    # Take top 2 promoters based on counts 
    gene_tss_df = None

    for gene in tss1kb_file['TargetGene']: 
        tss1kb_file_subset = tss1kb_file.loc[tss1kb_file['TargetGene'].str.contains(gene)]
        sorted_tss1kb_file_subset = tss1kb_file_subset.sort_values(by=['H3K27ac.RPM'], ascending=False)

        # ensure that distances between promoters are at least 200bp from each other
        top_two = filter_promoters_by_distance(sorted_tss1kb_file_subset)
        if gene_tss_df is None:
            gene_tss_df = top_two
        else:
            gene_tss_df = pd.concat([gene_tss_df, top_two])
        print(gene_tss_df)
    file_output = "H3K27ac."+os.path.basename(args.h3k27ac)+"_Expressed"
    gene_tss_df[['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand']].to_csv("ENSEMBL_GeneTSS.txt", sep="\t", index=False, header=False)
    gene_tss_df[['chr', 'start_Gene', 'end_Gene', 'TargetGene', 'score', 'strand']].to_csv("ENSEMBL_Genes.txt", sep="\t", index=False, header=False)

if __name__=="__main__":
    args = parseargs()
    process_genome_tss(args)

