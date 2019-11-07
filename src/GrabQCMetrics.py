#! bin/python3

import pandas as pd
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--preds', '-p',
                        help="Predictions File")
    parser.add_argument('--threshold', '-t',
                        default=0.022, type=float,
                        help="Threshold value for filtering enhancer gene links. Set as 0.022 default, following the Paper Threshold")
    parser.add_argument('--peakfile', '-k', default="",
			help="path to Peak file")
    return parser.parse_args()

# Grab QC Summary 
# Generates box-and-whisker plot for quantile-normalized and non-normalized DHS, H3K27ac data 
# Answers questions like: 
# Enhancers per chromosome
# Enhancers per gene
# Genes per enhancer
# E-G distance

def PeakFileQC(peakfile):
    if peakfile.endswith(".gz"):
        peaks = pd.read_csv(peakfile, compression="gzip", sep="\t", header=None)
    else:
        peaks = pd.read_csv(peakfile, sep="\t", header=None)
    peaks['dist'] = peaks[2]-peaks[1]
    peaks_array = list(peaks['dist'])
    PlotDistribution(peaks_array, WidthOfPeaks)
    
    with open("PeakFileQCSummary.txt","w") as f:
        f.write(str(peakfile))
        f.write("\n")
        f.write("Number of peaks: ")
        f.write(str(len(peaks['dist'])))
        f.write("\n")
        f.write("Max width of peak: ")
        f.write(max(peaks['dist']))
        f.write("\n")
        f.write("Mean and Stdev width of peaks: ")
        f.write(str(peaks['dist'].mean()))
        f.write("\t")
        f.write(str(peaks['dist'].std()))
        f.write("\n")
        f.close()

def PlotDistribution(array, title):
    ax = sns.distplot(array)
    ax.set_title(title)
    ax.set_ylabel('Estimated PDF of distribution')
    ax.set_xlabel('Counts')
    ax.get_figure()
    fig.savefig(str(title)+".png")

# Plotting Quantile Data
def PlotQuantileData(data):
    quantile_norm_df = pd.DataFrame()
    quantile_norm_df['readCount.quantile'] = thresholded_data['H3K27ac.ENCFF384ZZM_pooled.fc.bigWig.readCount.quantile']
    quantile_norm_df['dataset'] = 'H3K27ac.readCount.quantile'

    quantile_norm = pd.DataFrame()
    quantile_norm['readCount.quantile'] = thresholded_data['DHS.ENCFF801RMG.merged.nodup.pooled.fc.signal.bigwig.readCount.quantile']
    quantile_norm['dataset'] = "DHS.readCount.quantile"

def GrabQCMetrics(prediction_df):
#    data = pd.read_csv(preds, sep="\t")
#    genes = data['TargetGene'].drop_duplicates()
    # isolate dataframe to consider just genic or intergenic regions
#    subset_data = data.loc[data['class']=='genic' or data['class']=='intergenic']
    # set threshold for enhancer-gene links generated
    # setting this threshold to 0.22
#    thresholded_data = subset_data.loc[subset_data['ABC.Score']> args.threshold]
    # Grab thresholded_data summary
    # Grabs enhancers per gene information 
    GeneCounts = prediction_df.groupby(['TargetGene']).size()
    GeneCounts.to_csv("EnhancerPerGene.txt", sep="\t")
  
    GeneMean = prediction_df.groupby(['TargetGene']).size().mean()
    GeneStdev = prediction_df.groupby(['TargetGene']).size().std()
    # Grab Number of genes per enhancers 
    num_enhancers = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size()
    num_enhancers.to_csv("GenesPerEnhancer.txt", sep="\t")
    mean_genes_per_enhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size().mean()
    stdev_genes_per_enhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size().std()
    
    # Grab Number of Enhancer-Gene Pairs Per Chromsome
    enhancergeneperchrom = prediction_df.groupby(['chr']).size()
    enhancergeneprechrom.to_csv("EnhancerGenePairsPerChrom.txt", sep="\t")
    
    # Enhancer-Gene Distancee
    prediction_df['dist'] = prediction_df['end'] - prediction_df['start']
    distance = list(prediction_df['dist'])
    # Plot Distributions and save as png    
    PlotDistribution(num_enhancers, NumberOfGenesPerEnhancer)
    PlotDistribution(GeneCounts, NumberOfEnhancersPerGene)
    PlotDistribution(enhancergeneperchrom, EnhancersPerChromosome)
    PlotDistribution(distance, EnhancerGeneDistance)
     
    with open("QCSummary.txt", "w") as f:
        f.write("Average Number of Enhancers per Gene: ")
        f.write(str(GeneMean))
        f.write("\n")
        f.write("Standard Deviation of Enhancers per Gene:")
        f.write(str(GeneStdev))
        f.write("\n")
        f.write("Average Number of Genes linked to an Enhancer:")
        f.write(str(mean_genes_per_enhancer))
        f.write("\n")
        f.write("Standard Deviation of Genes linked to an Enhancer:")
        f.write(str(stdev_genes_per_enhancer))
        f.write("\n")
        f.write("Mean Enhancer-Gene Distance:")
        f.write(str(prediction_df['dist'].mean()))
        f.write("\n")
        f.write("Standard Deviation of Enhancer-Gene Distance:")
        f.write(str(prediction_df['dist'].std()))
        f.write("\n")
        f.close()   

def main():
    args = parse_args()
    GrabQCMetrics(args)
    if args.peakfile != "":
        PeakFileQC(args.peakfile)
