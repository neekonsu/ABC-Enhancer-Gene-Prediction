#! bin/python3

import pandas as pd
import numpy as np
import argparse


def parse_args(required_args=True):
    parser = argparse.ArgumentParser(description="Running QC metrics across different ABC maps")
    readable = argparse.FileType('r')
    parser.add_argument('--celltype_file', required=True, help="Path of different celltype files to compare")
    # replace textio wrapper returned by argparse with actual filename
    args = parser.parse_args()
    for name, val in vars(args).items():
        if hasattr(val, 'name'):
            setattr(args, name, val.name)
    print(args)
    return args

def RunQC(args):
    celltypes = pd.read(args.celltype_file, sep="\t", header=None)
    median_values_dict = {}
    cells = celltypes[0]
    peak_paths = celltypes[1]
    neighborhoods_paths = celltypes[2]
    index=0
    for k in cells:
        data = pd.read_csv(os.path.join(paths[index], "Summary.txt")), sep="\t", header=None)
        median_vals = read_file(peak_paths, neighborhoods_paths, index)
        median_values_dict[k] = median_vals
    
def read_file(peak_paths, neighborhood_paths, index):
    data = pd.read_csv(os.path.join(peak_paths[index], "PeakFileQCSummary.txt")), sep="\t", header=None)
    neighborhoods_data = pd.read_csv(os.path.join(neighborhood_paths[index], "QCSummary.txt")), sep="\t", header=None)
    median_vals = [str(i).split(':')[1] for i in data[0]]
    neighbor_median_vals = [str(i).split(':')[1] for i in neighborhoods_data[0]]
    return median_vals + neighbor_median_vals

def makeDataFrame(median_values_dict):
    epg=[]
    gpe=[]
    numpeaks=[]
    widthpeaks=[]
    dist=[]
    num_candidate=[]
    med_dist=[]
    10th=[]
    90th=[]
    celltype=[]
    num_enhancers_per_chrom=[]
    med_width_can=[]
    for i in median_values_dict.keys():
        celltype.append(i)
        peaks.append(dict_a[i][0])
        widthpeaks.append(dict_a[i][1])
        med_dist.append(dict_a[i][2])
        num_candidate.append(dict_a[i][3])
        med_width_can.append(dict_a[i][4])
        epg.append(dict_a[i][5])
        gpe.append(dict_a[i][6])
        dist.append(dict_a[i][7])
        num_enhancers_per_chrom(dict_a[i][8])
        10th.append(dict_a[i][9])
        90th.append(dict_a[i][10])

    df = pd.DataFrame()
    df['EnhancerPerGene'] = epg
    df['GenesPerEnhancer'] = gpe
    df['NumberOfPeaks'] = peaks
    df['WidthofPeaks'] = widthpeaks
    df['NumberOfCandidateRegions'] = num_candidate
    df['WidthOfCandidateRegions'] = med_width_can
    df['E-Gdistance'] = dist
    df['E-G10thquantile'] = 10th
    df['E-G90thquantile'] = 90th
    df['NumEnhancersPerChrom'] = num_enhancers_per_chrom
    df['celltype'] = celltype
    return df
