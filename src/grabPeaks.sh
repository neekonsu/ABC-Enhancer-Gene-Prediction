#! bin/bash 

#while read p;
#do
#	a=($p)
#	echo ${a[0]} ${a[1]}
#	macs2 callpeak -f BAM -g hs -p .1 --call-summits --outdir /srv/scratch/kmualim/ABC_data/PeakAndNeighborhoodFiles/Peaks_${a[0]} -t /srv/scratch/kmualim/ABC_data/ENCODEdata/cellline_files/${a[1]}
#done < jill_data_lookup.txt  


while read p;
do
	a=($p)
#	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak /srv/scratch/kmualim/ABC_data/PeakAndNeighborhoodFiles/Peaks_${a[0]}/NA_peaks.narrowPeak --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --regions_blacklist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed --regions_whitelist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 150000 --outDir /srv/scratch/kmualim/ABC_data/PeakAndNeighborhoodFiles/Peaks_${a[0]}/ --bam $2/${a[1]} --genome_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed	
#	python /users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py --candidate_enhancer_regions  /mnt/lab_data3/kmualim/PeakAndNeighborhoods/Peaks_MCF7/NA_peaks.narrowPeak.candidateRegions.bed  --genes /users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed --H3K27ac $2/${a[2]} --DHS $2/${a[1]} --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --outdir Neighborhoods_${a[0]}_qnorm --qnorm /users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/src/EnhancersQNormRef.K562.txt 
	python /users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/src/predict.py --scale_hic_using_powerlaw --threshold 0.02 --make_all_putative --hic_resolution 5000 --HiCdir /mnt/lab_data2/kmualim/data/send_to_Kristy/HiC/AverageHiC/raw --enhancers Neighborhoods_${a[0]}_qnorm/EnhancerList.txt --genes Neighborhoods_${a[0]}_qnorm/GeneList.txt --outdir Predictions_${a[0]}_qnorm --cellType ${a[0]} --chromosome ['chr22','chr21']
done < $1 

