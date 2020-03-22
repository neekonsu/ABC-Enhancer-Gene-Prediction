#! bin/python3

from metrics import *

if __name__=="__main__":
    macs_peaks = "/mnt/lab_data2/kmualim/data/send_to_Kristy/FilesToReproduceABC/PeakAndNeighborhoodFiles/Peaks_tmp/ENCFF030DCL.macs2_peaks.narrowPeak"
    genome_tss = "/users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
    outdir="."
    neighborhood_dir = "/mnt/lab_data2/kmualim/data/send_to_Kristy/FilesToReproduceABC/PeakAndNeighborhoodFiles//Neighborhoods"
    prediction="/mnt/lab_data2/kmualim/data/send_to_Kristy/FilesToReproduceABC/PredictionResults/Predictions_for_all_chr_no_threshold/EnhancerPredictionsAllPutative.txt.gz"
    grab_nearest_tss_from_peak(macs_peaks, genome_tss, outdir)
    prediction_df = pd.read_csv(prediction, sep="\t", compression="gzip")
    GrabQCMetrics(prediction_df, outdir)
    PeakFileQC(macs_peaks, outdir)
    NeighborhoodFileQC(neighborhood_dir, outdir)
    enhancers = "/mnt/lab_data2/kmualim/data/send_to_Kristy/FilesToReproduceABC/PeakAndNeighborhoodFiles//Neighborhoods/EnhancerList.txt"
    EnhancerList = pd.read_csv(enhancers, sep="\t")
    title="_QuantileNorm"
    PlotQuantilePlot(EnhancerList, title, outdir)
