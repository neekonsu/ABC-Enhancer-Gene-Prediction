version 1.0

workflow ABCpipeline {
    meta {
        description: "WDL version of the ABC pipeline"
    }
 
    input {
        File dnaseqbam
        File chrom_sizes
        File regions_blacklist
        File regions_whitelist
    }

    call makeCandidateRegions {
       input:
           bam = dnaseqbam,
           chrom_sizes = chrom_sizes,
           regions_blacklist = regions_blacklist,
           regions_whitelist = regions_whitelist
    }

    output {
       File candidateRegions = makeCandidateRegions.candidateRegions
    }
}



    task makeCandidateRegions {
        input {
            File bam
            File chrom_sizes
            File regions_blacklist
            File regions_whitelist
            Float pval_cutoff = 0.1
            Int peakExtendFromSummit = 250
            Int nStrongestPeaks = 3000
        }

        String docker_image = "quay.io/nbarkas/abc-general-container:latest"
        Int num_threads = 1
        String mem_size = "1 GB"
        

        command {
            set -euo pipeline

            mkdir outputs

            python src/makeCandidateRegions.py \
                --bam ~{bam} \
                --outDir outputs \
                --chrom_sizes ~{chrom_sizes} \
                --regions_blacklist ~{regions_blacklist} \
                --regions_whitelist ~{regions_whitelist} \
                --pval_cutoff ~{pval_cutoff} \
                --peakExtendFromSummit ~{peakExtendFromSummit} \
                --nStrongestPeaks ~{nStrongestPeaks}
        }
        output {
            # TODO: Add all the outputs
            File candidateRegions = "candidateRegions.bed"
        }
        runtime {
            docker: docker_image
            cpu: num_threads
            memory: mem_size
            disks: "local-disk" + ceil(size(bam, "GiB")) * 1.2
        }
    }


