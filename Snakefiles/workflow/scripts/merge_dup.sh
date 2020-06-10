#! bin/bash 

workingdir="/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/data/ENCODEPortal/workflow/scripts"
indir=$1
outdir=$2
input_file=$3
threads=$4

m=()
while read p;
do
	a=($p)
	echo ${a[1]} ${a[2]} ${a[3]}
	len_=$((${#a[@]}-1))
	for file in $(eval echo {1..$len_});
	do
		echo $file
		m+=("$indir/${a[$file]}")
	done

	samtools merge $outdir/${a[$len_]} $m &
#	then
#		samtools merge $outdir/${a[4]} $indir/${a[1]} $indir/${a[2]} $indir/${a[3]} --threads $threads &
#	else
#		samtools merge $outdir/${a[3]} $indir/${a[1]} $indir/${a[2]} --threads $threads &
#	fi
done < 	$input_file

