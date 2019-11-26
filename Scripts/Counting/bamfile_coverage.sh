#!/bin/bash

# 1 {input.ucsc2ncbi}
# 2 {output.counted}
# 3.. {params.sample_ids} {input.bam_files}

ucsc2ncbi=$1
output=$2
shift 2
# list of sample ids followed by their BAM files in the same order
samples=("$@")

# use tmp output to avoid snakemake from stopping before script is done
tmp=$(mktemp)
echo -e "sampleID\trecord_count" > $tmp # add head to file

# get chromosome names of different formats
chrNamesUCSC=$(cut -f1 $ucsc2ncbi | tr '\n' '|')
chrNamesNCBI=$(cut -f2 $ucsc2ncbi | tr '\n' '|')

# get number of samples to split args list into sample ID list and BAM file list
n_samples=$(expr ${#samples[@]} / 2)

# for each sample
for i in $(seq 1 $n_samples)
do
    sampleID=${samples[$i-1]}
    bam_file=${samples[$(expr $i + ${n_samples} - 1)]}

    # identify chromosome format
    bam_chr=$(samtools idxstats ${bam_file} | grep chr | wc -l)
    if [ ${bam_chr} -ne 0 ]
    then
        chrNames=$chrNamesUCSC
    else
        chrNames=$chrNamesNCBI
    fi

    # write coverage from idxstats into file
    count=$(samtools idxstats ${bam_file} | grep -E "^(${chrNames})" | cut -f3 | paste -sd+ - | bc)
    echo -e "${sampleID}\t${count}" >> $tmp
done

cp $tmp $output
rm $tmp

