#!/bin/bash

# 1 {input.ucsc2ncbi}
# 2 {output.counted}
# 3.. {params.sample_ids} {input.bam_files}

ucsc2ncbi=$1
output=$2
shift 2
samples=("$@")

tmp=$(mktemp)
echo -e "sampleID\trecord_count" > $tmp

chrNamesUCSC=$(cut -f1 $ucsc2ncbi | tr '\n' '|')
chrNamesNCBI=$(cut -f2 $ucsc2ncbi | tr '\n' '|')

n_samples=$(expr ${#samples[@]} / 2)

for i in $(seq 1 $n_samples)
do
    sampleID=${samples[$i-1]}
    bam_file=${samples[$(expr $i + ${n_samples} - 1)]}

    # get chromosomes
    bam_chr=$(samtools idxstats ${bam_file} | grep chr | wc -l)
    if [ ${bam_chr} -ne 0 ]
    then
        chrNames=$chrNamesUCSC
    else
        chrNames=$chrNamesNCBI
    fi

    count=$(samtools idxstats ${bam_file} | grep -E "^(${chrNames})" | cut -f3 | paste -sd+ - | bc)
    echo -e "${sampleID}\t${count}" >> $tmp
done

cp $tmp $output
rm $tmp

