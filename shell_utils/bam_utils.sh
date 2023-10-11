#!/usr/bin/bash

##### utils functions to fetch alignments, reads and other operation for bam files
# includes: 
# fetch_alignments: ./bam_utils.sh fetch_alignments bam_file_paths.txt sites_file.txt out_dir
# fetch_fastq: ./bam_utils.sh fetch_fastq bam_file_paths.txt sites_file.txt out_dir

#@ utils
function fetch_fastq() {
    # main: fetch fastq from bam files
    # usage: fetch_fastq bam_file_paths.txt sites_file.txt out_dir
    # result_name: fq_prefix_sites.fastq
    # eg: 17HanZZ0002_chr6_396535_397106.fastq
    file_paths=$1
    sites_file=$2
    out_dir=$3
    
    echo "Processing sites file: $sites_file"
    echo "Out dir: $out_dir"
    
    cat $sites_file | while read site
    do
        echo "Fetching alignments for $site"
        site_a=($site)
        site_s=$((${site_a[1]}-10000)) # 向左扩10kb
        site_e=$((${site_a[2]}+10000)) # 向右扩10kb
        f_site="${site_a[0]}:${site_s}-${site_e}"
        
        # continue
        num=0
        cat $file_paths | while read file
        do
            filename=${file##*/}
            basename=${filename%.*}
            result_file=${basename}_${site_a[0]}_${site_a[1]}_${site_a[2]}.fastq
            num=$((num+1))
            
            echo "result file: "$num" "$out_dir/$result_file
            cmd="samtools view -b $file $f_site | samtools sort -n | samtools fastq > $out_dir/$result_file"
            # echo "cmd: $cmd"
            $samtools view -b $file $f_site | $samtools sort -n | $samtools fastq > $out_dir/$result_file
        done
    done
}


