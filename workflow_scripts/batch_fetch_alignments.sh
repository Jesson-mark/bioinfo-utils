#!/usr/bin/bash

hifi_bam_dir=$out_dir/hpgp_bam/hifi_bam
cat inputs/hpgp_hifi_bam_paths.txt | while read bam
do
    filename=${bam##*/}
    basename=${filename%.*}

    out_bam=${hifi_bam_dir}/${basename}.bam
    samtools view -bh $bam --region-file $regions_flanked > $out_bam
    samtools index $out_bam
    echo "done for $bam"
done

# softwares
samtools="/public/home/fan_lab/wangjie/Programs/samtools-1.13/bin/samtools"

# main: fetch alignments from bam files
# usage: fetch_alignments bam_file_paths.txt sites_file.txt out_dir
# result_name: bam_prefix_sites.bam
# eg: 17HanZZ0002_chr6_396535_397106.bam
if [[ $# -eq 0 ]]; then
    echo "Usage: fetch_alignments bam_file_paths.txt sites_file.txt out_dir"
    echo '    Example:
    data_dir=view_same_motif_not_same_cnb
    sites_file=same_motif_not_same_cnb_10_loci.bed
    win_work_dir=D:/Desktop/Learning/GraduateStudent/vntr_analysis/Quartet_data/compare_twins

    # 生成路径
    mkdir -p ${data_dir}
    mkdir -p ${data_dir}/sub_bams
    mkdir -p ${data_dir}/igv_scripts
    mkdir -p ${data_dir}/igv_snapshots

    cp $sites_file ${data_dir}/sites.txt

    # make bam_paths.txt and sites.txt
    ls ~/Quartet_data/data/soft_linked_bam/LCL[5-6]_HiFi.bam > ${data_dir}/bam_paths.txt

    # 取出bam
    fetch_alignments.sh ${data_dir}/bam_paths.txt ${data_dir}/sites.txt ${data_dir}/sub_bams/
    # 生成igv脚本
    gen_igv_scripts.sh ${data_dir}/bam_paths.txt ${data_dir}/sites.txt ${data_dir}/igv_scripts ${win_work_dir}/${data_dir}/sub_bams ${win_work_dir}/${data_dir}/igv_snapshots 
    if [[ -f ${data_dir}/igv_scripts/merged.txt ]]; then
        rm -f ${data_dir}/igv_scripts/merged.txt
    fi
    cat ${data_dir}/igv_scripts/*txt > ${data_dir}/igv_scripts/merged.txt
    '
    exit 1
fi

file_paths=$1
sites_file=$2
out_dir=$3

num_sites=$(wc -l $sites_file | awk '{print $1}')
num=0

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
    cat $file_paths | while read file
    do
        filename=${file##*/}
        basename=${filename%.*}
        result_file=${basename}_${site_a[0]}_${site_a[1]}_${site_a[2]}.bam
        
        echo "result file: "$out_dir/$result_file
        # cmd="samtools view -b $file $f_site | samtools sort -n | samtools fastq > $out_dir/$result_file"
        $samtools view -b $file $f_site > $out_dir/$result_file
        $samtools index $out_dir/$result_file
    done
    num=$((num+1))
    echo "done for site $num / $num_sites"
done


