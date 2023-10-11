#!/usr/bin/bash
set -ue

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

# arguments by user
work_dir=$1
sample_id_file=$2
straglr_per_sample_dir=$3 # 
log_dir=$4
job_queue=$5

mkdir -p $log_dir

# fixed arguments
script=/public/home/fan_lab/wangjie/utils/workflow_scripts/call_vntr_scripts/scripts/compress_raw_vntr_results.pbs
threads=1

num=0
for sample_id in $(cat $sample_id_file)
do
    num=$((num+1))
    sample_out_dir=${straglr_per_sample_dir}/${sample_id}
    raw_bed=${sample_out_dir}/raw.bed
    raw_tsv=${sample_out_dir}/raw.tsv
    if [[ -f $raw_bed ]]; then
        raw_bed_size=$(ls -lh $raw_bed | cut -d' ' -f5)
        raw_tsv_size=$(ls -lh $raw_tsv | cut -d' ' -f5)
        echo -e "# $sample_id\t$num\t$raw_bed_size\t$raw_tsv_size"

        log="$log_dir/${sample_id}_"`date "+%Y-%m%d-%H%M"`".log"
        args=$(echo work_dir=$work_dir,straglr_results_dir=$sample_out_dir | tr -d ' ' )
        echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -o $log -v $args $script
    else
        echo "# $raw_bed do not exist! Will skip it"
    fi
done



