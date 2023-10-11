#!/usr/bin/bash
set -ue

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
source /public/home/fan_lab/wangjie/utils/workflow_scripts/job_manager_scripts/submit_job_utils.sh

work_dir=$1
batch_vntr_merge_done=$2
display_batch_size=$3

log_dir=$4
vntr_out_dir=$5 # per_sample_dir
running_samples=$6
has_run_samples=$7

# almost fixed args
clean_straglr_script=$8
ref_tr_bed=$9 # 包含全部unique TR的文件

# args to pbs
applied_memory=${10}

# ---------------------------BEGIN------------------------------ #
dir_exists $log_dir
myprint "work_dir is"$work_dir
echo "sample info is $batch_vntr_merge_done script is $clean_straglr_script"

# 寻找还未跑过程序的样本
curr_time=`date "+%Y-%m%d-%H%M"`
cut -f1 $batch_vntr_merge_done | grep -v "sample_id" > /tmp/tmp_sample_ids_${curr_time}.txt
cut -f1 $has_run_samples > /tmp/tmp_has_run_${curr_time}.txt
not_run_samples=$(setdiff /tmp/tmp_sample_ids_${curr_time}.txt /tmp/tmp_has_run_${curr_time}.txt)
not_run_samples_array=($not_run_samples)

rm /tmp/tmp_has_run_${curr_time}.txt /tmp/tmp_sample_ids_${curr_time}.txt

# 告知用户即将对下列样本跑程序
num_not_run_samples=${#not_run_samples_array[@]}
num=0
if [[ "$num_not_run_samples" != 0 ]]; then
    echo "Will run for these samples: "
    for sample_id in ${not_run_samples_array[*]}
    do
        num=$((num+1))
        sample_out_dir=${vntr_out_dir}/${sample_id}
        raw_bed=${sample_out_dir}/raw.bed
        raw_tsv=${sample_out_dir}/raw.tsv
        if [[ -f $raw_bed ]]; then
            raw_bed_size=$(ls -lh $raw_bed | cut -d' ' -f5)
            raw_tsv_size=$(ls -lh $raw_tsv | cut -d' ' -f5)
            echo -e "$sample_id\t$num\t$raw_bed_size\t$raw_tsv_size\t$batch_vntr_merge_done"

            if [[ $num == $display_batch_size ]]; then
                break
            fi
        else
            echo "$raw_bed do not exist! Will skip it"
        fi
    done
    # exit
else
    echo "All samples have been done."
    exit
fi

# 询问是否提交任务
if_submit

# 让用户输入 batch_size queue threads dry_run 
arguments=$(read_arguments)
args=($arguments)

job_queue=${args[0]}
threads=${args[1]}
bstart=${args[2]}
bend=${args[3]}
dry_run=${args[4]}

echo "will submit $bstart-$bend, job_queue is $job_queue, threads is $threads"

# 开始提交任务
# ---------------------------BEGIN------------------------------ #
num=0
for sample_id in ${not_run_samples_array[*]}
do
    num=$((num+1))
    if [[ $num -ge $bstart && $num -le $bend ]]; then # 只提交指定的样本
        # inputs and outputs
        sample_out_dir=${vntr_out_dir}/${sample_id}
        raw_bed=${sample_out_dir}/raw.bed
        raw_tsv=${sample_out_dir}/raw.tsv
        out_prefix=${sample_out_dir}/clean

        # submit jobs
        log="$log_dir/${sample_id}_"`date "+%Y-%m%d-%H%M"`".log"
        args=$(echo work_dir=$work_dir,sample_id=$sample_id,raw_bed=$raw_bed,raw_tsv=$raw_tsv, \
                out_prefix=$out_prefix,tr_region=$ref_tr_bed | tr -d ' ' )
        if [[ $dry_run == "y" ]]; then
            echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $clean_straglr_script
        else
            qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $clean_straglr_script
        fi

        # ----------------------------END------------------------------- #
        if [[ $dry_run == "y" ]]; then
            echo "dry_run mode is on. won't write sample_id into running_samples and has_run_samples"
            # echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -o $log -v $args $script
        else
            echo "$num: run for $sample_id, this sample has been added to $has_run_samples and $running_samples"
            echo -e "$sample_id\t$batch_vntr_merge_done\t$job_queue\t$threads\t`date "+%Y-%m%d-%H%M"`" >> $has_run_samples
            echo $sample_id >> $running_samples

            #sleep 1
        fi
    fi

done

echo -e "Successfully run program for those samples. \nBye. "
