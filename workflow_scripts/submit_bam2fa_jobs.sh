#!/usr/bin/bash
set -ue

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
source /public/home/fan_lab/wangjie/utils/workflow_scripts/job_manager_scripts/submit_job_utils.sh

work_dir=$(bash ../utils/get_workdir.sh)"/data"
myprint "work_dir is"$work_dir

sample_info=$1
display_batch_size=$2

# ---------------------------BEGIN------------------------------ #
bam_dir=/public/group_share_data/fan_lab/shali/90minority_30X/bam/1_results
file_suffix=".bam"
script=scripts/bam2fasta.pbs
out_dir=${work_dir}/ont_fasta
log_dir=logs/bam2fa
dir_exists $log_dir

applied_memory=30gb
job_type="bam2fa"
echo "sample info is $sample_info script is $script"

# sample ids
sample_id_dir="sample_ids"
running_samples="${sample_id_dir}/${job_type}_running.txt"
has_run_samples="${sample_id_dir}/${job_type}_has_run.txt"

# ----------------------------END------------------------------- #

# 寻找还未跑过程序的样本
curr_time=`date "+%Y-%m%d-%H%M"`
cut -f1 $sample_info | grep -v "sample_id" > /tmp/tmp_sample_ids_${curr_time}.txt
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
        input_bam=${bam_dir}/${sample_id}${file_suffix}
        if [[ -f $input_bam ]]; then
            bam_size=$(ls -lh $input_bam | cut -d' ' -f5)
            echo -e "$sample_id\t$num\t$bam_size\t$sample_info"
            if [[ $num == $display_batch_size ]]; then
                break
            fi
        else
            echo "$input_bam do not exist! Will skip it"
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
        bam=${bam_dir}/${sample_id}${file_suffix}
        out_fasta=${out_dir}/${sample_id}.fa.gz
    
        # submit jobs
        log="$log_dir/${sample_id}_"`date "+%Y-%m%d-%H%M"`".log"
        args=$(echo work_dir=$work_dir, \
            bam=$bam,out_fasta=$out_fasta, \
            threads=$threads | tr -d ' ' )

        if [[ $dry_run == "y" ]]; then
            echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script
        else
            qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script
        fi

        # ----------------------------END------------------------------- #
        if [[ $dry_run == "y" ]]; then
            echo "dry_run mode is on. won't write sample_id into running_samples and has_run_samples"
            # echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -o $log -v $args $script
        else
            echo "$num: run for $sample_id, this sample has been added to $has_run_samples and $running_samples"
            echo -e "$sample_id\t$sample_info\t$job_queue\t$threads\t`date "+%Y-%m%d-%H%M"`" >> $has_run_samples
            echo $sample_id >> $running_samples

            # sleep 1
        fi
    fi
done

echo -e "Successfully run program for those samples. \nBye. "
