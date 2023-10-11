#!/usr/bin/bash
set -ue
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
source /public/home/fan_lab/wangjie/utils/workflow_scripts/job_manager_scripts/submit_job_utils.sh

work_dir=$1
call_vntr_sample_batch_ids_file=$2
display_batch_size=$3
run_straglr_script=$4
ref_tr_prefix=$5

bam_dir=$6
bam_suffix=$7
vntr_out_dir=$8 # per_sample_dir
log_dir=$9
running_samples=${10}
has_run_samples=${11}

applied_memory=${12}
max_str_len=${13}
min_str_len=${14}
id_sep=${15}

# echo $log_dir
# echo $running_samples
# echo $has_run_samples
# echo $applied_memory
# echo $max_str_len" "$min_str_len
# exit 

# ---------------------------BEGIN------------------------------ #
myprint "work_dir is"$work_dir
echo "sample info is $call_vntr_sample_batch_ids_file script is $run_straglr_script"

# 寻找还未跑过程序的样本
curr_time=`date "+%Y-%m%d-%H%M"`
cut -f1 $call_vntr_sample_batch_ids_file | grep -v "sample_id" > /tmp/tmp_sample_ids_${curr_time}.txt
cut -f1 $has_run_samples > /tmp/tmp_has_run_${curr_time}.txt
not_run_samples=$(setdiff /tmp/tmp_sample_ids_${curr_time}.txt /tmp/tmp_has_run_${curr_time}.txt)
not_run_samples_array=($not_run_samples)

rm /tmp/tmp_has_run_${curr_time}.txt /tmp/tmp_sample_ids_${curr_time}.txt

# 告知用户即将对下列样本跑程序
num_not_run_samples=${#not_run_samples_array[@]}
num=0
if [[ "$num_not_run_samples" != 0 ]]; then
    echo "Will run for these samples: "
    for sid in ${not_run_samples_array[*]}
    do
        num=$((num+1))
        sample_id=$(echo $sid | cut -d ${id_sep} -f1)
        input_bam=${bam_dir}/${sample_id}${bam_suffix}
        if [[ -f $input_bam ]]; then
            bam_size=$(ls -lh $input_bam | cut -d' ' -f5)
            echo -e "$sid\t$num\t$bam_size\t$call_vntr_sample_batch_ids_file"
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
for sid in ${not_run_samples_array[*]}
do
    # sid是sample_batch_id
    
    num=$((num+1))
    if [[ $num -ge $bstart && $num -le $bend ]]; then # 只提交指定的样本
        # 获取sample_id与tr batch id
        sample_id=$(echo $sid | cut -d ${id_sep} -f1)
        batch=$(echo $sid | cut -d ${id_sep} -f2)

        # inputs and outputs
        input_bam=${bam_dir}/${sample_id}${bam_suffix}
        sample_out_dir=${work_dir}/${vntr_out_dir}/${sample_id}
        raw_dir=${sample_out_dir}/raw
        if [[ ! -f $input_bam ]]; then
            echo "$input_bam do not exist! Error!"
            exit 1
        fi

        dir_exists $sample_out_dir
        dir_exists $raw_dir

        # tr input
        tr_bed=${ref_tr_prefix}${batch}.tsv

        # 需要指定raw和clean的out_prefix
        out_prefix=${raw_dir}/${sid}

        # submit jobs
        log="$log_dir/${sid}_"`date "+%Y-%m%d-%H%M"`".log"
        args=$(echo work_dir=$sample_out_dir,tr_bed=$tr_bed, \
                bam=$input_bam,out_prefix=$out_prefix,threads=$threads, \
                max_str_len=$max_str_len,min_str_len=$min_str_len | tr -d ' ' )
        if [[ $dry_run == "y" ]]; then
            echo qsub -N ${sid} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $run_straglr_script
        else
            qsub -N ${sid} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $run_straglr_script
        fi

        # ----------------------------END------------------------------- #
        if [[ $dry_run == "y" ]]; then
            echo "dry_run mode is on. won't write sample_id into running_samples and has_run_samples"
        else
            echo "$num: run for $sid, this sample has been added to $has_run_samples and $running_samples"
            echo -e "$sid\t$call_vntr_sample_batch_ids_file\t$job_queue\t$threads\t`date "+%Y-%m%d-%H%M"`" >> $has_run_samples
            echo $sid >> $running_samples

            # sleep 1
        fi
    fi

done

echo -e "Successfully run program for those samples. \nBye. "
