#!/usr/bin/bash
set -ue

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

work_dir=/public/home/fan_lab/wangjie/Quartet_data/data
myprint "work_dir is"$work_dir

sample_info=$1
display_batch_size=$2

# ---------------------------BEGIN------------------------------ #
# 脚本
script=scripts/run_mosdepth.pbs
job_type="cal_depth"
echo "sample info is $sample_info script is $script"

# dirs and files
bam_dir=ONT_downsample/
file_suffix=".bam"
out_dir=mosdepth/per_sample_dir
log_dir=logs/cal_depth
dir_exists $log_dir
dir_exists $out_dir

# PBS arguments
applied_memory=10gb

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
    for sid in ${not_run_samples_array[*]}
    do
        num=$((num+1))
        input_bam=${bam_dir}/${sid}${file_suffix}
        if [[ -f $input_bam ]]; then
            bam_size=$(ls -lh $input_bam | cut -d' ' -f5)
            echo -e "$sid\t$num\t$bam_size\t$sample_info"
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
read -p "Do you want to run program for those samples? (y/n): " proceed

if [[ $proceed != 'y' ]]; then
    echo "Bye."
    exit
fi

# 请用户选择提交任务数目、队列、线程
read -p "Please input batch_size(n or start-end), jobqueue(dft, fat, fsh, mini), threads, (merge_cns, y|n): " batch_size queue threads
if [[ $queue == "dft" ]]; then
    job_queue="default"
elif [[ $queue == "fat" ]]; then
    job_queue="fat"
elif [[ $queue == "fsh" ]]; then
    job_queue="fsh_team"
elif [[ $queue == "mini" ]]; then
    job_queue="mini"
else
    echo "Wrong jobqueue! Must be one of dft, fat, fsh!"
    exit 1
fi

# 获取待提交的样本顺序范围
if [[ $batch_size =~ "-" ]]; then # 判断用户是指定了范围还是指定了提交任务数目
    bstart=$(echo $batch_size | cut -d- -f1) # bstart: batch_start
    bend=$(echo $batch_size | cut -d- -f2)
else
    bstart=1
    bend=$batch_size
fi

if [[ $bend -gt $display_batch_size ]]; then
    echo "batch_size is out of range!"
    exit
fi

echo "will submit $bstart-$bend, job_queue is $job_queue, threads is $threads"

# 开始提交任务
# ---------------------------BEGIN------------------------------ #
num=0
for sample_id in ${not_run_samples_array[*]}
do
    num=$((num+1))
    if [[ $num -ge $bstart && $num -le $bend ]]; then # 只提交指定的样本
        # inputs and outputs
        input_bam=${bam_dir}/${sample_id}${file_suffix}

        sample_out_dir=${out_dir}/${sample_id}
        dir_exists $sample_out_dir
        out_prefix=${sample_out_dir}/${sample_id}

        # submit jobs
        log="$log_dir/${sample_id}_"`date "+%Y-%m%d-%H%M"`".log"
        args=$(echo work_dir=$work_dir,input_bam=$input_bam, \
                out_prefix=$out_prefix, \
                threads=$threads | tr -d ' ' )
        qsub -N ${sample_id}d -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script

        # ----------------------------END------------------------------- #

        echo "$num: run for $sample_id, this sample has been added to $has_run_samples and $running_samples"
        echo -e "$sample_id\t$sample_info\t$job_queue\t$threads\t`date "+%Y-%m%d-%H%M"`" >> $has_run_samples
        echo $sample_id >> $running_samples

        sleep 1
    fi

done

echo -e "Successfully run program for those samples. \nBye. "
