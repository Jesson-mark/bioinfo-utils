#!/usr/bin/bash

# 将被取消的样本从running和has_run文件中去除

set -ue

job_config=$1
selected_job_type=$2

if [ $# -eq 0 ]; then
    echo "Usage: process_failed_sample.sh job_config.tsv job_type sample_id1 ... "
    echo "Available job_type"
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit
fi

found=False
while read aline
do
    line=($aline)
    shorter_job_type=${line[0]}
    job_type=${line[1]}

    if [[ $shorter_job_type == $selected_job_type ]]; then
        found=True

        log_dir=${line[2]}

        has_run_samples_ids=${line[3]}
        running_samples_ids=${line[4]}
        success_samples_ids=${line[5]}

        successful_log_dir=${log_dir}/success/
        break
    fi
done < $job_config

if [[ $found == False ]]; then
    echo "Wrong jobtype $selected_job_type"
    echo "Available job_type"
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit 1
fi

num_skip=2
num=0
for sample_id in "$@"; do
    # 跳过前两个参数
    num=$((num+1))
    if [[ $num -le $num_skip ]]; then
        continue
    fi

    echo "sample_id is "$sample_id

    # 从running_samples_ids中去除
    res=$(grep -n $sample_id $running_samples_ids)
    if [ $? == 1 ]; then
        echo "This sample $sample_id not exists in $running_samples_ids " 
        exit
    fi
    line_in_running_idx=$(echo $res | cut -d ':' -f1)

    echo "deleting $sample_id from $running_samples_ids"
    sed -i "${line_in_running_idx}d" $running_samples_ids

    # 从has_run_samples_ids中去除
    res=$(grep -n $sample_id $has_run_samples_ids)
    if [ $? == 1 ]; then
        echo "This sample $sample_id not exists in $has_run_samples_ids " 
        exit 1
    fi
    line_in_running_idx=$(echo $res | cut -d ':' -f1)

    echo "deleting $sample_id from $has_run_samples_ids "
    sed -i "${line_in_running_idx}d" $has_run_samples_ids

    echo ""
done

echo "Success!"
