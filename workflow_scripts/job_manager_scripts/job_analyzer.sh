#!/usr/bin/bash
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

function process_per_sample(){
    job_type=$1
    sample_id=$2
    log_dir=$3
    has_run_samples=$4

    job_time=$(get_log_time $sample_id $log_dir )
    job_info=$(grep $sample_id $has_run_samples | cut -f 3-5  )
    ji=($job_info) # jobqueue, threads, submit_time

    echo -e "$sample_id\t$job_time\t${ji[0]}\t${ji[1]}\t${ji[2]}"
}

job_config=job_config.tsv # 该文件需要存在于执行本脚本的目录下面

if [ $# -eq 0 ]; then
    echo -e "Usage: ./scripts/job_analyzer.sh job_type [info | sid] [sample_info | sample_id1, ... ] "
    echo -e "Available Commands: "
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit 1
fi

selected_job_type=$1
sample_info_or_sid=$2

echo "job_type is $selected_job_type, selection is $sample_info_or_sid" >&2 

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

# 输出header
echo -e "sample_id\tjob_time\tqueue\tthreads\tsubmit_time"

# 输出每个样本的信息
case $sample_info_or_sid in 
    "info")
        sample_info=$3

        cat $sample_info | while read line
        do
            l=($line)
            sample_id=${l[0]}
            process_per_sample $job_type $sample_id $successful_log_dir $has_run_samples_ids
        done
        ;;

    "sid")
        num_skip=2
        num=0

        for sample_id in "$@"
        do
            num=$((num+1))
            if [[ $num -le $num_skip ]]; then
                continue
            fi
            process_per_sample $job_type $sample_id $successful_log_dir $has_run_samples_ids
            
        done

        ;;

    *)
        echo "Wrong sample info or sid"
        exit 1
esac



