#!/usr/bin/bash

# set -ue

function find_log_file(){
    sample_id=$1
    log_dir=$2

    local curr_time=`date "+%Y-%m%d-%H%M"`
    tmp_file=/tmp/tmp_ls_err_${curr_time}
    log=$(ls $log_dir/*$sample_id* 2> $tmp_file)
    if [[ $? -eq 0 ]]; then
        echo $log
        rm $tmp_file
    else
        echo "Not found log of $sample_id in $log_dir"
        rm tmp_file
        return 1
    fi
}

job_config=job_config.tsv

if [ $# -eq 0 ]; then
    echo "Usage: find_job_log.sh job_type sample_id1 sample_id2 ... "
    echo -e "Available job types: "
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    echo "Example: find_job_log.sh n 17HanZZ0750 17HanZZ0771 17HanZZ0882"
    exit
fi

selected_job_type=$1
found=False
while read aline
do
    line=($aline)
    shorter_job_type=${line[0]}
    job_type=${line[1]}

    if [[ $shorter_job_type == $selected_job_type ]]; then
        found=True

        log_dir=${line[2]}
        break
    fi
done < $job_config

if [[ $found == False ]]; then
    echo "Wrong jobtype $selected_job_type"
    echo "Available job_type"
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit 1
fi

# echo ""
num_skip=1
num=0
for sample_id in "$@"
do
    num=$((num+1))
    if [[ $num -le $num_skip ]]; then
        continue
    fi

    find_log_file $sample_id $log_dir

done


