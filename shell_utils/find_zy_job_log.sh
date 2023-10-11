#!/usr/bin/bash

# set -ue

if [ $# -eq 0 ]; then
    echo "Usage: find_job_log.sh n sample_id1 sample_id2 ... "
    echo "Example: find_job_log.sh n 17HanZZ0750 17HanZZ0771 17HanZZ0882"
    exit
fi

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

job_type=$1
case $job_type in 
    "n")
        log_dir="logs/correct_ont_reads_by_nextdenovo"
        ;;
    
    "an")
        log_dir="logs/mm_align_nextdenovo_corrected_reads"
        ;;

    "tr")
        log_dir="logs/call_vntr/"
        ;;

    "tgtr")
        log_dir="logs/tg_call_vntr/"
        ;;

    "ct" | "clean_tr")
        log_dir="logs/clean_tr/"
        ;;
    *)
        echo "Wrong job type"
        exit
        ;;
esac

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


