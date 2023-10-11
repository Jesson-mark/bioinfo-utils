#!/usr/bin/bash
# set -ux

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
source util_func.sh

work_dir=$(get_workdir)
if [[ $? -ne 0 ]]; then
    myprint "Wrong username"
    exit
fi

myprint "work_dir is"$work_dir

username=${USER}
# cd $work_dir

if [[ $# -eq 0 ]]; then
    echo "Usage: check_vntr_progress.sh p | sid [sample_id ... ]"
    exit 1
fi


function check_per_sample(){
    sample_id=$1
    log=$2
    is_done=$(log_is_done $log)

    # 判断log是否有报错信息
    grep -v ONT_error_correction $log | grep -q -i -E 'error|fail|walltime'
    if [[ $? -eq 0 ]]; then
        has_error="has_error"
    else
        has_error="no_error"
    fi
    # 判断跑每一步的程序数目
    num_con_trf_input=$(grep 'construct trf input' $log | wc -l )
    num_perform_trf=$(grep 'perform trf' $log | wc -l )
    num_find_sim_pat=$(grep 'find_similar_long_patterns_gt' $log | wc -l )
    num_group_alleles=$(grep 'group alleles by locus' $log | wc -l )
    num_run_gmm=$(grep 'genotype vntrs by GMM' $log | wc -l )

    # 输出
    if [[ $is_done == "done" ]]; then
        echo -e "$sample_id\t$has_error\t\t$num_con_trf_input\t\t  $num_perform_trf\t\t\t$num_find_sim_pat\t\t\t$num_group_alleles\t  $num_run_gmm\t$is_done"
    else
        running_time=$(qstat -a |grep $username | grep $sample_id | awk '{print $1,$3,$7,$NF}')
        echo -e "$sample_id\t$has_error\t\t$num_con_trf_input\t\t  $num_perform_trf\t\t\t$num_find_sim_pat\t\t\t$num_group_alleles\t  $num_run_gmm\t$is_done\t$running_time"
    fi
}

log_dir='logs/call_vntr'
echo -e "sample_id\thas_error\tconstruct_trf_input\tperform_trf\tfind_sim_patterns\tgroup_alleles\trun_GMM\tis_done\trunning_time"

if [[ $1 == "p" ]]; then # 未指定样本
    logs=($(ls $log_dir/*log 2>/dev/null ))
    if [[ $? -eq 0 ]]; then
        for log in ${logs[*]}
        do
            sample_id=$(echo $log | cut -d/ -f3 | cut -d_ -f1-2)
            check_per_sample $sample_id $log
        done
    else
        echo "no sample is running"
    fi
elif [[ $1 == "sid" ]]; then # 指定多个样本
    num_skip=1
    num=0

    # 从命令行接收样本id
    for sample_id in "$@"
    do
        num=$((num+1))
        if [[ $num -le $num_skip ]]; then
            continue
        fi
        
        # sample_id=$1
        log=$(ls ${log_dir}/${sample_id}*log)
        if [[ $? -ne 0 ]]; then
            echo "Error! log of $sample_id not exists"
            exit 1
        fi
        check_per_sample $sample_id $log
    done
fi



