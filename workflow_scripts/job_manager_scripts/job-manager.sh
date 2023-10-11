#!/usr/bin/bash
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
source /public/home/fan_lab/wangjie/utils/shell_utils/envs.sh

function log_is_done(){
    log=$1
    is_done=$(check-logs.sh $log | head -1 | cut -f2)
    if [[ $is_done == "True" ]]; then
        echo "done"
    else
        echo "failed"
    fi
}

function check_jobs_from_id(){
    running_ids=$1
    log_dir=$2
    
    # 临时文件
    local curr_time=`date "+%Y-%m%d-%H%M-%N"`
    local tmp_ls_err=/tmp/tmp_ls_error_${curr_time}
    touch $tmp_ls_err
    local tmp_job_status=/tmp/tmp_job_status_${curr_time}

    # jobids=$(qstat -a |grep ${USER} |awk '{print $1}')
    # get-tsk-info.sh $jobids > $tmp_job_status
    # qstat -a | grep ${USER} > $tmp_job_status
    qmefull > $tmp_job_status
    
    num=0
    cat $running_ids | while read line
    do
        l=($line)
        sample_id=${l[0]}
        num=$((num+1))
        
        # 判断任务状态
        res=$(grep -w $sample_id $tmp_job_status)

        if [[ $? -eq 0 ]]; then
            # 判断是running, queued, canceled
            # 若是canceled，那么判断是successful还是failed
            info=($res)
            job_status=${info[9]}
            if [[ ${job_status} == "R" ]]; then
                status="running"
                # fa=${reads_dir}/${sample_id}.fa.gz
                # file_size=$(ls -lh $fa | awk '{print $5}')
                echo -e "$num\t$sample_id\t$status\t${info[0]}\t${info[2]}\t${info[6]}\t${info[10]}"
                
            elif [[ ${job_status} == "Q" ]]; then
                status="queued"
                echo -e "$num\t$sample_id\t$status\t${info[0]}\t${info[2]}\t${info[6]}\t${info[10]}"
                # echo -e "$num\t$sample_id\t$status"

            elif [[ ${job_status} == "C" ]]; then
                log=$(ls $log_dir/*$sample_id* 2>$tmp_ls_err) 
                if [[ $? -eq 0 ]]; then
                    status=$(log_is_done $log)
                    if [[ $status == "done" ]]; then
                        total_time=$(get_log_time $sample_id $log_dir )
                        # total_time=$(awk '$1!~"+"' $log | grep "Total time" | sed 's/Total time: //' )
                        echo -e "$num\t$sample_id\tdone\tTotal time: $total_time"
                    else
                        echo -e "$num\t$sample_id\tfailed"
                    fi
                else
                    status="canceled"
                    echo -e "$num\t$sample_id\t$status"
                fi

            else
                echo "Wrong status $res"
                exit
            fi

        else # 判断是successful还是failed
            log=$(ls $log_dir/*$sample_id* 2>$tmp_ls_err) 
            if [[ $? -eq 0 ]]; then
                status=$(log_is_done $log)
                if [[ $status == "done" ]]; then
                    total_time=$(get_log_time $sample_id $log_dir )
                    # total_time=$(awk '$1!~"+"' $log | sed 's/Total time: //' )
                    echo -e "$num\t$sample_id\tdone\tTotal time: $total_time"
                else
                    echo -e "$num\t$sample_id\tfailed"
                fi

            else
                # TODO: 如果没有搜到，那么说明这个任务被取消了，没有生成log文件
                # 因此这里需要打印canceled
                status="canceled"
                echo -e "$num\t$sample_id\t$status"

            fi
        fi
    done

    rm $tmp_job_status $tmp_ls_err
}

function check_different_jobs() {
    job_type=$1
    log_dir=$2
    running_ids=$3
    
    echo -e "sample_num\tsample_id\tjob_status\tJob_ID\tQueue\tTSK\tElapTime"

    # 临时文件
    curr_time_=`date "+%Y-%m%d-%H%M"`
    tmp_job_status=/tmp/job_status_tmptmptmp_${curr_time_}
    check_jobs_from_id $running_ids $log_dir > $tmp_job_status
    cat $tmp_job_status
	#cp $tmp_job_status tmp_job_status

    # 日志信息
    log_num=$(ls -l $log_dir/*log |wc -l)
    echo "running sample num: "$(wc -l $running_ids)
    echo "number of log files: $log_num"

    echo "Summary"
    cat $tmp_job_status | cut -f3,5 |sort |uniq -c

    echo "Successful jobs: "
    successful_ids=$(grep "done" $tmp_job_status | awk '{printf("%s ",$2)}')
    # successful_ids=$(check_logs.sh $log_expression | grep True | awk '{split($1,res,"/"); split(res[3], res2, "_"); printf("%s ",res2[1])}END{print ""}')
    echo "process_done_sample.sh job_config.tsv $job_type $successful_ids"

    echo "Failed jobs: " 
    canceled_job_ids=$(grep 'failed' $tmp_job_status | awk '{printf("%s ", $2)}END{print ""}')
    if [[ $canceled_job_ids != "" ]]; then
        echo "process_failed_sample.sh job_config.tsv $job_type $canceled_job_ids"
        echo "find_job_log.sh $job_type $canceled_job_ids"
        echo "$log_dir/failed"
    fi
    echo "Canceled jobs: " 
    canceled_job_ids=$(grep 'canceled' $tmp_job_status | awk '{printf("%s ", $2)}END{print ""}')
    if [[ $canceled_job_ids != "" ]]; then
        echo "process_failed_sample.sh job_config.tsv $job_type $canceled_job_ids"
    fi

    rm $tmp_job_status 
}

job_config=job_config.tsv

if [ $# -eq 0 ]; then
    echo -e "Usage: \n  ./job_manager.sh [command] \n"
    echo -e "Available Commands: "
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit 1
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

        has_run_samples_ids=${line[3]}
        running_samples_ids=${line[4]}
        success_samples_ids=${line[5]}

        check_different_jobs $shorter_job_type $log_dir $running_samples_ids
    fi
done < $job_config

if [[ $found == False ]]; then
    echo "Wrong jobtype $selected_job_type"
    echo "Available job_type"
    sed '1d' $job_config |cut -f1-2 | sed 's/\t/ (/' | awk '{print "  "$0")"}'
    exit 1
fi
