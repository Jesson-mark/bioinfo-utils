#!/usr/bin/bash
set -ue
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

function set_env(){
    call_vntr_config_file=$1
    call_vntr_fixed_config_file=$2
    eval $(parse_yaml $call_vntr_config_file "p_")

    cp $call_vntr_config_file $call_vntr_fixed_config_file

    work_dir=$p_work_dir
    raw_sample_id_file=$p_raw_sample_id_file
    project_dir=$p_project_dir
    call_vntr_shorter_job_type=$p_call_vntr_shorter_job_type
    call_vntr_job_type=$p_call_vntr_job_type
    clean_vntr_shorter_job_type=$p_clean_vntr_shorter_job_type
    clean_vntr_job_type=$p_clean_vntr_job_type
    submit_straglr_jobs_script=$p_submit_straglr_jobs_script
    submit_clean_straglr_jobs_script=$p_submit_clean_straglr_jobs_script
    run_straglr_script=$p_run_straglr_script
    tr_batch_id=$p_tr_batch_id
    ref_tr_prefix=$p_ref_tr_prefix
    id_sep=$p_id_sep

    # def dir and file paths
    sample_id_dir=${project_dir}/sample_ids
    vntr_out_dir=${project_dir}/per_sample_dir
    sample_id_file=${sample_id_dir}/sample_ids.txt

    # create dirs and sample id file
    mkdir -p ${project_dir} ${sample_id_dir} ${vntr_out_dir}
    cp $raw_sample_id_file $sample_id_file # 拷贝样本ID文件至项目目录

    # 生成 sample-batch-ids
    call_vntr_sample_batch_ids_file=${sample_id_dir}/call_vntr_sample_batch_ids.txt
    create_tr_sample_ids.sh $sample_id_file $tr_batch_id $call_vntr_sample_batch_ids_file $id_sep

    # 添加任务
    log_dir=${project_dir}/logs # 在project_dir里面创建
    mkdir -p $log_dir
    call_vntr_log_dir=${log_dir}/${call_vntr_job_type}/
    clean_vntr_log_dir=${log_dir}/${clean_vntr_job_type}/

    res=$(grep $call_vntr_job_type job_config.tsv | cat)
    if [[ $res != "" ]]; then
        echo "Error! $call_vntr_job_type already in job_config.tsv! Please use a new job_type"
        exit 1
    fi
    res=$(grep $clean_vntr_job_type job_config.tsv | cat)
    if [[ $res != "" ]]; then
        echo "Error! $clean_vntr_job_type already in job_config.tsv! Please use a new job_type"
        exit 1
    fi

    add_new_tsk.sh $call_vntr_shorter_job_type $call_vntr_job_type $sample_id_dir $call_vntr_log_dir
    add_new_tsk.sh $clean_vntr_shorter_job_type $clean_vntr_job_type $sample_id_dir $clean_vntr_log_dir

    compress_vntr_log_dir=${log_dir}/compress_vntr_files
    echo "call_vntr_sample_batch_ids_file: $call_vntr_sample_batch_ids_file" >> $call_vntr_fixed_config_file
    echo "vntr_out_dir: $vntr_out_dir" >> $call_vntr_fixed_config_file
    echo "call_vntr_job_type: $call_vntr_job_type" >> $call_vntr_fixed_config_file
    echo "clean_vntr_job_type: $clean_vntr_job_type" >> $call_vntr_fixed_config_file
    echo "log_dir: $log_dir" >> $call_vntr_fixed_config_file
    echo "compress_vntr_log_dir: $compress_vntr_log_dir" >> $call_vntr_fixed_config_file

    echo "Final config file is $call_vntr_fixed_config_file"
}

function submit_call_vntr(){
    eval $(parse_yaml $call_vntr_fixed_config_file "r_")
    display_batch_size=$1

    has_run_samples=$(grep -w $r_call_vntr_job_type job_config.tsv | cut -f4)
    running_samples=$(grep -w $r_call_vntr_job_type job_config.tsv | cut -f5)
    call_vntr_log_dir=$(grep -w $r_call_vntr_job_type job_config.tsv | cut -f3)

    echo "to submit"
    echo $call_vntr_log_dir
    echo $r_applied_memory $r_max_str_len $r_min_str_len
    bash $r_submit_straglr_jobs_script $r_work_dir $r_call_vntr_sample_batch_ids_file \
            $display_batch_size $r_run_straglr_script $r_ref_tr_prefix \
            $r_bam_dir $r_bam_suffix $r_vntr_out_dir $call_vntr_log_dir $running_samples $has_run_samples \
            $r_applied_memory $r_max_str_len $r_min_str_len ${r_id_sep}
}

function check_combine_batch_vntr(){
    eval $(parse_yaml $call_vntr_fixed_config_file "r_")
    sample_id_dir=${r_project_dir}/sample_ids
    vntr_out_dir=${r_vntr_out_dir}

    batch_vntr_check_done=${sample_id_dir}/batch_vntr_check_done.txt # 存放所有batch已经检测结束的样本id
    batch_vntr_merge_done=${sample_id_dir}/batch_vntr_merge_done.txt # 存放所有batch已经merge结束的样本id
    call_vntr_success=${sample_id_dir}/${r_call_vntr_job_type}_success.txt
    call_vntr_success_sampleids=${sample_id_dir}/${r_call_vntr_job_type}_success_sample_batch_ids.txt
    cut -f1 $call_vntr_success | cut -d $r_id_sep -f1 | sort |uniq > $call_vntr_success_sampleids

    myprint "Begin to check batch vntr"
    check_batch_vntr_res.sh $batch_vntr_check_done $r_tr_batch_id $call_vntr_success_sampleids $vntr_out_dir $r_id_sep

    myprint "Begin to combine batch vntr"
    #2 将所有batch检测结束的样本的raw.bed整合到一起
    combine_batch_vntr_res.sh $batch_vntr_check_done $batch_vntr_merge_done $vntr_out_dir $r_id_sep
}

function submit_clean_vntr(){
    display_batch_size=$1
    applied_memory=$2

    # read config
    eval $(parse_yaml $call_vntr_fixed_config_file "r_")
    vntr_out_dir=${r_vntr_out_dir}
    sample_id_dir=${r_project_dir}/sample_ids
    batch_vntr_merge_done=${sample_id_dir}/batch_vntr_merge_done.txt # 存放所有batch已经merge结束的样本id

    # read sample id file
    has_run_samples=$(grep -w $r_clean_vntr_job_type job_config.tsv | cut -f4)
    running_samples=$(grep -w $r_clean_vntr_job_type job_config.tsv | cut -f5)
    log_dir=$(grep -w $r_clean_vntr_job_type job_config.tsv | cut -f3)

    bash $r_submit_clean_straglr_jobs_script $r_work_dir $batch_vntr_merge_done $display_batch_size \
        $log_dir $vntr_out_dir $running_samples $has_run_samples $r_clean_straglr_script $r_ref_tr_bed $applied_memory
}

function submit_compress_vntr(){
    job_queue=$1

    eval $(parse_yaml $call_vntr_fixed_config_file "r_")
    # clean_vntr_success_samples=$(grep -w $r_clean_vntr_job_type job_config.tsv | cut -f6)
    sample_id_dir=${r_project_dir}/sample_ids
    batch_vntr_merge_done=${sample_id_dir}/batch_vntr_merge_done.txt # 存放所有batch已经merge结束的样本id

    sample_id_file=${batch_vntr_merge_done}
    log_dir=${r_compress_vntr_log_dir}

    bash $r_submit_compress_straglr_jobs_script $r_work_dir $sample_id_file $r_vntr_out_dir \
        $log_dir $job_queue $r_compress_raw_vntr_results_script

}


if [[ $# -eq 0 ]]; then
    echo "Usage: call_vntr_wrkf.sh cmd [display_batch_size]"
    echo "    avaliable cmds: set_env, submit_call_vntr, check_cmbn, submit_clean_vntr, compress"
    echo "    1. call_vntr_wrkf.sh set_env"
    echo "    2. call_vntr_wrkf.sh submit_call_vntr display_batch_size"
    echo "    3. call_vntr_wrkf.sh check_cmbn"
    echo "    4. call_vntr_wrkf.sh submit_clean_vntr display_batch_size applied_memory"
    echo "    5. call_vntr_wrkf.sh submit_compress_vntr job_queue"
    exit 
fi

# arguments by user
cmd=$1 # set_env, submit_call_vntr, check_cmbn, submit_clean_vntr, compress
call_vntr_config_file=call_vntr_config.yml
call_vntr_fixed_config_file=call_vntr_config.fixed.yml

##### core codes
case $cmd in 
    'set_env')
        if [[ ! -f $call_vntr_config_file ]]; then
            echo "$call_vntr_config_file not exists!"
            exit 1
        fi

        set_env $call_vntr_config_file $call_vntr_fixed_config_file

        ;;
    'submit_call_vntr')
        display_batch_size=$2
        submit_call_vntr $display_batch_size

        ;;
    'check_cmbn')
        check_combine_batch_vntr 
        ;;
    
    'submit_clean_vntr')
        display_batch_size=$2
        applied_memory=$3
        submit_clean_vntr $display_batch_size $applied_memory
        ;;
    
    'submit_compress_vntr')
        job_queue=$2
        submit_compress_vntr $job_queue
        ;;
    *)
        echo "Wrong cmd $1 !"
        exit 1
esac

# 提交任务脚本


