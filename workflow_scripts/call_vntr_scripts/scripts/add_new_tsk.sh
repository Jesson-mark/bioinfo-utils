if [[ $# -eq 0 ]]; then
    echo "bash scripts/add_new_tsk.sh $shorter_job_type $job_type $sample_id_dir $log_dir"
    exit 0
fi

set -ue
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

shorter_job_type=$1
job_type=$2
sample_id_dir=$3
log_dir=$4

if [[ ! -f job_config.tsv ]]; then
    echo -e "shorter_job_type\tjob_type\tlog_dir\thas_run_samples_ids\trunning_samples_ids\tsuccess_samples_ids" > job_config.tsv
fi

# generate has_run, running, success files
has_run_samples_ids="${sample_id_dir}/${job_type}_has_run.txt"
running_samples_ids="${sample_id_dir}/${job_type}_running.txt"
success_samples_ids="${sample_id_dir}/${job_type}_success.txt"

touch $success_samples_ids
touch $running_samples_ids
touch $has_run_samples_ids

dir_exists $log_dir
dir_exists ${log_dir}/success
dir_exists ${log_dir}/failed

echo -e "$shorter_job_type\t$job_type\t$log_dir\t$has_run_samples_ids\t$running_samples_ids\t$success_samples_ids" >> job_config.tsv
echo "Done for $job_type"
