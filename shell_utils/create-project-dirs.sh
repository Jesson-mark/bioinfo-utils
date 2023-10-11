#!/usr/bin/bash 

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"
if [ $# -eq 0 ]; then
    echo "Usage: "
    exit
fi

project_dir=$1
inputs_dir=${project_dir}/inputs 
outputs_dir=${project_dir}/outputs 
analyses_dir=${project_dir}/analyses
logs_dir=${project_dir}/logs

echo "Project_dir is ${project_dir}"
echo "Creating following dirs"
echo "${inputs_dir}" "${outputs_dir}" "${analyses_dir}" "${logs_dir}"

mkdir -p "${inputs_dir}" "${outputs_dir}" "${analyses_dir}" "${logs_dir}"


