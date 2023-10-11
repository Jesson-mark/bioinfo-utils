#!/usr/bin/bash 

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    exit
fi

########################## Parse arguments ##########################
project_dir=$1


scripts_dir=${project_dir}/scripts
logs_dir=${project_dir}/logs
workflow_script=${project_dir}/workflow.sh

mkdir -p "${scripts_dir}"
mkdir -p "${logs_dir}"
{
    echo '#!/usr/bin/bash '; echo ''
    echo 'set -ue'
    echo 'source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"'; echo ''
} > "$workflow_script"

echo "Generating done"

