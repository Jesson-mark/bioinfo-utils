#!/usr/bin/bash
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh
set -ue

# work_dir=/public/home/fan_lab/wangjie/CoreQueue/vntr

if [ $# -eq 0 ]; then
    echo "Usage: ./run.sh command "
    exit
fi

case $1 in
	"cmd")
		 ;;
	*)
		echo "Wrong"
esac
