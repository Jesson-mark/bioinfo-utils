#!/usr/bin/bash

set -ue
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

if [ $# -eq 0 ]; then
    echo "Usage: tar_dir dir"
    echo 'tar all sub_dirs, use: ls | while read dir; do if [[ -d $dir ]]; then tar_dir $dir & fi; done'
    exit 
fi

dir=$1
adir=${dir%%/} # 将路径末尾的/去掉
myprint "Dir is $adir"
starttime=$(date +'%Y-%m-%d %H:%M:%S')

if [ -d "$dir" ]; then
    myprint "taring $adir to ${adir}.tar"
    tar -cf "${adir}".tar "${adir}"
    myprint "removing $adir"
    rm -rf "${adir}"
else
    myprint "$adir is not a directory!"
fi

# 计时结束
endtime=$(date +'%Y-%m-%d %H:%M:%S')
start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数
end_seconds=$(date --date="$endtime" +%s);

echo "$start_seconds" "$end_seconds" | awk '{total_sec=$2-$1; min=int(total_sec/60); lef_sec=total_sec%60; printf("Total time: %d minutes, %d seconds.\n", min, lef_sec); }'

