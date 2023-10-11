#!/usr/bin/bash 

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    exit
fi

########################## Parse arguments ##########################
file1=$1
file2=$2


########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
myprint "Begin to analysis"

tmp_file1=${file1}.tmptmp
tmp_file2=${file2}.tmptmp

sort "$file1" > "$tmp_file1"
sort "$file2" > "$tmp_file2"

echo "Differences are"
diff "$tmp_file1" "$tmp_file2"

rm -f "$tmp_file1" "$tmp_file2"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
cal_used_time.sh "$starttime" "$endtime"
myprint "Program is done"
