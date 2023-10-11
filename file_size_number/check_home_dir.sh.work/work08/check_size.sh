#!/bin/sh
set -xve
hostname
cd /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work08
starttime=$(date +'%Y-%m-%d %H:%M:%S')

time du -sh /public/home/fan_lab/wangjie/gnomad_LD
touch /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work08/check_size.sh.done

endtime=$(date +'%Y-%m-%d %H:%M:%S')
start_seconds=$(date --date="$starttime" +%s)
end_seconds=$(date --date="$endtime" +%s)
echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\n", hours, left_mins);}'


