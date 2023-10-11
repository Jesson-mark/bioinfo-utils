#!/bin/bash

home_dir=/public/home/fan_lab/wangjie
group_share_dir=/public/group_share_data/fan_lab/wangjie
dwh_dir=/mnt/DWH/zone00/groups_data/fan_lab/wangjie/

echo -e "Begin to stats\n"
echo "home_dir is $home_dir"
echo "group_share_dir is $group_share_dir"
echo "dwh_dir is $dwh_dir"

starttime=`date +'%Y-%m-%d %H:%M:%S'`

# count 
size_of_home=$(du -sh $home_dir | awk '{print $1}')
num_files_home=$(find $home_dir | wc -l )
echo "size of home: $size_of_home"
echo "number of files of home: $num_files_home"

size_of_group_share=$(du -sh $group_share_dir | awk '{print $1}')
num_files_group_share=$(find $group_share_dir | wc -l )
echo "size of group_share: $size_of_group_share"
echo "number of files of group_share: $num_files_group_share"

size_of_dwh=$(du -sh $dwh_dir | awk '{print $1}')
echo "size of DWH: $size_of_dwh"

# done, print used time
endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数
end_seconds=$(date --date="$endtime" +%s);
echo -e "\nDone"
echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; min=int(total_sec/60); lef_sec=total_sec%60; printf("Total time: %d minutes, %d seconds.\n", min, lef_sec); }'

