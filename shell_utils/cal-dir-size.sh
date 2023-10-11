#!/usr/bin/bash
set -ue

work_dir=/public/home/fan_lab/wangjie/utils/file_size_number
cd $work_dir

curr_time=`date "+%Y-%m%d-%H%M"`
stats_dir_info.sh > dir_num_file_${curr_time}.txt

echo "out file is $work_dir/dir_num_file_${curr_time}.txt"


