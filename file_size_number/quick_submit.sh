#!/usr/bin/bash

curr_time=`date "+%Y-%m%d"`
nohup stats_dir_info.sh > file_storage/dir_num_file_${curr_time}.txt 2>&1 &

