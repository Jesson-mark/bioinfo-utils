#!/usr/bin/bash
source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

echo "number of files on home"
for dir in $(ls ~)
do
	count_files ~/$dir
done

echo "number of files on group_share home"
group_share_dir=/public/group_share_data/fan_lab/wangjie/

for dir in $(ls $group_share_dir)
do
	count_files $group_share_dir/$dir
done

