cd /public/home/fan_lab/wangjie/utils/file_size_number

while :
do
	curr_time=`date "+%Y-%m%d-%H%M"`
	out_file=file_storage/dir_num_file_${curr_time}.txt
	nohup stats_dir_info.sh > $out_file 2>&1 &
	echo "submitted job for $curr_time, out file is $out_file"
	sleep 1d
done

