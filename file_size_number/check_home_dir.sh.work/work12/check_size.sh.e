hostname
+ hostname
cd /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work12
+ cd /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work12
starttime=$(date +'%Y-%m-%d %H:%M:%S')
++ date '+%Y-%m-%d %H:%M:%S'
+ starttime='2023-06-24 11:11:50'

time du -sh /public/home/fan_lab/wangjie/miniconda3
+ du -sh /public/home/fan_lab/wangjie/miniconda3

real	8m47.759s
user	0m1.736s
sys	1m4.703s
touch /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work12/check_size.sh.done
+ touch /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work12/check_size.sh.done

endtime=$(date +'%Y-%m-%d %H:%M:%S')
++ date '+%Y-%m-%d %H:%M:%S'
+ endtime='2023-06-24 11:20:37'
start_seconds=$(date --date="$starttime" +%s)
++ date '--date=2023-06-24 11:11:50' +%s
+ start_seconds=1687576310
end_seconds=$(date --date="$endtime" +%s)
++ date '--date=2023-06-24 11:20:37' +%s
+ end_seconds=1687576837
echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\n", hours, left_mins);}'
+ echo 1687576310 1687576837
+ awk '{total_sec=$2-$1; mins=int(total_sec/60); hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\n", hours, left_mins);}'

