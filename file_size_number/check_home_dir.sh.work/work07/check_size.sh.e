hostname
+ hostname
cd /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work07
+ cd /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work07
starttime=$(date +'%Y-%m-%d %H:%M:%S')
++ date '+%Y-%m-%d %H:%M:%S'
+ starttime='2023-06-24 11:11:47'

time du -sh /public/home/fan_lab/wangjie/GIAB_trio
+ du -sh /public/home/fan_lab/wangjie/GIAB_trio

real	0m0.003s
user	0m0.000s
sys	0m0.003s
touch /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work07/check_size.sh.done
+ touch /public/home/fan_lab/wangjie/utils/file_size_number/check_home_dir.sh.work/work07/check_size.sh.done

endtime=$(date +'%Y-%m-%d %H:%M:%S')
++ date '+%Y-%m-%d %H:%M:%S'
+ endtime='2023-06-24 11:11:47'
start_seconds=$(date --date="$starttime" +%s)
++ date '--date=2023-06-24 11:11:47' +%s
+ start_seconds=1687576307
end_seconds=$(date --date="$endtime" +%s)
++ date '--date=2023-06-24 11:11:47' +%s
+ end_seconds=1687576307
echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\n", hours, left_mins);}'
+ echo 1687576307 1687576307
+ awk '{total_sec=$2-$1; mins=int(total_sec/60); hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\n", hours, left_mins);}'

