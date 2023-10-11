#!/usr/bin/bash

echo "Total: "
pestat | awk '{total_threads+=$5; used_threads+=$9}END{print "total: "total_threads", used: "used_threads}'

echo "default: "
pestat | awk '$1~/node.$/ || $1~/node..$/{total_threads+=$5; used_threads+=$9}END{print "total: "total_threads", used: "used_threads}'

echo "used by wangjie: $(qstat -a | grep wangjie | awk '$10=="R"{num+=$7}END{print num}')"
# echo "used by shfan: $(qstat -a | grep shfan | awk '$10=="R"{num+=$7}END{print num}')"
# echo "used by mppan: $(qstat -a | grep mppan | awk '$10=="R"{num+=$7}END{print num}')"

# 查看排队的任务        
num_queued_jobs=$(qstat -a default | awk '$10=="Q"' |wc -l)
echo "number of queued jobs on default: $num_queued_jobs"

num_queued_jobs=$(qstat -a fat | awk '$10=="Q"' |wc -l)
echo "number of queued jobs on fat: $num_queued_jobs"

num_queued_jobs=$(qstat -a fsh_team | awk '$10=="Q"' |wc -l)
echo "number of queued jobs on fsh_team: $num_queued_jobs"

num_queued_jobs=$(qstat -a mini | awk '$10=="Q"' |wc -l)
echo "number of queued jobs on mini: $num_queued_jobs"

echo "statistics of jobs queued on default: "
qstat -a default |awk '$10=="Q"' | awk '{print $2}' |sort |uniq -c

echo "Utils cmds:"
echo "  List IDs of all jobs:"
printf "    qme | awk '{printf("; echo -e '"%s ",$1)} END {print ""}\c'; printf "'\n"

echo "  List IDs of all queued jobs:"
printf "    qme | awk '%s{printf(" '$10=="Q"'; echo -e '"%s ",$1)} END {print ""}\c'; printf "'\n"

echo "  List IDs of all running jobs:"
printf "    qme | awk '%s{printf(" '$10=="R"'; echo -e '"%s ",$1)} END {print ""}\c'; printf "'\n"
