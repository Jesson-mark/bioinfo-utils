#!/usr/bin/bash

# set -ue

echo -e "logfile\tTime(in mins)\tTime(in hours)"
for logfile in "$@"; do
    grep -q "Program is done" $logfile
    if [ $? -eq 0 ]; then
        total_time_mins=$(grep -v '+' $logfile | grep "Total time" | sed 's/Total time: //' )
        total_time_hours=$(echo $total_time_mins | sed 's/ min.*//' | awk '{hours=$1/60; printf("%.1f hours\n", hours)}')
        echo -e "$logfile\t$total_time_mins\t$total_time_hours"
    else
        echo "Job of $logfile isn't completed"
    fi

done


