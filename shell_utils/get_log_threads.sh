#!/usr/bin/bash

# set -ue

echo -e "logfile\tthreads"
for logfile in "$@"; do
    grep -q "Program is done" $logfile
    if [ $? -eq 0 ]; then
        threads=$(grep -v '+' $logfile |grep 'Threads' | tail -n1 | sed 's/.*Threads is //' | sed 's/ ]//')
        echo -e "$logfile\t$threads"
    else
        echo "Job of $logfile isn't completed"
    fi

done


