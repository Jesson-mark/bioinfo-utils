#!/usr/bin/bash

# set -ue

if [ $# -eq 0 ]; then
    echo "Usage: get_tsk_mem.sh tsk1_id tsk2_id ... "
    echo "Example: get_tsk_mem.sh 10573774 10563360"
    exit
fi

# echo ""
for tsk_id in "$@"; do
    # 判断id是否包含admin1
    if [[ $tsk_id =~ "admin1" ]]; then
        ntsk_id=${tsk_id%.admin1}
        tsk_id=$ntsk_id
    fi

    tsk_info=$(qstat -a | grep $tsk_id )
    node_info=$(pestat | grep $tsk_id | awk '{print $1,$2}')
    mem=$(qstat -f $tsk_id | grep 'resources_used.mem' | tr -d ' ' | sed 's/resources_used.mem=//' | sed 's/kb//' | awk '{print "mem="$1 / 1000"Mb"}')
    # jobname=$(qstat -f $tsk_id | grep Job_Name | sed 's/.*= //')
    echo -e "$tsk_info\t$node_info\t$mem"
done


