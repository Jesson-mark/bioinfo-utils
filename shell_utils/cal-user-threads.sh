#!/usr/bin/bash

all_users=$(qstat -a | awk '$10=="R"{print $2}' | sort |uniq)
all_users_arr=($all_users)
for user in ${all_users_arr[*]}
do
    used_threads=$(qstat -a | grep $user | awk '$10=="R"{num+=$7}END{print num}')
    echo -e "$user\t$used_threads"
done > ~/used_threads

sort -k2,2rn ~/used_threads
rm ~/used_threads

