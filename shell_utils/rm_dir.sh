#!/usr/bin/bash
# usage:
# ./tar_dir.sh .
rootdir=$1

for dir in `ls $rootdir`
do
    if [ -d $dir ];then
        echo "removing $dir"
        rm -r $dir
    fi
done

