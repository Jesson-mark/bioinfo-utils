#!/usr/bin/bash 

if [ $# -eq 0 ]; then
    echo "Usage: check-files.sh file1 file2 ... "
    echo "Example1: check-files.sh a.txt b.txt c_dir"
    echo 'Example2: check-files.sh dir/AAA*'
    exit 0
fi

num_not_exists=0
num_exists=0

for afile in "$@"; 
do
    if [[ -f $afile || -d $afile ]]; then
        num_exists=$((num_exists+1))
    else
        num_not_exists=$((num_not_exists+1))
        echo "File or dir '$afile' not exists!"
    fi
done
