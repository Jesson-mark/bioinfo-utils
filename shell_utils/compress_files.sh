#!/usr/bin/bash 

set -ue

if [ $# -eq 0 ]; then
    echo "Usage: compress_files.sh file1 file2 ... "
    echo "Function: Compress files by gzip"
    echo "Example1: compress_files.sh logs/A.log logs/B.log"
    echo 'Example2: compress_files.sh logs/AAA*'
    exit 0
fi

for file in "$@"; 
do
    echo "compressing $file"
    gzip $file

done




