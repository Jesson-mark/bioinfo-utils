#!/usr/bin/bash

# function left_join() {}
# set -uex
if [[ $# -eq 0 ]]; then
    echo "Usage: left_join.sh query_file target_file [qry_idx, tgt_idx]"
    echo "  If qry_idx and tgt_idx not provided, this script will use the first column of both file to left join. "
    exit
fi

if [[ $# -eq 2 ]]; then
    qry_file=$1
    tgt_file=$2
    qry_col_idx=1
    tgt_col_idx=1

elif [[ $# -eq 4 ]]; then
    qry_file=$1
    tgt_file=$2
    qry_col_idx=$3
    tgt_col_idx=$4

else
    echo "Wrong number of arguments! Please see the usage carefully! "
    exit
fi

# temporary files
curr_time=`date "+%Y-%m%d-%H%M"`

# TODO 若输入文件不是tab分隔，得到的结果可能列数对不上，需要检查一下

awk -v OFS="\t" -F "\t" '
    NR==FNR{
        tgt_col_idx="'$tgt_col_idx'";
        arr[tgt_col_idx]=$0; 
        next 
    }
    NR>FNR{
        qry_col_idx="'$qry_col_idx'"
        if(qry_col_idx in arr){
            split(arr[qry_col_idx], res, "\t");
            
            printf("%s\t",$0); 
            for(i=1;i<=length(res);i++){
                if(i != tgt_col_idx){
                    if(i == NF){
                        printf("%s",res[i]); 
                    } else {
                        printf("%s\t",res[i]);
                    }
                }
            }
            print "";
        }
    }' $tgt_file $qry_file


