#!/usr/bin/bash

if [[ $# -eq 0 ]]; then
    echo "Usage: extract_coverage.sh mosdepth_per_sample_dir success_samples_ids"
    exit 1
fi

function extract_mosdepth_coverage(){
    mosdepth_summary=$1
    mean_cov=$(grep -w "total" $mosdepth_summary | cut -f4 )
    echo $mean_cov
}

out_dir=$1
success_samples_ids=$2

num=0
echo -e "sample_id\tmean_depth"
cut -f1 $success_samples_ids | while read sid
do
    mosdepth_summary=${out_dir}/${sid}/${sid}.mosdepth.summary.txt
    if [[ -f $mosdepth_summary ]]; then
        cov=$(extract_mosdepth_coverage $mosdepth_summary)
        echo -e "$sid\t$cov"
    else
        echo "Error! File $mosdepth_summary do not exists"
        exit 1
    fi
done
