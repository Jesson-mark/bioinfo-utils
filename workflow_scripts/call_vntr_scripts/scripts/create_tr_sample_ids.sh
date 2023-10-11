#!/usr/bin/bash

if [[ $# -eq 0 ]]; then
    echo "Usage: bash scripts/create_call_tr_sample_info.sh sample_id_file out_id_file"
    echo "sample_id_file: first column must be sample id, no header included"
    exit 1
fi

sample_ids=$1
tr_batch_id=$2 # VNTR的batch ID
out_id_file=$3
id_sep=$4 # 这个用来分割样本ID和batch ID，样本ID中不可以包含这个ID

# 判断样本ID是否包含 id_sep
res=$(grep ${id_sep} $sample_ids)
if [[ $res != "" ]]; then
    echo "Error! sample ID($sample_ids) cannot contain id_sep($id_sep)"
    exit 1
fi

cut -f1 $sample_ids | while read sid
do
    for bid in $(cat $tr_batch_id)
    do
        echo "${sid}${id_sep}${bid}"
    done
done > ${out_id_file}

echo "done, total sample batchs is $(wc -l ${out_id_file})"

