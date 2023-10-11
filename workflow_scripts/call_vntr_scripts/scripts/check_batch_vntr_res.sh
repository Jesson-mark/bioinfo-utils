#!/usr/bin/bash

# 目的
# 判断哪些样本跑完了所有的batch，并将该样本记录下来
# 合并该样本所有batch的bed文件

if [[ $# -eq 0 ]]; then
    echo "Usage:"
	echo '    check_batch_vntr_res.sh $batch_vntr_check_done $tr_batch_id $call_vntr_success_sampleids $vntr_out_dir'
    exit 1
fi

batch_vntr_check_done=$1
tr_batch_id=$2 # VNTR的batch ID
call_vntr_success_sampleids=$3
vntr_out_dir=$4
id_sep=$5

# 检查哪些样本的全部batch跑完了
if [[ ! -f $batch_vntr_check_done ]]; then
    touch $batch_vntr_check_done
fi

num=0
num_batches=$(wc -l $tr_batch_id | cut -d' ' -f1)
cat $call_vntr_success_sampleids | while read sample_id
do
    res=$(grep $sample_id $batch_vntr_check_done | cat)
    num=$((num+1))
    if [[ $res == "" ]]; then # 该样本并没有全部跑完

        sample_out_dir=${vntr_out_dir}/${sample_id}
        raw_dir=${sample_out_dir}/raw
        
        num_done_batch=0
        for bid in $(cat $tr_batch_id)
        do
            raw_bed=${raw_dir}/${sample_id}${id_sep}${bid}.bed
            if [[ -f $raw_bed ]]; then # 该batch已完成
                num_done_batch=$((num_done_batch+1))
            fi
        done

        echo -e "$sample_id\t$num_done_batch/${num_batches}"
        if [[ $num_done_batch == ${num_batches} ]]; then
            echo $sample_id >> $batch_vntr_check_done
        fi
    else
        echo "$num all batches are done for $sample_id "
    fi
done


