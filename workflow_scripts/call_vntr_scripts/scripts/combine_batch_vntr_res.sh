#!/usr/bin/bash

# 将30个batch的straglr结果合并到一起，并将原来的文件压缩起来

set -ue

source /public/home/fan_lab/wangjie/utils/shell_utils/utils.sh

if [[ $# -eq 0 ]]; then
	echo "Usage:"
    echo '    combine_batch_vntr_res.sh $batch_vntr_check_done $batch_vntr_merge_done $vntr_out_dir'
    exit 1
fi

batch_vntr_check_done=$1
batch_vntr_merge_done=$2
vntr_out_dir=$3
id_sep=$4

if [[ ! -f $batch_vntr_merge_done ]]; then
    touch $batch_vntr_merge_done 
fi

num=0
cat $batch_vntr_check_done | while read sample_id
do
    num=$((num+1))

    # 判断该样本是否已经combine并压缩
    res=$(grep $sample_id $batch_vntr_merge_done | cat)
    if [[ $res == "" ]]; then # 该样本的bed文件并没有被combine
        # dir path
        sample_out_dir=${vntr_out_dir}/${sample_id}
        raw_dir=${sample_out_dir}/raw

        # 合并raw bed文件
        head -1 ${raw_dir}/${sample_id}${id_sep}b00.bed > ${sample_out_dir}/raw_bed_header
        sed -n '2p' ${raw_dir}/${sample_id}${id_sep}b00.tsv > ${sample_out_dir}/raw_tsv_header
        cat ${raw_dir}/*bed | grep -v "#chrom" | cat ${sample_out_dir}/raw_bed_header - > ${sample_out_dir}/raw.bed
        cat ${raw_dir}/*tsv | awk '$0!~"#"' | cat ${sample_out_dir}/raw_tsv_header - > ${sample_out_dir}/raw.tsv

        echo "$num done for $sample_id "
        echo $sample_id >> $batch_vntr_merge_done
    else
        echo "$num already done for $sample_id "
    fi
done
echo "All done. Bye"
