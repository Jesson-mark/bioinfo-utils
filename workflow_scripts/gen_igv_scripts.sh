#!/usr/bin/bash
set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

# 功能: 以locus为单位，对多个bam在指定locus的reads比对情况截图，png文件命名格式: {locus_name}_{sample_names}.png
# 每个locus的igv脚本存在一个文件中

# sampleid_bampaths.txt: 第一列是sample id，第二列是要去查看的bam路径, 不包含header
# 每一个bam都要包含sites_file.txt里面指定的locus的reads
# bams_per_snapshot: 可以指定每个截图里面加载的bam文件数目
# igv_scripts_dir: linux dir
# bam_dir: windows dir
# snapshot_dir: windows dir

##### update history
# 2022-10-26: copied from /public/home/fan_lab/wangjie/others/yechao/call_vntr/scripts/gen_igv_scripts_v2.sh

if [[ $# -eq 0 ]]; then
    echo "usage:  "
	echo "  gen_igv_scripts.sh $sampleid_bampaths $bams_per_snapshot $sites_file $igv_scripts_dir $bam_dir $snapshot_dir"
    echo "  sampleid_bampaths.txt: 第一列是sample id，第二列是要去查看的bam路径, 不包含header, 每一个bam都要包含sites_file.txt里面指定的locus的reads"
    echo "  igv_scripts_dir: linux dir"
    echo "  bam_dir: windows dir"
    echo "  snapshot_dir: windows dir"
	echo "  Example: "
	echo "    create dirs"
	echo "    view_dir="
	echo '    mkdir ${view_dir}'
	echo '    mkdir ${view_dir}/sub_bams'
	echo '    mkdir ${view_dir}/igv_scripts'
	echo '    mkdir ${view_dir}/snapshots'
    exit 1
fi

file_paths=$1
bams_per_snapshot=$2
sites_file=$3
igv_scripts_dir=$4
bam_dir=$5 # win dir
snapshot_dir=$6 # win dir

genome=hg38

echo "Processing sites file: $sites_file"
echo "bam_dir: $bam_dir"
echo "snapshot_dir: $snapshot_dir"

cut -f1-3 $sites_file | while read site
do
    echo "Generating igv scripts for $site"
    read -r chr start end <<< $site
    echo "chr is $chr, start is $start, end is $end"
    
    locus="${chr}:$((${start}-20))-$((${end}+20))" # 左右各扩20bp
    chr_name=${chr}_${start}_${end}
    igv_scripts_file=${igv_scripts_dir}/${chr_name}.txt

    # generate igv script
    echo "new" > $igv_scripts_file
    echo "genome $genome" >> $igv_scripts_file

    total_num=$(wc -l $file_paths | awk '{print $1}')
    num=0
    sample_names=""
    cat $file_paths | while read line
    do
        l=($line)
        sample_id=${l[0]}
        if [[ $sample_id == 'sample_id' ]]; then
            continue
        fi
        
        bampath=${l[1]}
        filename=$(get_filename $bampath)
        bam_path=${bam_dir}/${filename}

        echo "load $bam_path" >> $igv_scripts_file
        num=$((num+1))
        sample_names=$sample_names"_"$sample_id
        # echo "$total_num $num "

        mod=$(awk 'BEGIN{print '"$num"' % '"$bams_per_snapshot"'}')
        if [[ $num -eq $total_num ]]; then
            load_done=True
            last_snapshot=True
        elif [[ $mod -eq 0 ]]; then  
            load_done=True
            last_snapshot=False
        else
            load_done=False
            last_snapshot=False
        fi

        if [[ $load_done == True ]]; then
            snapshot_filename=${chr_name}${sample_names}.png
            
            echo "snapshotDirectory $snapshot_dir" >> $igv_scripts_file
            echo "goto $locus" >> $igv_scripts_file
            echo "region $chr $start $end" >> $igv_scripts_file
            echo "sort base" >> $igv_scripts_file
            echo "snapshot $snapshot_filename" >> $igv_scripts_file
            echo "" >> $igv_scripts_file

            sample_names=""

            if [[ $last_snapshot == False ]]; then
                echo "new" >> $igv_scripts_file
                echo "genome $genome" >> $igv_scripts_file
            fi
        fi
    done
    echo "done, igv_scripts_file is $igv_scripts_file"
done



