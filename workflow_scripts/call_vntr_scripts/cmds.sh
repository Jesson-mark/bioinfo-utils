############# 这里是使用straglr检测vntr的示例脚本
##### 制作示例数据
mkdir ref_vntr # 这里仅存放少量vntr位点，以保证用例的正常运行
mkdir bams # 仅存放HG002的部分bam
mkdir sample_ids

echo "HG002" > sample_ids/test_sample_ids.txt

# 制作ref vntr
straglr_simpleRepeat_bed=/public/home/fan_lab/wangjie/STR/ref_vntrs/straglr_input/straglr_simpleRepeat_unique.bed
head -60 $straglr_simpleRepeat_bed > ref_vntr/60_straglr_ref_vntrs.bed
split -l 12 -d --additional-suffix=.tsv ref_vntr/60_straglr_ref_vntrs.bed ref_vntr/splitted/ref_vntr_b # 切成5个batch
ls ref_vntr/splitted/ |cut -d_ -f3 | cut -d. -f1 > ref_vntr/batch_vntr_ids.txt # 获取vntr的batch id

# 获取示例bam
hg002_hifi_bam=/public/home/fan_lab/wangjie/GIAB_trio/GIAB_bams/bams/HG002_hifi.bam
cut -f1-3 ref_vntr/60_straglr_ref_vntrs.bed > ref_vntr/60_straglr_ref_vntrs.pos.bed
sambamba slice -L ref_vntr/60_straglr_ref_vntrs.pos.bed $hg002_hifi_bam -o bams/HG002_part_HiFi.bam
samtools sort bams/HG002_part_HiFi.bam -o bams/HG002_part_HiFi.sort.bam
samtools index bams/HG002_part_HiFi.sort.bam
rm -f bams/HG002_part_HiFi.bam

##### 根据样本ID制作sample id，每个样本后面加入要跑的batch vntr的ID
# 例如对HG002检测VNTR，共有3个batch(b00, b01, b02)
# 经过处理，新的ID为: HG002_b00, HG002_b01, HG002_b02
sample_id_file=sample_ids/test_sample_ids.txt
tr_batch_id=ref_vntr/batch_vntr_ids.txt
out_id_file=sample_ids/call_vntr_sample_batch_ids.txt
./scripts/create_tr_sample_ids.sh $sample_id_file $tr_batch_id $out_id_file

##### 生成任务config文件
# 检测vntr
shorter_job_type="vntr"
job_type="call_vntr"
sample_id_dir="sample_ids"
log_dir="logs/${job_type}/"
add_new_tsk.sh $shorter_job_type $job_type $sample_id_dir $log_dir

# 进行clean
shorter_job_type="cln-tr"
job_type="clean_vntr"
sample_id_dir="sample_ids"
log_dir="logs/${job_type}/"
add_new_tsk.sh $shorter_job_type $job_type $sample_id_dir $log_dir

##### 提交call vntr任务
# 提交任务
bash scripts/submit_straglr_jobs.sh sample_ids/call_vntr_sample_batch_ids.txt 10

##### 整合所有batch的结果并clean成vcf文件
#1 检查所有batch都检测结束的样本
batch_vntr_check_done=sample_ids/batch_vntr_check_done.txt # 存放所有batch已经检测结束的样本id
batch_vntr_merge_done=sample_ids/batch_vntr_merge_done.txt # 存放所有batch已经merge结束的样本id
call_vntr_success=sample_ids/call_vntr_success.txt
call_vntr_success_sampleids=sample_ids/call_vntr_success_sample_batch_ids.txt
cut -f1 $call_vntr_success | cut -d_ -f1 | sort |uniq > $call_vntr_success_sampleids
vntr_out_dir=vntr_results/per_sample_dir

tr_batch_id=ref_vntr/batch_vntr_ids.txt
check_batch_vntr_res.sh $batch_vntr_check_done $tr_batch_id $call_vntr_success_sampleids $vntr_out_dir

#2 将所有batch检测结束的样本的raw.bed整合到一起
combine_batch_vntr_res.sh $batch_vntr_check_done $batch_vntr_merge_done $vntr_out_dir

#3 对bed文件进行clean，得到vcf文件
bash scripts/submit_clean_straglr_jobs.sh sample_ids/batch_vntr_merge_done.txt 10

#4 将raw.bed和raw.tsv压缩起来，并删除raw目录
work_dir=/public/home/fan_lab/wangjie/utils/workflow_scripts/call_vntr_scripts
sample_id_file=sample_ids/batch_vntr_merge_done.txt
straglr_per_sample_dir=vntr_results/per_sample_dir # 
log_dir=logs/compress_vntr_files
job_queue=mini

submit_compress_straglr_jobs.sh $work_dir $sample_id_file $straglr_per_sample_dir $log_dir $job_queue > scripts/cmds_to_compress_vntr_files.sh


##### 流程
call_vntr_wrkf.sh set_env
call_vntr_wrkf.sh submit_call_vntr 10

# 目前开发完了提交检测vntr任务，clean，compress任务对应的脚本还未开发


