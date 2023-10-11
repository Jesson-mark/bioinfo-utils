#!/usr/bin/env python
#%%
import os
import sys
import yaml
import time
import shutil

if len(sys.argv) <= 1:
    print("Usage: python gen_scripts.py config.yml")
    sys.exit(1)

# config_file = "hg003_config.yml"
config_file = sys.argv[1]

with open(config_file, 'r') as ifl:
    config = yaml.safe_load(ifl)

sample_id = config["sample_id"] # args
work_dir = os.path.abspath(config["work_dir"]) # args
input_files = config["input_files"] # args
job_queue = config['job_queue']
total_threads = config["total_threads"] # args
parallel_jobs = config["parallel_jobs"] # args
applied_memory = config["applied_memory"] # args
sort_options_mem = config["sort_options_mem"] # args
minimap2_I_mem = config["minimap2_I_mem"] # args
genome_size = config["genome_size"] # args
enable_soft_clipping = config["enable_soft_clipping"]
cal_depth_threads = config['cal_depth_threads']

if not os.path.exists(work_dir):
    os.makedirs(work_dir)
print("work dir is", work_dir)
threads_each_job = total_threads // parallel_jobs

# 打印并保存生成脚本的命令
gen_cmd = "python %s"%(" ".join(sys.argv))
print("cmd is: %s"%(gen_cmd))
with open("%s/generate_cmd.sh"%(work_dir), 'w') as ofl:
    ofl.write(gen_cmd + "\n")

# 将 config_file 拷贝到工作路径下面
shutil.copy(config_file, work_dir)

# 创建目录
correction_dir = '%s/01_correction'%(work_dir)
alignment_dir = '%s/02_alignment'%(work_dir)
mean_depth_dir = '%s/03_mean_depth'%(work_dir)
scripts_dir = '%s/scripts'%(work_dir)
logs_dir = os.path.abspath('%s/logs'%(work_dir))
if not os.path.exists(correction_dir):
    print('making dir %s'%(correction_dir))
    os.makedirs(correction_dir)
if not os.path.exists(alignment_dir):
    print('making dir %s'%(alignment_dir))
    os.makedirs(alignment_dir)
if not os.path.exists(mean_depth_dir):
    print('making dir %s'%(mean_depth_dir))
    os.makedirs(mean_depth_dir)
if not os.path.exists(scripts_dir):
    print('making dir %s'%(scripts_dir))
    os.makedirs(scripts_dir)
if not os.path.exists(logs_dir):
    print('making dir %s'%(logs_dir))
    os.makedirs(logs_dir)

curr_time = time.strftime("%Y_%m_%d-%H_%M", time.localtime())

#1 correction
input_fofn = "%s/input.fofn"%(correction_dir)
run_cfg = "%s/run.cfg"%(correction_dir)
corrected_reads = os.path.abspath("%s/corrected_reads.fa.gz"%(correction_dir))
correct_log = "%s/correct_%s.log"%(logs_dir, curr_time)

#2 alignment
corrected_bam_prefix = "%s/corrected"%(alignment_dir)
corrected_bam = "%s.q20.sort.bam"%(corrected_bam_prefix)
align_log = "%s/align_%s.log"%(logs_dir, curr_time)

#3 mean_depth
mean_depth_prefix = "%s/mean_depth"%(mean_depth_dir)
cal_depth_log = "%s/cal_depth_%s.log"%(logs_dir, curr_time)

# 脚本
run_nextdenovo_script = "%s/run_nextdenovo.pbs"%(scripts_dir)
run_minimap2_script = "%s/run_minimap2.pbs"%(scripts_dir)
run_mosdepth_script = "%s/run_mosdepth.pbs"%(scripts_dir)

# 一些实用字符串
linux_date_str = "`date +'%Y-%m-%d %H:%M:%S'`"
cal_total_time = "echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); " + \
                    'hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\\n", hours, left_mins);' + \
                    "}'"

# %%
def gen_nextdenovo_script():
    # 1. 生成input.fofn, run.cfg文件
    # 2. 生成 run_nextdenovo.pbs 脚本，并将上述两个文件写进去
    # 3. 脚本里面跑完以后就建立一个done文件

    #1 generate input.fofn
    with open(input_fofn, 'w') as ofl:
        for file in config['input_files']:
            # todo 判断文件是否存在，若不存在就报错

            ofl.write("%s\n"%(os.path.abspath(file)))
    print('write %s done'%(input_fofn))

    #2 generate run.cfg
    with open(run_cfg, 'w') as ofl:
        ofl.write('[General]\n')
        ofl.write('job_type = local\njob_prefix = %s\n'%(sample_id))
        ofl.write('task = correct\nrewrite = yes\nparallel_jobs = %s\n'%(parallel_jobs))
        ofl.write('input_type = raw\nread_type = ont\n')
        ofl.write('input_fofn = input.fofn\n\n')

        ofl.write('[correct_option]\nread_cutoff = 1k\n')
        ofl.write('genome_size = %s\npa_correction = %s\n'%(genome_size, parallel_jobs))
        ofl.write('sort_options = -m %s -t %s\n'%(sort_options_mem, threads_each_job))
        ofl.write('minimap2_options_raw = -I %s -t %s\n'%(minimap2_I_mem, threads_each_job))
        ofl.write('correction_options = -p %s -b\n'%(threads_each_job))
    print('write %s done'%(run_cfg))

    #3 generate run_nextdenovo_script
    with open(run_nextdenovo_script, 'w') as ofl:
        ofl.write('#PBS -j oe\n#PBS -N %s_correct\n#PBS -l nodes=1:ppn=%s\n#PBS -o %s\n'%(
                        sample_id, total_threads, correct_log))
        ofl.write('#PBS -q %s\n#PBS -l mem=%s\n\n'%(job_queue, applied_memory))

        ofl.write('########################## Presets ##########################\n')
        ofl.write('set -ue\nsource "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"\n\n')

        ofl.write('########################## Arguments ##########################\n')
        ofl.write('# work_dir sample_out_dir sample_id corrected_reads threads\n')
        ofl.write('#!!! Note: This script will enter sample_out_dir, so corrected_reads must be absolute path\n\n')
        # write arguments
        ofl.write('work_dir=%s\nsample_out_dir=01_correction\n'%(os.path.abspath(work_dir)))
        ofl.write('sample_id=%s\ncorrected_reads=%s\nthreads=%s\n\n'%(sample_id, 
                                                            corrected_reads, total_threads))

        ofl.write('########################## softwares ##########################\n')
        ofl.write('seqkit=/public/home/fan_lab/wangjie//Programs/bin/seqkit\n')
        ofl.write('export PATH=/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin:$PATH \n')
        ofl.write('export PATH=/public/home/fan_lab/wangjie/Programs/NextDenovo:$PATH\n\n')

        ofl.write('########################## Change workdir ##########################\n')
        ofl.write('# change workdir\n')
        ofl.write('sample_work_dir=${work_dir}/${sample_out_dir}\n')
        ofl.write('myprint "Working dir is: "${sample_work_dir}\n')
        ofl.write('cd ${sample_work_dir}\n\n')

        ofl.write('########################## Begin ##########################\n')
        ofl.write('myprint "Begin to analysis"\n')
        ofl.write("starttime=`date +'%Y-%m-%d %H:%M:%S'`\n")
        ofl.write('nextDenovo run.cfg\n')
        ofl.write('myprint "Correct done for sample $sample_id"\n\n')
        ofl.write('myprint "Merging cns fasta ..."\n')
        ofl.write('concat_threads=4\n')
        ofl.write('time $seqkit scat -i fasta -j ${concat_threads} -f 02.cns_align/01.seed_cns.sh.work/ -o $corrected_reads \n')
        
        ofl.write('set +e # 遇到报错也不退出shell\n\n')
        ofl.write('myprint "removing db_split result of $sample_id"\n')
        ofl.write('myprint "cmd is: rm 01.raw_align/*.2bit, rm 01.raw_align/input.seed*, rm 01.raw_align/.input.seed*"\n')
        ofl.write('rm 01.raw_align/*.2bit\n')
        ofl.write('rm 01.raw_align/.input.seed*\n\n')
        ofl.write('myprint "removing sort_align result of $sample_id"\n')
        ofl.write('myprint "cmd is: rm 01.raw_align/04.sort_align.sh.work/sort_align*/*ovl, rm 01.raw_align/04.sort_align.sh.work/sort_align*/*ovl.bl"\n')
        ofl.write('rm 01.raw_align/04.sort_align.sh.work/sort_align*/*ovl\n')
        ofl.write('rm 01.raw_align/04.sort_align.sh.work/sort_align*/*ovl.bl\n\n')
        ofl.write('myprint "removing seed_cns result of $sample_id"\n')
        ofl.write('myprint "cmd is: rm 02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta"\n')
        ofl.write('rm 02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta\n')
        ofl.write('rm 02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta.idx\n\n')
        ofl.write('# 打包压缩日志文件\n')
        ofl.write('myprint "will pack and compress results to job_history"\n')
        ofl.write('if [[ ! -d "job_history" ]]; then\n')
        ofl.write('\tmyprint "dir job_history not exists, will make it "\n')
        ofl.write('\tmkdir job_history\nfi\n')
        ofl.write('curr_time=`date "+%Y-%m%d-%H%M"`\n')
        ofl.write('tar_result=%s_$curr_time.tar.gz\n'%(sample_id))
        ofl.write('myprint "cmd is: tar -zvcf $tar_result 01.raw_align 02.cns_align 03.ctg_graph"\n')
        ofl.write('tar -zvcf job_history/$tar_result 01.raw_align 02.cns_align 03.ctg_graph\n\n')
        ofl.write('myprint "done for packing and compressing, result is $tar_result, will remove those three dirs"\n')
        ofl.write('rm -rf 01.raw_align 02.cns_align 03.ctg_graph\n\n')
        ofl.write('touch correction.done\n\n')

        ofl.write('########################## Done ##########################\n')
        ofl.write("endtime=`date +'%Y-%m-%d %H:%M:%S'`\n")
        ofl.write('start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数\n')
        ofl.write('end_seconds=$(date --date="$endtime" +%s);\n')
        ofl.write("echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); ")
        ofl.write('hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\\n", hours, left_mins);')
        ofl.write("}'\n\n")
        ofl.write('myprint "Threads is $threads"\n')
        ofl.write('myprint "Program is done"\n')

    print('write %s done'%(run_nextdenovo_script))

gen_nextdenovo_script()

# %%
def gen_align_script():
    with open(run_minimap2_script, 'w') as ofl:
        ofl.write('#PBS -j oe\n#PBS -N %s_align\n#PBS -l nodes=1:ppn=%s\n#PBS -o %s\n'%(
                        sample_id, total_threads, align_log))
        ofl.write('#PBS -q %s\n\n'%(job_queue))
        ofl.write('########################## Presets ##########################\n')
        ofl.write('set -ue\n')
        ofl.write('source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"\n')
        ofl.write('########################## Arguments ##########################\n')
        ofl.write('#arguments: work_dir, data, preset out_prefix, threads, enable_soft_clipping\n')
        
        # write arguments
        ofl.write('work_dir=%s\ndata=%s\n'%(os.path.abspath(work_dir), corrected_reads))
        ofl.write('preset=ont\nout_prefix=%s\nthreads=%s\n'%(corrected_bam_prefix, total_threads))
        ofl.write('enable_soft_clipping=%s\n\n'%(enable_soft_clipping))

        ofl.write('########################## softwares ##########################\n')
        ofl.write('samtools="/public/home/fan_lab/wangjie/Programs/samtools-1.13/bin/samtools"\n')
        ofl.write('minimap2="/public/home/fan_lab/wangjie/Programs/minimap2_2.24/minimap2"\n')
        ofl.write('genome="/public/home/fan_lab/wangjie/genome/hg38/hg38_22_XYM.fa"\n\n')
        ofl.write('########################## Change workdir ##########################\n')
        ofl.write('myprint "Working dir is: "${work_dir}\n')
        ofl.write('cd ${work_dir}\n\n')
        ofl.write('########################## Begin ##########################\n')
        ofl.write('# scripts\n')
        ofl.write("starttime=%s\n"%(linux_date_str))
        ofl.write('myprint "Begin to analysis"\n\n')
        ofl.write('out_raw_bam=${out_prefix}.q20.bam\n')
        ofl.write('out_sort_bam=${out_prefix}.q20.sort.bam\n\n')
        ofl.write('if [ $preset == "hifi" ]; then\n')
        ofl.write('\tif [[ $enable_soft_clipping == "y" ]]; then\n')
        ofl.write('\t\t$minimap2 -t ${threads} -ax map-hifi -Y $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\telse\n')
        ofl.write('\t\t$minimap2 -t ${threads} -ax map-hifi $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\tfi\n')
        ofl.write('elif [ $preset == "ont" ]; then\n')
        ofl.write('\tif [[ $enable_soft_clipping == "y" ]]; then\n')
        ofl.write('\t\t$minimap2 -t ${threads} -ax map-ont -Y $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\telse\n')
        ofl.write('\t\t$minimap2 -t ${threads} -ax map-ont $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\tfi\n')
        ofl.write('elif [ $preset == "fasta" ]; then\n')
        ofl.write('\tif [[ $enable_soft_clipping == "y" ]]; then\n')
        ofl.write('\t\t$minimap2 -t ${threads} -a -Y $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\telse\n')
        ofl.write('\t\t$minimap2 -t ${threads} -a $genome ${data} | $samtools view -bS -q 20 -@ ${threads} > ${out_raw_bam}\n')
        ofl.write('\tfi\n')
        ofl.write('fi\n\n')
        ofl.write('myprint "Converting sam to bam, sorting bam..."\n')
        ofl.write('$samtools sort -l 5 -@ ${threads} -o ${out_sort_bam} ${out_raw_bam}\n')
        ofl.write('$samtools index -@ ${threads} ${out_sort_bam}\n\n')
        ofl.write('myprint "raw bam size are as follows"\n')
        ofl.write('ls -lh ${out_raw_bam}\n\n')
        ofl.write('myprint "Cleaning unsorted bam..."\n')
        ofl.write('rm ${out_raw_bam}\n\n')
        ofl.write('########################## Done ##########################\n')
        ofl.write('# 计时结束\n')
        ofl.write("endtime=`date +'%Y-%m-%d %H:%M:%S'`\n")
        ofl.write('start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数\n')
        ofl.write('end_seconds=$(date --date="$endtime" +%s);\n\n')
        ofl.write('%s\n\n'%(cal_total_time))
        # ofl.write("echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); ")
        # ofl.write('hours=int(mins/60); left_mins=mins%60; printf("Total time: %d hours, %d minutes.\\n", hours, left_mins);')
        # ofl.write("}'\n\n")
        ofl.write('myprint "Threads is $threads"\n')
        ofl.write('myprint "Program is done"\n')
    print('write %s done'%(run_minimap2_script))


gen_align_script()

#%%
def gen_cal_depth_script():
    with open(run_mosdepth_script, 'w') as ofl:
        ofl.write('#PBS -j oe\n#PBS -N %s_cal_depth\n#PBS -l nodes=1:ppn=%s\n#PBS -o %s\n'%(
                        sample_id, cal_depth_threads, cal_depth_log))
        ofl.write('#PBS -q %s\n\n'%(job_queue))

        ofl.write('########################## Presets ##########################\n')
        ofl.write('set -ue\n')
        ofl.write('source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"\n\n')
        ofl.write('########################## Arguments ##########################\n')
        ofl.write('#arguments: work_dir input_bam out_prefix threads \n')
        
        # write arguments
        ofl.write('work_dir=%s\ninput_bam=%s\n'%(os.path.abspath(work_dir), corrected_bam))
        ofl.write('out_prefix=%s\nthreads=%s\n\n'%(mean_depth_prefix, cal_depth_threads))
        
        ofl.write('########################## softwares ##########################\n')
        ofl.write('mosdepth="/public/home/fan_lab/wangjie/Programs/bin/mosdepth"\n\n')
        ofl.write('########################## Change workdir ##########################\n')
        ofl.write('myprint "Working dir is: "${work_dir}\n')
        ofl.write('cd ${work_dir}\n\n')
        ofl.write('########################## Begin ##########################\n')
        ofl.write('# scripts\n')
        ofl.write('starttime=%s\n'%(linux_date_str))
        ofl.write('myprint "Begin to analysis"\n')
        ofl.write('myprint "out_prefix is $out_prefix"\n')
        ofl.write('myprint "input bam is $input_bam"\n')
        ofl.write('myprint "Threads is $threads"\n\n')
        ofl.write('$mosdepth -n --fast-mode --by 500 ${out_prefix} ${input_bam} -t ${threads}\n\n')
        ofl.write('########################## Done ##########################\n')
        ofl.write('# 计时结束\n')
        ofl.write('endtime=%s\n'%(linux_date_str))
        ofl.write('start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数\n')
        ofl.write('end_seconds=$(date --date="$endtime" +%s);\n\n')
        ofl.write('%s\n\n'%(cal_total_time))
        ofl.write("echo $start_seconds $end_seconds | awk '{total_sec=$2-$1; mins=int(total_sec/60); " + \
                    'left_secs=mins%60; printf("Total time: %d minutes, %d seconds.\\n", mins, left_secs);' + \
                    "}'\n"
                )
        ofl.write('myprint "Program is done"\n')
    print('write %s done'%(run_mosdepth_script))

gen_cal_depth_script()
