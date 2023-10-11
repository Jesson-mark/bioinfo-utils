#PBS -j oe

### other parameters
#PBS -l nodes=1:ppn=1
#PBS -q default
#PBS -o test.log
#PBS -N job_name

########################## Presets ##########################
set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

########################## Arguments ##########################
# work_dir sample_out_dir sample_id corrected_reads

########################## softwares ##########################
seqkit=/public/home/fan_lab/wangjie//Programs/bin/seqkit
export PATH=/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin:$PATH # 
export PATH=/public/home/fan_lab/wangjie/Programs/NextDenovo:$PATH

########################## Change workdir ##########################
myprint "Working dir is: ""${work_dir}"
cd "${work_dir}"

########################## Begin ##########################
# scripts
starttime=$(date +'%Y-%m-%d %H:%M:%S')
myprint "Begin to analysis"
ls ~


########################## Done ##########################

# 计时结束
endtime=$(date +'%Y-%m-%d %H:%M:%S')
cal_used_time.sh "$starttime" "$endtime"

myprint "Threads is $threads"
myprint "Program is done"

