# run TRiCoLOR

```bash
mkdir vntr_results vntr_results/tricolor vntr_results/tricolor/sensor_outputs vntr_results/tricolor/refer_outputs

sensor_script=/public/home/fan_lab/wangjie/utils/workflow_scripts/tricolor_sensor.pbs
refer_script=/public/home/fan_lab/wangjie/utils/workflow_scripts/tricolor_refer.pbs
mmidir=/public/group_share_data/fan_lab/wangjie/genome/hg38/tricolor_mmi_index

work_dir=
sample_id_file=
haplotagged_bam_dir=
log_dir=
sensor_out_dir=vntr_results/tricolor/sensor_outputs
refer_out_dir=vntr_results/tricolor/refer_outputs
sensor_threads=1
refer_threads=20
for sample_id in $(cat $sample_id_file)
do
    echo -e "\n######### For sample "$sample_id
    haplotagged_bam=${haplotagged_bam_dir}/${sample_id}.bam
    sample_sensor_out_dir=${sensor_out_dir}/${sample_id}
    sample_refer_out_dir=${refer_out_dir}/${sample_id}

    # sensor
    echo "#1 sensor script"
    jobname=SN_${sample_id}
    log=${log_dir}/${jobname}.log
    args=$(echo work_dir=$work_dir,input_bam=$haplotagged_bam, \
            out_dir=$sample_sensor_out_dir | tr -d ' ' )
    echo qsub -N ${jobname} -q fat -o $log -l nodes=1:ppn=${sensor_threads} -v $args $sensor_script

    # refer
    echo "#2 refer script"
    bedfile=${sample_sensor_out_dir}/TRiCoLOR.srt.bed.gz
    samplename=${sample_id}
    jobname=RF_${sample_id}
    log=${log_dir}/${jobname}.log
    args=$(echo work_dir=$work_dir,samplename=$samplename,input_bam=$haplotagged_bam, \
            bedfile=$bedfile,out_dir=$sample_refer_out_dir,threads=$refer_threads,mmidir=$mmidir | tr -d ' ' )
    echo qsub -N ${jobname} -q fat -o $log -l nodes=1:ppn=${refer_threads} -v $args $refer_script
done > scripts/cmds_to_tricolor_uncor_ont_vntr.sh

```

# Template

将该文件作为提交任务的模板，在该文件的基础上进行修改即可 ~/utils/workflow_scripts/submit_bam2fa_jobs.sh

