job_queue=fsh_team
applied_memory=10gb
script=../../ont_correction/scripts/minimap2_align.pbs

# args to script
work_dir=/public/home/fan_lab/wangjie/Quartet_data/vntr/analysis
data=sub_bams/hifi_LCL5_minimap2_s_chr1_10628_10800.fq.gz
preset=hifi
out_prefix=hifi_LCL5_chr1_10628_10800
out_dir=realign
threads=2
enable_soft_clipping=y
args=$(echo work_dir=$work_dir,preset=$preset, \
            data=$data,out_prefix=$out_prefix,out_dir=$out_dir, \
            threads=$threads, \
            enable_soft_clipping=$enable_soft_clipping | tr -d ' ' )
log=realign_${out_prefix}.log

qsub -N align -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script
