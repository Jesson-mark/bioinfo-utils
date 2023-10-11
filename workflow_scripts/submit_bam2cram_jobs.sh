ref_genome=/public/home/fan_lab/wangjie/genome/hg38/hg38_22_XYM.fa
threads=20
work_dir=/public/home/fan_lab/wangjie/GIAB_trio/GIAB_bams
sample_ids=sample_ids/sample_ids.txt
bam_dir=bams
cram_dir=crams
script=scripts/convert_bam2cram.pbs

for sample_id in $(cat $sample_ids)
do
    bam=${bam_dir}/${sample_id}.bam
    cram=${cram_dir}/${sample_id}.cram
    log=logs/bam2cram/${sample_id}.log

    args=$(echo work_dir=$work_dir,raw_bam=$bam,out_cram=$cram,ref_genome=$ref_genome,threads=$threads)
    echo qsub -o $log -l nodes=1:ppn=$threads -q fat -N ${sample_id} -v $args $script
done


