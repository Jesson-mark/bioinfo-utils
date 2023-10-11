# create dirs
view_dir=$1
mkdir ${view_dir}
mkdir ${view_dir}/sub_bams
mkdir ${view_dir}/igv_scripts
mkdir ${view_dir}/snapshots

# get region
bcftools query -f '%CHROM\t%POS\t%END\n' 05_vntr/merged_vntrs/ecorv_set/batch1.ecorv.vcf.gz |he -4 > $view_dir/region.bed

# get bam path
sample_ids=(17HanZZ0003 17HanZZ0010 17HanZZ0014 17HanZZ0015 17HanZZ0018 17HanZZ0019 17HanZZ0020)
raw_bam_dir=/mnt/DWH/zone00/groups_data/fan_lab/wangjie/ZYQueue/corrected_reads/bams
[[ -f ${view_dir}/bam_paths.txt ]] && rm -f ${view_dir}/bam_paths.txt
for sample_id in ${sample_ids[*]}
do
    bam=${raw_bam_dir}/${sample_id}.q20.sort.bam
    echo $bam >> ${view_dir}/bam_paths.txt
done 

# get sub bam
fetch_alignments.sh ${view_dir}/bam_paths.txt ${view_dir}/region.bed ${view_dir}/sub_bams

# gen bam path
sample_ids=(17HanZZ0003 17HanZZ0010 17HanZZ0014 17HanZZ0015 17HanZZ0018 17HanZZ0019 17HanZZ0020)
[[ -f ${view_dir}/sampleid_bampaths.txt ]] && rm -f ${view_dir}/sampleid_bampaths.txt
for sample_id in ${sample_ids[*]}
do
    bam=${view_dir}/sub_bams/${sample_id}.q20.sort.bam
    echo -e "$sample_id\t"$bam >> ${view_dir}/sampleid_bampaths.txt
done 

# generate igv scripts
sampleid_bampaths=${view_dir}/sampleid_bampaths.txt
bams_per_snapshot=4
sites_file=${view_dir}/region.bed
igv_scripts_dir=${view_dir}/igv_scripts
win_workdir=D:/Desktop/Learning/GraduateStudent/ZYQueue/ONT_error_correction/05_vntr/merged_vntrs/ecorv_set/view_four_vntrs
bam_dir=${win_workdir}/sub_bams
snapshot_dir=${win_workdir}/snapshots
gen_igv_scripts.sh $sampleid_bampaths $bams_per_snapshot $sites_file $igv_scripts_dir $bam_dir $snapshot_dir

cat ${igv_scripts_dir}/chr*txt > ${igv_scripts_dir}/merged.txt

exit 
# create dirs
view_dir=$1
region=$2
raw_bam_path=$3
view_multi_copy_vntrs

mkdir ${view_dir}
mkdir ${view_dir}/sub_bams
mkdir ${view_dir}/igv_scripts
mkdir ${view_dir}/snapshots

# get region
cp $region ${view_dir}

# get bam path
sample_ids=(17HanZZ0003 17HanZZ0010 17HanZZ0014 17HanZZ0015 17HanZZ0018 17HanZZ0019 17HanZZ0020)
raw_bam_dir=/mnt/DWH/zone00/groups_data/fan_lab/wangjie/ZYQueue/corrected_reads/bams
[[ -f ${view_dir}/bam_paths.txt ]] && rm -f ${view_dir}/bam_paths.txt
for sample_id in ${sample_ids[*]}
do
    bam=${raw_bam_dir}/${sample_id}.q20.sort.bam
    echo $bam >> ${view_dir}/bam_paths.txt
done 

# get sub bam
fetch_alignments.sh ${view_dir}/bam_paths.txt ${view_dir}/region.bed ${view_dir}/sub_bams

# gen bam path
sample_ids=(17HanZZ0003 17HanZZ0010 17HanZZ0014 17HanZZ0015 17HanZZ0018 17HanZZ0019 17HanZZ0020)
[[ -f ${view_dir}/sampleid_bampaths.txt ]] && rm -f ${view_dir}/sampleid_bampaths.txt
for sample_id in ${sample_ids[*]}
do
    bam=${view_dir}/sub_bams/${sample_id}.q20.sort.bam
    echo -e "$sample_id\t"$bam >> ${view_dir}/sampleid_bampaths.txt
done 

# generate igv scripts
sampleid_bampaths=${view_dir}/sampleid_bampaths.txt
bams_per_snapshot=4
sites_file=${view_dir}/region.bed
igv_scripts_dir=${view_dir}/igv_scripts
win_workdir=D:/Desktop/Learning/GraduateStudent/ZYQueue/ONT_error_correction/05_vntr/merged_vntrs/ecorv_set/view_four_vntrs
bam_dir=${win_workdir}/sub_bams
snapshot_dir=${win_workdir}/snapshots
gen_igv_scripts.sh $sampleid_bampaths $bams_per_snapshot $sites_file $igv_scripts_dir $bam_dir $snapshot_dir

cat ${igv_scripts_dir}/chr*txt > ${igv_scripts_dir}/merged.txt


