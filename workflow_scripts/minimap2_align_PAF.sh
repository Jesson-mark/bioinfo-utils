#PBS -j oe
#arguments: work_dir, data, name, threads

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

# variables
cd $work_dir
myprint "Working dir is: "$work_dir

# scripts
minimap2=$HOME"/Programs/bin/minimap2"
samtools=$HOME"/Programs/samtools-1.13/bin/samtools"
genome=$HOME"/genome/hg38/hg38_22_XYM.fa"

myprint "Begin align to genome"
$minimap2 -x asm5 -t $threads $genome $data > ${name}.paf

myprint "Program is done"

