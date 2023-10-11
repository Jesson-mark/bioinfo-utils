#/usr/bin/bash
set -ue

if [[ $# -eq 0 ]]; then
    echo "Usage: liftover.sh raw_bed new_bed unmapped direction"
    echo "Avaliable direction:"
    echo "    hg19_to_hg38, hg38_to_hg19"
    exit
fi

raw_bed=$1
new_bed=$2
unmapped=$3
direction=$4

echo "raw_bed is $raw_bed"
echo "new_bed is $new_bed"
echo "unmapped is $unmapped"
echo "direction is $direction"

liftOver=/public/home/fan_lab/wangjie/Programs/bin/liftOver
chainfile=~/genome/liftover/hg19ToHg38.over.chain.gz

if [[ $direction == 'hg19_to_hg38' ]]; then
    chainfile=~/genome/liftover/hg19ToHg38.over.chain.gz
elif [[ $direction == 'hg38_to_hg19' ]]; then
    chainfile=~/genome/liftover/hg38ToHg19.over.chain.gz
else
    echo "Wrong direction"
    exit 1
fi

$liftOver $raw_bed $chainfile $new_bed $unmapped
