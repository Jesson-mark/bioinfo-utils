#!/usr/bin/bash 

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

echo "----- hg19 to hg38"
echo 'hg19_to_hg38_chain_file=/public/home/fan_lab/wangjie/genome/liftover/hg19ToHg38.over.chain.gz'
echo 'liftOver "$hg19_pos" $hg19_to_hg38_chain_file "$hg38_pos" "$unmapped"'

echo ""
echo "----- hg38 to hg19"
echo 'hg38_to_hg19_chain_file=/public/home/fan_lab/wangjie/genome/liftover/hg38ToHg19.over.chain.gz'
echo 'liftOver "$hg38_pos" $hg38_to_hg19_chain_file "$hg19_pos" "$unmapped"'


