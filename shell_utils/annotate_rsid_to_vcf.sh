#!/usr/bin/bash

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

# softwares
table_annovar=$HOME"/Programs/annovar/table_annovar.pl"

if [[ $# -eq 0 ]]; then
    echo "Usage: scripts/utils/annotate_rsid_to_vcf.sh vcf_file out_prefix genome_version threads "
    exit 
fi

vcf_file=$1
out_prefix=$2
genome_version=$3 # hg19 or hg38
threads=$4

# def vars
out_vcf=${out_prefix}.RSID.vcf.gz
anno_prefix=${out_prefix}.annovar
avinput=${anno_prefix}.avinput
rsid_file=${anno_prefix}.rsid.txt.gz

# begin
starttime=$(date +'%Y-%m-%d %H:%M:%S')

myprint "Converting VCF to avinput"
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\n' "$vcf_file" > "${avinput}"

myprint "Annotating avinput"
humandb=$HOME/AnnovarHumandb/humandb/
$table_annovar -out "${anno_prefix}" -buildver "${genome_version}" "${avinput}" "$humandb" \
    --protocol avsnp150 -operation f -nastring . --thread "$threads"

myprint "Fetching rsid"
sed '1d' "${anno_prefix}"."${genome_version}"_multianno.txt |cut -f1-2,6 | bgzip > "$rsid_file"
tabix -s1 -b2 -e2 "${rsid_file}"

myprint "Annotating rsid to $vcf_file"
bcftools annotate -a "${rsid_file}" -c CHROM,POS,ID "${vcf_file}" -o "${out_vcf}" --threads "${threads}" 
tabix -p vcf "${out_vcf}"

myprint "Removing files of ANNOVAR"
myprint "Removing annovar files"
ls -lh "${anno_prefix}"* > "${out_prefix}".snapshot.annovar.files
rm -f "${anno_prefix}"*

# 计时结束
endtime=$(date +'%Y-%m-%d %H:%M:%S')
start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数
end_seconds=$(date --date="$endtime" +%s);

echo "$start_seconds" "$end_seconds" | awk '{total_sec=$2-$1; min=int(total_sec/60); lef_sec=total_sec%60; printf("Total time: %d minutes, %d seconds.\n", min, lef_sec); }'

myprint "Program is done"
