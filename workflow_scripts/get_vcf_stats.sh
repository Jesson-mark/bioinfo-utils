#!/usr/bin/bash 
set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

if [[ $# -eq 0 ]]; then
    echo "Usage: get_vcf_stats.sh vcf_file out_dir threads"
    echo "Function: "
    echo "    1. get SNP ID"
    echo "    2. get SNP POS and ID(bed)"
    echo "    3. get sample names"
    echo "    4. cal SNP freq(by plink)"
    exit 
fi

vcf_file=$1
out_dir=$2
threads=$3
mkdir -p "$out_dir"

myprint "fetching SNP ID"
bcftools query -f "%ID\n" "${vcf_file}" > "${out_dir}"/snp_ids.txt

myprint "fetching SNP POS and ID"
bcftools query -f "%CHROM\t%POS\t%ID\n" "${vcf_file}" | awk '{printf("%s\t%s\t%s\t%s\n", $1,$2-1,$2,$3)}' | bgzip > "${out_dir}"/snp_pos.bed.gz
tabix -p bed "${out_dir}"/snp_pos.bed.gz

myprint "fetching sample names"
# get_vcf_samples=/public/home/fan_lab/wangjie/utils/shell_utils/get_vcf_samples.sh
bcftools query -l "${vcf_file}" > "${out_dir}"/sample_names.txt

myprint "calculating frequency"
out_prefix="${out_dir}"/freq
plink2 --vcf "$vcf_file" --freq cols=+pos --out "${out_prefix}" --threads "$threads"
bgzip -f "${out_prefix}".afreq
tabix -f -S1 -s1 -b2 -e2 "${out_prefix}".afreq.gz


