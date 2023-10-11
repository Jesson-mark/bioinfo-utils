#/usr/bin/bash

if [[ $# -eq 0 ]]; then
    echo "Usage: get_vcf_samples vcf_file > results.txt"
    exit 0
fi

vcf_file=$1
bcftools view -h $vcf_file | tail -1 | cut -f10- | sed 's/\t/\n/g' 
