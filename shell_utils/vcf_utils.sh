
function get_vcf_samples() {
    vcf=$1
    bcftools view -h $vcf | tail -1 | cut -f10- | sed 's/\t/\n/g' 
}

