process BurgessFromVCF {
    label 'burgess'

    cpus 1
    memory { 1.GB * task.attemp }
    time { 1.hour * task.attemp }

    tag "${sample}"

    input:
        tuple val(sample), path(vcf), path(tabix)
        path(burgess_vcf)
        path(burgess_summary)
    
    output:
        tuple val(sample), path(profile)
    
    script:
        profile = "${sample}_burgess"

        """
        set -eo pipefail

        ${params.tools.bcftools} annotate -a ${burgess_vcf} -c ID -o ./${vcf/.vcf/_annotated.vcf} ${vcf};

        ${params.tools.plink} --vcf ./${vcf/.vcf/_annotated.vcf} \
				--score ${burgess_summary} 3 4 9 header \
				--out ${profile} --make-bed;
        """
}
