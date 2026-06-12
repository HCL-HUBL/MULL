process BurgessScore {
    label 'burgess'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

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

        ${params.tools.plink} --vcf ./${vcf/.vcf/_annotated.vcf} \
            --score ${burgess_summary} 3 4 9 header \
            --out ${profile} --make-bed;
        """
}
