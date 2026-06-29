// Computes the Burgess PRS score for one sample (VCF):
process BurgessScore {
    label 'burgess'
    tag "${sample}"

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    input:
        tuple val(sample), path(vcf), path(tbi)
        val(burgess_summary)

    output:
        path(profile)
    
    script:
        outname = "${sample}_burgess"
        profile = "${sample}_burgess.profile"

        """
        set -eo pipefail

        ${params.tools.plink} --vcf ./${vcf} \
            --score ${burgess_summary} 3 4 9 header \
            --out ${outname} --make-bed;
        """
}