process BurgessBamToVcf {
    label 'burgess'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(burgess_vcf), path(burgess_tbi)

    output:
        tuple val(sample), path(ann_vcf)

    script:
        ann_vcf = "${sample}_annotated.vcf.gz"

        """
        set -eo pipefail

        ${params.tools.gatk} HaplotypeCaller \
            -R ${params.ref} \
            --max-reads-per-alignment-start 0 --force-call-filtered-alleles true \
            --output-mode EMIT_ALL_ACTIVE_SITES -A AlleleFraction \
            --intervals ${burgess_vcf} --alleles ${burgess_vcf} \
            -I ${bam} -O ${sample}.vcf;
        
		${params.tools.bgzip} -f ${sample}.vcf;
		${params.tools.tabix} -f ${sample}.vcf.gz;

        ${params.tools.bcftools} annotate -a ${burgess_vcf} -c ID -o ${ann_vcf} ${sample}.vcf.gz;
        """
}
