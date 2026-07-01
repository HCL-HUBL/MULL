// This process takes a bam as input and genotypes the SNPs present in "burgess_vcf"
// Then it annotates the SNP IDs with the ids present in "burgess_vcf"
process BurgessBamToVcf {
    label 'burgess'
    tag "${sample}"

    cpus 3
    memory { 1.GB * task.attempt }
    time { 10.minute * task.attempt }

    input:
        tuple val(sample), path(bam), path(bai)
        val(burgess_vcf)

    output:
        tuple val(sample), path(ann_vcf), path(ann_tabix)

    script:
        ann_vcf = "${sample}_annotated.vcf.gz"
        ann_tabix = "${ann_vcf}.tbi"

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
        ${params.tools.tabix} -f ${ann_vcf};
        """
}
