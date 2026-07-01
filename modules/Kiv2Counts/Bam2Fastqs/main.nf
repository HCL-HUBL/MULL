// Usually KILDA starts from Fastqs, for the Lp(a) model it is easier to start 
// with BAMs for genotyping and Burgess computations, hence this process:
process Bam2Fastqs {
    label 'kiv2Counts'
    tag "${sample}"

    cpus 1
    memory { 100.MB * task.attempt }
    time { 15.minute * task.attempt }

    input:
        tuple val(sample), path(bam), path(bai)

    output:
        tuple val(sample), path(fastq)
    
    script:
        fastq = "${sample}.fastq"

        """
        set -eo pipefail

        ${params.tools.samtools} fastq \
            -o ${fastq} ${bam}
        """
}