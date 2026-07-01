// Counting all the kmers in the fastq files
process CountKmers {
    label 'kiv2Counts'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 15.minute * task.attempt }

    input:
        tuple val(sample), path(fastqs)
        val(norm_kiv2_fasta)
    
    output:
        path(jellyfish_kmers)
    
    script:
        jellyfish_kmers = "${sample}.kmers"

        """
        set -eo pipefail
        
        ${params.tools.jellyfish} count -t ${task.cpus} \
            -m ${params.kmer_size} -s 100M -C \
            -o ${jellyfish_kmers} --if=${norm_kiv2_fasta} <(zcat -f ${fastqs.join(" ")});
        """
}