// Querying the relevant kmers (KIV2 + normalisation regions)
process DumpKmers {
    publishDir "${params.outdir}/counts/", mode: 'copy'

    label 'kiv2Counts'

    cpus 1
    memory { 100.MB * task.attempt }
    time { 10.minute * task.attempt }

    input:
        path(jellyfish_kmers)
        
    output:
        path(counts)
    
    script:
        counts = "${jellyfish_kmers}".replaceAll(/.kmers/, ".counts")
        
        """
        set -eo pipefail
        ${params.tools.jellyfish} dump ${jellyfish_kmers} -c -t > ${counts}
        """
}