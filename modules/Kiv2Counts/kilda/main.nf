// Counting the number of kiv2, from all the counts files
process kilda {
    publishDir "${params.outdir}/", mode: 'copy'
    
    label 'kiv2Counts'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    input:
       path(kiv2_kmers)
       path(norm_kmers)
       path(counts_list)
       path(quantiles)
    
    output:
        path(kilda_CNs)
    
    script:
        kilda_out = "kilda_out/"
        kilda_CNs = "${kilda_out}/kilda_kiv2.CNs"
        
        rsid_param = ""
        if("${params.rsids_list}" != "") {
            rsid_param = " -r ${params.rsids_list}"
        }

        q_params = ""
        if(quantiles) q_params = "-q ${quantiles}"

        """
        set -eo pipefail

        kilda.py \
            -c ${counts_list} -o ${kilda_out} -v \
            -k ${kiv2_kmers} -l ${norm_kmers} \
            ${rsid_param} ${q_params}
        """
}