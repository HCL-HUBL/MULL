process BurgessScoresConcat {
    publishDir "${params.outdir}/", mode: 'copy'
    label 'burgess'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    input:
        path(profiles)
    
    output:
        path(profiles_concat)

    script:
        profiles_concat = "${params.outname}.profiles"

        """
        set -eo pipefail

        # From multiple space delimited file to tsv:
        cat ${profiles} | grep -v -w "FID" | sed 's/ \\+/\\t/g' | cut -f 2,3,4,5,6,7 > ${profiles_concat}_raw
        # Some profiles do not have a FID (in this case, we print the IID twice):
        awk '{if(\$2 == -9) {print  \$1"\\t"\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5} else {print}}' ${profiles_concat}_raw > ${profiles_concat}
        """
}
