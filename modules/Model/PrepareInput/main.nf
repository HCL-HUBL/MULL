process PrepareInput {
    label 'Model'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.hour * task.attempt }

    input:
        path(profiles_concat)
        path(kilda_CNs)
    
    output:
        path(merged_data)
    
    script:
        merged_data = "${params.outname}_merged.tsv"
        
        """
        #!/usr/bin/env Rscript

        profiles <- read.table(file = "${profiles_concat}", header = F, sep = "\t")
        kilda <- read.table(file = "${kilda_CNs}", header = T, sep = "\t")
        """
}