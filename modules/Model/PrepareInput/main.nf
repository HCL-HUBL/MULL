process PrepareInput {
    publishDir "${params.outdir}/", mode: 'copy'
    
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

        kilda <- read.table(file = "${kilda_CNs}", header = T, sep = "\\t")

        # Computing if variants are present (if alt kmers > 0.25 total, then yes): 
        seuil <- 0.25
        for(c in seq(from = 3, to = ncol(kilda)-1, by = 2)) {
            ratio <- kilda[c+1] / (kilda[c+1] + kilda[c])
            names(ratio) <- gsub(pattern = "_alt", replacement = "", x = names(ratio))

            kilda[[names(ratio)]] <- ifelse(test = (ratio > seuil & !is.na(ratio)),
                                            yes = 1, no = 0)
        }

        # Adding the burgess PRS score:
        profiles <- read.table(file = "${profiles_concat}", header = F, sep = "\\t")
        colnames(profiles) <- c("FID", "IID", "PHENO", "CNT", "CNT2", "burgess_score")

        kiv2_burgess <- merge(x = kilda, y = profiles,
                            by.x = "ID", by.y = "IID",
                            all.x = F, all.y = F)

        kiv2_burgess\$KIV2_groups <- -1
        kiv2_burgess\$KIV2_groups[kiv2_burgess\$quantile <= 2] <- 0
        kiv2_burgess\$KIV2_groups[kiv2_burgess\$quantile > 2] <- 1

        snps <- c("rs41272114", "rs10455872", "rs3798220", "rs41259144", "rs140570886", 
                  "rs41267809", "rs41270998", "rs41267807", "rs1853021", "rs1800769", 
                  "rs3124784", "rs41272112", "rs11591147")
        
        model_input <- kiv2_burgess[,c("ID", "KIV2_groups", "burgess_score", snps)]

        write.table(x = model_input, file = "${merged_data}", quote = F, sep = "\\t", row.names = F)
        """
}