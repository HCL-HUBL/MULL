process PredictLpa {
    publishDir "${params.outdir}/", mode: 'copy'
    
    label 'Model'

    cpus 1
    memory { 1.GB * task.attempt }
    time { 10.minute * task.attempt }

    input:
        path(merged_data)
        path(lpa_model)
    
    output:
        path(predictions)
    
    script:
        predictions = "${params.outname}_predictions.tsv"
        
        """
        #!/usr/bin/env Rscript

        load("${lpa_model}")

        model_input <- read.table(file = "${merged_data}", header = T, sep = "\\t")
        
        snps <- c("rs41272114", "rs10455872", "rs3798220", "rs41259144", "rs140570886", 
                  "rs41267809", "rs41270998", "rs41267807", "rs1853021", "rs1800769", 
                  "rs3124784", "rs41272112", "rs11591147")

        for(snp in snps) {
            model_input[[snp]] <- as.factor(model_input[[snp]])
        }

        predictions <- data.frame(SVM_lpa_g_L = round(predict(svm_poly_fit1, model_input),2),
                                  model_input)

        write.table(x = predictions, file = "${predictions}", quote = F, sep = "\\t", row.names = F)
        """
}