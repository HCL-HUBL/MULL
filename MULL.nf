#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.moduleBinaries = true

include { BurgessBamToVcf }     from './modules/Burgess/BurgessBamToVcf'
include { BurgessScore }        from './modules/Burgess/BurgessScore'
include { BurgessScoresConcat } from './modules/Burgess/BurgessScoresConcat'

include { Bam2Fastqs }          from './modules/Kiv2Counts/Bam2Fastqs'
include { CountKmers }          from './modules/Kiv2Counts/CountKmers'
include { CreateFastaKmers }    from './modules/Kiv2Counts/CreateFastaKmers'
include { CreateSampleMap }     from './modules/Kiv2Counts/CreateSampleMap'
include { DumpKmers }           from './modules/Kiv2Counts/DumpKmers'
include { kilda }               from './modules/Kiv2Counts/kilda'

include { PrepareInput }        from './modules/Model/PrepareInput'
include { PredictLpa }          from './modules/Model/PredictLpa'


// Default parameters:
params.samplesheet      = ""
params.ref              = ""
params.outname          = "mull"
params.outdir           = "./mull_out/"
params.quantile         = ""

// Fixed parameters:
params.burgess_vcf     = "${projectDir}/data/burgess_43_snps_sorted.vcf.gz"
params.burgess_summary = "${projectDir}/data/burgess_43_snps.summary"
params.kiv2_kmers      = "${projectDir}/data/KIV2_hg38_kmers_6copies_specific.fasta"
params.norm_kmers      = "${projectDir}/data/LPA_hg38_kmers_1copies_specific.fasta"
params.model           = "${projectDir}/data/SVM_lpa_model.RData"
params.kmer_size        = 31
params.rsids_list       = "${projectDir}/data/lpa_model_rsids.tsv"

// Checking the input parameters:
if(params.samplesheet == "")        error("\nERROR: Samplesheet is missing (see config: 'samplesheet')")
if(params.ref == "")                error("\nERROR: Path to the reference fasta file is missing (see config: 'ref')")

// Generating Channels:
samplesheet_ch  = Channel.fromPath("${params.samplesheet}")
kiv2_kmers_ch   = Channel.fromPath("${params.kiv2_kmers}")
norm_kmers_ch   = Channel.fromPath("${params.norm_kmers}")
rsids_list_ch   = Channel.fromPath("${params.rsids_list}")
model_ch        = Channel.fromPath("${params.model}")

workflow BURGESS {
    take:
        input_bams

    main:
        BurgessBamToVcf(input_bams, 
                        "${params.burgess_vcf}")
       
        BurgessScore(BurgessBamToVcf.out,
                     "${params.burgess_summary}")

        BurgessScoresConcat(BurgessScore.out.collect())
    
    emit:
        BurgessScoresConcat.out
}


workflow KILDA {
    take:
        input_bams
        kiv2_kmers
        norm_kmers
        rsids
        quantiles

    main:
        norm_kiv2_fasta = CreateFastaKmers(kiv2_kmers, norm_kmers, rsids)

        Bam2Fastqs(input_bams)        
        jelly_kmers_ch = CountKmers(Bam2Fastqs.out, norm_kiv2_fasta.first())
        counts_ch = DumpKmers(jelly_kmers_ch).collect()
        counts_list_ch = CreateSampleMap(counts_ch)
        
        kilda(kiv2_kmers, norm_kmers, counts_list_ch, quantiles)

    emit:
        kilda.out
}


workflow MODEL {
    take:
        profiles
        kilda
        model
    
    main:
        PrepareInput(profiles, kilda)
        PredictLpa(PrepareInput.out, model)
}


workflow {
    // The samplesheet should contain the sample id, the bam and the bai:
    input_bams_ch = samplesheet_ch.splitCsv(sep: '\t').map{ [ it[0], it[1], it[2] ] }

    BURGESS(input_bams_ch)
    burgess_profiles = BURGESS.out


    if(params.quantiles == "") { quantiles_ch = [] } 
    else { quantiles_ch = Channel.fromPath(params.quantiles) }

    KILDA(input_bams_ch,
          kiv2_kmers_ch,
          norm_kmers_ch,
          rsids_list_ch,
          quantiles_ch)
    kilda_cn = KILDA.out


    MODEL(burgess_profiles,
          kilda_cn,
          model_ch)
}