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

// Default parameters:
params.bams_dir         = ""
params.ref              = ""
params.outname          = "mull"
params.outdir           = "./mull_out/"

params.burgess_vcf      = "${projectDir}/data/burgess_43_snps_sorted.vcf.gz"
params.burgess_summary  = "${projectDir}/data/burgess_43_snps.summary"

params.kmer_size        = 31
params.kiv2_kmers       = "${projectDir}/data/KIV2_hg38_kmers_6copies_specific.fasta"
params.norm_kmers       = "${projectDir}/data/LPA_hg38_kmers_1copies_specific.fasta"
params.rsids_list       = "${projectDir}/data/lpa_model_rsids.tsv"
params.quantiles        = ""


// Checking the input parameters:
if(params.bams_dir == "")           error("\nERROR: Path to BAMs folder is missing (see config: 'bams_dir')")
if(params.ref == "")                error("\nERROR: Path to the reference fasta file is missing (see config: 'ref')")

if(params.burgess_vcf == "")        error("\nERROR: Path to the burgess VCF file is missing (see config: 'burgess_vcf')")
if(params.burgess_summary == "")    error("\nERROR: Path to the burgess Summary file is missing (see config: 'burgess_summary')")

if(params.kmer_size <= 0)            error("\nERROR: The kmer size must be positive (see config: 'kmer_size')")
if(params.kiv2_kmers == "")         error("\nERROR: Path to the KIV-2 specific kmers fasta file is missing (see config: 'kiv2_kmers')")
if(params.norm_kmers == "")         error("\nERROR: Path to the normalisation specific kmers fasta file is missing (see config: 'norm_kmers')")
if(params.rsids_list == "")         error("\nERROR: Path to the rsids list for the Lp(a) model is missing (see config: 'rsids_list')")


// Generating Channels:
input_bams_ch = Channel.fromPath("${params.bams_dir}/*.bam").map { bam_file ->
                                                                        def sample_id = bam_file.baseName
                                                                        def bai_file = file("${params.bams_dir}/${sample_id}.bam.bai")
                                                                        tuple(sample_id, bam_file, bai_file) }

kiv2_kmers_ch = Channel.fromPath("${params.kiv2_kmers}")
norm_kmers_ch = Channel.fromPath("${params.norm_kmers}")
rsids_list_ch = Channel.fromPath("${params.rsids_list}")


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


workflow {

    BURGESS(input_bams_ch)
    burgess_profiles = BURGESS.out

    if(params.quantiles == "") {
        quantiles_ch = []
    } else {
        quantiles_ch = Channel.fromPath(params.quantiles)
    }

    KILDA(input_bams_ch,
          kiv2_kmers_ch,
          norm_kmers_ch,
          rsids_list_ch,
          quantiles_ch) 
}