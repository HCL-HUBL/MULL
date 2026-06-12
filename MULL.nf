#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.moduleBinaries = true

include { BurgessBamToVcf }        from './modules/BurgessBamToVcf'


// Default parameters:
params.bams_dir         = ""
params.ref              = ""
params.burgess_vcf      = "${projectDir}/data/burgess_43_snps_sorted.vcf.gz"
params.burgess_summary  = "${projectDir}/data/burgess_43_snps.summary"


// Checking the input parameters:
if(params.bams_dir == "")           error("\nERROR: Path to BAMs folder is missing (see config: 'bams_dir')")
if(params.ref == "")                error("\nERROR: Path to the reference fasta file is missing (see config: 'ref')")
if(params.burgess_vcf == "")        error("\nERROR: Path to the burgess VCF file is missing (see config: 'burgess_vcf')")
if(params.burgess_summary == "")    error("\nERROR: Path to the burgess Summary file is missing (see config: 'burgess_summary')")


// Generating Channels:
input_bams_ch = Channel.fromPath("${params.bams_dir}/*.bam").map { bam_file ->
                                                                        def sample_id = bam_file.baseName
                                                                        def bai_file = file("${params.bams_dir}/${sample_id}.bam.bai")
                                                                        tuple(sample_id, bam_file, bai_file) }

burgess_vcf_ch = tuple(Channel.fromPath("${params.burgess_vcf}"),
                       Channel.fromPath("${params.burgess_vcf}.tbi"))


workflow BURGESS {
    take:
        input_bams_ch
        burgess_vcf_ch

    main:
        BurgessBamToVcf(input_bams_ch, 
                        burgess_vcf_ch)
}

workflow {
    BURGESS(input_bams_ch,
            burgess_vcf_ch)
}