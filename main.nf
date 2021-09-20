#!/usr/bin/env nextflow

nextflow.enable.dsl=2

 
// Set input parameters

// Set input gene/snp id file
params.snp_file  = "$PWD/testdata/plink_data/variants_ids_v2.tsv"
snp_file_ch = channel.fromPath(file(params.snp_file))

// Set input plink directory
params.plink_input_dir = "$PWD/testdata/plink_data/lof_hc"

// Set number of groups to be generated
params.n_sets = 10

// Run SKAT at variant level (default: gene level)
params.run_variants = false


/*
 * Check snp input file and split into gene/SNP sets 
 */
process generateSNPsets {
    
    input:
      path f
    
    output:
      path '*.set_*.tsv', emit: ch_split_snps

    """
    generate_snp_sets.R -i $f --n_sets ${params.n_sets}
    """
}


/*
 * Run SKAT tets per gene chunks
 */
process runSKAT {

    input:
      path sid

    output:
      path '*.skat_results.tsv', emit: ch_results
    
    script:
      if( params.run_variants )
        """
        skat_pipeline.R -i ${params.plink_input_dir} -s ${sid} --run_variant
        """
      else
        """
        skat_pipeline.R -i ${params.plink_input_dir} -s ${sid}
        """
}


/*
 * Main workflow
 */
workflow {
    
    // 1. split gene/snp ids into chuncks
    split_snp_ch = generateSNPsets(snp_file_ch)
    
    // 2. run skat-o test per chuncks
    skat_ch = runSKAT(split_snp_ch.flatten())
    
    // 3. combine results and write to disk
    results = skat_ch
              .collectFile(name: 'merged_skat_results.csv', skip: 1, keepHeader: true)
              .subscribe { println( "Results from SKAT test wrote to ${it}.") }

}

