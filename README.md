# wes_skat

### Nextflow workflow to run the SNP-set Kernel Association Test (SKAT) on whole-exome data (WES)


SKAT [1] is an emerging statistical method to test for association between variants/genes and continuous/binary phenotypes. The original SKAT R-package is single-thread based, meaning that the runtime can be long when analysing large dataset (thousands of genes/features) and multiple conditions (variant types and phenotypes).

Here, it is offered a work-around solution based on Nextflow (https://www.nextflow.io), following the split-apply-combine paradigm. The version presented here has been used for analysing large-scale exome-sequencing case-control cohorts, using the SKAT-O method.


#### Install/Dependences

Nextflow installation:
`curl -s https://get.nextflow.io | bash`

See more info about Nexflow [here](https://www.nextflow.io).

R-dependesnces:

`R >= 3.5.2, R-package ‘SKAT’ v2.2.2, R-package ‘optparse’ v1.6.6`


#### Useage

Run SKAT to generate variant stats as showed below.

Command:

```
./nextflow run main.nf --snp_file path/to/snp_id/file \
                       --plink_input_file path/to/plink/files \
                       --n_sets 10 \
                       --run_variants
```

Parameters:

 `snp_file`: Path to SNP id file. Expected a TSV file with two-columns, no-header.
               First column refers to genes/features, second to SNP IDs.
               
 `plink_input_file`: Path to PLINK-formatted file(s) basename (without the extensions bim, bam and fam).
 
 `n_sets`: Number of sets/chuncks to be generated.
 
 `run_variants`: If set, run SKAT at variant level. Otherwise, run it at gene level (default).

Terminal output:

```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [awesome_archimedes] - revision: 1586fead8e
executor >  local (11)
[db/c534db] process > generateSNPsets (1) [100%] 1 of 1 ✔
[b5/034e96] process > runSKAT (6)         [100%] 10 of 10 ✔
Results from SKAT test wrote to .../6d/merged_skat_results.csv
```


#### Results

Example output SKAT gene stats [here](https://github.com/enriquea/wes_skat/blob/master/data/skat.gene.results.tsv).

Example output SKAT variant stats [here](https://github.com/enriquea/wes_skat/blob/master/data/skat.variant.results.tsv).



#### References 

[1] Zhao Z, Bi W, Zhou W, VandeHaar P, Fritsche LG, Lee S. UK Biobank Whole-Exome Sequence Binary Phenome Analysis with Robust Region-Based Rare-Variant Test. Am J Hum Genet. 2020 Jan 2. doi: 10.1016/j.ajhg.2019.11.012. PMID: 31866045.

