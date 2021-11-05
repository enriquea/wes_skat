#!/usr/bin/env Rscript

library(SKAT)
library(data.table)
library(optparse)

# arguments
option_list = list(
  make_option(c("-i", "--path_plink"), type="character", default=NULL, 
              help="Path to bim, fam and bed PLINK-like format files", metavar="character"),
  make_option(c("-s", "--sid_file"), type="character", default=NULL, 
              help="Path to gene/variant set id file (two-columns, no header)", metavar="character"),
  make_option(c("--covariates_file"), type="character", default=NULL, 
              help="Path to sample covariates file (tab-separated, first column with sample ids)", metavar="character"),
  # make_option(c("-o", "--output_file"), type="character", default=NULL, 
  #            help="Path to output file with SKAT test results", metavar="character"),
  make_option(c("--run_variant"), type="logical", default=FALSE, action = 'store_true', 
              help="Run SKAT at variant level", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_dir = opt$path_plink
setid_file = opt$sid_file
covs_file = opt$covariates_file

# input_dir = "./testdata/plink_data/lof_hc"
# setid_file = "./testdata/plink_data/variants_ids_v2.tsv"
# covs_file = "./testdata/plink_data/covs.tsv"

# check dir with plink file paths
if (is.null(input_dir)){
   print_help(opt_parser)
  stop("The path to directory with PLINK files must be specified...")
}

# check covs file
if (is.null(covs_file)){
  # print_help(opt_parser)
  warning("Path to covariates file not specified. Running SKAT without covariates")
}

parse_covs_file <- function(covs_file){
   covs <- read.delim(covs_file,
                      stringsAsFactors = F)
   return(covs)
}


#' get_phe_from_fam
#' 
#' Re-code phenotype field (column 6) from PLINK-formatted fam file
#'
#' @param fam_file PLINK-like fam dataframe
#'
#' @return Formatted fam dataframe
#' @export
#'
#' @examples
recode_phe_from_fam <- function(fam_file){
  # recode (PLINK) binary phenotype 
  fam <- read.delim(fam_file, 
                    header = F, 
                    stringsAsFactors = F)
  
  fam$V6 <- ifelse(fam$V6==2, 1, 0)
  return(fam)
}


#' run_skat_variant
#' 
#' Run skat-o robust per variants
#'
#' @param ssd.info SSD.info file
#' @param obj output object from SKAT_Null_Model
#' @param index gene/feature index
#'
#' @return SKAT stats per variants
#' @export
#'
#' @examples
run_skat_variant <- function(ssd.info, obj, index){
  
  z <- Get_Genotypes_SSD(ssd.info, Set_Index = index)
  res <- SKATBinary_Robust(z,
                           obj,
                           method="SKATO",
                           weights.beta=c(1,25), 
                           weights = NULL, 
                           r.corr=1, 
                           impute.method = "bestguess", 
                           is_check_genotype=TRUE,
                           is_dosage = FALSE, 
                           missing_cutoff=0.15, 
                           max_maf=0.001, 
                           estimate_MAF=1)
  
  
  # getting variant stats per gene
  gene =            as.character(ssd.info$SetInfo$SetID[[index]])
  p.value_gene =    as.double(res$p.value)
  Q =               as.double(res$Q)
  mac.gene =        as.numeric(res$mac)
  p.value_variant = as.double(res$p.value_singlevariant)
  test.snp.mac =    as.numeric(res$test.snp.mac)
  
  # retrive SNP ids
  # TODO: if only one SNP per gene, the current SKAT implementation fail keeping the SNP id.
  # TODO: temporary, name it as singleton
  snp_id = names(res$test.snp.mac)
  if(length(snp_id) <= 1){
     snp_id=ifelse(is.null(snp_id), 'single_variant', snp_id)
  }
  
  
  variant_stats <- data.frame(gene=gene,
                              p.value_gene=p.value_gene,
                              Q=Q,
                              mac.gene=mac.gene,
                              snp_id = snp_id,
                              p.value_variant=p.value_variant,
                              test.snp.mac=test.snp.mac)
  
  return(variant_stats)
}


#' run_skat_gene
#' 
#' Run skat-o robust per genes/features
#'
#' @param ssd.info SSD.info file
#' @param obj output object from SKAT_Null_Model
#'
#' @return SKAT stats per genes/features
#' @export
#'
#' @examples 
run_skat_gene <- function(ssd.info, obj){
  
  SKATSSDall_rare <- SKATBinary_Robust.SSD.All(ssd.info,
                                               obj,
                                               method="SKATO",
                                               weights.beta=c(1,25), 
                                               weights = NULL, 
                                               r.corr=1, 
                                               impute.method = "bestguess", 
                                               is_check_genotype=TRUE,
                                               is_dosage = FALSE, 
                                               missing_cutoff=0.15, 
                                               max_maf=0.001, 
                                               estimate_MAF=1)
  
  res <- SKATSSDall_rare$results
  return(res)
}


# set bed file
file_bed <- paste0(input_dir, '.bed')

# set bim file 
file_bim <- paste0(input_dir, '.bim')

# set fam file
file_fam <- paste0(input_dir, '.fam')

# set file sid
file_setid <- setid_file

# set ssd file
file_ssd <- paste0(setid_file, '.ssd')

# set info file
file_info <- paste0(setid_file, '.info')


# Generate SSD object
Generate_SSD_SetID(File.Bed=file_bed,
                   File.Bim=file_bim,
                   File.Fam=file_fam,
                   File.SetID=file_setid,
                   File.SSD=file_ssd,
                   File.Info=file_info)


SSD.INFO <- Open_SSD(file_ssd, file_info)


## building NULL model
fam = recode_phe_from_fam(file_fam)

# build NULL model with covs (if a covs file exists)
if(!is.null(covs_file)){
  covs <- parse_covs_file(covs_file = covs_file)
  fam_covs <- merge.data.frame(fam,
                               covs,
                               by.x = 'V2',
                               by.y = 's',
                               all.x = T)
  
   print(head(fam_covs))

  obj <- SKAT_Null_Model(fam_covs$V6 ~ 1 + fam_covs$PC1 + fam_covs$PC2 + fam_covs$PC3 + fam_covs$PC4 + fam_covs$PC5,
                         out_type="D") # without covariantes
} else {
  obj <- SKAT_Null_Model(fam$V6 ~ 1, out_type="D") # without covariantes
}


# Run Skat test either per genes or variants
if(!opt$run_variant){
  
  #  Run SKATBinary_Robust per genes
  res <- run_skat_gene(ssd.info = SSD.INFO, obj = obj)
  
} else {
  
  # Run SKATBinary_Robust per variants 
  res.list <- list()
  for(index in 1:SSD.INFO$nSets){
    tmp <- run_skat_variant(ssd.info = SSD.INFO, obj = obj, index = index)
    res.list[[index]] <- tmp
  }
  
  res <- as.data.frame(data.table::rbindlist(res.list))
}

# close SSD file
Close_SSD()

# export results
output_file_path = paste0(setid_file, '.skat_results.tsv')
write.table(res,
            file = output_file_path,
            row.names = F,
            quote = F,
            sep = '\t')