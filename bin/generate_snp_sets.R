#!/usr/bin/env Rscript

## author: eam
## date: 2021-09-17

library(optparse)

# arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to gene/variant id file (two-columns, no header)", metavar="character"),
  # make_option(c("-o", "--output_path"), type="character", default='./snp_sets', 
  #            help="Path to output file(s) with sets of variant ids", metavar="character"),
  make_option(c("-n", "--n_sets"), type="integer", default=10, 
              help="Number of sets to be generated", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_file <- opt$input
output_dir <- opt$output_path
n_sets <- opt$n_sets

# input_file = '/Users/enrique/projects/local/wes_skat/testdata/plink_data/variants_ids_v2.tsv'
# output_dir = './snp_sets'

# check args
if (is.null(input_file)){
  print_help(opt_parser)
  stop("Files with gene/snp ids must be specified.")
}


df <- read.delim(input_file, header = F)
names(df) <- c('genes', 'snp_id')

gene_list <- unique(as.character(df[,'genes']))

size_group = length(gene_list)/n_sets

chunks <- split(
  gene_list, ceiling(seq_along(gene_list)/size_group)
)

# check if output dir exist or create it, and write results
# if(!dir.exists(output_dir)){
#  output_dir = dir.create(output_dir, showWarnings = T)
# } 

for(index in 1:length(chunks)){
  set <- subset(df, genes %in% chunks[[index]])
  names(set) <- NULL # remove header
  file_path = paste0(input_file, '.set_', index, '.tsv')
  write.table(set,
              file = file_path,
              sep = '\t',
              row.names = F,
              quote = F)
}
  




