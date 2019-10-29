#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(argparse)
})

args = commandArgs(TRUE)
if (length(args) == 0) {
    message('Run convert-legacy-rdata.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Convert .Rdata file from older versions of facets-suite to the new .rds format.')
parser$add_argument('-i', '--input-file', required = TRUE,
                    help = 'Input file')
parser$add_argument('-o', '--output-prefix', required = FALSE,
                    help = 'Name of output file [default input basename]')
args = parser$parse_args()

output = ifelse(is.null(args$output_prefix),
                paste0(gsub('.Rdata', '', basename(args$input_file)), '.rds'),
                paste0(args$output_prefix, '.rds'))

# Convert input to new format -------------------------------------------------------------------------------------

load(args$input_file)

new_format = list(
    snps = out$jointseg,
    segs = fit$cncf,
    purity = fit$purity,
    ploidy = fit$ploidy,
    diplogr = fit$dipLogR,
    alballogr = out$alBalLogR,
    flags = out$flags,
    em_flags = fit$emflags,
    loglik = fit$loglik
)

saveRDS(new_format, output)
