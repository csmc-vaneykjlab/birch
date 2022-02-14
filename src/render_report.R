.libPaths( c("C:/Users/BhatA/Documents/R/R-3.6.3/library/" , .libPaths() ) )

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(assertthat)
  library(readr)
  #library(purrr)
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  #library(patchwork)
  library(DEP)
  library(rmarkdown)
  library(proBatch)
  library(gtable)
  library(grid)
  library(gridExtra)
  library(sqldf)
  library(pcaMethods)
  library(ComplexHeatmap)
  library(SummarizedExperiment)
  library(tools)
  library(missRanger)
  library(cowplot)
  library(wrMisc)
  library(reshape)
  
})

#----------------Process script parameters-----------------
# params duplicated here and in .Rmd in order to preserve --help msg

options(warn=-1)

### Define CLI
option_list <- list(
  make_option(c("-N", "--input_norm"), action="store", default=NA, type="character",
              help="Input sample TIC normalized data file (required)"),
  make_option(c("-n", "--input_unnorm"), action="store", default=NA, type="character",
              help="Input sample TIC UNnormalized data file (required)"),
  make_option(c("-m", "--metadata_annotation"), action="store", default=NA, type="character",
              help="Input annotation file for metadata. (required)"),
  #make_option(c("-r", "--iRT_annotation"), action="store", default=NA, type="character",
  #            help="Input file with two columns 'Protein' and 'peptide_group_label' identifying iRT fragments (required)"),
  make_option(c("-o", "--output_dir"), action="store", default=NA, type='character',
              help="Specify the output directory. Use absolute path. (required) "),
  make_option(c("-v", "--outfile_prefix"), action="store", default="", type="character",
              help="Filename prefix for the output files, default none"),
  make_option(c("-b", "--batch_column"), action="store", default="Digestion_batch", type="character",
              help="Which annotation file column will batch correction be performed on? Either 'Digestion_batch' or 'MS_batch'"),
  make_option(c("-c", "--cols_of_interest"), action="store", default="Digestion_batch,MS_batch", type='character',
              help="Comma separated list of columns to consider."),
  make_option(c("-s", "--sample_threshold"), action="store", default=0.7, type="double",
              help="Keep only the samples with at least X% of data available across samples, default 0.7 (70%)"),
  make_option(c("-e", "--expgroup_threshold"), action="store", default=0.5, type="double",
              help="Keep only features with at least X% of data available across any experimental group, default 0.5 (50%)"),
  make_option(c("-t", "--batch_threshold"), action="store", default=0.7, type="double",
              help="Keep only features with at least X% of data available across all batches, default 0.7 (70%)"),
  make_option(c("-p", "--samples_for_correlation"), action="store", default="DR:Digestion_batch,TR:Technical_batch", type="character",
              help="Digestion and Technical Rep Sample keywords from attribute_experimental group along with column to plot"),
  make_option(c("-i", "--imputation_method"), action="store", default="ranger", type="character",
              help='Specify one of 
                         "zero","minimum","colmedian","rowmedian",
                         "knnmethod","seqknn","bpca","svdmethod",
                         "lls","mle","qrilc","mindet","minprob",
                         "impseq","impseqrob",
                         "mice-norm","mice-cart","trknn",
                         "rf","pi","grr","gms", or "halfminimum", 
                    default "ranger"'),
  make_option(c("--iRT_protein_name"), action="store", default="irt_protein", type="character",
              help="Name of the iRT peptide, default 'irt_protein'"),
  make_option(c("-q", "--quantile_norm"), action="store_true", default=FALSE, type="logical",
              help="When flag is set, perform quantile normalization instead of using norm input from ProEpic"),
  make_option(c("-B", "--batch_correction_first"), action="store_true", default=FALSE, type="logical",
              help="When flag is set, perform batch correction before imputation. By default, performs imputation before batch correction."),
  make_option(c("-w", "--pdf_width"), action="store", default=12, type="double",
              help="Set the width of the PDF output in inches, default 12"),
  make_option(c("-l", "--pdf_length"), action="store", default=7, type="double",
              help="Set the height of the PDF output in inches, default 7")
)

args = parse_args(OptionParser(option_list=option_list))
#--------------------Pass args to .Rmd---------------------

get_script_directory <- function(){
  commandArgs() %>% 
    tibble::enframe(name=NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value) %>% 
    dirname()
}
outpath <- args$output_dir

TEMP_OUTDIR <- file.path(
  outpath, paste0("temp_","probatch_report","_files"), 
  fsep = "\\"
)
dir.create(TEMP_OUTDIR)

tryCatch({
  rmarkdown::render(
    file.path('C:\\Users\\BhatA\\Box\\Batch_correction\\Batch-Correction-tool\\src\\Probatch_Report.Rmd'),
    output_format = "html_document",
    output_file = paste0("probatch_report",".html"),
    output_dir = TEMP_OUTDIR,
    intermediates_dir = TEMP_OUTDIR,
    envir = new.env(),
    quiet = TRUE,
    params = list(
      input_norm = args$input_norm,
      input_unnorm = args$input_unnorm,
      metadata_annotation = args$metadata_annotation,
      output_dir = args$output_dir,
      outfile_prefix = args$outfile_prefix,
      batch_column = args$batch_column,
      cols_of_interest = args$cols_of_interest,
      sample_threshold = args$sample_threshold,
      expgroup_threshold = args$expgroup_threshold,
      batch_threshold = args$batch_threshold,
      samples_for_correlation = args$samples_for_correlation,
      imputation_method = args$imputation_method,
      iRT_protein_name = args$iRT_protein_name,
      quantile_norm = args$quantile_norm,
      batch_correction_first = args$batch_correction_first,
      pdf_width = args$pdf_width,
      pdf_length = args$pdf_length
    )
  )
  
  output_files <- list.files(TEMP_OUTDIR, full.names=T)
  copy_success <- file.copy(output_files, to=outpath, overwrite=T, recursive=T, copy.date=T)
  
  if(!all(copy_success)) {
    files_failed <- output_files[!copy_success]
    stop(paste(
      paste0("FAILED TO COPY THE FOLLOWING TEMP FILES TO ", outpath, ":"),
      files_failed,
      sep = "\n"
    ))
  }
}, error = function(e) {
  write(paste("Failed to render HTML report.", e, sep="\n"), stderr())
}, finally = {
  unlink(TEMP_OUTDIR, recursive=T)
})
