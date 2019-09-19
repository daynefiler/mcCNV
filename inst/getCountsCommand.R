#!/usr/bin/env Rscript

"Usage: getCountsCommand.R -b <bamfile> -i <intfile> -o <outfile> 

Options:
  -h --help
  -b <bamfile> --bamfile <bamfile>  the path to the bamfile for processing
  -i <intfile> --intfile <intfile>  the path to the interval object
  -o <outfile> --outfile <outfile>  the output file path & name
" -> doc

require(cnvR, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

cnvGetCounts(bamfile = opt[["--bamfile"]],
             int = readRDS(opt[["--intfile"]]),
             outfile = opt[["--outfile"]],
             results = FALSE, verbose = TRUE)
