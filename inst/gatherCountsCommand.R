#!/usr/bin/env Rscript

"Usage: gatherCountsCommand.R -o <outfile> <counts>... 

Options:
  -h --help
  -o <outfile> --outfile <outfile>  the output file path & name
  <counts>  the path(s) to the count files for combining
" -> doc

require(cnvR, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

saveRDS(opt, "testOpt.RDS")

cnvGatherCounts(input = opt$counts,
                outfile = opt$outfile,
                results = FALSE)
