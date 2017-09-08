#!/usr/bin/env Rscript

"Usage: getCounts.Rscript -b <bamfile> -i <intfile> -r <refname> -o <outfile> 

Options: 
  -h --help
  -b <bamfile> --bamfile <bamfile>  the path to the bamfile for processing
  -i <intfile> --intfile <intfile>  the path to the interval file
  -r <refname> --refname <refname>  the reference to process, e.g. the 
                                    chromosome 
  -o <outfile> --outfile <outfile>  the output file path & name
" -> doc

require(cnvR)
require(docopt); opt <- docopt::docopt(doc)

cnvGetCounts(bamfile = opt[["--bamfile"]], intfile = opt[["--intfile"]],
             refname = opt[["--refname"]], outfile = opt[["--outfile"]],
             results = FALSE, verbose = TRUE)
