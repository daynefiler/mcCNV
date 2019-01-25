library(data.table)
library(GenomicRanges)
library(ape)
library(rtracklayer)

vol5 <- "/Users/dayne/Mounted_Drives/vol5"
wdir <- file.path(vol5, "dfiler", "CNV")
bait <- fread(file.path(wdir, "otoscope", "otoscopev4.final.txt"), 
              select = 1:4)
setnames(bait, c("chrom", "start", "end", "bait"))

u <- file.path("ftp:/", "ftp.ncbi.nih.gov", "genomes", "refseq", 
               "vertebrate_mammalian", "Homo_sapiens", "all_assembly_versions",
               "GCF_000001405.25_GRCh37.p13", 
               "GCF_000001405.25_GRCh37.p13_genomic.gff.gz")
trck <- readGFF(u, 
                columns = c("seqid", "type", "start", "end", "strand"),
                tags = c("ID", "Dbxref", "Parent", "gene", "genome", 
                         "chromosome"))
trck <- as.data.table(trck)
chrMap <- trck[genome %in% c("chromosome", "mitochondrion"), 
               list(seqid, chromosome)]
## Force mito chromosome to M
chrMap[seqid == "NC_012920.1", chromosome := "M"]
chrMap[ , chrom := paste0("chr", chromosome)]
setkey(chrMap, seqid)
setkey(trck, seqid)
trck <- chrMap[ , list(seqid, chrom)][trck]

exon <- as.data.table(with(trck[type == "exon" & !is.na(chrom)],
                           reduce(GRanges(seqnames = gene,
                                          ranges = IRanges(start, end),
                                          strand = strand))))
setkey(exon, seqnames)
setkey(trck, gene)

geneMap <- unique(trck[!is.na(gene) & !is.na(chrom), list(gene, chrom, seqid)])
setkey(geneMap, gene)
geneMap <- geneMap[J(exon[ , unique(seqnames)])]
## Remove chrY from the map
geneMap <- geneMap[chrom != "chrY"]
exon <- exon[geneMap]
setorder(exon, chrom, start)
exon[ , exonName := sprintf("%s:%s", seqnames, LETTERS[1:.N]), by = seqnames]
erng <- with(exon, GRanges(seqid, IRanges(start, end, names = exonName)))

setkey(bait, chrom)
setkey(chrMap, chrom)
bait <- chrMap[bait]
brng <- reduce(with(bait, GRanges(seqid, IRanges(start, end))))
hits <- findOverlaps(brng, erng)
sbst <- erng[unique(subjectHits(hits))]

saveRDS(sbst, file = file.path(wdir, "realData37", "otoIntervals.RDS"))

## Old approach using UCSC
# library(DBI)
# library(RMySQL)
# tmpConf <- tempfile()
# cat(c("[ucsc]", 
#       "user=genome", 
#       "password=", 
#       "host=genome-mysql.soe.ucsc.edu", 
#       "port=3306", 
#       "no-auto-rehash"), 
#     file = tmpConf,
#     sep = "\n")
# con <- dbConnect(drv = MySQL(), groups = "ucsc", 
#                  dbname = "hg19", default.file = tmpConf)
# geneQ <- "SELECT name, chrom, exonCount, exonStarts, exonEnds FROM knownGene;"
# gene <- as.data.table(dbGetQuery(con, geneQ))
# symbQ <- "SELECT kgID as name, geneSymbol as gene FROM kgXref;"
# symb <- as.data.table(dbGetQuery(con, symbQ))
# chromAlias <- as.data.table(dbGetQuery(con, "SELECT * FROM chromAlias;"))
# dbDisconnect(con)
# 
# setkey(gene, name)
# setkey(symb, name)
# gene <- symb[gene]
# 
# pc <- function(x) {
#   eval(parse(text = sprintf("c(%sNULL)", paste(x, collapse = ""))))
# }
# setnames(gene, c("exonStarts", "exonEnds"), c("start", "end"))
# exon <- gene[ , list(start = pc(start), end = pc(end)), by = list(chrom, gene)]
# exon <- unique(exon)
# setorder(exon, chrom, start, end)
# exon[ , cg := paste(chrom, gc, sep = ":")]
# erng <- with(exon, GRanges(chrom, IRanges(start, end)))
# erng <- reduce(erng)
# erng <- as.data.table(erng)
# setkey(erng, seqnames, start)
# setkey(exon, chrom, start)
# exon <- exon[erng[ , list(seqnames, start)]]
# erng <- with(exon, GRanges(chrom, IRanges(start, end, names = gene)))
# 
# brng <- with(bait, GRanges(chrom, IRanges(start, end, names = bait)))
# 
# hits <- findOverlaps(brng, erng)
# sbst <- erng[subjectHits(hits)]
# trgt <- as.data.table(sbst)
# trgt[ , gene := names(sbst)]
# setkey(trgt, seqnames)
# setkey(chromAlias, chrom)
# trgt <- chromAlias[source == "refseq", list(alias, chrom)][trgt]


