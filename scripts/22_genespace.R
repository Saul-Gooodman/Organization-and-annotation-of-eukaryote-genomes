suppressPackageStartupMessages({
  library(GENESPACE)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide working directory (wd) as arg[1].")
wd <- args[1]
pepDir <- file.path(wd, "peptide")
bedDir <- file.path(wd, "bed")

pepFiles <- list.files(pepDir, pattern="\\.fa(sta)?$", full.names=TRUE)
if (length(pepFiles) < 2) stop("Need >=2 genomes (peptide FASTA) to run GENESPACE.")

genomeIDs <- sub("\\.fa(sta)?$", "", basename(pepFiles))
df <- data.frame(
  genome = genomeIDs,
  pep    = file.path(pepDir, paste0(genomeIDs, ".fa")),
  bed    = file.path(bedDir,  paste0(genomeIDs, ".bed")),
  stringsAsFactors = FALSE
)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/usr/local/bin/MCScanX"
)

out <- run_genespace(
  gpar,
  overwrite = TRUE
)

refGenome <- if ("TAIR10" %in% df$genome) "TAIR10" else df$genome[1]

pan <- query_pangenes(
  out,
  bed = NULL,
  refGenome = refGenome,
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf
)
saveRDS(pan, file = file.path(wd, "pangenome_matrix.rds"))
cat("[OK] GENESPACE finished.\n")
cat("[OK] Pangenome matrix: ", file.path(wd, "pangenome_matrix.rds"), "\n", sep="")
