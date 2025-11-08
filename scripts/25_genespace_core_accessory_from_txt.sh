#!/usr/bin/env bash
#SBATCH --job-name=genespace_core_from_txt
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

set -euo pipefail

BASE="/data/users/yliu2/Organization_and_annotation"
WD="$BASE/genespace"
SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/genespace_latest.sif"

IN_FILE="$WD/pangenes/TAIR10_pangenes.txt.gz"
OUT_SUM="$WD/pangenes/core_accessory_summary.txt"
OUT_PAM="$WD/pangenes/pangene_presence_absence.csv"

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"
echo "[INFO] Using pangenes file: $IN_FILE"

if [ ! -s "$IN_FILE" ]; then
  echo "[FATAL] Missing: $IN_FILE" >&2
  exit 2
fi


apptainer exec --cleanenv \
  --env LANG=C.UTF-8,LC_ALL=C.UTF-8 \
  --bind "$BASE":/work \
  "$SIF" bash -c "
Rscript - << 'RS'
suppressPackageStartupMessages({
  library(data.table)
})

inFile <- '/work/genespace/pangenes/TAIR10_pangenes.txt.gz'
outSum <- '/work/genespace/pangenes/core_accessory_summary.txt'
outPAM <- '/work/genespace/pangenes/pangene_presence_absence.csv'

cat('[INFO] Reading:', inFile, '\\n')
dt <- suppressWarnings(fread(inFile))
nm <- tolower(names(dt))


col_og     <- names(dt)[match('og', nm)]
col_genome <- names(dt)[match('genome', nm)]
if (is.na(col_og) || is.na(col_genome)) {
  stop('Could not find required columns og/genome. Got: ', paste(names(dt), collapse=', '))
}


gene_candidates <- c('gene','id','pgrepid','pgid','ofid')
col_gene <- NA_character_
for (cand in gene_candidates) {
  hit <- which(nm == cand)
  if (length(hit) == 1L) { col_gene <- names(dt)[hit]; break }
}
if (is.na(col_gene)) {
  
  warning('No explicit gene ID column found; gene-level counts will be skipped.')
}


setnames(dt, c(col_og, col_genome), c('og','genome'))
if (!is.na(col_gene)) setnames(dt, col_gene, 'gene')


dt[, pres := 1L]
pam <- dcast(unique(dt[, .(og, genome, pres)]), og ~ genome, value.var='pres', fill=0)


req <- c('TAIR10','Yuwei')
miss <- setdiff(req, names(pam))
if (length(miss) > 0) {
  stop('Missing genomes in matrix: ', paste(miss, collapse=', '),
       '. Available: ', paste(names(pam), collapse=', '))
}

core_ogs       <- pam[TAIR10==1 & Yuwei==1, og]
tair_only_ogs  <- pam[TAIR10==1 & Yuwei==0, og]
yuwei_only_ogs <- pam[TAIR10==0 & Yuwei==1, og]


n_core_genes <- n_tair_only_genes <- n_yuwei_only_genes <- NA_integer_
if ('gene' %in% names(dt)) {
  core_genes       <- unique(dt[og %in% core_ogs])
  tair_only_genes  <- unique(dt[og %in% tair_only_ogs  & genome=='TAIR10'])
  yuwei_only_genes <- unique(dt[og %in% yuwei_only_ogs & genome=='Yuwei'])
  n_core_genes        <- nrow(core_genes)
  n_tair_only_genes   <- nrow(tair_only_genes)
  n_yuwei_only_genes  <- nrow(yuwei_only_genes)
}


con <- file(outSum, open='wt')
writeLines('GENESPACE core/accessory summary (TAIR10 vs Yuwei)', con)
writeLines(sprintf('Date: %s', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), con)
writeLines('', con)
writeLines(sprintf('Core orthogroups: %d', length(core_ogs)), con)
writeLines(sprintf('TAIR10-unique orthogroups: %d', length(tair_only_ogs)), con)
writeLines(sprintf('Yuwei-unique orthogroups: %d', length(yuwei_only_ogs)), con)
writeLines('', con)
if (!is.na(n_core_genes)) {
  writeLines('Gene-level counts:', con)
  writeLines(sprintf('Core genes (both): %d', n_core_genes), con)
  writeLines(sprintf('TAIR10-only genes: %d', n_tair_only_genes), con)
  writeLines(sprintf('Yuwei-only genes: %d', n_yuwei_only_genes), con)
}
close(con)


fwrite(pam, outPAM)
cat('[OK] Summary  -> ', outSum, '\\n', sep='')
cat('[OK] Matrix   -> ', outPAM, '\\n', sep='')
RS
"

echo "[INFO] End=$(date -Iseconds)"
echo "[DONE] Core/accessory summary generated."