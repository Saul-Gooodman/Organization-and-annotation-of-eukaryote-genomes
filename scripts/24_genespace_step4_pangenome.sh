#!/usr/bin/env bash
#SBATCH --job-name=genespace_summary
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

set -euo pipefail

BASE="/data/users/yliu2/Organization_and_annotation"
WD="$BASE/genespace"
SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/genespace_latest.sif"

IN_RDS="$WD/pangenes/pangenome_matrix.rds"
OUT_TXT="$WD/pangenes/core_accessory_summary.txt"

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"
echo "[INFO] IN_RDS=$IN_RDS"
echo "[INFO] OUT_TXT=$OUT_TXT"

apptainer exec --cleanenv \
  --env LANG=C.UTF-8,LC_ALL=C.UTF-8 \
  --env-unset LMOD_CMD \
  --env-unset MODULEPATH \
  --env-unset MODULESHOME \
  --bind "$BASE":/work \
  "$SIF" bash -c "
    Rscript - << 'RS'
    suppressPackageStartupMessages(library(GENESPACE))

    pan <- readRDS('/work/genespace/pangenes/pangenome_matrix.rds')
    nm <- names(pan)
    if ('pangeneMat' %in% nm) {
      pam <- pan\$pangeneMat
    } else if ('ogByGenome' %in% nm) {
      pam <- pan\$ogByGenome
    } else {
      stop('Cannot find pangene matrix. Available names:', paste(nm, collapse=', '))
    }

    req <- c('TAIR10', 'Yuwei')
    if (!all(req %in% colnames(pam))) {
      stop('Missing TAIR10 or Yuwei in matrix. Have:', paste(colnames(pam), collapse=', '))
    }

    core_ogs       <- rownames(pam)[pam[,'TAIR10']==1 & pam[,'Yuwei']==1]
    tair_only_ogs  <- rownames(pam)[pam[,'TAIR10']==1 & pam[,'Yuwei']==0]
    yuwei_only_ogs <- rownames(pam)[pam[,'TAIR10']==0 & pam[,'Yuwei']==1]

    og2genes <- NULL
    if ('og2Gene' %in% names(pan)) og2genes <- pan\$og2Gene
    if (is.null(og2genes)) {
      maybe <- c('og2gene', 'ogToGene', 'ogMembers')
      for (m in maybe) if (m %in% names(pan)) og2genes <- pan[[m]]
    }

    out <- file('/work/genespace/pangenes/core_accessory_summary.txt', 'wt')
    writeLines('GENESPACE core/accessory summary (TAIR10 vs Yuwei)', out)
    writeLines(sprintf('Date: %s', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), out)
    writeLines('', out)
    writeLines(sprintf('Core orthogroups: %d', length(core_ogs)), out)
    writeLines(sprintf('TAIR10-unique orthogroups: %d', length(tair_only_ogs)), out)
    writeLines(sprintf('Yuwei-unique orthogroups: %d', length(yuwei_only_ogs)), out)
    writeLines('', out)

    if (!is.null(og2genes)) {
      toCount <- function(ogs, genome=NULL) {
        x <- og2genes[og2genes\$og %in% ogs, , drop=FALSE]
        if (!is.null(genome)) x <- x[x\$genome==genome,,drop=FALSE]
        nrow(unique(x))
      }
      writeLines('Gene-level counts:', out)
      writeLines(sprintf('Core genes (both): %d', toCount(core_ogs)), out)
      writeLines(sprintf('TAIR10-only genes: %d', toCount(tair_only_ogs,'TAIR10')), out)
      writeLines(sprintf('Yuwei-only genes: %d', toCount(yuwei_only_ogs,'Yuwei')), out)
    } else {
      writeLines('No gene-level mapping found.', out)
    }

    close(out)
    cat('[OK] Summary saved to', '/work/genespace/pangenes/core_accessory_summary.txt', '\\n')
RS
  "

echo "[INFO] End=$(date -Iseconds)"
echo "[DONE] Summary generation complete."