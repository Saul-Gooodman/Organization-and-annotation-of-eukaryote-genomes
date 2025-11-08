#!/bin/bash
#SBATCH --job-name=circlize_TE
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail

module load R/4.3.2-foss-2021a
export R_LIBS_USER="/data/users/yliu2/Rlibs/4.3"
mkdir -p "$R_LIBS_USER"

Rscript -e 'repos <- c(CRAN="https://cloud.r-project.org");
            pkgs <- c("optparse","data.table","dplyr","stringr","circlize");
            inst <- rownames(installed.packages(lib.loc=Sys.getenv("R_LIBS_USER")));
            need <- setdiff(pkgs, inst);
            if(length(need)) install.packages(need, repos=repos, lib=Sys.getenv("R_LIBS_USER"), Ncpus=2)'

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
GFF="$OUTDIR/assembly.p_ctg.fa.mod.EDTA.TEanno.gff3"
FAI="$OUTDIR/assembly.p_ctg.fa.fai"
CIRCDIR="$OUTDIR/circos"
mkdir -p "$CIRCDIR"


if [[ ! -s "$FAI" ]]; then
  module spider SAMtools >/dev/null 2>&1 || true
  module load SAMtools || true
  if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$OUTDIR/assembly.p_ctg.fa"
  else
    echo "WARNING: samtools not available; please create $FAI manually."
    exit 1
  fi
fi

echo ">>> GFF: $GFF"
echo ">>> FAI: $FAI"
echo ">>> OUT: $CIRCDIR"

Rscript "$WORKDIR/scripts/05_circlize_TE_density.R" \
  --gff "$GFF" \
  --fai "$FAI" \
  --outdir "$CIRCDIR" \
  --prefix "assembly_p_ctg" \
  --scaf_top 10 \
  --win 100000 \
  --super "LTR/Gypsy,LTR/Copia,DNA/TIR,LINE,SINE"

echo ">>> Done."