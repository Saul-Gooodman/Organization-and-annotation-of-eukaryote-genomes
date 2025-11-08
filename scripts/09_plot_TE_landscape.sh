#!/bin/bash
#SBATCH --job-name=TE_lands
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail

module load R/4.3.2-foss-2021a

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
LAND="$OUTDIR/landscape"
PLOTDIR="$OUTDIR/plots"
mkdir -p "$PLOTDIR"

BASE="$LAND/genome.mod.out"   # base name, not a directory

if [[ ! -s "${BASE}.divsum" || ! -s "${BASE}.summary" ]]; then
  echo "Missing ${BASE}.divsum or ${BASE}.summary" >&2
  exit 2
fi

Rscript /data/courses/assembly-annotation-course/CDS_annotation/scripts/plot_div.R \
  -i "$BASE" -o "$PLOTDIR"