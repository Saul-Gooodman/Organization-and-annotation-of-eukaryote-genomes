#!/bin/bash
#SBATCH --job-name=TE_landscape
#SBATCH --partition=pshort_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

set -euo pipefail

module load R/4.3.2-foss-2021a

BASE="/data/users/yliu2/Organization_and_annotation/results/EDTA_annotation"
LAND="$BASE/landscape"
PLOTS="$BASE/plots"
mkdir -p "$LAND/Plots"


ACT="genome.mod.out"
EXP="assembly.fasta.mod.out"


REAL_ANNO_DIR="$BASE/assembly.p_ctg.fa.mod.EDTA.anno"
EXP_ANNO_DIR="$LAND/assembly.fasta.mod.EDTA.anno"


mkdir -p "$LAND"
if [ -d "$REAL_ANNO_DIR" ]; then
  ln -sfn "$REAL_ANNO_DIR" "$EXP_ANNO_DIR"
fi


for suf in "Rclass.tab" "Rfam.tab" "Rname.tab"; do
  src="$LAND/${ACT}.landscape.Div.${suf}"
  dst="$EXP_ANNO_DIR/${EXP}.landscape.Div.${suf}"
  if [ -f "$src" ] && [ ! -f "$dst" ]; then
    ln -s "$src" "$dst"
  fi
done


[ -f "$LAND/${ACT}.landscape.Div.Rclass.tab" ] && ln -sfn "$LAND/${ACT}.landscape.Div.Rclass.tab" "$LAND/${ACT}.divsum"
[ -f "$LAND/${ACT}.landscape.Div.Rfam.tab" ]   && ln -sfn "$LAND/${ACT}.landscape.Div.Rfam.tab"   "$LAND/${ACT}.summary"


cd "$LAND"


Rscript /data/courses/assembly-annotation-course/CDS_annotation/scripts/06-plot_div.R


# Rscript /data/courses/assembly-annotation-course/CDS_annotation/scripts/06-plot_div.R \
#   -i "$LAND/$ACT" \
#   -o "$PLOTS"

echo "Done. Check plots in: $PLOTS and $LAND"