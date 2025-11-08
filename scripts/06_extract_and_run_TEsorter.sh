#!/bin/bash
#SBATCH --job-name=TEsorter_CG
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
TELIB="$OUTDIR/assembly.p_ctg.fa.mod.EDTA.TElib.fa"
CONTAINER_TS="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"

module load SeqKit/2.3.1 || true

mkdir -p "$OUTDIR/te_clades"
cd "$OUTDIR/te_clades"


if [ ! -s Copia_sequences.fa ]; then
  seqkit grep -r -p "Copia" "$TELIB" > Copia_sequences.fa || true
fi
if [ ! -s Gypsy_sequences.fa ]; then
  seqkit grep -r -p "Gypsy" "$TELIB" > Gypsy_sequences.fa || true
fi


if [ ! -s Copia_sequences.fa ] && [ ! -s Gypsy_sequences.fa ]; then
  echo "No Copia or Gypsy sequences found in $TELIB"
  exit 0
fi

if [ -s Copia_sequences.fa ]; then
  apptainer exec --bind "$WORKDIR" "$CONTAINER_TS" \
    TEsorter Copia_sequences.fa -db rexdb-plant
fi
if [ -s Gypsy_sequences.fa ]; then
  apptainer exec --bind "$WORKDIR" "$CONTAINER_TS" \
    TEsorter Gypsy_sequences.fa -db rexdb-plant
fi

echo "TEsorter finished. Key outputs: *.rexdb-plant.cls.tsv and *.rexdb-plant.dom.faa"