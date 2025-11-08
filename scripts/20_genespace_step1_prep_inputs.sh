#!/usr/bin/env bash
#SBATCH --job-name=genespace_prep
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"

#------------------------------------------------------------
# Paths and variables
#------------------------------------------------------------
BASE="/data/users/yliu2/Organization_and_annotation"
FINAL="$BASE/gene_annotation/final"
OUTDIR="$BASE/genespace"
ACC="Yuwei"

mkdir -p "$OUTDIR" "$BASE/logs" "$OUTDIR/bed" "$OUTDIR/peptide"

GFF="$FINAL/filtered.genes.renamed.gff3"
LONGEST_PROT="$FINAL/assembly_p_ctg.longest.proteins.fasta"

#------------------------------------------------------------
# Check inputs
#------------------------------------------------------------
if [ ! -s "$GFF" ]; then
  echo "[FATAL] GFF3 not found: $GFF" >&2
  exit 1
fi
if [ ! -s "$LONGEST_PROT" ]; then
  echo "[FATAL] Protein FASTA not found: $LONGEST_PROT" >&2
  exit 1
fi

echo "[INFO] Using GFF3: $GFF"
echo "[INFO] Using Proteins: $LONGEST_PROT"

#------------------------------------------------------------
# 1) Generate BED (0-based start)
#------------------------------------------------------------
grep -P "\tgene\t" "$GFF" > "$OUTDIR/temp_genes.gff3"

awk 'BEGIN{OFS="\t"} {
  split($9,a,";"); split(a[1],b,"=");
  print $1, $4-1, $5, b[2]
}' "$OUTDIR/temp_genes.gff3" > "$OUTDIR/bed/${ACC}.bed"

#------------------------------------------------------------
# 2) Copy longest protein FASTA
#------------------------------------------------------------
cp "$LONGEST_PROT" "$OUTDIR/peptide/${ACC}.fa"

#------------------------------------------------------------
# Cleanup and report
#------------------------------------------------------------
rm -f "$OUTDIR/temp_genes.gff3"

echo "[OK] Created:"
echo "  $OUTDIR/bed/${ACC}.bed"
echo "  $OUTDIR/peptide/${ACC}.fa"

echo "[INFO] End=$(date -Iseconds)"