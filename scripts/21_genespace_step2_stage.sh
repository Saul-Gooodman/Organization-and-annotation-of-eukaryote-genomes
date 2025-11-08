#!/usr/bin/env bash
#SBATCH --job-name=gs_fix_yuwei_v2
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail

BASE="/data/users/yliu2/Organization_and_annotation"
FINAL="$BASE/gene_annotation/final"
WD="$BASE/genespace"

PROT_SRC="$FINAL/assembly_p_ctg.longest.proteins.fasta"
[ -s "$PROT_SRC" ] || PROT_SRC="$FINAL/assembly_p_ctg.all.maker.proteins.renamed.filtered.fasta"

FAA_OUT="$WD/peptide/Yuwei.fa"

echo "[INFO] Source protein file: $PROT_SRC"
echo "[INFO] Output file: $FAA_OUT"

# 1) Clean headers (remove isoform suffixes, replace ':' with '_')
awk '
  BEGIN{RS=">"; ORS=""}
  NR>1{
    split($1, head, "\n")
    id=head[1]
    sub(/-R[A-Z]+$/, "", id)
    gsub(":", "_", id)
    seq=""
    for (i=2; i<=length(head); i++) seq=seq head[i]
    seq=gensub(/\n/,"","g",seq)
    if (seq != "") print ">" id "\n" seq "\n"
  }
' "$PROT_SRC" > "$FAA_OUT"

# 2) Report stats
records=$(grep -c '^>' "$FAA_OUT")
echo "[INFO] Wrote $records protein records to $FAA_OUT"

# 3) Sanity check for amino acid characters
grep -v '^>' "$FAA_OUT" | head -5