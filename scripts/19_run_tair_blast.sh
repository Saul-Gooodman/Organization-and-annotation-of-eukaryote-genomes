#!/usr/bin/env bash
#SBATCH --job-name=tair_blast
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"

# 1) Paths
BASE="/data/users/yliu2/Organization_and_annotation"
FINAL="$BASE/gene_annotation/final"

# Query protein set (filtered + renamed from previous step)
QUERY="$FINAL/assembly_p_ctg.all.maker.proteins.renamed.filtered.fasta"

# Use the pre-built TAIR10 representative peptide BLAST DB (no extension!)
TAIR_DB_BASE="/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_pep_20110103_representative_gene_model"

echo "[INFO] QUERY=$QUERY"
echo "[INFO] TAIR_DB_BASE=$TAIR_DB_BASE"

# 2) Modules
module load BLAST+/2.15.0-gompi-2021a

# 3) Sanity checks
if [ ! -s "$QUERY" ]; then
  echo "[FATAL] Missing or empty query FASTA: $QUERY" >&2
  exit 2
fi
# BLAST DB presence check (.pin/.psq files)
for ext in pin psq phr; do
  test -s "${TAIR_DB_BASE}.${ext}" || { echo "[FATAL] Missing BLAST DB file: ${TAIR_DB_BASE}.${ext}" >&2; exit 3; }
done

# 4) Run BLASTP against TAIR10
OUT_BASENAME="$FINAL/tair_blastp"
BLAST_RAW="${OUT_BASENAME}.out"
BESTHITS="${OUT_BASENAME}.out.besthits"
SUMMARY="${OUT_BASENAME}_summary.txt"
FLC_HITS="${OUT_BASENAME}_FLC_hits.txt"

echo "[INFO] Running BLASTP ..."
blastp \
  -query "$QUERY" \
  -db "$TAIR_DB_BASE" \
  -evalue 1e-5 \
  -num_threads "${SLURM_CPUS_PER_TASK:-12}" \
  -max_target_seqs 10 \
  -outfmt 6 \
  -out "$BLAST_RAW"

# 5) Keep best hit per query (sort by qseqid then evalue asc; evalue=col 11 for outfmt 6 std)
sort -k1,1 -k11,11g "$BLAST_RAW" | sort -u -k1,1 --merge > "$BESTHITS"

# 6) Quick stats
total_prot=$(grep -c '^>' "$QUERY" || echo 0)
with_hits=$(wc -l < "$BESTHITS" || echo 0)
pct=$(awk -v a="$with_hits" -v b="$total_prot" 'BEGIN{if(b>0) printf("%.2f",100*a/b); else print "0.00"}')

# 7) Find FLC (AT5G10140) best hits
grep -i 'AT5G10140' "$BESTHITS" > "$FLC_HITS" || true

# 8) Summary
{
  echo "TAIR10 BLAST summary"
  echo "Timestamp: $(date -Iseconds)"
  echo "Query proteins: $QUERY"
  echo "TAIR DB base:  $TAIR_DB_BASE"
  echo "-------------------------------------------"
  echo "Total proteins: $total_prot"
  echo "With best hits: $with_hits (${pct}%)"
  echo "Besthits file:  $BESTHITS"
  echo
  echo "FLC (AT5G10140) hits:"
  if [ -s "$FLC_HITS" ]; then
    cat "$FLC_HITS"
  else
    echo "None found in besthits."
  fi
} | tee "$SUMMARY"

echo "[INFO] Done. Summary at: $SUMMARY"
echo "[INFO] End=$(date -Iseconds)"