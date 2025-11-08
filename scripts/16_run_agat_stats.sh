#!/bin/bash
#SBATCH --job-name=agat_stats
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/16_agat_stats_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/16_agat_stats_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch



set -euo pipefail

WRK="/data/users/yliu2/Organization_and_annotation/gene_annotation/final"
GFF="$WRK/filtered.genes.renamed.gff3"
OUT_TXT="$WRK/annotation_statistics.txt"
OUT_TSV="$WRK/annotation_statistics.tsv"

echo "[INFO] Node: $(hostname)"
echo "[INFO] GFF: $GFF"
test -s "$GFF" || { echo "[ERROR] Missing or empty GFF file: $GFF" >&2; exit 2; }


# Try to locate the AGAT Singularity image automatically
CAND_DIR="/data/courses/assembly-annotation-course/CDS_annotation/containers"
SIF=$(ls -1 "$CAND_DIR"/agat*.sif 2>/dev/null | head -n1 || true)

# ------------------------------------------------------------
# Function A: Run AGAT if available
# ------------------------------------------------------------
run_agat() {
  echo "[INFO] Using AGAT image: $SIF"
  apptainer exec --bind /data/users,/data/courses "$SIF" \
    agat_sp_statistics.pl -i "$GFF" -o "$OUT_TXT"

  
  grep -E 'Number of genes|Number of mRNA|Number of exons|Number of CDS' "$OUT_TXT" \
    | sed -E 's/^[[:space:]]+//; s/[[:space:]]+/ /g; s/: /\t/' > "$OUT_TSV" || true

  echo "[OK] AGAT statistics written to $OUT_TXT"
}

# ------------------------------------------------------------
# Function B: Fallback if AGAT is missing
# ------------------------------------------------------------
run_fallback() {
  echo "[WARN] No AGAT image found. Running fallback statistics..."

  genes=$(awk -F'\t' '$3=="gene"{c++} END{print c+0}' "$GFF")
  mrna=$(awk -F'\t' '$3=="mRNA"{c++} END{print c+0}' "$GFF")
  exons=$(awk -F'\t' '$3=="exon"{c++} END{print c+0}' "$GFF")
  cds=$(awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$GFF")

  {
    echo "=================================================="
    echo "Annotation Statistics (Fallback Mode)"
    echo "Timestamp: $(date -Iseconds)"
    echo "=================================================="
    echo
    echo "Counts:"
    printf "  Genes:  %d\n" "$genes"
    printf "  mRNA:   %d\n" "$mrna"
    printf "  Exons:  %d\n" "$exons"
    printf "  CDS:    %d\n" "$cds"
    echo
    echo "Top 10 scaffolds by mRNA count:"
    awk -F'\t' '$3=="mRNA"{c[$1]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' "$GFF" \
      | sort -k2,2nr | head -10 | awk '{printf "  %-12s %d\n",$1,$2}'
    echo
    echo "=================================================="
  } | tee "$OUT_TXT"

  {
    echo -e "metric\tvalue"
    echo -e "genes\t$genes"
    echo -e "mRNA\t$mrna"
    echo -e "exons\t$exons"
    echo -e "CDS\t$cds"
  } > "$OUT_TSV"

  echo "[OK] Fallback stats written to $OUT_TXT"
}

# ------------------------------------------------------------
# Execution logic
# ------------------------------------------------------------
if [ -n "${SIF:-}" ] && [ -s "$SIF" ]; then
  run_agat
else
  run_fallback
fi

echo "[DONE] Annotation statistics completed successfully."