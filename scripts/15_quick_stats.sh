#!/usr/bin/env bash
set -euo pipefail

FINAL="${FINAL:-/data/users/yliu2/Organization_and_annotation/gene_annotation/final}"
OUTTXT="${OUTTXT:-$FINAL/annotation_quick_summary.txt}"

GFF="$FINAL/filtered.genes.renamed.gff3"
PROT="$FINAL/assembly_p_ctg.all.maker.proteins.renamed.filtered.fasta"
TRAN="$FINAL/assembly_p_ctg.all.maker.transcripts.renamed.filtered.fasta"

# BUSCO short summary (proteins)
BUSCO_DIR=$(ls -d "$FINAL"/busco_proteins_* 2>/dev/null | sort | tail -n1 || true)
BUSCO_SUMMARY="${BUSCO_DIR}/run_brassicales_odb10/short_summary.txt"

for f in "$GFF" "$PROT" "$TRAN" "$BUSCO_SUMMARY"; do
  if [ ! -s "$f" ]; then
    echo "[ERROR] Missing or empty file: $f" >&2
    exit 13
  fi
done

# Counts
genes=$(awk -F'\t' '$3=="gene"{c++} END{print c+0}' "$GFF")
mrna=$(awk -F'\t' '$3=="mRNA"{c++} END{print c+0}' "$GFF")
cds=$(awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$GFF")
exon=$(awk -F'\t' '$3=="exon"{c++} END{print c+0}' "$GFF")
utr5=$(awk -F'\t' '$3=="five_prime_UTR"{c++} END{print c+0}' "$GFF")
utr3=$(awk -F'\t' '$3=="three_prime_UTR"{c++} END{print c+0}' "$GFF")

# AED CDF quick peek (first & last few lines if present)
AEDTXT="$FINAL/assembly_p_ctg.all.maker.noseq.renamed.gff.AED.txt"
aed_head="N/A"; aed_tail="N/A"
if [ -s "$AEDTXT" ]; then
  aed_head=$(head -n 6 "$AEDTXT")
  aed_tail=$(tail -n 6 "$AEDTXT")
fi

{
  echo "MAKER merged outputs summary (quick)"
  date -Iseconds
  echo "Final dir: $FINAL"
  echo
  echo "Files:"
  printf "  %-50s %s\n" "$(basename "$GFF")" "$(wc -l < "$GFF") lines"
  printf "  %-50s %s\n" "$(basename "$PROT")" "$(grep -c '^>' "$PROT") records"
  printf "  %-50s %s\n" "$(basename "$TRAN")" "$(grep -c '^>' "$TRAN") records"
  echo
  echo "Feature counts from GFF:"
  printf "  gene\t%d\n" "$genes"
  printf "  mRNA\t%d\n" "$mrna"
  printf "  CDS\t%d\n" "$cds"
  printf "  exon\t%d\n" "$exon"
  printf "  five_prime_UTR\t%d\n" "$utr5"
  printf "  three_prime_UTR\t%d\n" "$utr3"
  echo
  echo "BUSCO short summary:"
  echo "  file: $BUSCO_SUMMARY"
  sed -n '1,200p' "$BUSCO_SUMMARY"
  echo
  echo "AED CDF (head):"
  echo "$aed_head"
  echo "AED CDF (tail):"
  echo "$aed_tail"
} | tee "$OUTTXT"

echo "[OK] Wrote $OUTTXT"
