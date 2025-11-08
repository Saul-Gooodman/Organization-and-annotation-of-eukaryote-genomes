#!/bin/bash
#SBATCH --job-name=filter_refine
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ----- paths -----
WORKDIR="/data/users/yliu2/Organization_and_annotation"
ANNDIR="$WORKDIR/gene_annotation"
FINAL="$ANNDIR/final"
LOGDIR="$WORKDIR/logs"
mkdir -p "$FINAL" "$LOGDIR"

# MAKER container and InterProScan container
MAKER_SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"
IPR_SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/interproscan_latest.sif"
IPR_DATA="/data/courses/assembly-annotation-course/CDS_annotation/data/interproscan-5.70-102.0/data"

# Outputs produced in step 5
PREFIX="assembly_p_ctg"
GFF_NOSEQ="$ANNDIR/${PREFIX}.all.maker.noseq.gff"
PROT_FASTA="$ANNDIR/${PREFIX}.all.maker.proteins.fasta"
TRAN_FASTA="$ANNDIR/${PREFIX}.all.maker.transcripts.fasta"

# UCSC utils for faSomeRecords
module load UCSC-Utils/448-foss-2021a || true

# A helper to run a command inside MAKER container with necessary binds
run_maker() {
  apptainer exec --bind /data/users,/data/courses "$MAKER_SIF" "$@"
}

# ----- sanity checks -----
for f in "$GFF_NOSEQ" "$PROT_FASTA" "$TRAN_FASTA"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: required file not found or empty: $f" >&2
    exit 2
  fi
done

echo "[1] Prepare files and rename IDs ..."
cp -f "$GFF_NOSEQ"   "$FINAL/${PREFIX}.all.maker.noseq.renamed.gff"
cp -f "$PROT_FASTA"  "$FINAL/${PREFIX}.all.maker.proteins.renamed.fasta"
cp -f "$TRAN_FASTA"  "$FINAL/${PREFIX}.all.maker.transcripts.renamed.fasta"

cd "$FINAL"

# Choose a clean 3–4 letter prefix for IDs
ACC_PREFIX="YWL"   # change if you prefer another accession code

# Make ID map and apply to GFF/FASTA
run_maker maker_map_ids --prefix "$ACC_PREFIX" --justify 7 \
  "${PREFIX}.all.maker.noseq.renamed.gff" > id.map

run_maker map_gff_ids   id.map "${PREFIX}.all.maker.noseq.renamed.gff"
run_maker map_fasta_ids id.map "${PREFIX}.all.maker.proteins.renamed.fasta"
run_maker map_fasta_ids id.map "${PREFIX}.all.maker.transcripts.renamed.fasta"

echo "[2] Run InterProScan (Pfam) ..."
# InterProScan output
IPR_OUT="output.iprscan.tsv"

apptainer exec \
  --bind "$IPR_DATA":/opt/interproscan/data \
  --bind "$WORKDIR" \
  --bind /data/courses \
  --bind "$SCRATCH":/temp \
  "$IPR_SIF" \
  /opt/interproscan/interproscan.sh \
  -appl pfam --disable-precalc -f TSV \
  --goterms --iprlookup --seqtype p \
  -i "${PREFIX}.all.maker.proteins.renamed.fasta" \
  -o "$IPR_OUT"

echo "[3] Update GFF with InterProScan results ..."
run_maker ipr_update_gff \
  "${PREFIX}.all.maker.noseq.renamed.gff" "$IPR_OUT" \
  > "${PREFIX}.all.maker.noseq.renamed.iprscan.gff"

echo "[4] Compute AED distribution ..."
run_maker perl "$(run_maker which AED_cdf_generator.pl)" -b 0.025 \
  "${PREFIX}.all.maker.noseq.renamed.gff" \
  > "${PREFIX}.all.maker.renamed.gff.AED.txt"

# Quick AED summary
awk 'NF==2{n++; s+=$2; if($2<=0.5) c05++} END{printf("AED lines:%d  mean=%.3f  <=0.5:%d (%.1f%%)\n", n, (n?n? s/n:0:0), c05+0, (n?100*c05/n:0))}' \
  "${PREFIX}.all.maker.renamed.gff.AED.txt" \
  | tee -a "${PREFIX}.filter_refine.log"

echo "[5] Quality filtering (AED<1 and/or Pfam present) ..."
run_maker perl "$(run_maker which quality_filter.pl)" \
  -s "${PREFIX}.all.maker.noseq.renamed.iprscan.gff" \
  > "${PREFIX}.all.maker.noseq.renamed_iprscan_quality_filtered.gff"

echo "[6] Keep core gene features (gene/mRNA/exon/CDS/UTRs) ..."
grep -P "\t(gene|mRNA|exon|CDS|five_prime_UTR|three_prime_UTR)\t" \
  "${PREFIX}.all.maker.noseq.renamed_iprscan_quality_filtered.gff" \
  > "${PREFIX}.filtered.genes.renamed.gff3"

# Sanity check of feature types
cut -f3 "${PREFIX}.filtered.genes.renamed.gff3" | sort | uniq -c | sort -nr \
  | tee "${PREFIX}.filtered.genes.feature_counts.txt"

echo "[7] Extract surviving mRNA IDs and filter FASTAs ..."
grep -P "\tmRNA\t" "${PREFIX}.filtered.genes.renamed.gff3" \
  | awk -F'\t' '{print $9}' | sed -E 's/.*ID=([^;]+).*/\1/' > list.txt

# transcripts
faSomeRecords "${PREFIX}.all.maker.transcripts.renamed.fasta" list.txt \
  > "${PREFIX}.all.maker.transcripts.renamed.filtered.fasta"

# proteins
faSomeRecords "${PREFIX}.all.maker.proteins.renamed.fasta"   list.txt \
  > "${PREFIX}.all.maker.proteins.renamed.filtered.fasta"

# Brief stats
echo "[Summary] Counts after filtering:" | tee -a "${PREFIX}.filter_refine.log"
echo "  mRNA kept: $(wc -l < list.txt)"     | tee -a "${PREFIX}.filter_refine.log"
echo "  Transcripts FASTA: $(grep -c '^>' ${PREFIX}.all.maker.transcripts.renamed.filtered.fasta)" | tee -a "${PREFIX}.filter_refine.log"
echo "  Proteins   FASTA: $(grep -c '^>' ${PREFIX}.all.maker.proteins.renamed.filtered.fasta)"     | tee -a "${PREFIX}.filter_refine.log"

echo "Done. Outputs in: $FINAL"