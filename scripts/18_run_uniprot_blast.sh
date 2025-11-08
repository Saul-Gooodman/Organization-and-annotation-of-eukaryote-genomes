#!/bin/bash
#SBATCH --job-name=uniprot_blast_map
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=10
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail

# ------------------ paths ------------------
BASE="/data/users/yliu2/Organization_and_annotation"
FINAL="$BASE/gene_annotation/final"
LOGS="$BASE/logs"
mkdir -p "$LOGS"
cd "$FINAL"

# Inputs produced in previous steps
PROT_IN="$FINAL/assembly_p_ctg.all.maker.proteins.renamed.filtered.fasta"
GFF_IN="$FINAL/filtered.genes.renamed.gff3"

# UniProt (reviewed) DB (already formatted on the course share)
DB_UNIPROT="/data/courses/assembly-annotation-course/CDS_annotation/data/uniprot/uniprot_viridiplantae_reviewed.fa"

# Outputs
BLAST_RAW="$FINAL/uniprot_blastp.out"
BLAST_BEST="$FINAL/uniprot_blastp.out.besthits"
PROT_OUT="$FINAL/maker_proteins.filtered.fasta.Uniprot"
GFF_OUT="$FINAL/filtered.genes.renamed.gff3.Uniprot.gff3"

# MAKER container for functional mapping
SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"
RUN="apptainer exec --bind /data/users,/data/courses $SIF"

echo "[INFO] node=$(hostname)"
echo "[INFO] FINAL=$FINAL"
echo "[INFO] PROT_IN=$PROT_IN"
echo "[INFO] GFF_IN=$GFF_IN"
echo "[INFO] DB_UNIPROT=$DB_UNIPROT"
echo "[INFO] OUT (BLAST) = $BLAST_RAW ; $BLAST_BEST"
echo "[INFO] OUT (mapped) = $PROT_OUT ; $GFF_OUT"

# ------------------ sanity checks ------------------
for f in "$PROT_IN" "$GFF_IN" "$DB_UNIPROT"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing or empty: $f" >&2; exit 2; }
done

# ------------------ BLASTP vs UniProt ------------------
module load BLAST+/2.15.0-gompi-2021a

echo "[STEP] Running BLASTP against UniProt reviewed (Viridiplantae) ..."
blastp \
  -query "$PROT_IN" \
  -db "$DB_UNIPROT" \
  -num_threads "${SLURM_CPUS_PER_TASK:-10}" \
  -outfmt 6 \
  -evalue 1e-5 \
  -max_target_seqs 10 \
  -out "$BLAST_RAW"

echo "[STEP] Selecting best hit per query ..."
# outfmt 6 'std' columns: ... 11=evalue, 12=bitscore
# Sort by: qseqid, evalue asc, bitscore desc
sort -k1,1 -k11,11g -k12,12nr "$BLAST_RAW" \
| awk ' !seen[$1]++ ' > "$BLAST_BEST"

echo "[INFO] BLAST raw lines: $(wc -l < "$BLAST_RAW")"
echo "[INFO] Unique queries with hits: $(awk '{print $1}' "$BLAST_BEST" | sort -u | wc -l)"

# ------------------ Map UniProt functions back to FASTA/GFF ------------------
echo "[STEP] Mapping UniProt functional annotations into FASTA/GFF via MAKER tools ..."
$RUN which maker_functional_fasta
$RUN which maker_functional_gff

$RUN maker_functional_fasta "$DB_UNIPROT" "$BLAST_BEST" "$PROT_IN" > "$PROT_OUT"
$RUN maker_functional_gff   "$DB_UNIPROT" "$BLAST_BEST" "$GFF_IN"  > "$GFF_OUT"

# ------------------ Quick summary ------------------
n_prot=$(grep -c '^>' "$PROT_IN" || echo 0)
n_hit=$(awk '{print $1}' "$BLAST_BEST" | sort -u | wc -l)
pct=$(python3 - <<PY
n=$n_prot
h=$n_hit
print(f"{(100*h/n):.2f}" if n>0 else "0.00")
PY
)

echo "==========================================="
echo "UniProt BLAST summary"
echo "Total proteins:  $n_prot"
echo "With best hits:  $n_hit (${pct}%)"
echo "Besthits file:   $BLAST_BEST"
echo "Mapped FASTA:    $PROT_OUT"
echo "Mapped GFF3:     $GFF_OUT"
echo "==========================================="

echo "[DONE] UniProt functional mapping completed."