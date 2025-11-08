#!/usr/bin/env bash
#SBATCH --job-name=summarize_maker
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ------------------- Paths -------------------
WORKDIR="/data/users/yliu2/Organization_and_annotation"
ANNDIR="${WORKDIR}/gene_annotation"
OUT_PREFIX="assembly_p_ctg"
DATASTORE_INDEX="${ANNDIR}/${OUT_PREFIX}.maker.output/${OUT_PREFIX}_master_datastore_index.log"

# Maker helper binaries
MAKERBIN="/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin"

# Outputs
GFF_ALL="${ANNDIR}/${OUT_PREFIX}.all.maker.gff"
GFF_NOSEQ="${ANNDIR}/${OUT_PREFIX}.all.maker.noseq.gff"
FA_TXT="${ANNDIR}/${OUT_PREFIX}.all.maker.transcripts.fasta"
FA_PRO="${ANNDIR}/${OUT_PREFIX}.all.maker.proteins.fasta"
SUMMARY="${ANNDIR}/${OUT_PREFIX}.maker.summary.txt"

mkdir -p "${WORKDIR}/logs" "${ANNDIR}"

have_cmd() { command -v "$1" >/dev/null 2>&1; }
die() { echo "ERROR: $*" >&2; exit 2; }
nonempty_or_die() { [[ -s "$1" ]] || die "Expected non-empty file: $1"; }

[[ -f "$DATASTORE_INDEX" ]] || die "Datastore index not found: $DATASTORE_INDEX"

# Resolve gff3_merge/fasta_merge
if [[ -x "${MAKERBIN}/gff3_merge" && -x "${MAKERBIN}/fasta_merge" ]]; then
  GFF3_MERGE="${MAKERBIN}/gff3_merge"
  FASTA_MERGE="${MAKERBIN}/fasta_merge"
elif have_cmd gff3_merge && have_cmd fasta_merge; then
  GFF3_MERGE="gff3_merge"
  FASTA_MERGE="fasta_merge"
else
  die "Cannot find gff3_merge/fasta_merge. Check MAKERBIN or PATH."
fi

echo "Using:"
echo "  GFF3_MERGE = ${GFF3_MERGE}"
echo "  FASTA_MERGE = ${FASTA_MERGE}"
echo "  INDEX = ${DATASTORE_INDEX}"
echo

# 1) Merge GFF with embedded sequences
echo "[1/5] Merging GFF (with sequences) -> ${GFF_ALL}"
TMP="${GFF_ALL}.tmp"
"${GFF3_MERGE}" -s -d "${DATASTORE_INDEX}" > "${TMP}"
nonempty_or_die "${TMP}"
mv -f "${TMP}" "${GFF_ALL}"

# 2) Merge GFF without sequences
echo "[2/5] Merging GFF (no sequences) -> ${GFF_NOSEQ}"
TMP="${GFF_NOSEQ}.tmp"
"${GFF3_MERGE}" -n -s -d "${DATASTORE_INDEX}" > "${TMP}"
nonempty_or_die "${TMP}"
mv -f "${TMP}" "${GFF_NOSEQ}"

# 3) Merge FASTA (transcripts/proteins)
echo "[3/5] Merging FASTA -> ${FA_TXT} / ${FA_PRO}"
"${FASTA_MERGE}" -d "${DATASTORE_INDEX}" -o "${OUT_PREFIX}"
nonempty_or_die "${FA_TXT}"
nonempty_or_die "${FA_PRO}"

# 4) Quick stats
echo "[4/5] Computing quick stats -> ${SUMMARY}"
GENE_MODELS=$(awk -F'\t' '$3=="mRNA"' "${GFF_ALL}" | wc -l | awk '{print $1}')
TOP_SCAF=$(awk -F'\t' '$3=="mRNA"{c[$1]++} END{for(k in c) print c[k]"\t"k}' "${GFF_ALL}" | sort -nr | head -20 || true)
FEAT_COUNTS=$(awk -F'\t' '($0 !~ /^#/ && NF>=3){c[$3]++} END{for(k in c) printf "%s\t%d\n", k, c[k]}' "${GFF_ALL}" | sort -k2,2nr || true)

{
  echo "MAKER merged outputs summary"
  echo "Timestamp: $(date -Is)"
  echo "Output prefix: ${OUT_PREFIX}"
  echo "Datastore index: ${DATASTORE_INDEX}"
  echo
  echo "Files:"
  printf "  %-45s %s\n" "$(basename "${GFF_ALL}")" "$(wc -l < "${GFF_ALL}" 2>/dev/null || echo 0) lines"
  printf "  %-45s %s\n" "$(basename "${GFF_NOSEQ}")" "$(wc -l < "${GFF_NOSEQ}" 2>/dev/null || echo 0) lines"
  printf "  %-45s %s\n" "$(basename "${FA_TXT}")" "$(grep -c '^>' "${FA_TXT}" 2>/dev/null || echo 0) records"
  printf "  %-45s %s\n" "$(basename "${FA_PRO}")" "$(grep -c '^>' "${FA_PRO}" 2>/dev/null || echo 0) records"
  echo
  echo "Approx. number of gene models (mRNA count): ${GENE_MODELS}"
  echo
  echo "Top 20 scaffolds by mRNA count:"
  echo -e "mRNA\tScaffold"
  [[ -n "${TOP_SCAF}" ]] && echo "${TOP_SCAF}" || echo "(none)"
  echo
  echo "Feature type breakdown (from merged GFF):"
  echo -e "Feature\tCount"
  [[ -n "${FEAT_COUNTS}" ]] && echo "${FEAT_COUNTS}" || echo "(none)"
} > "${SUMMARY}"
nonempty_or_die "${SUMMARY}"

# 5) Done
echo "[5/5] Done."
echo "GFF (with seq):  ${GFF_ALL}"
echo "GFF (no seq):    ${GFF_NOSEQ}"
echo "Transcripts:     ${FA_TXT}"
echo "Proteins:        ${FA_PRO}"
echo "Summary:         ${SUMMARY}"