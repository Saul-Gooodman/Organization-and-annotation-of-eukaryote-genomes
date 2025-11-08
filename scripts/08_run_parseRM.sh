#!/bin/bash
#SBATCH --job-name=parseRM
#SBATCH --partition=pshort_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
ANNODIR="$OUTDIR/assembly.p_ctg.fa.mod.EDTA.anno"
LAND="$OUTDIR/landscape"
mkdir -p "$LAND"

RMOUT="$ANNODIR/assembly.p_ctg.fa.mod.out"
if [[ ! -f "$RMOUT" && -f "$ANNODIR/assembly.p_ctg.fa.mod.out.gz" ]]; then
  RMOUT="$ANNODIR/assembly.p_ctg.fa.mod.out.gz"
fi
if [[ ! -f "$RMOUT" ]]; then
  echo "ERROR: RepeatMasker .out(.gz) not found under $ANNODIR" >&2
  exit 2
fi

module load BioPerl/1.7.8-GCCcore-10.3.0

cd "$LAND"
ln -sf "$RMOUT" ./genome.mod.out

if [[ ! -f ./parseRM.pl ]]; then
  curl -L -o parseRM.pl https://raw.githubusercontent.com/4ureliek/Parsing-RepeatMasker-Outputs/master/parseRM.pl
  chmod +x parseRM.pl
fi

perl ./parseRM.pl -i genome.mod.out -l 50,1 -v > parseRM.log 2>&1 || {
  echo "parseRM.pl failed. See $LAND/parseRM.log" >&2
  exit 3
}

# Normalize output names
# Common parseRM outputs include files with .summary and .divsum suffixes
DIVSUM_CAND=$(ls -1 *.divsum 2>/dev/null | head -n1 || true)
if [[ -z "${DIVSUM_CAND}" ]]; then
  echo "WARNING: No *.divsum produced by parseRM.pl; will let R fallback on .summary" >&2
else
  cp -f "$DIVSUM_CAND" "$LAND/assembly_p_ctg.parseRM.divsum.tsv"
fi

SUMMARY_CAND=$(ls -1 *.summary 2>/dev/null | head -n1 || true)
if [[ -n "${SUMMARY_CAND}" ]]; then
  cp -f "$SUMMARY_CAND" "$LAND/assembly_p_ctg.parseRM.summary.tsv"
fi

echo "Done parseRM. Outputs in: $LAND"