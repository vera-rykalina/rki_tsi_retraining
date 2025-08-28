#!/bin/bash

set -euo pipefail

usage() {
  echo "Usage: $0 -i <id_list.txt> -b <base_dir> [-o <output_subdir_name>]"
  echo
  echo "  -i, --id-list        Text file with sample IDs (one per line)"
  echo "  -b, --base-dir       Base directory containing GAG_Rohdaten and PRRT-INT_Rohdaten"
  echo "  -o, --output-subdir  (Optional) Name of output folder inside base_dir (default: merged_validation)"
  echo
  echo "Example:"
  echo "  $0 -i ids.txt -b FG18_HIV_Pipelines/HIV-phyloTSI/Fastq_von_KM"
  exit 1
}

OUTPUT_SUBDIR="merged_validation"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--id-list)
      ID_LIST="$2"
      shift 2
      ;;
    -b|--base-dir)
      BASE_DIR="$2"
      shift 2
      ;;
    -o|--output-subdir)
      OUTPUT_SUBDIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Check required args
if [[ -z "${ID_LIST:-}" || -z "${BASE_DIR:-}" ]]; then
  echo "Error: Missing required arguments."
  usage
fi

# Mount projekte if not mounted
echo "Checking if projekte is already mounted..."
if ! mountpoint -q "$TMPDIR/projekte"; then
  echo "Requesting HPC access and mounting projekte..."
  hpc-ticket
  hpc-mount projekte
fi

if ! mountpoint -q "$TMPDIR/projekte"; then
  echo "Error: Failed to mount projekte."
  exit 2
fi

# Construct full paths to folders on mounted projekte
BASE_MOUNTED="$TMPDIR/projekte/$BASE_DIR"
GAG_DIR="$BASE_MOUNTED/GAG_Rohdaten"
PRRT_DIR="$BASE_MOUNTED/PRRT-INT_Rohdaten"
OUTPUT_DIR="$BASE_MOUNTED/$OUTPUT_SUBDIR"

# Check folders exist
if [[ ! -d "$GAG_DIR" ]]; then
  echo "Error: GAG_Rohdaten folder not found at: $GAG_DIR"
  exit 3
fi

if [[ ! -d "$PRRT_DIR" ]]; then
  echo "Error: PRRT-INT_Rohdaten folder not found at: $PRRT_DIR"
  exit 3
fi

mkdir -p "$OUTPUT_DIR"

DATE_PREFIX=$(date +%Y%m%d)

echo "Starting merging of fastq files..."

while IFS= read -r ID || [[ -n "$ID" ]]; do
  echo "Processing sample ID: $ID"

  # Find GAG files
  GAG_R1=$(find "$GAG_DIR" -type f -name "*${ID}*_R1_*.fastq.gz" | head -n 1)
  GAG_R2=$(find "$GAG_DIR" -type f -name "*${ID}*_R2_*.fastq.gz" | head -n 1)

  # Find PRRT files (search recursively)
  PRRT_R1=$(find "$PRRT_DIR" -type f -name "*${ID}*_R1_*.fastq.gz" | head -n 1)
  PRRT_R2=$(find "$PRRT_DIR" -type f -name "*${ID}*_R2_*.fastq.gz" | head -n 1)

  # Check if all files found
  if [[ -z "$GAG_R1" || -z "$GAG_R2" || -z "$PRRT_R1" || -z "$PRRT_R2" ]]; then
    echo "Warning: Missing files for ID $ID, skipping..."
    [[ -z "$GAG_R1" ]] && echo "  Missing GAG R1"
    [[ -z "$GAG_R2" ]] && echo "  Missing GAG R2"
    [[ -z "$PRRT_R1" ]] && echo "  Missing PRRT R1"
    [[ -z "$PRRT_R2" ]] && echo "  Missing PRRT R2"
    continue
  fi

  OUT_R1="${OUTPUT_DIR}/${DATE_PREFIX}_${ID}_GAG-PRRT-INT_validation_R1.fastq.gz"
  OUT_R2="${OUTPUT_DIR}/${DATE_PREFIX}_${ID}_GAG-PRRT-INT_validation_R2.fastq.gz"

  echo " Merging:"
  echo "   $GAG_R1 + $PRRT_R1 -> $OUT_R1"
  echo "   $GAG_R2 + $PRRT_R2 -> $OUT_R2"

  {
    gzip -dc "$GAG_R1"
    gzip -dc "$PRRT_R1"
  } | gzip > "$OUT_R1"

  {
    gzip -dc "$GAG_R2"
    gzip -dc "$PRRT_R2"
  } | gzip > "$OUT_R2"

  echo "  Done merging $ID."

done < "$ID_LIST"

echo "All done. Merged files saved in $OUTPUT_DIR"

echo "Unmounting projekte..."
hpc-umount projekte