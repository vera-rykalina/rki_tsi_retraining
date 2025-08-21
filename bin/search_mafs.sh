#!/bin/bash

usage() {
  echo "Usage: $0 -i <id_list.txt> -s <search_root_dir> -o <output_dir>"
  echo
  echo "  -i, --id-list      Text file with sample IDs (one per line)"
  echo "  -s, --search-dir   Relative path under projekte (e.g. FG18_HIV_Pipelines/HIV-phyloTSI/HIVtime_single_full_length_samples_v2/)"
  echo "  -o, --output-dir   Directory where matched .csv files will be copied"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--id-list)
      ID_LIST="$2"
      shift 2
      ;;
    -s|--search-dir)
      SEARCH_DIR="$2"
      shift 2
      ;;
    -o|--output-dir)
      OUTPUT_DIR="$2"
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

if [[ -z "$ID_LIST" || -z "$SEARCH_DIR" || -z "$OUTPUT_DIR" ]]; then
  echo "Error: Missing required arguments."
  usage
fi

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

MOUNTED_DIR="$TMPDIR/projekte/$SEARCH_DIR"

# Prepare output directory (no extra subfolder)
mkdir -p "$OUTPUT_DIR"

# Optional: copy all *_maf.csv files to a temp folder to speed searching
TMP_ALL="$TMPDIR/all_csv"
mkdir -p "$TMP_ALL"
find "$MOUNTED_DIR" -name "*_maf.csv" -exec cp {} "$TMP_ALL" \;

count=0

while IFS= read -r ID || [[ -n "$ID" ]]; do
  match="$TMP_ALL/${ID}_maf.csv"
  if [[ -f "$match" ]]; then
    cp "$match" "$OUTPUT_DIR/"
    echo "Copied: $match"
    ((count++))
  else
    echo "Warning: File not found for ID '$ID'" >&2
  fi
done < "$ID_LIST"

echo "Total matching maf CSVs copied: $count"

echo "Unmounting projekte..."
hpc-umount projekte