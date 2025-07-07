#!/bin/bash

# Created: 03.07.2025
# Author: Vera Rykalina
# Purpurse: Get raw fastq sutible for HIVtime files by IDs from BEEHIVE data. 
# Usage: ./find_fastq_beehive.sh /mnt/wissdaten/rykalinav/FG18_HIV_NGS_Rohdaten/02_Serokonverter/BEEHIVE/reads_combined_20005
# where to copy -> /mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_20005
# Comment: Do not use more than 10 files per run (no space on the disk)


RANGE="20005"
# Get input data (tsi IDs)
VAR1=$(cat /scratch/rykalinav/rki_tsi_retraining/data/beehive_${RANGE}.txt)

if [ -d "/mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_${RANGE}" ]; then rm -Rf /mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_${RANGE}; fi

# Create dir for output files
mkdir -p /mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_${RANGE}


# This is a loop for copying the target files
COUNTER=0

for i in ${VAR1}
    do
    echo "Looking for fastq files with ID ${i}."
    find "$1" -type f -name "*${i}_[0-2].fastq.gz" -exec cp "{}" "/mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_${RANGE}" \;
    let COUNTER++
    done

echo "-------------------------------------------------------------------"
echo "${COUNTER} IDs are found."
echo "-------------------------------------------------------------------"

FILE_NUM=$(find /mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_raw_${RANGE} -type f | wc -l)
echo "${FILE_NUM} files have been copied."
echo "-------------------------------------------------------------------"
