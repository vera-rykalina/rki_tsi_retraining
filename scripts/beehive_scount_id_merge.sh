#!/bin/bash
#SBATCH --time=0-20:00:00
#SBATCH --partition main
#SBATCH --mem=210GB
 

# Created: 04.07.2025
# Author: Vera Rykalina
# Purpurse: Renaming BEEHIVE data so that a file name is combined: seq_id plus scount (RKI ID). 
# Usage: srun --pty bash -i beehive_scount_id_merge.sh 


# get ticket
hpc-ticket
 
# mount projecte
hpc-mount projekte

# define a variable for a folder
INIT="19960"

# create a folder
if [ -d "$TMPDIR/mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI//TSI_beehive_renamed_${INIT}" ]; 
then rm -Rf $TMPDIR/mnt/projekte/rykalinav/FG18_HIV_Pipelines/HIV-phyloTSI//TSI_beehive_renamed_${INIT}; fi

mkdir -p $TMPDIR/projekte/FG18_HIV_Pipelines/HIV-phyloTSI/TSI_beehive_renamed_${INIT}

# cd to HIV-phyloTSI
cd $TMPDIR/projekte/FG18_HIV_Pipelines/HIV-phyloTSI

# copy fastq.gz files
cp TSI_beehive_raw_${INIT}/*.fastq.gz TSI_beehive_renamed_${INIT}
cp TSI_beehive_raw_${INIT}/*_${INIT}.txt TSI_beehive_renamed_${INIT}

# cd to the folder with files
cd TSI_beehive_renamed_${INIT}

# rename file to add scount 
while IFS= read -r line; do
    oldfile1=`echo $line | awk '{print $1}'`
    newfile1=`echo $line | awk '{print $3}'`
    oldfile2=`echo $line | awk '{print $2}'`
    newfile2=`echo $line | awk '{print $4}'`
    echo "Moving $oldfile1 to $newfile1"
    echo "Moving $oldfile2 to $newfile2"
    mv "$oldfile1" "$newfile1"
    mv "$oldfile2" "$newfile2"
  
done < beehive_scount_${INIT}.txt


# unmount projekte
hpc-umount projekte
