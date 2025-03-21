#!/bin/bash
#SBATCH --time=0-20:00:00
#SBATCH --partition main
#SBATCH --mem=210GB
 
# get ticket
hpc-ticket
 
# mount wissdaten
hpc-mount wissdaten
 
# create a folder
mkdir -p $TMPDIR/md5checks

# copy fastq.gz files
cp $TMPDIR/wissdaten/FG18_HIV_NGS_Rohdaten/02_Serokonverter/BEEHIVE/reads_19960/*.fastq.gz $TMPDIR/md5checks
cp $TMPDIR/wissdaten/FG18_HIV_NGS_Rohdaten/02_Serokonverter/BEEHIVE/reads_20004/*.fastq.gz $TMPDIR/md5checks
cp $TMPDIR/wissdaten/FG18_HIV_NGS_Rohdaten/02_Serokonverter/BEEHIVE/reads_20005/*.fastq.gz $TMPDIR/md5checks

# copy previously created md5sums
cp $TMPDIR/wissdaten/FG18_HIV_NGS_Rohdaten/02_Serokonverter/BEEHIVE/cs-transfer/checksums/*.txt $TMPDIR/md5checks


# unmount wissdaten
hpc-umount wissdaten
 
# check md5 sum
cd $TMPDIR/md5checks

for i in *.fastq.gz; do md5sum $i | awk '{print $1}' > ${i%.fastq.gz}_md5_local.txt; done

rm $TMPDIR/md5checks/19960_3_153_1_md5_local.txt
rm $TMPDIR/md5checks/19960_3_153_2_md5_local.txt

rm *.fastq.gz

# compare md5sums
for i in *_local.txt; do diff -y $i ${i%_local.txt}_sanger.txt >> md5_sums.txt; done
for i in *_local.txt; do diff -s $i ${i%_local.txt}_sanger.txt >> md5_diff.txt; done

# move results to Projekte
hpc-ticket

hpc-mount projekte

cp $TMPDIR/md5checks/*.txt $TMPDIR/projekte/FG18_HIV_Pipelines/HIV-phyloTSI/BEEHIVE_SUM_RESULTS

hpc-umount projekte
