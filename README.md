# Model Generator Pipeline
This repo contains a fully automated Nextflow pipeline for remodeling to be able to generate models for different lab settings (partial or full _gag_ and _pol_ genomic regions). 

## Installation
The pipeline is written in Nextflow, which can be easily installed using conda.

- Install conda if not installed:

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

- Install mamba if not installed (recommended):

```sh
conda install mamba -n base -c conda-forge
```

- Clone the repo

```sh
git clone https://github.com/vera-rykalina/rki_tsi_retraining.git  
```

- Install nextflow, using **nextflow.yml** file

```sh
conda env create -n nextflow -f env/nextflow.yml
```

All other software tools with their dependencies get installed automatically within the pipeline via conda directives. 


## Prerequisites
The pipeline retrieves the data from the project folder outside of the HPC environment. As this process is automated, the .keytab file should be created.

- 1. Check if keytab already exists

```bash
cd  
ls -al
```

Look for a folder called `.keytab`.  
If not present, proceed to create it.


- 2. Create Keytab (Login Node)

```bash
hpc-keytab
```

Follow the instructions:

```text
>> RKI HPC Keytab Creator
Please be sure to read: https://confluence.rki.local/x/Xa3ZB
continue with creation? [y/N] y
```


- 3. Test Keytab Setup

#### Remove existing tickets:<p>

```bash
kdestroy
```

- 4. create new one:

```bash
hpc-ticket
klist
```

Expected output:
```
Default principal: username@RKI.LOCAL
Valid starting... Expires...
```



## Usage (for more information use **--help**)

- Activate the *nextflow* environment:
```sh
conda activate nextflow
```

- Run: 

```sh
```bash
nextflow retraining.nf -c retraining.config \
--dataset inputs/tsi_seroconverter_beehive.xlsx \
--primers inputs/primers_sk_validation.fasta \
--modelname SK \
--outdir Results \
-profile rki_slurm,rki_conda \
-resume
```


## Pipeline workflow
![Plot](/images/retraining_500_dpi.png)