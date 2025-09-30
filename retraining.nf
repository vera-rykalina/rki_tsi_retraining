#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


/*
nextflow retraining.nf \
         -c retraining.config \
         -profile rki_slurm,rki_conda \
         --outdir Results \
         --dataset inputs/tsi_seroconverter_beehive.xlsx \
         --primers inputs/primers_sk_validation.fasta \
         --modelname SK \
         -resume 
*/


// help message
params.help = false

if (params.help) { exit 0, helpMSG() }


// Check required parameters
if (!params.outdir) {
  error "Missing output directory, use [--outdir]"
}

if (!params.dataset) {
  error "Missing input dataset file, use [--dataset]"
}

if (!params.primers) {
  error "Missing input primers file, use [--primers]"
}

if (!params.modelname) {
  error "Missing model name (e.g. IGS, SK, etc), use [--modelname]"
}


params.profile = null
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }



def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\u001B[31m";
    c_dim = "\033[2m";
    log.info """
  

    ${c_blue}TSI Remodling Pipeline${c_blue}
    ====================================================
    Author: Vera Rykalina
    ${c_blue}Affiliation: Robert Koch Institute${c_blue}
    Created: 22 August 2025
    ====================================================
  

    ${c_yellow}Usage examples:${c_reset}
    nextflow retraining.nf \
                --dataset inputs/tsi_seroconverter_beehive.xlsx \
                --primers inputs/primers_sk_validation.fasta
                --modelname SK \
                --oudir MyModel \
                -profile rki_slurm,rki_conda \
                -c retraining.config \ 
                -resume
   
    
    ${c_green}Required settings:${c_reset}  
    
    --dataset                         Path to a .xlsx annotated FL data [No change needed: use inputs/tsi_seroconverter_beehive.xlsx!]
    
    --outdir                          Name for an output directory, e.g. MyModel, Results etc. [string]

    --primers                         Path of a FASTA file containing the primer sequences (see examples in the input folder)

    --modelname                       Name for a model folder, e.g. SK or IGS [default: IGS]

    ${c_green}Optional input settings:${c_reset}
     
    --sheet                           Sheet in the .xlsx dataset [default: 1]

    --mafs_patstats_search_dir        Path to the folder where FL analysed data are [default: FG18_HIV_Pipelines/HIV-phyloTSI/HIVtime_single_full_length_samples_v2/]  
    
    -profile                          Sets a profile [rki_slurm,rki_conda]

    """
}




params.sheet = 1
params.mafs_patstats_search_dir = "FG18_HIV_Pipelines/HIV-phyloTSI/HIVtime_single_full_length_samples_v2/" // change accordingly; a folder where are all single scount_maf.csv.
params.reference = "${projectDir}/inputs/ref_3455.fasta"
params.refdata_tsi_hxb2 = "${projectDir}/inputs/HXB2_refdata.csv"
modelname = params.modelname



process GET_REGIONS_FOR_MASKING {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/01_primer_product_regions", mode: "copy", overwrite: true
  debug true

  input:
    path ref_fasta
    path primer_fasta
    
    
  output:
    path "regions.csv"
    
  script:
    """
    primer_finder.py \
     --ref ${ref_fasta} \
     --primers ${primer_fasta} \
     --output regions.csv
    
    """
  }


process SELECT_SAMPLES {
  label "low"
  conda "${projectDir}/envs/retrainingR.yml"
  publishDir "${params.outdir}/02_filtered_dataset", mode: "copy", overwrite: true
  
  input:
    path excel_file
    val sheet

  output:
    path "metadata.csv", emit: Csv
    path "metadata.txt", emit: Txt

  script:
  """
  select_samples.R \
    --input ${excel_file} \
    --output-csv metadata.csv \
    --output-txt metadata.txt \
    --sheet ${sheet}
  """
}

process COPY_SELECT_MAFS {
    label "low"
    //publishDir "${params.outdir}/selected_mafs", mode: 'copy', overwrite: true
    //debug true

    input:
      path id_list_txt
      val search_dir
 

    output:
      path "*.csv"

    script:
    """

    search_mafs.sh \
      --id-list ${id_list_txt} \
      --search-dir ${search_dir} \
      --output-dir .

    """
}

process CONCATINATE_MAFS {
  label "low"
  publishDir "${params.outdir}/03_concatinated_maf", mode: "copy", overwrite: true
  debug true

  input:
    path maf_csvs
    
  output:
    path "*.csv"
    
  script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' ${maf_csvs} > concatinated_maf.csv
    
    """
  }

process MAF_MASK_POSITIONS {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/04_maf_postions_masked", mode: "copy", overwrite: true
  //debug true

  input:
    path concat_maf
    path regions_to_mask
    
  output:
    path "maf_positions_masked.csv"
    
  script:
    """
    maf_primer_masker.py \
     --maf ${concat_maf} \
     --regions ${regions_to_mask} \
     --output maf_positions_masked.csv
    
    """
  }


process COPY_SELECT_PATSTATS {
    label "low"
    //publishDir "${params.outdir}/selected_patstats", mode: 'copy', overwrite: true
    //debug true

    input:
      path id_list_txt
      val search_dir
 

    output:
      path "*.csv"

    script:
    """

    search_patstats.sh \
      --id-list ${id_list_txt} \
      --search-dir ${search_dir} \
      --output-dir .

    """
}


process REMOVE_PROP_GP {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  //publishDir "${params.outdir}/cleaned_patstats", mode: "copy", overwrite: true
  //debug true

  input:
    path patstats_csv

  output:
    path "*.csv"

  script:
    // Extract ID from filename without extension
    def id = patstats_csv.simpleName.split('_')[0]

    """
    remove_prop_gp_cols.py \
      --input ${patstats_csv} \
      --output ${id}_no_prop_patStats.csv
    """
}


process CONCATINATE_PATSTATS {
  label "low"
  publishDir "${params.outdir}/05_concatinated_patstats", mode: "copy", overwrite: true
  //debug true

  input:
    path patstats_csvs
    
  output:
    path "*.csv"
    
  script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' ${patstats_csvs} > concatinated_patstats.csv
    
    """
  }

process PATSTATS_MASK_POSITIONS {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/06_patstats_postions_masked", mode: "copy", overwrite: true
  //debug true

  input:
    path concat_patstats
    path regions_to_mask
    
  output:
    path "patstats_positions_masked.csv"
    
  script:
    """
    patstats_primer_masker.py \
     --patstats ${concat_patstats} \
     --regions ${regions_to_mask} \
     --output patstats_positions_masked.csv
    
    """
  }

process CALCULATE_MEANS {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/07_all_features_means", mode: "copy", overwrite: true
  //debug true

  input:
    path refdata_tsi_hxb2
    path maf_positions_masked
    path patstats_positions_masked
    
  output:
    path "all_features_means.csv"
    
  script:
    """
    calculate_means.py \
     --ref ${refdata_tsi_hxb2} \
     --maf ${maf_positions_masked} \
     --patstats ${patstats_positions_masked} \
     --output all_features_means.csv
    
    """
  }


process GET_RETRAINING_DF {
  label "low"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/08_metadata_features", mode: "copy", overwrite: true
  //debug true

  input:
    path metadata_df
    path features_means_df
  
    
  output:
    path "retraining_df.csv"
    
  script:
    """
    merge_metadata_features.py \
     --metadata ${metadata_df} \
     --features ${features_means_df} \
     --output retraining_df.csv
    
    """
  }


process FEATURE_SELECTION_REPORTS_BASE {
  label "medium"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/09_features_reports", mode: "copy", overwrite: true
  //debug true

  input:
    path retraining_df

  output:
    path "report_base/top10_base.csv"
    path "report_base/best_features_base.txt", emit: Txt

    
  script:
    """
    feature_selection.py \
     --input ${retraining_df} \
     --output report_base \
     --amplicons
    
    """
  }


process FEATURE_SELECTION_REPORTS_BASE_AMP {
  label "high"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/09_features_reports", mode: "copy", overwrite: true
  //debug true

  input:
    path retraining_df

  output:
    path "report_base_amp/top10_base_amp.csv"
    path "report_base_amp/best_features_base_amp.txt", emit: Txt


    
  script:
    """
    feature_selection.py \
     --input ${retraining_df} \
     --output report_base_amp \
     --amplicons \
     --include_is_amp
    
    """
  }



  process FEATURE_SELECTION_REPORTS_BASE_VL {
  label "high"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/09_features_reports", mode: "copy", overwrite: true
  //debug true

  input:
    path retraining_df

  output:
    path "report_base_vl/top10_base_vl.csv"
    path "report_base_vl/best_features_base_vl.txt", emit: Txt

    
  script:
    """
    feature_selection.py \
     --input ${retraining_df} \
     --output report_base_vl \
     --amplicons \
     --include_viral_load
    
    """
  }

process TRAINING {
  label "high"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/10_training_outputs", mode: "copy", overwrite: true
  //debug true

  input:
    path refdata_tsi_hxb2
    path regions_csv
    path retraining_df
    path features_list
    val modelname

  output:
    path "model_${modelname}"
    path "reports"

  script:
    """
    training.py \
     --input ${retraining_df} \
     --features ${features_list} \
     --modeldir model_${modelname} \
     --report reports \
     --modelname ${modelname} \
     --amplicons    

    mv ${refdata_tsi_hxb2} model_${modelname}
    mv ${regions_csv} model_${modelname}
    """
  }

workflow {
  // inputs
  ch_dataset = Channel.fromPath ( params.dataset, checkIfExists: true )
  ch_primers = Channel.fromPath ( params.primers, checkIfExists: true )
  ch_refdata = Channel.fromPath ( params.refdata_tsi_hxb2, checkIfExists: true )
  
  // general processes
  ch_masking_regions = GET_REGIONS_FOR_MASKING ( params.reference, ch_primers )
  ch_samples = SELECT_SAMPLES ( ch_dataset, params.sheet )

  // processes maf
  ch_mafs = COPY_SELECT_MAFS ( ch_samples.Txt, params.mafs_patstats_search_dir )
  ch_concat_maf = CONCATINATE_MAFS ( ch_mafs.collect() )
  ch_maf_pos_masked = MAF_MASK_POSITIONS ( ch_concat_maf, ch_masking_regions )
  
  // processes patstats
  ch_patstats = COPY_SELECT_PATSTATS ( ch_samples.Txt, params.mafs_patstats_search_dir )
  ch_patstats_no_prop_gp = REMOVE_PROP_GP ( ch_patstats.flatten() )
  ch_concat_patstat = CONCATINATE_PATSTATS ( ch_patstats_no_prop_gp.collect() )
  ch_patstats_pos_masked = PATSTATS_MASK_POSITIONS ( ch_concat_patstat, ch_masking_regions )  

  // remodelling
  ch_all_features_means = CALCULATE_MEANS ( params.refdata_tsi_hxb2, ch_maf_pos_masked, ch_patstats_pos_masked )
  ch_retraining_df = GET_RETRAINING_DF ( ch_samples.Csv, ch_all_features_means )
  ch_base_amp = FEATURE_SELECTION_REPORTS_BASE_AMP (ch_retraining_df )
  ch_base_vl = FEATURE_SELECTION_REPORTS_BASE_VL (ch_retraining_df )
  ch_base = FEATURE_SELECTION_REPORTS_BASE ( ch_retraining_df )
  ch_model = TRAINING ( ch_masking_regions, ch_refdata, ch_retraining_df, ch_base_amp.Txt, params.modelname )
}



