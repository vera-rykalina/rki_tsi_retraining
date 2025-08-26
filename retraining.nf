#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Check required parameters
if (!params.outdir) {
  error "Missing output directory, use [--outdir]"
}
if (!params.dataset) {
  error "Missing input dataset file, use [--dataset]"
}


params.sheet = 1
params.mafs_patstats_search_dir = "FG18_HIV_Pipelines/HIV-phyloTSI/HIVtime_single_full_length_samples_v2/" // change accordingly; a folder where are all single scount_maf.csv.
params.reference = "${projectDir}/inputs/ref_3455.fasta"
params.refdata_tsi_hxb2 = "${projectDir}/inputs/HXB2_refdata.csv"



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
    publishDir "${params.outdir}/03_selected_mafs", mode: 'copy', overwrite: true
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
  publishDir "${params.outdir}/04_concatinated_maf", mode: "copy", overwrite: true
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
  publishDir "${params.outdir}/05_maf_postions_masked", mode: "copy", overwrite: true
  debug true

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
    publishDir "${params.outdir}/06_selected_patstats", mode: 'copy', overwrite: true
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
  publishDir "${params.outdir}/07_cleaned_patstats", mode: "copy", overwrite: true
  debug true

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
  publishDir "${params.outdir}/08_concatinated_patstats", mode: "copy", overwrite: true
  debug true

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
  publishDir "${params.outdir}/09_patstats_postions_masked", mode: "copy", overwrite: true
  debug true

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
  publishDir "${params.outdir}/10_all_features_means", mode: "copy", overwrite: true
  debug true

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
  publishDir "${params.outdir}/11_metadata_features", mode: "copy", overwrite: true
  debug true

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


process FEATURE_SELECTION_REPORTS {
  label "medium"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/12_features_reports", mode: "copy", overwrite: true
  //debug true

  input:
    path retraining_df

  output:
    path "features_without_vl.csv"
    path "features_with_vl.csv"
    path "features.txt", emit: Txt

    
  script:
    """
    feature_selection.py \
     --input ${retraining_df} \
     --output-csv-with-vl features_with_vl.csv \
     --output-csv-without-vl features_without_vl.csv \
     --output-txt features.txt \
     --amplicons
    
    """
  }

process TRAINING {
  label "medium"
  conda "${projectDir}/envs/phylo_tsi.yml"
  publishDir "${params.outdir}/13_training_outputs", mode: "copy", overwrite: true
  debug true

  input:
    path retraining_df
    path features_list

  output:
    path "model_igs"
    path "metric.csv"

    
  script:
    """
    training.py \
     --input ${retraining_df} \
     --features ${features_list} \
     --modeldir model_igs \
     --report metric.csv \
     --amplicons    
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
  ch_features = FEATURE_SELECTION_REPORTS ( ch_retraining_df )
  ch_model_igs = TRAINING ( ch_retraining_df, ch_features.Txt )
}



