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


process SELECT_SAMPLES {
  conda "${projectDir}/envs/retrainingR.yml"
  publishDir "${params.outdir}/01_filtered_dataset", mode: "copy", overwrite: true
  
  input:
    path excel_file
    val sheet

  output:
    path "retraining_samples.csv", emit: Csv
    path "retraining_samples.txt", emit: Txt

  script:
  """
  select_samples.R \
    --input ${excel_file} \
    --output-csv retraining_samples.csv \
    --output-txt retraining_samples.txt \
    --sheet ${sheet}
  """
}

process COPY_SELECT_MAFS {
    label "select_mafs_patstats"
    publishDir "${params.outdir}/02_selected_mafs", mode: 'copy', overwrite: true
    debug true

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


process COPY_SELECT_PATSTATS {
    label "select_mafs_patstats"
    publishDir "${params.outdir}/04_selected_patstats", mode: 'copy', overwrite: true
    debug true

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



workflow {
  ch_dataset = Channel.fromPath(params.dataset, checkIfExists: true)
  ch_samples = SELECT_SAMPLES(ch_dataset, params.sheet)
  ch_mafs = COPY_SELECT_MAFS ( ch_samples.Txt, params.mafs_patstats_search_dir )
  ch_concat_maf = CONCATINATE_MAFS ( ch_mafs.collect() )
  ch_patstats = COPY_SELECT_PATSTATS ( ch_samples.Txt, params.mafs_patstats_search_dir )

}



