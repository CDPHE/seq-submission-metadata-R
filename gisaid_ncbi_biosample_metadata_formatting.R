#!/usr/bin/Rscript

# Mandy Waters
# July 2021
# Script to read in COVID sequencing output metadata and write out a GISAID and NCBI biosample formatted metadata sheet for samples with >50% coverage
# Usage (short) = ./gisaid_ncbi_biosample_metadata_formatting.R -m "path/to/metadata_sheets" -f fasta_file.fa -s submitter_name -t "sequencing type"
# Usage (long) = ./gisaid_ncbi_biosample_metadata_formatting.R --metadata "path/to/metadata_sheets"  --fasta_name fasta_file.fa --submitter_name submitter_name --seq_type "sequencing type"
# Usage notes: only strings with spaces need quotes
# Usage example = ./gisaid_ncbi_biosample_metadata_formatting.R -m /home/mandy_waters -f GISAID_upload_sequences.fa -s mandy_waters -t "Oxford Nanopore GridION"

# load packages, but don't print the messages to terminal
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))

# set argument / options
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-m", "--metadata"), action="store", default=NA, type='character',
              help="path to metadata sheets"),
  make_option(c("-f", "--fasta_name "), action="store", default=NA, type='character',
              help="fasta file ename to upload"),
  make_option(c("-s", "--submitter_name"), action="store", default=NA, type='character',
              help="submitter name"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Make the program not be verbose."),
  make_option(c("-t", "--seq_type"), action="store", default=NA, type='character',
              help="sequencing type as a string in quotes")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$v) {
  # print the variables to the terminal
  cat("\nmetadata path:\n")
  cat(opt$metadata)
  cat("\n\nname of fasta file to submit:\n")
  cat(opt$fasta_name)
  cat("\n\nsubmitter name:\n")
  cat(opt$submitter_name)
  cat("\n\nsequencing type:\n")
  cat(opt$seq_type)
  cat("\n\n")
}

# if all options are filled in then fill in the metadata otherwise error out
if(!is.na(opt$metadata) & !is.na(opt$submitter_name) & !is.na(opt$seq_type) & !is.na(opt$fasta_name)){
  if(opt$seq_type == "Illumina MiSeq" | opt$seq_type == "Illumina NextSeq" | opt$seq_type == "Oxford Nanopore GridION"){

# set working directory
setwd(opt$metadata)
date = Sys.Date()

# Reading in all metadata sheets, subsetting out the columns needed, pulling samples with >= 50% coverage
metadata_readin = ldply(list.files(pattern = "results_.*.csv"), read.csv, header=TRUE, na.strings = "", check.names = FALSE)
cols_to_keep = as.data.frame(c("accession_id", "percent_non_ambigous_bases", "collection_date", "fasta_header", "seq_run"))
colnames(cols_to_keep) = "col_names"
metadata_subset_cols = metadata_readin[,colnames(metadata_readin) %in% cols_to_keep$col_names]
metadata_drop_no_assembly = metadata_subset_cols[metadata_subset_cols$percent_non_ambigous_bases != 'sample failed assembly',]
metadata_50_cov = metadata_drop_no_assembly[metadata_drop_no_assembly$percent_non_ambigous_bases >= 50.0,]

# Creating an eepty matrix to fill in with metadata in GISAID format
gisad_metadata = as.data.frame(matrix("",ncol=29,nrow=nrow(metadata_50_cov)))
colnames(gisad_metadata) = c("submitter","fn","covv_virus_name","covv_type","covv_passage","covv_collection_date","covv_location","covv_add_location","covv_host","covv_add_host_info","covv_gender","covv_patient_age","covv_patient_status","covv_specimen","covv_outbreak","covv_last_vaccinated","covv_treatment","covv_seq_technology","covv_assembly_method","covv_coverage","covv_orig_lab","covv_orig_lab_addr","covv_provider_sample_id","covv_subm_lab","covv_subm_lab_addr","covv_subm_sample_id","covv_authors","covv_comment","comment_type")

rename_fasta = as.data.frame(matrix("",ncol=2,nrow=nrow(metadata_50_cov)))
colnames(rename_fasta) = c("accesion_id", "gisaid_name")

# Filling in columns whose values do not change or will be put in with options
gisad_metadata$submitter = opt$submitter_name
gisad_metadata$fn	= opt$fasta_name
gisad_metadata$covv_type = "betacoronavirus"
gisad_metadata$covv_passage = "Original"
gisad_metadata$covv_location = "North America / USA / Colorado"
gisad_metadata$covv_host = "Human"
gisad_metadata$covv_gender = "unknown"
gisad_metadata$covv_patient_age	= "unknown"
gisad_metadata$covv_patient_status = "unknown"
gisad_metadata$covv_seq_technology = opt$seq_type #Oxford Nanopore GridION, Illumina MiSeq, Illumina NextSeq
gisad_metadata$covv_orig_lab = "Colorado Department of Public Health and Environment"
gisad_metadata$covv_orig_lab_addr = "8100 E. Lowry Blvd. Denver CO, 80230"
gisad_metadata$covv_subm_lab = "Colorado Department of Public Health and Environment"
gisad_metadata$covv_subm_lab_addr = "8100 E. Lowry Blvd. Denver CO, 80230"
gisad_metadata$covv_authors = "Laura Bankers, Molly C. Hetherington-Rauth, Diana Ir, Alexandria Rossheim, Michael Martin, Mandy Waters, Shannon R. Matzinger, Sarah Elizabeth Totten, Emily A. Travanty"

# Filling in data that has unique info for each row
for (i in 1:nrow(gisad_metadata)){
  gisad_metadata[i,6] = metadata_50_cov[i,3] # covv_collection_date
  gisad_metadata[i,3] = paste("hCoV-19/USA/CO-CDPHE-", metadata_50_cov[i,1], "/", str_sub(gisad_metadata[i,6],start=1,end=4),sep="") # covv_virus_name
  rename_fasta[i,1] = paste("CO-CDPHE-", metadata_50_cov[i,1], sep = "")
  rename_fasta[i,2] = paste("hCoV-19/USA/CO-CDPHE-", metadata_50_cov[i,1], "/", str_sub(gisad_metadata[i,6],start=1,end=4),sep="")
}

write.csv(gisad_metadata, file = paste("gisaid_submission", date, "metadata.csv",sep="_"), row.names = FALSE, quote = TRUE)
write.csv(rename_fasta, "fasta_rename_accession_to_gisaid_id.csv", row.names = FALSE, quote = FALSE)

# Creating the NCBI biosample metadata table
ncbi_biosample_metadata = as.data.frame(matrix("",ncol=48,nrow=nrow(metadata_50_cov)))
colnames(ncbi_biosample_metadata)= c("sample_name","sample_title","bioproject_accession","organism","collected_by","collection_date","geo_loc_name","host","host_disease","isolate","isolation_source","antiviral_treatment_agent","collection_device","collection_method","date_of_prior_antiviral_treat","date_of_prior_sars_cov_2_infection","date_of_sars_cov_2_vaccination","exposure_event","geo_loc_exposure","gisaid_accession","gisaid_virus_name","host_age","host_anatomical_material","host_anatomical_part","host_body_product","host_disease_outcome","host_health_state","host_recent_travel_loc","host_recent_travel_return_date","host_sex","host_specimen_voucher","host_subject_id","lat_lon","passage_method","passage_number","prior_sars_cov_2_antiviral_treat","prior_sars_cov_2_infection","prior_sars_cov_2_vaccination","purpose_of_sampling","purpose_of_sequencing","sars_cov_2_diag_gene_name_1","sars_cov_2_diag_gene_name_2","sars_cov_2_diag_pcr_ct_value_1","sars_cov_2_diag_pcr_ct_value_2","sequenced_by","vaccine_received","virus_isolate_of_prior_infection","description")

ncbi_biosample_metadata$sample_title = "PCR tiled amplicon WGS of SARS-CoV-2"
ncbi_biosample_metadata$bioproject_accession = "PRJNA686984"
ncbi_biosample_metadata$organism = "Severe acute respiratory syndrome coronavirus 2"
ncbi_biosample_metadata$collected_by = "CDPHE"
ncbi_biosample_metadata$geo_loc_name = "USA: Colorado"
ncbi_biosample_metadata$host = "Homo sapiens"
ncbi_biosample_metadata$host_disease = "severe accute respiratory syndrome"
ncbi_biosample_metadata$isolation_source = "patient isolate"

for (i in 1:nrow(gisad_metadata)){
ncbi_biosample_metadata[i,1] = paste("CO-CDPHE-", metadata_50_cov[i,1],sep="") # sample_name
ncbi_biosample_metadata[i,10] = paste("CO-CDPHE-", metadata_50_cov[i,1],sep="") # isolate
ncbi_biosample_metadata[i,6] = metadata_50_cov[i,3] # collection_date
}
} else {
  cat("ERROR: check the spelling of your sequence type input: sequence type should match 'Illumina MiSeq', 'Illumina NextSeq', or 'Oxford Nanopore GridION'\n\n")
  }
} else {
  cat("ERROR: you didn't specify all varialbles. Please supply all values\n\n") # print error messages
}

write.table(ncbi_biosample_metadata, file=paste("ncbi_biosample_submission", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = FALSE, sep = '\t')
