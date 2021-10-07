# Mandy Waters
# July 2021
# Script to read in COVID sequencing output metadata and write out a GISAID and NCBI biosample formatted metadata sheet for samples with >50% coverage
# Usage (short) = Rscript gisaid_ncbi_biosample_metadata_formatting.R -m "path/to/metadata_sheets" -f fasta_file.fa -s submitter_name -t "sequencing type" -r "rerun_file"
# Usage (long) = Rscript gisaid_ncbi_biosample_metadata_formatting.R --metadata "path/to/metadata_sheets"  --fasta_name fasta_file.fa --submitter_name submitter_name --seq_type "sequencing type" --rerun_file /path/to/rerun_filename.tsv
# Usage notes: only strings with spaces need quotes
# Usage example = Rscript gisaid_ncbi_biosample_metadata_formatting.R -m /home/mandy_waters/metadata -f GISAID_upload_sequences.fa -s mandy_waters -t "Oxford Nanopore GridION" -r /home/mandy_waters/rerun_samples_12Sep2021.tsv

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
              help="REQUIRED: path to metadata sheets"),
  make_option(c("-f", "--fasta_name"), action="store", default=NA, type='character',
              help="REQUIRED: fasta file name to upload to GISAID"),
  make_option(c("-s", "--submitter_name"), action="store", default=NA, type='character',
              help="REQUIRED: submitter name for GISAID"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="DEFAULT: should the program print extra stuff out? [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="OPTIONAL: Make the program not be verbose"),
  make_option(c("-r", "--rerun_file"), action="store", default=NA, type='character',
              help="REQUIRED: path and file name for text file with rerun information as downloaded from Google Drive"),
  make_option(c("-c", "--completed_file"), action="store", default=NA, type='character',
              help="REQUIRED: path and file name for text file with completed submission information as downloaded from Google Drive"),
  make_option(c("-t", "--seq_type"), action="store", default=NA, type='character',
              help="REQUIRED: sequencing type as a string in quotes that has to match 'Illumina MiSeq', 'Illumina NextSeq', or 'Oxford Nanopore GridION'")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$v) {
  # print the variables to the terminal
  cat("\nmetadata path:\n")
  cat(opt$metadata)
  cat("\n\nname of fasta file to submit to GISAID:\n")
  cat(opt$fasta_name)
  cat("\n\nsubmitter name for GISAID:\n")
  cat(opt$submitter_name)
  cat("\n\nsequencing type:\n")
  cat(opt$seq_type)
  cat("\n\nname of rerun file:\n")
  cat(opt$rerun_file)
  cat("\n\nname of completed submissions file:\n")
  cat(opt$completed_file)
  cat("\n\n")
}

# if all options are filled in then fill in the metadata otherwise error out
if(!is.na(opt$metadata) & !is.na(opt$submitter_name) & !is.na(opt$seq_type) & !is.na(opt$fasta_name) & !is.na(opt$rerun_file) & !is.na(opt$completed_file)){
  if(opt$seq_type == "Illumina MiSeq" | opt$seq_type == "Illumina NextSeq" | opt$seq_type == "Oxford Nanopore GridION"){

# set working directory
setwd(opt$metadata)
date = Sys.Date()

# Reading in all metadata sheets, subsetting out the columns needed, pulling samples with >= 50% coverage
metadata_readin = ldply(list.files(pattern = "results_.*.tsv"), read.delim, header=TRUE, na.strings = c("","sample failed assembly"), check.names = FALSE)
cols_to_keep = as.data.frame(c("accession_id", "percent_non_ambigous_bases", "collection_date", "fasta_header", "seq_run"))
colnames(cols_to_keep) = "col_names"
metadata_subset_cols = metadata_readin[,colnames(metadata_readin) %in% cols_to_keep$col_names]
patterns <- c("Blank", "NC", "NTC", "ExtractionPositive", "DiaplexPositive", "POS")
metadata_no_blank_nc = filter(metadata_subset_cols, !grepl(paste(patterns, collapse="|"), accession_id))
metadata_50_cov = metadata_no_blank_nc[metadata_no_blank_nc$percent_non_ambigous_bases >= 50.0,]
dups = as.data.frame(metadata_50_cov[duplicated(metadata_50_cov$accession_id),c(1,5)])
colnames(dups) = c("duplicate_accessions", "seq_run")

# print warning if accession is in batch more than once
if(length(dups>0)){
  cat("\nWARNING: the follwoing samples were duplicated in the batch. The sample with higher coverage will be retained\n")
  print(dups, row.names = FALSE)
  cat('\n')
}

# Droopping the sample with lower coverage
metadata_50_cov = metadata_50_cov[order(metadata_50_cov[,"accession_id"],-metadata_50_cov[,"percent_non_ambigous_bases"]),]
metadata_50_cov = metadata_50_cov[!duplicated(metadata_50_cov$accession_id),]

# pulling samples that have missing or iincorrect (priror to 2020) collection date
missing_collection_date = as.data.frame(metadata_50_cov[is.na(metadata_50_cov$collection_date),c(1,2,3,5)])
missing_collection_date_removed = as.data.frame(metadata_50_cov[!is.na(metadata_50_cov$collection_date),c(1,2,3,5)])
collection_date_wrong = as.data.frame(missing_collection_date_removed[str_sub(missing_collection_date_removed$collection_date,start=1,end=4) < 2020, c(1,2,3,4)])
project_list = as.data.frame(unique(metadata_readin$seq_run))
colnames(project_list) = "projects"

# Read in the Re-run file and strip off info in parens in the first run column
rerun_readin = read.delim(opt$rerun_file, header = TRUE)
rerun_readin$first_run = gsub("\\s*\\([^\\)]+\\)","",as.character(rerun_readin$first_run))

# Read in the completed tab 
completed_readin = read.delim(opt$completed_file, header = TRUE)
completed_accessions = completed_readin[completed_readin$accession_id %in% metadata_50_cov$accession_id,]

# Find if samples are part of the rerun set
rerun_check = rerun_readin[rerun_readin$accession_id %in% metadata_readin$accession_id,c(1,2,3,4)]

# Printing warnings and files for checks
# missing colsetlection date
if(length(missing_collection_date >0)){
cat("\nWARNING: the samples below are missing a collection date:\n")
print(missing_collection_date, row.names = FALSE)
cat('\n')
write.table(missing_collection_date, paste("samples_missing_collection_date_", date, ".tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = '\t')
}

# wrong collection date
if(length(collection_date_wrong > 0)){
  cat("\nWARNING: the samples below have a collection date before 2020:\n")
  print(collection_date_wrong, row.names = FALSE)
  cat('\n')
  write.table(collection_date_wrong, paste("samples_wrong_collection_date_", date, ".tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = '\t')
}

cat("Checking and deleting metadata for samples that have been re-run after projects in this batch\n\n")
# rerun samples
if(length(rerun_check > 0)){
  cat("\nWARNING: the samples below have been submitted under multiple projects.  If you are processing the project in 'first_run' these samples will be removed from the dataset:\n")
  print(rerun_check[,c(1,2,3)], row.names = FALSE)
  cat('\n')
  write.table(rerun_check, paste("samples_rerun_", date, ".tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = '\t')
  
  # Removing samples if this is the first round of sequencing
  for (i in 1:nrow(rerun_check)){
    metadata_50_cov = metadata_50_cov[!(metadata_50_cov$accession_id == rerun_check[i,1] & metadata_50_cov$seq_run != rerun_check[i,3]),] 
  }
}

# already submitted
if(length(completed_accessions > 0)){
  cat("\nWARNING: the samples below have already been submitted and are being deleted from the dataset:\n")
  completed_accessions_print_out = completed_accessions[,c(1,2,6,12,13,14,15,16)]
  print(completed_accessions_print_out, row.names = FALSE)
  cat('\n')
  # Removing samples if this is the first round of sequencing
  for (i in 1:nrow(completed_accessions)){
    metadata_50_cov = metadata_50_cov[!(metadata_50_cov$accession_id == completed_accessions[i,2]),] 
  }
  write.table(completed_accessions, paste("samples_already_submitted_", date, ".tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = '\t')
}


# Creating an empty matrix to fill in with metadata in GISAID format
gisad_metadata = as.data.frame(matrix("",ncol=29,nrow=nrow(metadata_50_cov)))
colnames(gisad_metadata) = c("submitter","fn","covv_virus_name","covv_type","covv_passage","covv_collection_date","covv_location","covv_add_location","covv_host","covv_add_host_info","covv_gender","covv_patient_age","covv_patient_status","covv_specimen","covv_outbreak","covv_last_vaccinated","covv_treatment","covv_seq_technology","covv_assembly_method","covv_coverage","covv_orig_lab","covv_orig_lab_addr","covv_provider_sample_id","covv_subm_lab","covv_subm_lab_addr","covv_subm_sample_id","covv_authors","covv_comment","comment_type")

rename_fasta = as.data.frame(matrix("",ncol=3,nrow=nrow(metadata_50_cov)))
colnames(rename_fasta) = c("accesion_id", "gisaid_name", "proj_num")

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
  rename_fasta[i,3] = metadata_50_cov[i,5]
}

write.csv(gisad_metadata, file = paste("gisaid_submission", date, "metadata.csv",sep="_"), row.names = FALSE, quote = TRUE)
write.csv(rename_fasta, "fasta_rename_accession_to_gisaid_id.csv", row.names = FALSE, quote = FALSE)
write.table(rename_fasta, "fasta_rename_accession_to_gisaid_id.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

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
  cat("ERROR: you didn't specify all varialbles. Please supply all required options. Type Rscript gisaid_ncbi_biosample_metadata_formatting.R -h for list of required options\n\n") # print error messages
}

write.table(ncbi_biosample_metadata, file=paste("ncbi_biosample_submission", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = FALSE, sep = '\t')
