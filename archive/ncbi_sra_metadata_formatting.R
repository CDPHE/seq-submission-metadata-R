# Mandy Waters
# July 2021
# Script to read in biosample data from NCBI and write out SRA metadata format
# Usage (short) = ./ncbi_sra_metadata_formatting.R -p /path/to/working/dir -t "sequencing platform"
# Usage (long) = ./ncbi_sra_metadata_formatting.R -path /path/to/working/dir -seq_type "sequencing platform"
# Usage notes: only strings with spaces need quotes.  It will throw an error if you don't supply the path and sequence type name or if the sequence type is misspelled.
# Usage example = ./ncbi_sra_metadata_formatting.R -p /mandy_waters/biosample_info -t "Oxford Nanopore GridION"

# load packages, but don't print the messages to terminal
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))

# set argument / options
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-p", "--path"), action="store", default=NA, type='character',
              help="path to set working directory"),
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
  cat("\nworking directory:\n")
  cat(opt$path)
  cat("\n\nsequencing type:\n")
  cat(opt$seq_type)
  cat("\n\n\n\n")
}

# if all options are filled in and sequencing type is spelled correctly then fill in the metadata otherwise error out
if(!is.na(opt$path) & !is.na(opt$seq_type)){
  if(opt$seq_type == "Illumina MiSeq" | opt$seq_type == "Illumina NextSeq" | opt$seq_type == "Oxford Nanopore GridION"){

# set working directory and get date
setwd(opt$path)
date = Sys.Date()

# Reading in biosample data and renaming headers
sra = read.delim("BioSampleObjects.txt", header = TRUE)
colnames(sra) = c("Accession","Sample_Name","SPUID","Organism","Tax_ID","Isolate","BioProject")
sra = sra %>% dplyr::mutate(new_id = str_extract(Sample_Name, "[^-]+$"))

# pulling primer information
if (file.exists("filtered_results_subset_metadata.tsv")) {
  filtered_results = read.delim("filtered_results_subset_metadata.tsv", header = TRUE)
    if("primer_set" %in% colnames(filtered_results)){
      primer_name = unique(filtered_results$primer_set)
    }else{
     primer_name = "Artic V3"
    }
} else {
  primer_name = "Artic V3"
}

# Creating an empty matrix to fill in with metadata 
sra_metadata = as.data.frame(matrix("",ncol=14,nrow=nrow(sra)))
colnames(sra_metadata) = c("bioproject_accession","biosample_accession","library_ID","title","library_strategy","library_source","library_selection","library_layout","platform","instrument_model","design_description","filetype","filename","filename2")

# Filling in columns whose values do not change or will be put in with options
sra_metadata$bioproject_accession = sra$BioProject
sra_metadata$biosample_accession = sra$Accession
sra_metadata$title = "PCR tiled amplicon WGS of SARS-CoV-2"
sra_metadata$library_strategy = "AMPLICON"
sra_metadata$library_source = "VIRAL RNA"
sra_metadata$library_selection = "PCR"
sra_metadata$design_description = paste("Whole genome sequencing of SARS-CoV-2 using the", primer_name, "protocol", sep = " ")
sra_metadata$filetype = "fastq"

primer = 1

# Filling in data that has unique info for each row
for (i in 1:nrow(sra_metadata)){
  sra_metadata[i,3] = primer # library_ID
  
  if (opt$seq_type == "Illumina MiSeq"){
    sra_metadata[i,8] = "paired" # library_layout
    sra_metadata[i,9] = "ILLUMINA" # platform
    sra_metadata[i,10] = "Illumina MiSeq" # instrument model
    sra_metadata[i,13] = paste(sra[i,8], "R1_001.fastq.gz", sep = "_") # filename
    sra_metadata[i,14] = paste(sra[i,8,], "R2_001.fastq.gz", sep = "_") # filename2
  }else if (opt$seq_type == "Illumina NextSeq"){
    sra_metadata[i,8] = "single" # library_layout
    sra_metadata[i,9] = "ILLUMINA" # platform
    sra_metadata[i,10] = "NextSeq 550" # instrument model
    sra_metadata[i,13] =  paste(sra[i,8],"fastq.gz", sep =".") # filename
    sra_metadata[i,14] = ""
  }else if(opt$seq_type == "Oxford Nanopore GridION"){
    sra_metadata[i,8] = "single" # library_layout
    sra_metadata[i,9] = "OXFORD_NANOPORE" # platform
    sra_metadata[i,10] = "GridION" # instrument model
    sra_metadata[i,13] =  paste(sra[i,2],"fastq.gz", sep =".") # filename
    sra_metadata[i,14] = ""
  }
  primer = primer + 1
}

  } else {
    cat("ERROR: check the spelling of your sequence type input: sequence type should match 'Illumina MiSeq', 'Illumina NextSeq', or 'Oxford Nanopore GridION'\n\n")
    }
} else {
  cat("ERROR: you didn't specify all varialbles. Please supply all values\n\n") # print error messages
}

# write out SRA metadata file
write.table(sra_metadata, file=paste("ncbi_sra_submission", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = TRUE, sep = '\t')

