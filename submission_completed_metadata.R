# Mandy Waters
# July 2021
# Script to pull in information needed to add content to COMPLETED_submissions_metadata tab in SUBMISSIONS_and_compiled_metadata.xlsx
# Usage (short) = ./submission_completed_metadata.R -p /path/to/metadata/files -s submitter_name
# Usage (long) = ./submission_completed_metadata.R --path /path/to/metadata/files --submitter_name submitter_name
# Usage notes: only strings with spaces need quotes.  It will throw an error if you don't supply the path and sequence type name or if the sequence type is misspelled.
# Usage example = ./submission_completed_metadata.R -p /mandy_waters/metadata -t "Oxford Nanopore GridION"

# load packages, but don't print the messages to terminal
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(tibble))

# set argument / options
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-p", "--path"), action="store", default=NA, type='character',
              help="path to set working directory with metadata files"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Make the program not be verbose."),
  make_option(c("-s", "--submitter_name"), action="store", default=NA, type='character',
              help="submitter name")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$v) {
  # print the variables to the terminal
  cat("\nworking directory:\n")
  cat(opt$path)
  cat("\n\nsubmitter name:\n")
  cat(opt$submitter_name)
  cat("\n\n")
}

# if all options are filled in then fill in the metadata otherwise error out
if(!is.na(opt$path) & !is.na(opt$submitter_name)){
  
  # set working directory
  setwd(opt$path)
  date = Sys.Date()
    
    # Reading in necessary metadata
    ncbi_in = ldply(list.files(pattern = "metadata-.*-processed-ok.tsv"), read.delim, header=TRUE)
    ncbi_split = ncbi_in %>% separate(filename, c("accession_short",NA,NA), sep = "[.]", remove = FALSE)
    gisaid_in = ldply(list.files(pattern = "gisaid_hcov-19_.*.tsv"), read.delim, header = TRUE)
      if(ncbi_split[1,12] == "OXFORD_NANOPORE"){
      gisaid_split = gisaid_in %>% separate(Virus.name, c(NA,NA,"accession_short", "year"), sep = "/", remove = FALSE)
      }else{gisaid_split = gisaid_in %>% separate(Virus.name, c(NA,NA,"accession_short", "year"), sep = "/", remove = FALSE)
      gisaid_split$accession_short = gsub('CO-CDPHE-','',gisaid_split$accession_short)
      }
    gisaid_filtered = gisaid_split[gisaid_split$accession_short %in% ncbi_split$accession_short,]
    ncbi_gisaid_merge = merge(x=ncbi_split, y=gisaid_filtered[,c(1:5,15:16)], by = "accession_short", all.x=TRUE)
# renaming the file if NCBI sends it as accessions.txt    
    if (file.exists("accessions.txt")) {
      file.rename("accessions.txt", "seqids.txt")
    } else {
      cat("GenBank file has expected file name\n")
    }
    genbank = read.delim("seqids.txt", header=FALSE)
    genbank_split = genbank %>% separate(V1, c("sub","Virus.name"), sep = " ", remove = FALSE)
    ncbi_gisaid_genbank_merge = merge(x=ncbi_gisaid_merge, y=genbank_split[,3:4], by = "Virus.name", all.x = TRUE)
    if(ncbi_split[1,12] == "OXFORD_NANOPORE"){
      ncbi_gisaid_genbank_merge_split = ncbi_gisaid_genbank_merge %>% separate(accession_short, c(NA,NA,"cdphe_accesssion"), sep = "-", remove = FALSE)
    }else{ncbi_gisaid_genbank_merge_split = ncbi_gisaid_genbank_merge
    ncbi_gisaid_genbank_merge_split$cdphe_accesssion = ncbi_gisaid_genbank_merge$accession_short
    ncbi_gisaid_genbank_merge_split$accession_short = paste("CO-CDPHE-", ncbi_gisaid_genbank_merge$accession_short, sep = "")
    }

    # Read in results file to get seq_run number
    metadata_readin = ldply(list.files(pattern = "results_.*.tsv"), read.delim, header=TRUE, na.strings = c("","sample failed assembly"), check.names = FALSE)
    cols_to_keep = as.data.frame(c("accession_id","seq_run"))
    colnames(cols_to_keep) = "col_names"
    metadata_subset_cols = metadata_readin[,colnames(metadata_readin) %in% cols_to_keep$col_names]
    
    col_order <- c("cdphe_accesssion", "Virus.name", "accession_short", "title", "platform", "Collection.date", "Lineage", "Clade", "bioproject_accession", "biosample_accession", "accession", "V2", "Accession.ID")
    final_metadata <- ncbi_gisaid_genbank_merge_split[, col_order]
    final_metadata = final_metadata %>% add_column(submitter = opt$s, .before="cdphe_accesssion")
    final_metadata = final_metadata %>% add_column(isolate = "patient isolate", .before="Collection.date")
    final_metadata = final_metadata %>% add_column(submission_data = date, .after="Accession.ID")
    colnames(final_metadata) = c("submitter","accession_id","Virus name","isolate/sample_name","sample_title","instrument_model","Isolation Source","Collection date","Lineage","Clade","BioProject","BioSample","SRA","GenBank","GISAID","Submission Date")
    
    missing_genbank = final_metadata[is.na(final_metadata$GenBank),]
    missing_genbank_merged = merge(x=missing_genbank, y=metadata_subset_cols, by="accession_id", all.x = TRUE)
    final_metadata = final_metadata[!is.na(final_metadata$GenBank),]
    
    if(length(missing_genbank >0)){
      cat('\nSome of the samples were rejected by GenBank.  Check ONGOING_genabnk_missing_date_metadata.tsv for list. \n')
      missing_genbank_merged = subset(missing_genbank_merged, select=-`Submission Date`)
      missing_genbank_merged = missing_genbank_merged %>% relocate(GISAID, .after = Clade)
      missing_genbank_merged = missing_genbank_merged %>% relocate(seq_run, .after = submitter)
      write.table(missing_genbank_merged, paste("ONGOING_genbank_missing", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = TRUE, sep = '\t')
    }
    
} else {
  cat("ERROR: you didn't specify all varialbles. Please supply all values\n\n") # print error messages
}

write.table(final_metadata, file=paste("COMPLETED", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = TRUE, sep = '\t')

